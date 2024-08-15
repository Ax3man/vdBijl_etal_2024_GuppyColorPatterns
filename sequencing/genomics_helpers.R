prep_gwas_table <- function(
    trait, produce_plots = FALSE, pval_column, fdr_level = 0.05, calc_q = TRUE, reference = 'female',
    # whether to apply the Fraser remapping of LG12 coordinates, don't use if e.g. matching to vcf
    adjust_chr12 = TRUE
) {
  require(qvalue)
  if (reference == 'female') {
    gwas <- data.table::fread(glue::glue('sequencing/gwas/gemma_output/{trait}.assoc.txt'), data.table = FALSE)
    gwas$chr2 <- fct_relabel(
      factor(gwas$chr),
      \(l) case_when(l == 'NC_024238.1' ~ 'MT',
                     str_starts(l, 'NC') ~ (str_sub(l, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                     TRUE ~ 'Un'
      ))
    scaff_sizes <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      dplyr::select(-V4, -V5)
  } else {
    if (reference != 'male') stop('`reference` should be "male" or "female".')
    source('sequencing/male_reference_gwas/update_chr12_liftover.R')

    gwas <- data.table::fread(glue::glue('sequencing/male_reference_gwas/gemma_output/{trait}.assoc.txt'), data.table = FALSE)
    gwas$chr2 <- fct_relabel(factor(gwas$chr), \(l) ifelse(str_starts(l, '0'), 'Un', l))
    # Note that something in the pipeline (plink?) has assumed that chr23 is the X (assumed human).
    gwas$chr2 <- factor(gwas$chr2, levels = c(1:22, 'X', 'Un'), labels = c(1:23, 'Un'))
    gwas$chr <- ifelse(gwas$chr == 'X', '23', gwas$chr)

    gwas <- mutate(gwas, chr = fct_relevel(factor(chr), as.character(1:23)))

    # for the male reference, we update the chr12 coordinates, following Fraser et al. and relevel the chr
    if (adjust_chr12) {
      gwas <- mutate(gwas, ps = ifelse(chr == '12', update_chr12(ps), ps))
    }

    scaff_sizes <- data.table::fread('sequencing/male_reference_gwas/male_reference.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      dplyr::select(-V4, -V5) |>
      mutate(
        chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr),
        chr = fct_relevel(factor(chr), as.character(1:23))
      ) |>
      arrange(chr) |>
      mutate(cum_scaff_start = cumsum(lag(scaff_len, default = 0)))
  }

  # make qqplot
  if (produce_plots) {
    opa <- par(mfrow = c(2, 1))
    hist(gwas[[pval_column]])
    GWASTools::qqPlot(gwas[[pval_column]], thinThreshold = 4)
    par(opa)
  }

  if (calc_q) {
    q <- qvalue(gwas[[pval_column]], fdr.level = fdr_level)
    summary(q)
    gwas$qvalue <- q$qvalues
    gwas$significant <- q$significant
  }

  gwas %>%
    left_join(scaff_sizes, join_by(chr)) %>%
    arrange(chr2) |>
    mutate(cum_ps = ps + cum_scaff_start)
}

read_kmer_blast <- function(trait, parse_pvals = FALSE) {
  file <- glue::glue('sequencing/kmer_gwas/kmer_blast2ref/{trait}.tsv')
  d <- data.table::fread(file) %>% as_tibble()
  colnames(d) <- c(
    'kmer', 'chr', 'perc_ident', 'align_length', 'nr_mismatch', 'nr_gaps', 'kmer_start', 'kmer_end',
    'start', 'end', 'e_value', 'bit_score'
  )
  if (parse_pvals) {
    d <- d %>%
      separate(kmer, into = c('kmer', 'p_value'), sep = '_') %>%
      mutate(p_value = parse_number(p_value))
  }

  scaff_sizes <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
    dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
    dplyr::select(-V4, -V5) %>%
    mutate(chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                            str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                            TRUE ~ 'Un'))

  d %>%
    left_join(scaff_sizes, join_by(chr)) %>%
    mutate(cum_start = start + cum_scaff_start)
}

scaff_labs <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
  dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
  mutate(
    scaff_mid = cum_scaff_start + 0.5 * scaff_len,
    chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                     str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                     TRUE ~ 'Un')
  ) %>%
  filter(!(chr2 %in% c('Un', 'MT'))) %>%
  dplyr::select(-V4, -V5)

male_scaff_labs <- data.table::fread('sequencing/male_reference_gwas/male_reference.fna.fai') %>%
  dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
  dplyr::select(-V4, -V5) |>
  mutate(
    chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr),
    chr = fct_relevel(factor(chr), as.character(1:23))
  ) |>
  arrange(chr) |>
  mutate(
    cum_scaff_start = cumsum(lag(scaff_len, default = 0)),
    scaff_mid = cum_scaff_start + 0.5 * scaff_len,
    chr2 = chr
  ) |>
  filter(!str_starts(chr, '0'))

reduce_peaks <- function(d, min_dist = 1e4, pval_column = 'p_lrt') {
  # sort SNPs by significance, and assign windows around them
  d <- arrange(d, !!sym(pval_column)) %>%
    mutate(.start = ps - min_dist, .end = ps + min_dist)

  out <- dplyr::slice(d, 1)
  while (TRUE) {
    a <- anti_join(dplyr::slice(d, -1), out, join_by(chr, between(x$ps, y$.start, y$.end)))
    if (nrow(a) == nrow(d)) break
    d <- a
    out <- bind_rows(out, dplyr::slice(d, 1))
  }
  out <- dplyr::select(out, -.start, -.end)
  return(out)
}
