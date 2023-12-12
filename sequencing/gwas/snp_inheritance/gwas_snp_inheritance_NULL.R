suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(tidyverse)
})
vcf <- 'sequencing/gwas/filtered.vcf.gz'

## Define helper functions
prep_gwas_table <- function(trait, produce_plots = TRUE, pval_column, fdr.level = 0.05) {
  require(qvalue)
  gwas <- data.table::fread(glue::glue('sequencing/gwas/gemma_output/{trait}.assoc.txt'), data.table = FALSE)

  # make qqplot
  if (produce_plots) {
    opa <- par(mfrow = c(2, 1))
    hist(gwas[[pval_column]])
    GWASTools::qqPlot(gwas[[pval_column]], thinThreshold = 4)
    par(opa)
  }

  q <- qvalue(gwas[[pval_column]], fdr.level = fdr.level)
  summary(q)
  gwas$qvalue <- q$qvalues
  gwas$significant <- q$significant
  gwas$chr2 <- fct_relabel(
    factor(gwas$chr),
    \(l) case_when(l == 'NC_024238.1' ~ 'MT',
                   str_starts(l, 'NC') ~ (str_sub(l, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                   TRUE ~ 'Un'
    ))
  scaff_sizes <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
    dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
    dplyr::select(-V4, -V5)

  gwas %>%
    left_join(scaff_sizes, join_by(chr)) %>%
    mutate(cum_ps = ps + cum_scaff_start)
}

source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('sequencing/gwas/peak_viz/panel_functions/make_linkage_panel.R')

get_linkage <- function(geno, workers = 4, min_MAF = 0.1) {
  require(future.apply)

  sampling <- get_sampling_structure()

  suppressMessages(source('quant_gen/prepare_pedigrees.R'))
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  # #source('selection_decisions/compile_decisions.R')
  # groups <- selection %>%
  #   dplyr::select(fish_id, replicate, selection) %>%
  #   mutate(fish_id = tolower(fish_id)) %>%
  #   mutate(across(everything(), as.factor))

  dat <- geno %>%
    as.data.frame() %>%
    mutate(MAF = ifelse(AF1 < 0.5, AF1, 1 - AF1)) %>%
    filter(MAF >= min_MAF) %>%
    dplyr::select(chr = seqnames, start, end, alt, sampleNames, GT) %>%
    inner_join(sampling, join_by(sampleNames == sample_name)) %>%
    #inner_join(groups, join_by(fish_id)) %>%
    add_patriline(ped_df) %>%
    mutate(
      dose = case_when(GT == '0/0' ~ 0, GT == '0/1' ~ 1, GT == '1/0' ~ 1, GT == '1/1' ~ 2, TRUE ~ NA_integer_),
      fish_idX = fish_id
    )

  #plan(multisession, workers = workers)
  link <- dat %>%
    group_by(chr, start, end, alt) %>%
    group_nest()
  X <- link$data
  link$out <- future.apply::future_lapply(
    X,
    \(x) safe_snp_sex_linkage_AIC(x, A_small, X_small),
    future.seed = TRUE
  )

  link <- link %>%
    unnest(cols = out) %>%
    dplyr::select(-data)
  plan(sequential)

  return(link)
}

# Load GWAS results
g <- prep_gwas_table('car_PIE', FALSE, 'p_SHet', fdr.level = 0.01) %>%
  filter(!significant, !(chr2 %in% c(12, 'Un', 'MT'))) %>%
  slice_sample(n = 1000)

regions <- GRanges(
  seqnames = g$chr,
  ranges = IRanges(start = g$ps, end = g$ps)
)
geno <- readVcfAsVRanges(vcf, param = ScanVcfParam(which = regions)) %>%
  as.data.frame() %>%
  filter(
    !(sampleNames %in% c('NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355'))
  )
cat('\nLoaded genotypes...\n')
link <- get_linkage(geno, workers = 4, min_MAF = 0.1)

ggplot(link, aes(AIC_weight)) + geom_boxplot() + facet_grid(rows = vars(source))
link %>% group_by(chr, start, alt) %>% slice_max(AIC_weight) %>% ungroup() %>% dplyr::count(source)
  #ggplot(aes(source)) + geom_bar()

data.table::fwrite(link, glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance_NULL.csv'))


# Load GWAS results
g12 <- prep_gwas_table('car_PIE', FALSE, 'p_SHet', fdr.level = 0.01) %>%
  filter(!significant, chr2 == 12) %>%
  slice_sample(n = 1000)

regions12 <- GRanges(
  seqnames = g12$chr,
  ranges = IRanges(start = g12$ps, end = g12$ps)
)
geno12 <- readVcfAsVRanges(vcf, param = ScanVcfParam(which = regions12)) %>%
  as.data.frame() %>%
  filter(
    !(sampleNames %in% c('NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355'))
  )
cat('\nLoaded genotypes...\n')
link12 <- get_linkage(geno12, workers = 4, min_MAF = 0.1)

ggplot(link12, aes(AIC_weight)) + geom_histogram() + facet_grid(rows = vars(source))
link12 %>% group_by(chr, start, alt) %>% slice_max(AIC_weight) %>% ungroup() %>% dplyr::count(source)
#ggplot(aes(source)) + geom_bar()

data.table::fwrite(link12, glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance_NULL_12.csv'))
