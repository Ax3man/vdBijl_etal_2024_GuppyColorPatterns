snp_sex_linkage_AIC <- function(dat, Amat, Xmat) {
  suppressPackageStartupMessages(require(sommer))
  h2 <- function(m) {
    vc <- summary(m)$varcomp
    vc[1, 1] / sum(vc[, 1])
  }

  mA <- mmer(dose ~ 1, ~ vsr(fish_id, Gu = Amat), data = dat, verbose = FALSE, dateWarning = FALSE)
  mX <- mmer(dose ~ 1, ~ vsr(fish_idX, Gu = Xmat), data = dat, verbose = FALSE, dateWarning = FALSE)
  mY <- mmer(dose ~ 1, ~ patriline, data = dat, verbose = FALSE, dateWarning = FALSE)
  AICs <- c(mA$AIC, mX$AIC, mY$AIC)
  h2s <- c(h2(mA), h2(mX), h2(mY))
  data.frame(
    source = c('auto', 'X', 'Y'),
    AIC_weight = MuMIn::Weights(AICs) %>% as.numeric(),
    h2 = h2s
  )
}
safe_snp_sex_linkage_AIC <- \(d, a, x) {
  purrr::safely(snp_sex_linkage_AIC, otherwise = data.frame())(d, a, x)$result
}

get_linkage <- function(geno, trait, workers = 25, min_MAF = 0.1) {
  require(furrr)

  workers <- pmin(parallelly::availableCores(logical = FALSE) - 1, workers)

  sampling <- get_sampling_structure()
  source('quant_gen/prepare_pedigrees.R')
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  # load previous results (only significant variants)
  prev <- data.table::fread(glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv'))

  dat <- geno %>%
    as.data.frame() %>%
    mutate(MAF = ifelse(AF1 < 0.5, AF1, 1 - AF1)) %>%
    filter(MAF >= min_MAF) %>%
    dplyr::select(chr = seqnames, start, end, alt, sampleNames, GT) %>%
    inner_join(sampling, join_by(sampleNames == sample_name)) %>%
    add_patriline(ped_df) %>%
    mutate(
      dose = case_when(GT == '0/0' ~ 0, GT == '0/1' ~ 1, GT == '1/0' ~ 1, GT == '1/1' ~ 2, TRUE ~ NA_integer_),
      fish_idX = fish_id
    )

  done <- semi_join(prev, dat, join_by(chr, start, end, alt))
  to_do <- anti_join(dat, prev, join_by(chr, start, end, alt))

  #plan(multisession, workers = workers)
  link <- to_do %>%
    group_by(chr, start, end, alt) %>%
    group_nest()
  # link$out <- pbapply::pblapply(
  #   link$data,
  #   \(x) safe_snp_sex_linkage_AIC(x, A_small, X_small)
  # )
  link$out <- future_map(
      link$data,
      \(x) safe_snp_sex_linkage_AIC(x, A_small, X_small),
      .options = furrr_options(seed = TRUE)
  )
  #plan(sequential)

  link <- link %>%
    unnest(cols = out) %>%
    dplyr::select(-data) %>%
    bind_rows(done)

  return(link)
}

make_linkage_panel <- function(link, viz_range) {
  suppressWarnings(
    ggplot(link, aes(start, AIC_weight, color = source)) +
      geom_point() +
      geom_smooth(method = 'gam', se = FALSE, method.args = list(family = 'binomial')) +
      scale_color_manual(values = c('auto' = 1, 'X' = 'firebrick', 'Y' = 'blue')) +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      labs() +
      coord_cartesian(ylim = c(0, 1))
  )
}


