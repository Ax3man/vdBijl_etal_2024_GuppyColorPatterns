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

get_linkage <- function(geno, trait, workers = 25, min_MAF = 0.1, reference = 'female', use_previous = TRUE) {
  suppressPackageStartupMessages(require(furrr))

  workers <- pmin(parallelly::availableCores(logical = FALSE) - 1, workers)

  cat('Preparing data...\n')

  source('quant_gen/prepare_pedigrees.R')
  sampling <- get_sampling_structure() |> add_patriline(ped_df)
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  # load previous results (only significant variants)
  if (use_previous) {
    if (reference == 'female') {
      prev_file <- glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')
    } else {
      prev_file <- glue::glue('sequencing/male_reference_gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')
    }
    if (file.exists(prev_file)) {
      prev <- data.table::fread(prev_file)
    } else {
      prev <- data.frame(chr = character(), start = integer(), end = integer(), alt = character())
    }
  }

  dat <- geno %>%
    as.data.frame() %>%
    mutate(MAF = ifelse(AF1 < 0.5, AF1, 1 - AF1)) %>%
    filter(MAF >= min_MAF) %>%
    dplyr::select(chr = seqnames, start, end, alt, sampleNames, GT) %>%
    inner_join(sampling, join_by(sampleNames == sample_name)) %>%
    mutate(
      dose = case_when(GT == '0/0' ~ 0, GT == '0/1' ~ 1, GT == '1/0' ~ 1, GT == '1/1' ~ 2, TRUE ~ NA_integer_),
      fish_idX = fish_id
    )

  if (use_previous) {
    done <- semi_join(prev, dat, join_by(chr, start, end, alt))
    to_do <- anti_join(dat, prev, join_by(chr, start, end, alt))
  } else {
    to_do <- dat
  }

  cat('Fitting models...\n')
  if (workers > 1) plan(multisession, workers = workers)
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
      .options = furrr_options(
        seed = TRUE,
        globals = c('A_small', 'X_small', 'safe_snp_sex_linkage_AIC', 'snp_sex_linkage_AIC')
      )
  )
  if (workers > 1) plan(sequential)

  cat('Post-processing...\n')
  link <- link %>%
    unnest(cols = out) %>%
    dplyr::select(-data)

  if (use_previous) {
    link <- bind_rows(link, done)
  }
  cat('Done!\n')
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


