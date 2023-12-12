combine_gwas_sHet <- function(trait, workers = 42, overwrite = FALSE) {
  require(tidyverse)
  require(future)
  source('sequencing/gwas/CPASSOC/FunctionSet.R')

  outfile <- paste0('sequencing/gwas/gemma_output/', trait, '.assoc.txt')
  if (!overwrite && file.exists(outfile)) return(invisible(NULL))

  # Find the correct files
  if (str_detect(trait, 'Yhap', negate = TRUE)) {
    pattern <- paste0(trait, '_\\d+\\.assoc\\.txt$')
  } else {
    pattern <- paste0(str_remove(trait, '_Yhap'), '_\\d+_Yhap\\.assoc\\.txt$')
  }
  files <- list.files('sequencing/gwas/gemma_output', pattern, full.names = TRUE)
  assoc <- combine_assoc(files)

  # Obtain the genetic correlations
  cor_mat <- assoc %>%
    # filter out significant SNPs as recommended by CPASSOC
    filter(if_all(starts_with('Z'), \(x) abs(x) < 1.96)) %>%
    # remove MT, sex chromosomes and unplaced scaffolds as LD may be very high
    filter(!chr %in% c('NC_024238.1', 'NC_024342.1'), str_starts(chr, 'NW', negate = TRUE)) %>%
    group_by(chr) %>%
    # Reduce LD by taking every 100th SNP
    dplyr::slice(seq(1, n(), 100)) %>%
    ungroup() %>%
    dplyr::select(starts_with('Z')) %>%
    as.matrix() %>%
    cor()

  print(cor_mat)

  if (length(files) != nrow(cor_mat)) {
    stop('Files and correlations have mismatched dimensions.')
  }

  Z <- dplyr::select(assoc, starts_with('Z'))
  message('NOTE: Assuming a sample size of 297 individuals.')
  S <- rep(297, ncol(Z))
  B <- cor_mat

  cat('Estimating gamma...\n')
  para <- EstimateGamma(N = 1e6, SampleSize = S, CorrMatrix = B, correct = 1, isAllpossible = T )

  cat('Calculating SHet statistics...\n')
  old_plan <- plan(multisession, workers = workers)
  x <- SHet(X = Z, SampleSize = S, CorrMatrix = B, correct = 1, isAllpossible = T)
  plan(old_plan)

  cat('Calculating p-values...\n')
  assoc$p_SHet <- pgamma(q = x - para[3], shape = para[1], scale = para[2], lower.tail = F)
  assoc$log_p_SHet <- pgamma(q = x - para[3], shape = para[1], scale = para[2], lower.tail = F, log.p = TRUE)

  data.table::fwrite(assoc, outfile)

  return(invisible(NULL))
}
purrr::walk(
  c('car_PIE', 'mel_PIE'),
  combine_gwas_sHet, overwrite = TRUE
)

