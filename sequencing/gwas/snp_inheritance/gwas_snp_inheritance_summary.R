library(tidyverse)

make_summary <- function(file) {
  snp_inheritance <- data.table::fread(file) %>%
    group_by(chr, start, end, alt) %>%
    slice_max(AIC_weight, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(source = if_else(AIC_weight < 0.8, 'Not called', source))

  cat('Whole genome:\n')
  snp_inheritance %>%
    dplyr::count(source) %>%
    mutate(total = sum(n), frac = n / total) %>%
    print()

  cat('Autosomes only:\n')
  snp_inheritance %>%
    filter(chr != 'NC_024342.1', str_starts(chr, 'NW', negate = TRUE)) %>%
    dplyr::count(source) %>%
    mutate(total = sum(n), frac = n / total) %>%
    print()
  return(invisible(NULL))
}

make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/car_PIE.csv')
# Autosomes: (66 + 1523) = 1589 out of 4841 called variants or, 32.8%

make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/mel_PIE.csv')
# Autosomes: (98 + 1223) = 1321 out of 3691 called variants or, 35.8%

# compare to the NULL set
make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance_NULL.csv')
# Autosomes: (6 + 3) = 9 out of 947 called variants or, 0.95%

# and the random set from LG12:
make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance_NULL_12.csv')
# LG12: (305 + 418) = 723 out of 862 called variants or, 83.9%
