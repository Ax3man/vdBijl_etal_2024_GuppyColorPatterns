library(tidyverse)

make_summary <- function(file = NULL, df = NULL) {
  if (!is.null(file)) df <- data.table::fread(file)

  snp_inheritance <- df %>%
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
# Autosomes: (66 + 1523) = 1589 out of 4841 called variants are sex-linked or, 32.8%
#            1523 out of 4841 called variants are Y-linked or, 31.4%
#.             66 out of 4841 called variants are X-linked or, 1.4%

make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/mel_PIE.csv')
# Autosomes: (98 + 1223) = 1321 out of 3691 called variants are sex-linked or, 35.8%
#            1223 out of 3691 called variants are Y-linked or, 33.1%
#              98 out of 3691 called variants are X-linked or, 2.7%

# compare to the NULL set
make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance_NULL.csv')
# (6 + 3) = 9 out of 947 called variants are sex-linked or, 0.95%
#           3 out of 947 called variants are Y-linked or, 0.32%


# and the random set from LG12:
make_summary('sequencing/gwas/snp_inheritance/gwas_snp_inheritance_NULL_12.csv')
# LG12: (305 + 418) = 723 out of 862 called variants are sex-linked or, 83.9%
#       418 out of 723 called variants are Y-linked or, 57.8%
