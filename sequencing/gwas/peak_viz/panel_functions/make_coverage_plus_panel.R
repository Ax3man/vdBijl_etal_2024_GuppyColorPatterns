#trait <- 'pa_car_1'
#chr <- 'NC_024333.1'
#pos <- 593580

coverage_plus <- function(chr, pos, coverage_width = 1e3) {

  sampling <- get_sampling_structure()
  source('quant_gen/prepare_pedigrees.R')
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  ind_cov_10kb <- get_region_coverage(chr, pos, c(pos - coverage_width/2, pos + coverage_width/2)) %>%
    group_by(sample_name) %>%
    summarise(rel_coverage = mean(rel_coverage))

  md <- ind_cov_10kb %>%
    inner_join(get_phenotypes(trait), join_by(sample_name))

  # Load the SNP GRM from GEMMA
  GRM <- get_GRM()

  m <- mmer(
    phenotype_value ~ rel_coverage,
    random = ~vsr(sample_name, Gu = GRM),
    data = md,
  )
  betas <- summary(m)$betas %>% mutate(p_approx = pt(abs(t.value), df = 250, lower.tail = FALSE) * 2)
  h2 <- summary(m)$varcomp[1, 'VarComp'] / sum(summary(m)$varcomp[, 'VarComp'])

  p1 <- ggplot(md, aes(rel_coverage, factor(phenotype_value), fill = factor(phenotype_value))) +
    geom_violin(show.legend = FALSE) +
    ggforce::geom_sina(show.legend = FALSE) +
    scale_fill_manual(values = c('grey80', 'darkorange')) +
    labs(
      y = 'phenotype', x = 'relative coverage',
      subtitle = glue::glue('p = {round(betas[2, "p_approx"], 3)}; h2 = {round(h2, 2)}')
    )

  p2 <- md %>%
    mutate(dose = rel_coverage, fish_idX = fish_id) %>%
    add_patriline(ped_df) %>%
    snp_sex_linkage_AIC() %>%
    ggplot(aes(AIC_weight, 'a', fill = source)) +
    geom_col() +
    scale_fill_manual(values = c('grey40', 'firebrick', 'navy')) +
    coord_cartesian(expand = FALSE) +
    labs(x = 'AIC weight', y = NULL, subtitle = 'Inheritance of coverage') +
    theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())

  p1 + p2
}


coverage_linkage <- function(region_cov, sig_range) {

  sampling <- get_sampling_structure()
  source('quant_gen/prepare_pedigrees.R')
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  sig_range_coverage <- semi_join(
    region_cov,
    data.frame(start = sig_range[1], end = sig_range[2]),
    join_by(overlaps(window_start, window_end, start, end))
  ) %>%
    group_by(sample_name) %>%
    summarise(mean_rel_coverage = mean(rel_coverage)) %>%
    inner_join(sampling, join_by(sample_name))

  sig_range_coverage %>%
    mutate(dose = mean_rel_coverage, fish_idX = fish_id) %>%
    add_patriline(ped_df) %>%
    snp_sex_linkage_AIC(A_small, X_small) %>%
    ggplot(aes('a', AIC_weight, fill = source)) +
    geom_col() +
    scale_fill_manual(values = c('grey40', 'firebrick', 'navy')) +
    coord_cartesian(expand = FALSE) +
    labs(y = 'AIC weight', x = NULL, subtitle = 'Inheritance of coverage') +
    theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())
}
