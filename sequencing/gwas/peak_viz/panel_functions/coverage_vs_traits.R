coverage_vs_traits <- function(region_cov, sig_range) {

  sig_range_coverage <- semi_join(
    region_cov,
    data.frame(start = sig_range[1], end = sig_range[2]),
    join_by(overlaps(window_start, window_end, start, end))
  ) %>%
    group_by(sample_name) %>%
    summarise(mean_rel_coverage = mean(rel_coverage))

  GRM <- get_GRM()

  traits <- c('pa_car_1', 'pa_car_2', 'pa_car_3', 'size_car_4', 'pa_car_5', 'pa_car_6', 'pa_car_7')

  data_list <- map(
    traits,
    \(p) get_phenotypes(p) %>%
      inner_join(sig_range_coverage, join_by(sample_name)) %>%
      drop_na()
  ) %>% setNames(traits)

  p_vals <- map_dbl(
    data_list,
    \(d) {
      m1 <- mmer(phenotype_value ~ 1 + mean_rel_coverage, random = ~vsr(sample_name, Gu = GRM), data = d)
      m2 <- mmer(phenotype_value ~ 1, random = ~vsr(sample_name, Gu = GRM), data = d)
      a <- anova(m1, m2)
      # anova.mmer only gives the p-value as a string (sigh), so calculate p manually:

      pchisq(abs(-2 * diff(a$loLik)), df = -diff(a$Df), lower.tail = FALSE)
    }
  ) %>% stack() %>% setNames(c('p', 'trait'))

  p1 <- ggplot(
    bind_rows(data_list, .id = 'trait') %>% filter(str_starts(trait, 'pa_')),
    aes(factor(phenotype_value), mean_rel_coverage)
  ) +
    geom_violin(col = 'grey20', fill = 'darkorange', alpha = 0.5, linewidth = 0.3) +
    stat_summary(fun.data = 'mean_cl_boot') +
    geom_text(
      aes(y = Inf, x = 1.5, label = format.pval(p, digits = 1)),
      p_vals %>% filter(str_starts(trait, 'pa_')),
      vjust = 1.5
    ) +
    facet_grid(~trait) +
    scale_fill_manual(values = c('0' = NA, '1' = 'darkorange'), guide = 'none', na.value = NA) +
    labs(y = 'mean relative coverage', x = 'trait present') +
    theme(axis.line = element_blank(), panel.grid.major.y = element_line())

  center_x <- mean(range(sig_range_coverage$mean_rel_coverage, na.rm = TRUE))
  p2 <- ggplot(
    bind_rows(data_list, .id = 'trait') %>% filter(str_starts(trait, 'size_')),
    aes(mean_rel_coverage, phenotype_value)
  ) +
    geom_point() + geom_smooth(method = 'lm', se = FALSE, color = 'darkorange') +
    geom_text(
      aes(y = Inf, x = center_x, label = format.pval(p, digits = 1)),
      p_vals %>% filter(str_starts(trait, 'size_')),
      vjust = 1.5
    ) +
    facet_grid(~trait) +
    labs(y = 'ornament size (z-score)', x = 'mean relative coverage') +
    theme(axis.line = element_blank())

  p1 + p2 + plot_layout(widths = c(6, 1))
}
