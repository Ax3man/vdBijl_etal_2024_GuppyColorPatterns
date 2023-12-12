SNP_vs_traits <- function(
    geno, pos,
    traits = c('pa_car_1', 'pa_car_2', 'pa_car_3', 'size_car_4', 'pa_car_5', 'pa_car_6', 'pa_car_7')
) {
  snp <- filter(geno, start == pos)
  if (n_distinct(snp$end) > 1) snp <- filter(snp, end == min(end))

  GRM <- get_GRM()

  data_list <- map(
    traits,
    \(p) get_phenotypes(p) %>%
      inner_join(snp, join_by(sample_name == sampleNames)) %>%
      mutate(dose = case_when(GT == '0/0' ~ 0, GT == '0/1' ~ 1, GT == '1/1' ~ 2, TRUE ~ NA)) %>%
      dplyr::select(fish_id, sample_name, phenotype_value, dose, GT) %>%
      drop_na()
  ) %>% setNames(traits)

  p_vals <- map_dbl(
    data_list,
    \(d) {
      m1 <- mmer(phenotype_value ~ dose, random = ~vsr(sample_name, Gu = GRM), data = d)
      m2 <- mmer(phenotype_value ~ 1, random = ~vsr(sample_name, Gu = GRM), data = d)
      a <- anova(m1, m2)
      # anova.mmer only gives the p-value as a string (sigh), so calculate p manually:
      pchisq(abs(-2 * diff(a$loLik)), df = -diff(a$Df), lower.tail = FALSE)
    }
  ) %>% stack() %>% setNames(c('p', 'trait'))

  if (any(str_starts(traits, 'pa_'))) {
    p1 <- ggplot() +
      geom_bar(
        aes(GT, fill = factor(phenotype_value)),
        bind_rows(data_list, .id = 'trait') %>% filter(str_starts(trait, 'pa_')),
        position = 'fill', width = 0.95
      ) +
      geom_text(
        aes(y = Inf, x = 1.5, label = format.pval(p, digits = 1)),
        p_vals %>% filter(str_starts(trait, 'pa_')),
        vjust = 1.5
      ) +
      facet_grid(~trait) +
      #coord_cartesian(expand = FALSE) +
      scale_fill_manual(values = c('0' = 'grey', '1' = 'darkorange'), guide = 'none', na.value = NA) +
      scale_y_continuous(expand = expansion(c(0, 0.1))) +
      labs(y = 'proportion with trait', x = 'genotype at top SNP') +
      theme(axis.line = element_blank(), panel.grid.major.y = element_line())
  } else { p1 <- NULL}

    if (any(str_starts(traits, 'size_'))) {
      p2 <- ggplot(
        bind_rows(data_list, .id = 'trait') %>% filter(str_starts(trait, 'size_')),
        aes(GT, phenotype_value)
      ) +
        geom_violin(col = NA, fill = 'darkorange', alpha = 0.5) +
        stat_summary(fun.data = 'mean_cl_boot') +
        geom_text(
          aes(y = Inf, x = 1.5, label = format.pval(p, digits = 1)),
          p_vals %>% filter(str_starts(trait, 'size_')),
          vjust = 1.5
        ) +
        facet_grid(~trait) +
        labs(y = 'ornament size (z-score)', x = NULL) +
        theme(axis.line = element_blank())
    } else { p2 <- NULL}

  return(list(p1, p2))
}

make_ornament_panel <- function(ornament, geno, chr, pos) {
  orn <- get_ornaments()

  geno %>%
    dplyr::filter(.data$seqnames == .env$chr, .data$start == .env$pos) %>%
    full_join(orn, join_by(sampleNames == sample_name)) %>%
    ggplot(aes(GT, fill = .data[[ornament]] == 1)) +
    geom_bar() +
    scale_fill_manual(values = c('TRUE' = 'darkorange', 'FALSE' = 'grey'), name = NULL,
                      labels = c('TRUE' = 'present', 'FALSE' = 'absent')) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'genotype', y = 'count')
}
