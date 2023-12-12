make_snp_depth_panel <- function(geno, viz_range, pheno = NULL) {
  mean_cov <- get_mean_coverage()

  d <- inner_join(
    as.data.frame(geno),
    mean_cov,
    join_by(sampleNames == sample_name)
  ) %>%
    mutate(rel_depth = totalDepth / mean_coverage)

  if (is.null(pheno)) {
    ggplot(d, aes(start, rel_depth)) +
      geom_hline(yintercept = 1, lty = 2) +
      stat_summary(geom = 'point') +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      labs(y = 'median SNP cov')
  } else {
    inner_join(
      d,
      pheno,
      join_by(sampleNames == sample_name)
    ) %>%
      ggplot(aes(start, rel_depth, color = factor(phenotype_value))) +
      geom_hline(yintercept = 1, lty = 2) +
      stat_summary(geom = 'line', fun = mean) +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      scale_color_manual(
        values = c('0' = 'black', '1' = 'darkorange'),
        name = 'phenotype',
        labels = c('absent', 'present')
      ) +
      labs(y = 'median SNP cov')
  }
}

make_region_depth_panel <- function(
    region_cov, viz_range, pheno = NULL, sig_range = NULL,
    alpha = 0.1, linewidth = 0.5
) {
  if (is.null(pheno)) {
    p <- ggplot(region_cov, aes((window_start + window_end) / 2, rel_coverage, group = sample_name)) +
      geom_line(alpha = alpha, linewidth = linewidth) +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      labs(y = 'relative coverage')
  } else {
    p <- region_cov %>%
      inner_join(pheno, join_by(sample_name)) %>%
      ggplot(aes(
        (window_start + window_end) / 2, rel_coverage,
        group = sample_name, color = factor(phenotype_value))) +
      geom_line(alpha = alpha, linewidth = linewidth) +
      stat_summary(aes(group = phenotype_value), geom = 'line', linewidth = 1, fun = median) +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      scale_color_manual(
        values = c('0' = 'black', '1' = 'darkorange'),
        name = 'phenotype',
        labels = c('absent', 'present')
      ) +
      labs(y = 'relative coverage')
  }
  if (!is.null(sig_range)) {
    p <- p + annotate(
      geom = 'segment', x = sig_range[1], xend = sig_range[2], y = Inf, yend = Inf,
      color = 'darkorange', linewidth = 2
    )
  }
  return(p)
}

make_allelic_depth_panel <- function(geno, chr, pos) {
  geno %>%
    dplyr::filter(.data$seqnames == .env$chr, .data$start == .env$pos) %>%
    pivot_longer(refDepth:altDepth, names_to = 'allele', values_to = 'allelic_depth') %>%
    left_join(get_mean_coverage(), join_by(sampleNames == sample_name)) %>%
    mutate(allelic_depth = allelic_depth / mean_coverage) %>%
    ggplot(aes(y = allelic_depth * 2, x = GT, fill = factor(allele, levels = c('altDepth', 'refDepth')))) +
    geom_bar(stat = 'summary', fun = 'median', position = 'stack') +
    scale_y_continuous(
      name = 'est. allele count',
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_manual(values = c('lightblue', 'grey20'), labels = c('alt', 'ref'), name = NULL) +
    labs(x = 'genotype')
}
