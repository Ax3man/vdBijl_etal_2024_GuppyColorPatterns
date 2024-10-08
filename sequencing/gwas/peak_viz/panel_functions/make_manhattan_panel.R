make_manhattan_panel <- function(gwas, viz_range) {
  if (n_distinct(gwas$chr) > 1) stop('This function only plots zoomed-in single chr manhattans')

  ggplot(gwas, aes(ps, -log10(p_lrt))) +
    geom_point() +
    scale_x_continuous(
      limits = viz_range,
      name = 'position (kb)',
      label = \(x) scales::label_comma()(x / 1e3)
    ) +
    ylab(expression(-log[10](italic(p[assoc]))))
}

make_manhattan_linkage_panel <- function(gwas, link, viz_range, pval_column, AIC_limit = 0.8) {
  link <- link %>% group_by(start, end) %>% filter(AIC_weight > AIC_limit)
  link$chr <- ifelse(str_starts(link$chr, 'chr'), as.character(parse_number(link$chr)), link$chr)

  left_join(
    gwas,
    link,
    join_by(chr, ps == start, allele1 == alt)
  ) %>%
    mutate(
      source = coalesce(source, 'not called'),
      source = factor(source, c('auto', 'X', 'Y', 'not called'))
    ) %>%
    ggplot(aes(ps, -log10(!!sym(pval_column)), color = source)) +
    geom_point(show.legend = TRUE, size = 0.1) +
    scale_color_manual(
      name = glue::glue('Inferred pattern\nof inheritance\n(>{x}% support)', x = AIC_limit * 100),
      limits = c('auto', 'X', 'Y', 'not called'),
      labels = c('autosomal', 'X-linked', 'Y-linked', 'not called'),
      values = c('auto' = 1, 'X' = 'firebrick', 'Y' = 'blue', 'not called' = 'grey60'),
    ) +
    scale_x_continuous(
      limits = viz_range,
      name = 'position (kb)',
      label = \(x) scales::label_comma()(x / 1e3)
    ) +
    ylab(expression(-log[10](italic(p))))
}

make_kmer_manhattan_panel <- function(trait, chr, viz_range) {
  read_kmer_blast(trait) %>%
    filter(
      .data$chr == .env$chr,
      end > viz_range[1], start < viz_range[2],
    ) %>%
    slice_min(e_value, by = kmer) %>%
    slice_max(align_length - nr_mismatch - nr_gaps, by = kmer) %>%
    ggplot(aes((start + end) / 2)) +
    geom_histogram(binwidth = 1e3, fill = 'black') +
    scale_x_continuous(
      limits = viz_range,
      name = 'Position (kb)',
      label = \(x) scales::label_comma()(x / 1e3)
    ) +
    scale_y_continuous(expand = c(0, 0), name = 'Significant 31-mers')
}
