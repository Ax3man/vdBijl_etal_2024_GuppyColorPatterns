make_excess_het_panel <- function(geno, viz_range) {

  calc_excess_het <- function(n00, n01, n11) {
    n <- n00 + n01 + n11
    n0 <- 2 * n00 + n01
    n1 <- 2 * n11 + n01
    p <- n0 / (n0 + n1)
    q <- 1 - p
    # only test for *excess* heterozygosity, so test one-sided
    binom.test(n01, n, 2 * p * q, alternative = 'greater')[['p.value']]
  }
  excess_het <- as.data.frame(geno) %>%
    dplyr::count(start, GT) %>%
    complete(start, GT, fill = list(n = 0)) %>%
    group_by(start) %>%
    summarise(p_excess_het = calc_excess_het(n[GT == '0/0'], n[GT == '0/1'], n[GT == '1/1']))

  ggplot(excess_het, aes(start, -log10(p_excess_het))) +
    geom_point() +
    scale_x_continuous(
      limits = viz_range,
      name = 'position (kb)',
      label = \(x) scales::label_comma()(x / 1e3)
    ) +
    ylab(expression(-log[10](italic(p[excess~het]))))
}
