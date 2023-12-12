make_QC_panel <- function(geno, viz_range, type = 'site') {
  if (type == 'site') {
    site_qual <- as.data.frame(geno) %>%
      dplyr::select(start, end, alt, QUAL, MQ) %>%
      distinct()

    ggplot() +
      geom_hline(aes(color = 'lower limit (Phred 30)', yintercept = 30), show.legend = TRUE) +
      geom_line(aes(start, QUAL, color = 'QUAL'), site_qual) +
      geom_line(aes(start, MQ, color = 'MQ'), site_qual) +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      ) +
      labs(y = 'Phred score', color = NULL) +
      coord_cartesian(ylim = c(NA, 60))
  } else {
    as.data.frame(geno) %>%
      group_by(start, end, alt) %>%
      summarise(pct_poor_calls = mean(GQ < 10), .groups = 'drop') %>%
      ggplot(aes(start, pct_poor_calls)) +
      geom_point() +
      scale_y_continuous(labels = scales::percent, name = '% GQ < 10') +
      scale_x_continuous(
        limits = viz_range,
        name = 'position (kb)',
        label = \(x) scales::label_comma()(x / 1e3)
      )
    #ggplot(as.data.frame(geno), aes(start, GQ)) +
      #geom_point(size = 0.3, alpha = 0.3)
      #geom_line(aes(group = sampleNames))
  }

}
