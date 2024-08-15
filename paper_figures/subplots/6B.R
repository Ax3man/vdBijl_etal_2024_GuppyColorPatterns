library(patchwork)
library(ggtext)
source('sequencing/gwas/peak_viz/peak_viz_continuous.R')
source('paper_figures/theme.R')

if (!file.exists('paper_figures/subplots/6B_unadjusted.rds')) {
  panel_B <- peak_viz_continuous(
    trait = 'car_PIE',
    ornament = 'car_1',
    chr = 'NC_024333.1',
    pos = 593580,

    include_kmers = TRUE,
    file = NULL,
    pval_column = 'p_SHet',
    fdr.level = 0.05,
    viz_range = c(593580 - 2e4, 593580 + 2e4)
  )
  write_rds(panel_B, 'paper_figures/subplots/6B_unadjusted.rds')
} else {
  panel_B <- read_rds('paper_figures/subplots/6B_unadjusted.rds')
}

eb <- element_blank()
new_x <- scale_x_continuous(
  name = 'Position (kb)',
  label = \(x) scales::label_comma()(x / 1e3),
  expand = c(0, 0)
)
coord <- coord_cartesian(xlim = c(593580 - 1.4e4, 593580 + 2e4), clip = 'on')

panel_B[[1]]$layers[[1]]$geom_params$arrow_body_height <- unit(2.5, 'mm')
panel_B[[1]]$layers[[1]]$geom_params$arrowhead_width <- unit(2, 'mm')
panel_B[[1]]$layers[[1]]$geom_params$arrowhead_height <- unit(3.5, 'mm')
panel_B[[1]]$layers[[2]]$geom_params$arrow_body_height <- unit(2.5, 'mm')
panel_B[[1]]$layers[[2]]$geom_params$arrowhead_width <- unit(2, 'mm')
panel_B[[1]]$layers[[2]]$geom_params$arrowhead_height <- unit(3.5, 'mm')
panel_B[[1]]$layers[[3]]$aes_params$size <- 2

panel_B[[2]]$layers[[1]]$aes_params$linewidth <- 0.2
panel_B[[2]]$layers[[1]]$aes_params$alpha <- 0.05
panel_B[[2]]$layers[[2]]$aes_params$linewidth <- 0.4

panel_B[[5]][[1]]$layers[[2]]$aes_params$size <- 2

Fig6B <- wrap_plots(
  panel_B[[1]] + new_x + coord,
  panel_B[[2]] + new_x + coord + labs(y = 'Relative coverage') +
    theme(axis.text.x = eb, axis.line.x = eb, axis.ticks.x = eb, axis.title.x = eb),
  panel_B[[3]] + new_x + coord +
    geom_point(size = 0.6, show.legend = TRUE) +
    annotate(
      x = 593580, y = -log10(panel_B[[3]]$data$p_SHet[panel_B[[3]]$data$ps == 593580]),
      geom = 'point', shape = 1, color = 'red', size = 2, stroke = 0.6
    ) +
    theme(axis.text.x = eb, axis.line.x = eb, axis.ticks.x = eb, axis.title.x = eb),
  panel_B[[4]] + new_x + coord,
  panel_B[[5]] + plot_layout(guides = 'collect'),
  snp_effects[['NC_024333.1/593580/593580/T']],
  ncol = 1, heights = c(2, 3, 3, 3, 3 ,3)
)


