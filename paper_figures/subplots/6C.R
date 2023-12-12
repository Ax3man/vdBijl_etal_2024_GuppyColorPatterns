library(patchwork)
library(ggtext)
source('sequencing/gwas/peak_viz/peak_viz_continuous.R')
source('paper_figures/theme.R')

if (!file.exists('paper_figures/subplots/6C_unadjusted.rds')) {
  panel_C <- peak_viz_continuous(
    trait = 'car_PIE',
    ornament = 'car_6',
    chr = 'NC_024333.1',
    pos = 20103960,

    include_kmers = TRUE,
    file = NULL,
    pval_column = 'p_SHet',
    fdr.level = 0.05,
    viz_range = c(20103960 - 1.5e4, 20103960 + 2.5e4)
  )
  write_rds(panel_C, 'paper_figures/subplots/6C_unadjusted.rds')
} else {
  panel_C <- read_rds('paper_figures/subplots/6C_unadjusted.rds')
}

eb <- element_blank()
new_x <- scale_x_continuous(
  name = 'Position (kb)',
  label = \(x) scales::label_comma()(x / 1e3),
  expand = c(0, 0)
)
coord <- coord_cartesian(xlim = c(20103960 - 1e4, 20103960 + 2e4), clip = 'on')

panel_C[[1]]$layers[[1]]$geom_params$arrow_body_height <- unit(2.5, 'mm')
panel_C[[1]]$layers[[1]]$geom_params$arrowhead_width <- unit(2, 'mm')
panel_C[[1]]$layers[[1]]$geom_params$arrowhead_height <- unit(3.5, 'mm')
panel_C[[1]]$layers[[2]]$geom_params$arrow_body_height <- unit(2.5, 'mm')
panel_C[[1]]$layers[[2]]$geom_params$arrowhead_width <- unit(2, 'mm')
panel_C[[1]]$layers[[2]]$geom_params$arrowhead_height <- unit(3.5, 'mm')
panel_C[[1]]$layers[[3]]$aes_params$size <- 2

panel_C[[2]]$layers[[1]] <- NULL
panel_C[[2]]$scales$scales[[1]]$name <- c('Inferred pattern\nof inheritance\n(>80% support)')

panel_C[[3]]$scales$scales[[2]]$name <- 'Significant 31-mers'

panel_C[[4]]$layers[[1]]$aes_params$linewidth <- 0.2
panel_C[[4]]$layers[[1]]$aes_params$alpha <- 0.05
panel_C[[4]]$layers[[2]]$aes_params$linewidth <- 0.4

panel_C[[6]][[1]] <- panel_C[[6]][[1]] +
  facet_wrap(vars("<img src='ornament_analysis/ornament_figure_labels/car_6.png' width='40' />")) +
  theme(strip.text = element_markdown(), strip.background = eb)
tmp <- panel_C[[6]][[1]]
panel_C[[6]][[1]] <- panel_C[[6]][[2]] + labs(x = 'Genotype at\ntop variant')
panel_C[[6]][[1]]$scales$scales[[1]]$name <- 'Est. allele count'
panel_C[[6]][[1]]$scales$scales[[2]]$name <- 'Allele'

panel_C[[6]][[2]] <- tmp + ylab('Count') + labs(x = 'Genotype at\ntop variant')
panel_C[[6]][[2]]$scales$scales[[1]]$name <- 'Ornament'

Fig6C <- wrap_plots(
  panel_C[[1]] + new_x + coord,
  panel_C[[4]] + new_x + coord + labs(y = 'Relative coverage') +
    theme(axis.text.x = eb, axis.line.x = eb, axis.ticks.x = eb, axis.title.x = eb),
  #panel_C[[3]] + new_x + coord,
  panel_C[[2]] + new_x + coord +
    geom_point(size = 0.6) +
    annotate(
      x = 20103960, y = -log10(panel_C[[2]]$data$p_SHet[panel_C[[2]]$data$ps == 20103960]),
      geom = 'point', shape = 1, color = 'red', size = 2, stroke = 0.6
    ),
  panel_C[[6]] + plot_layout(guides = 'collect'),
  snp_effects[[1]] + theme(text = element_text(size = 7)) +
    guides(fill = guide_colorbar(barwidth = unit(3, 'mm'), barheight = unit(15, 'mm'))),
  ncol = 1, heights = c(2, 3, 3, 3, 3)
)
