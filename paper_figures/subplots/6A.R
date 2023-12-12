library(patchwork)
library(ggtext)
source('sequencing/gwas/peak_viz/peak_viz_continuous.R')
source('paper_figures/theme.R')

if (!file.exists('paper_figures/subplots/6A_unadjusted.rds')) {
  panel_A <- peak_viz_continuous(
    trait = 'car_PIE',
    ornament = 'car_6',
    chr = 'NC_024344.1',
    pos = 26587331,

    include_kmers = TRUE,
    file = NULL,
    pval_column = 'p_SHet',
    fdr.level = 0.05,
    viz_range = c(26550000, 26620000)
  )
  write_rds(panel_A, 'paper_figures/subplots/6A_unadjusted.rds')
} else {
  panel_A <- read_rds('paper_figures/subplots/6A_unadjusted.rds')
}

snp_effects <- read_rds('paper_figures/subplots/6_snp_effects.rds')

eb <- element_blank()
new_x <- scale_x_continuous(
  name = 'Position (kb)',
  label = \(x) scales::label_comma()(x / 1e3),
  expand = c(0, 0)
)
coord <- coord_cartesian(xlim = c(26553000, 26595000), clip = 'on')

panel_A[[1]]$layers[[1]]$geom_params$arrow_body_height <- unit(2.5, 'mm')
panel_A[[1]]$layers[[1]]$geom_params$arrowhead_width <- unit(2, 'mm')
panel_A[[1]]$layers[[1]]$geom_params$arrowhead_height <- unit(3.5, 'mm')
panel_A[[1]]$layers[[2]]$geom_params$arrow_body_height <- unit(2.5, 'mm')
panel_A[[1]]$layers[[2]]$geom_params$arrowhead_width <- unit(2, 'mm')
panel_A[[1]]$layers[[2]]$geom_params$arrowhead_height <- unit(3.5, 'mm')
panel_A[[1]]$layers[[3]]$aes_params$size <- 2

panel_A[[2]]$layers[[1]] <- NULL

panel_A[[3]]$scales$scales[[2]]$name <- 'Significant 31-mers'

panel_A[[4]]$layers[[1]]$aes_params$linewidth <- 0.2
panel_A[[4]]$layers[[1]]$aes_params$alpha <- 0.05
panel_A[[4]]$layers[[2]]$aes_params$linewidth <- 0.4

panel_A[[6]][[1]] <- panel_A[[6]][[1]] +
  facet_wrap(vars("<img src='ornament_analysis/ornament_figure_labels/car_6.png' width='40' />")) +
  theme(strip.text = element_markdown(), strip.background = eb)
tmp <- panel_A[[6]][[1]]
panel_A[[6]][[1]] <- panel_A[[6]][[2]] + labs(x = 'Genotype at\ntop variant')#, tag = 'D')
panel_A[[6]][[1]]$scales$scales[[1]]$name <- 'Est. allele count'
panel_A[[6]][[1]]$scales$scales[[2]]$name <- 'Allele'

panel_A[[6]][[2]] <- tmp + labs(x = 'Genotype at\ntop variant', y = 'Count')#, tag = 'E')
panel_A[[6]][[2]]$scales$scales[[1]]$name <- 'Ornament'


Fig6A <- wrap_plots(
  panel_A[[1]] + new_x + coord,# + labs(tag = 'A') + theme(text = element_text(size = 7)),
  panel_A[[4]] + new_x + coord + labs(y = 'Relative coverage') + #, tag = 'B') +
    theme(axis.text.x = eb, axis.line.x = eb, axis.ticks.x = eb, axis.title.x = eb),
  #panel_A[[3]] + new_x + coord,
  panel_A[[2]] + new_x + coord +
    geom_point(size = 0.6) +
    annotate(
      x = 26587331, y = -log10(panel_A[[2]]$data$p_SHet[panel_A[[2]]$data$ps == 26587331]),
      geom = 'point', shape = 1, color = 'red', size = 2, stroke = 0.6
    ),# +
    #labs(tag = 'C'),
  panel_A[[6]],
  snp_effects[[3]],
  ncol = 1, heights = c(2, 3, 3, 3, 3)
)


