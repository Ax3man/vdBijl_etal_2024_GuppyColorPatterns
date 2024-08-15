library(patchwork)

# Paper: The width of figures, when printed, will usually be 5.7 cm (2.24 inches or 1 column), 12.1
# cm (4.76 inches or 2 columns), or 18.4 cm (7.24 inches or 3 columns). Bar graphs, simple line
# graphs, and gels may be reduced to a smaller width. Symbols and lettering should be large enough
# to be legible after reduction [a reduced size of about 7 points (2.5 mm) high, and not smaller
# than 5 points]. Avoid wide variation in type size within a single figure. In laying out
# information in a figure, the objective is to maximize the space given to presentation of the data.
# Avoid wasted white space and clutter.

source('paper_figures/subplots/1_phenotyping_pipeline.R')
source('paper_figures/subplots/1_random_patrilines.R')
source('paper_figures/subplots/1_population_incidence_heatmaps.R')
fig1A <- wrap_elements(full = phenotyping_pipeline + labs(tag = 'A'))
fig1B <- ((example_patrilines + labs(tag = 'B')) | pop_incidence_heatmaps) + plot_layout(widths = c(4, 3))
fig1 <- (fig1A / fig1B) + plot_layout(heights = c(1, 2.2))
ggsave('paper_figures/Fig1.png', fig1, units = 'cm', width = 18.4, height = 9)
ggsave('paper_figures/Fig1.pdf', fig1, units = 'cm', width = 18.4, height = 9, device = cairo_pdf)

source('paper_figures/subplots/car_pedigree.R')
source('paper_figures/subplots/car_selection_maps.R')
source('paper_figures/subplots/embedding_selection.R')

Fig2 <- wrap_elements(full = P_car_pedigree) +
  (wrap_elements(full = P_car_selmap) /
     wrap_elements(embedding_selection) +
                  plot_layout(heights = c(1.6, 1))) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A')
ggsave('paper_figures/Fig2.png', Fig2, units = "cm", width = 18.4, height = 10)
ggsave('paper_figures/Fig2.pdf', Fig2, units = "cm", width = 18.4, height = 10, device = cairo_pdf)

source('paper_figures/subplots/car_pixel_h2.R')
source('paper_figures/subplots/mel_pixel_h2.R')

Fig3 <- (
  (P_car_pixel_h2) +
  (P_mel_pixel_h2 + theme(strip.text.y.left = element_blank())) &
  theme(legend.position = 'bottom')
) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect')
ggsave('paper_figures/Fig3.png', Fig3, units = "cm", width = 12.1, height = 7)
ggsave('paper_figures/Fig3.pdf', Fig3, units = "cm", width = 12.1, height = 7, device = cairo_pdf)

source('paper_figures/subplots/car_ornament_img.R')
source('paper_figures/subplots/car_ornament_presence_size.R')
source('paper_figures/subplots/car_ornament_h2.R')
#source('paper_figures/subplots/car_ornament_embedding.R')

Fig4 <- (P_car_orn_img) /
  (P_car_orn_h2 + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  (P_car_orn_pres + theme(legend.position = 'right', legend.justification = c(0, 0))) /
  P_car_orn_size /
  #(P_car_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(0.4, 1, 1, 1))
ggsave('paper_figures/Fig4.png', Fig4, units = "cm", width = 18.4, height = 7.5)
ggsave('paper_figures/Fig4.pdf', Fig4, units = "cm", width = 18.4, height = 7.5, device = cairo_pdf)

source('paper_figures/Fig5.R')
ggsave('paper_figures/Fig5.png', Fig5, width = 18.4, height = 12, units = 'cm', dpi = 600)
#ggsave('paper_figures/Fig5.pdf', Fig5, width = 18.4, height = 12, units = 'cm', dpi = 600, device = cairo_pdf)


source('paper_figures/subplots/6A.R')
source('paper_figures/subplots/6B.R')
source('paper_figures/subplots/6C.R')
Fig6 <- (Fig6A + plot_annotation(tag_levels = c('A')) & theme(legend.position = 'none')) |
  (Fig6B & theme(legend.position = 'none')) |
  (Fig6C)

ggsave('paper_figures/Fig6_unlabeled.png', Fig6, units = 'cm', width = 18, height = 15)
ggsave('paper_figures/Fig6_unlabeled.pdf', Fig6, units = 'cm', width = 18, height = 15, device = cairo_pdf)

source('paper_figures/subplots/7_embedding_space.R')
ggsave('paper_figures/Fig7.png', Fig7, units = 'cm', width = 12.1, height = 6.5)
