library(tidyverse)

source('paper_figures/subplots/car_ornament_img.R')
source('paper_figures/subplots/car_ornament_embedding.R')

source('paper_figures_supplement/subplots/mel_ornament_img.R')
source('paper_figures_supplement/subplots/mel_ornament_embedding.R')

car <- (P_car_orn_img) /
  (P_car_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_layout(heights = c(0.6, 2))

mel <- (P_mel_orn_img) /
  (P_mel_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_layout(heights = c(0.6, 2))

combined <- wrap_plots(car, mel, nrow = 2) + plot_annotation(tag_levels = 'A')
ggsave(
  'paper_figures_supplement/ornament_embeddings.png',
  combined, units = "cm", width = 18.4, height = 8
)
