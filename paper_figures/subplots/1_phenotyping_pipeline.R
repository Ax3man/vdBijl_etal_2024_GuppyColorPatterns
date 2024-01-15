library(tidyverse)
source('paper_figures/theme.R')

im <- imager::load.image('paper_figures/subplots/guppy colour pipeline_WvdBedit.png') %>%
  as.data.frame(wide = 'c')
phenotyping_pipeline <- ggplot(im, aes(x, y, fill = rgb(c.1, c.2, c.3, c.4))) +
  geom_raster() +
  scale_fill_identity() +
  scale_y_reverse() +
  coord_equal() +
  theme_void(base_size = 7, base_family = 'Arial')
