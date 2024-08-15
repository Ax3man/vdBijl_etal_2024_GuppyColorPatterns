library(tidyverse)
source('paper_figures/theme.R')

model_results <- read_rds('visualization/per_pixel_models/car_selection_model_results.rds')
contrasts <- map_dfr(model_results, 'contrast')
contrasts$generation_label <- with(contrasts, case_when(
  generation == 'F1' ~ 'F[1]', generation == 'F2' ~ 'F[2]', generation == 'F3' ~ 'F[3]',
))

complete_fish <- imager::load.image('data/extracted_fish_warped/replicate_3/gen_2/20210602_IMG_2097.png') %>%
  as.data.frame(wide = 'c') %>% filter(c.4 > 0.1) %>% dplyr::select(x, y, alpha = c.4)
eb <- element_blank()

# note that estimates are one the logg odds ratio (not response) scale
P_car_selmap <- ggplot(contrasts, aes(x, y, fill = exp(-estimate))) +
  geom_raster(
    aes(x = x, y = y, alpha = alpha),
    fill = 'grey60', data = complete_fish, show.legend = FALSE, inherit.aes = FALSE
  ) +
  geom_raster() +
  scale_alpha_identity(guide = 'none') +
  scale_y_reverse() +
  scico::scale_fill_scico(
    palette = 'vik',
    trans = 'log',
    limits = c(1/31, 31),
    oob = scales::squish,
    breaks = c(1/30, 1/10, 1/3, 1, 3/1, 10, 30),
    labels = c('1/30', '1/10', '1/3', '1', '3', '10', '30'),
    guide = guide_colorbar(title.position = 'top')
  ) +
  facet_grid(generation_label ~ ., switch = 'y', labeller = label_parsed) +
  coord_fixed(expand = FALSE) +
  labs(fill = 'Incidence\nodds ratio') +
  theme(
    legend.position = 'right',
    legend.key.height = grid::unit(30, "points"), legend.key.width = unit(0.4, 'lines'),
    legend.title.align = 0.5,
    strip.background = eb, strip.text.y.left = element_text(angle = 0),
    axis.text = eb, axis.title = eb, axis.ticks = eb, axis.line = eb,
  )

# old version with horizontal color bar
# P_car_selmap <- ggplot(contrasts, aes(x, y, fill = exp(-estimate))) +
#   geom_raster(
#     aes(x = x, y = y, alpha = alpha),
#     fill = 'grey60', data = complete_fish, show.legend = FALSE, inherit.aes = FALSE
#   ) +
#   geom_raster() +
#   scale_alpha_identity(guide = 'none') +
#   scale_y_reverse() +
#   scico::scale_fill_scico(
#     palette = 'vik',
#     trans = 'log',
#     limits = c(1/31, 31),
#     oob = scales::squish,
#     breaks = c(1/30, 1/10, 1/3, 1, 3/1, 10, 30),
#     labels = c(
#       expression(paste('< ', over(1, 30))),
#       expression(over(1, 10), over(1, 3), over(1, 1), over(3, 1), over(10, 1)),
#       expression(paste('> ', over(30, 1)))),
#     guide = guide_colorbar(title.position = 'top')
#   ) +
#   facet_grid(generation_label ~ ., switch = 'y', labeller = label_parsed) +
#   coord_fixed(expand = FALSE) +
#   labs(fill = 'down > up        Incidence odds ratio        up > down') +
#   theme(
#     legend.position = 'bottom',
#     legend.key.width = grid::unit(35, "points"), legend.key.height = unit(0.4, 'lines'),
#     legend.title.align = 0.5,
#     strip.background = eb, strip.text.y.left = element_text(angle = 0),
#     axis.text = eb, axis.title = eb, axis.ticks = eb, axis.line = eb,
#   )


# calculate rough average effect sizes for different body areas
ggplot2::`%+%`(P_car_selmap, filter(contrasts, x < 250))
filter(contrasts, x < 250, generation == 'F3') %>% pull(estimate) %>% mean() %>% {exp(-.)}

ggplot2::`%+%`(P_car_selmap, filter(contrasts, x > 310, y > 75))
filter(contrasts, x > 310, y > 75, generation == 'F3') %>% pull(estimate) %>% mean() %>% {exp(-.)}

ggplot2::`%+%`(P_car_selmap, filter(contrasts, x > 330, y < 65))
filter(contrasts, x > 330, y < 65, generation == 'F3') %>% pull(estimate) %>% mean() %>% {exp(.)}
