library(tidyverse)
source('paper_figures/theme.R')

model_results <- read_rds('visualization/per_pixel_models/mel_selection_model_results.rds')
contrasts <- map_dfr(model_results, 'contrast')
contrasts$generation_label <- with(contrasts, case_when(
  generation == 'F1' ~ 'F[1]', generation == 'F2' ~ 'F[2]', generation == 'F3' ~ 'F[3]',
))

complete_fish <- imager::load.image('data/extracted_fish_warped/replicate_3/gen_2/20210602_IMG_2097.png') %>%
  as.data.frame(wide = 'c') %>% filter(c.4 > 0.1) %>% dplyr::select(x, y, alpha = c.4)
eb <- element_blank()

# note that estimates are one the logg odds ratio (not response) scale
P_mel_selmap <- ggplot(contrasts, aes(x, y, fill = exp(-estimate))) +
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
    labels = c(
      expression(paste('< ', over(1, 30))),
      expression(over(1, 10), over(1, 3), over(1, 1), over(3, 1), over(10, 1)),
      expression(paste('> ', over(30, 1)))),
    guide = guide_colorbar(title.position = 'top')
  ) +
  facet_grid(generation_label ~ ., switch = 'y', labeller = label_parsed) +
  coord_fixed(expand = FALSE) +
  labs(fill = 'Incidence odds ratio\ndown > up                                                      up > down') +
  theme(
    legend.position = 'bottom',
    legend.key.width = grid::unit(35, "points"), legend.key.height = unit(0.8, 'lines'),
    legend.title.align = 0.5,
    strip.background = eb, strip.text.y.left = element_text(angle = 0),
    axis.text = eb, axis.title = eb, axis.ticks = eb, axis.line = eb,
  )

library(tidyverse)
source('paper_figures/theme.R')
source('selection_decisions/compile_decisions.R')

ped_df <- data.table::fread('data/pedigree.csv') %>%
  as_tibble() %>%
  mutate(
    sex = str_sub(animal, 1, 1),
    grow_tank = ifelse(is.na(sire), 'source_pop', paste(sire, dam, date_of_birth, sep = '_'))
  ) %>%
  mutate(across(animal:dam, tolower))

mel <- data.table::fread('photo_database.csv') %>%
  group_by(replicate, generation, fish_id, facing_direction) %>%
  summarise(mel_perc_v2 = mean(mel_perc_v2), .groups = 'drop_last') %>%
  summarise(mel_perc_v2 = mean(mel_perc_v2), .groups = 'drop') %>%
  mutate(
    replicate = factor(replicate, paste0('replicate_', 1:3), paste('Replicate', 1:3)),
    generation = factor(
      generation, c('parental_gen_1', 'gen_2', 'gen_3', 'gen_4'), c('P', 'F1', 'F2', 'F3'))
  )

pd <- mel %>%
  left_join(dplyr::select(selection, fish_id, selection) %>% mutate(fish_id = tolower(fish_id))) %>%
  left_join(ped_df, by = c('fish_id' = 'animal')) %>%
  left_join(
    dplyr::select(mel, fish_id, sire_mel_perc = mel_perc_v2, sire_generation = generation),
    by = c('sire' = 'fish_id')
  )

# Expect 300 missing values (P generation has no sires)
P_mel_pedigree <- ggplot(pd, aes(x = mel_perc_v2, y = generation, color = selection)) +
  geom_hline(yintercept = 1:4, col = 'black', linewidth = 0.2) +
  ggridges::geom_density_ridges(
    data = filter(pd, generation == 'P'),
    fill = 1, alpha = 0.5, color = 'black', panel_scaling = FALSE, scale = 1, size = 0.25,
    quantile_lines = TRUE, quantiles = 0.5
  ) +
  ggridges::geom_density_ridges(
    aes(fill = selection, group = interaction(selection, generation)), filter(pd, generation != 'P'),
    alpha = 0.5, color = 'black', panel_scaling = FALSE, scale = .6, size = 0.25,
    quantile_lines = TRUE, quantiles = 0.5
  ) +
  geom_segment(aes(xend = sire_mel_perc, yend = sire_generation), alpha = 0.1, linewidth = 0.15) +
  #geom_point(alpha = .6, size = 1, shape = '|') +
  #facet_wrap( ~ replicate) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(
    limits = rev(c('P', 'F1', 'F2', 'F3')),
    labels = rev(c('P', expression(F[1]), expression(F[2]), expression(F[3]))),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_color_manual(
    values = c('navy', 'black', 'firebrick'),
    labels = c('Down-selected', 'Up-selected'),
    guide = 'none'
  ) +
  scale_fill_manual(
    values = c('navy', 'firebrick'),
    labels = c('Down-selected', 'Up-selected'),
    guide = guide_legend(override.aes = list(alpha = 1), )
  ) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = 'Black coloration (% body area)', y = 'Generation', fill = NULL) +
  theme(legend.position = 'top')

P_mel_pedigree | P_mel_selmap

wrap_elements(full = P_mel_pedigree) +
  wrap_elements(full = P_mel_selmap) +
  plot_annotation(tag_levels = 'A')

ggsave('paper_figures_supplement/mel_selection_response.png', units = "cm", width = 18.4, height = 8)
