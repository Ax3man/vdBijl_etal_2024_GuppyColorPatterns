library(tidyverse)
source('paper_figures/theme.R')
source('selection_decisions/compile_decisions.R')

ped_df <- data.table::fread('data/pedigree.csv') %>%
  as_tibble() %>%
  mutate(
    sex = str_sub(animal, 1, 1),
    grow_tank = ifelse(is.na(sire), 'source_pop', paste(sire, dam, date_of_birth, sep = '_'))
  )

pd <- selection %>%
  dplyr::select(-sire, -dam) %>%
  left_join(ped_df, by = c('fish_id' = 'animal')) %>%
  left_join(
    dplyr::select(selection, fish_id, sire_car_perc = car_perc, sire_generation = generation),
    by = c('sire' = 'fish_id')
  ) %>%
  mutate(replicate = factor(replicate, paste0('replicate_', 1:3), paste('Replicate', 1:3)))

# Expect 300 missing values (P generation has no sires)
P_car_pedigree <- ggplot(pd, aes(x = car_perc, y = generation, color = selection)) +
    geom_hline(yintercept = 1:4, col = 'black', linewidth = 0.2) +
    ggridges::geom_density_ridges(
      data = filter(pd, generation == 'P'),
      fill = 1, alpha = 0.5, color = 'black', panel_scaling = FALSE, scale = 3, linewidth = 0.25,
      quantile_lines = TRUE, quantiles = 0.5
    ) +
    ggridges::geom_density_ridges(
      aes(fill = selection, group = interaction(selection, generation)), filter(pd, generation != 'P'),
      alpha = 0.5, color = 'black', panel_scaling = FALSE, scale = .6, linewidth = 0.25,
      quantile_lines = TRUE, quantiles = 0.5
    ) +
    geom_segment(aes(xend = sire_car_perc, yend = sire_generation), alpha = 0.1, linewidth = 0.15) +
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
    coord_cartesian(xlim = c(0, 30)) +
    labs(x = 'Orange coloration (% body area)', y = 'Generation', fill = NULL) +
    theme(legend.position = 'top')
#P_car_pedigree
