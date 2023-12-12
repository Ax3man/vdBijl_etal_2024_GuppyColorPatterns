# Selection decision script for Replicate 1, third selection (on the second offspring, F2).
library(tidyverse); theme_set(theme_bw())

photo_database <- data.table::fread('photo_database.csv') %>%
  as_tibble() %>%
  filter(replicate == 'replicate_1', generation == 'gen_3')
photo_database %>% group_by(fish_id) %>% filter(n() < 3)

pedigree <- data.table::fread('data/pedigree.csv') %>%
  mutate(
    across(everything(), tolower),
    across(c(dam, sire), ~na_if(.x, ""))
  )
fathers <- data.table::fread('selection_decisions/finished/rep1_gen2_previous_decisions.csv') %>%
  select(fish_id, selection)
previous_decisions <- data.table::fread('selection_decisions/rep1_gen3_previous_decisions.csv')

ggplot(photo_database, aes(car_perc, fct_reorder(fish_id, car_perc), color = facing_direction)) +
  geom_point(shape = '|', size = 2) + expand_limits(x = 0) +
  scale_color_manual(values = c('firebrick',  'navy')) +
  theme_minimal() +
  labs(y = NULL, x = '% carotenoid coloration', color = 'side')
ggplot(photo_database, aes(mel_perc, fct_reorder(fish_id, mel_perc), color = facing_direction)) +
  geom_point(shape = '|', size = 2) + expand_limits(x = 0) +
  scale_color_manual(values = c('firebrick',  'navy')) +
  theme_minimal() +
  labs(y = NULL, x = '% melanic coloration', color = 'side')


# Repeatabilities
library(lmerTest); library(performance)
icc(lmer(car_perc ~ 1 + (1 | fish_id), photo_database), by_group = FALSE)
icc(lmer(car_perc ~ 1 + (1 | fish_id/facing_direction), photo_database), by_group = FALSE)
icc(lmer(mel_perc ~ 1 + (1 | fish_id), photo_database), by_group = FALSE)
icc(lmer(mel_perc ~ 1 + (1 | fish_id/facing_direction), photo_database), by_group = FALSE)

pics_summ <- photo_database %>%
  group_by(date, fish_id, facing_direction) %>%
  summarise(car_perc = mean(car_perc), .groups = 'drop_last') %>%
  summarise(car_perc = mean(car_perc), .groups = 'drop') %>%
  left_join(pedigree, c('fish_id' = 'animal')) %>%
  left_join(fathers, c('sire' = 'fish_id')) %>%
  drop_na(sire, dam)
stopifnot(all(!is.na(pics_summ$selection)))
ggplot(pics_summ, aes(car_perc)) + geom_density() +
  geom_point(aes(y = -0.03), shape = '|', size = 7, alpha = 0.3) +
  stat_summary(aes(y = 0), orientation = 'y', fun.data = 'mean_cl_boot') +
  geom_boxplot(aes(y = -0.1), width = 0.1) + facet_grid(selection ~ .) + expand_limits(x = 0)

up <- filter(pics_summ, selection == 'up_selected') %>%
  # select top three per sire
  group_by(sire) %>%
  filter(rank(-car_perc) %in% 1:3) %>%
  ungroup() %>%
  # select top 30
  filter(rank(-car_perc) %in% 1:30)
down <- filter(pics_summ, selection == 'down_selected') %>%
  # We lost m3-165, so need to remove him here
  filter(fish_id != 'm3-165') %>%
  # select bottom three per sire
  group_by(sire) %>%
  filter(rank(car_perc) %in% 1:3) %>%
  ungroup() %>%
  # select bottom 30
  filter(rank(car_perc) %in% 1:30)
decisions <- pics_summ %>%
  arrange(-date, desc(fish_id)) %>%
  transmute(
    date, fish_id, sire, dam, selection, car_perc,
    decisions = case_when(
      fish_id %in% up$fish_id ~ 'move to up',
      fish_id %in% down$fish_id ~ 'move to down',
      TRUE ~ 'sacrifice'
    )
  )
latest_decisions <- anti_join(
  decisions, previous_decisions,
  c("date", "fish_id", "decisions")
)

data.table::fwrite(up, 'selection_decisions/rep1_gen3_up.csv')
data.table::fwrite(down, 'selection_decisions/rep1_gen3_down.csv')
data.table::fwrite(decisions, 'selection_decisions/rep1_gen3_previous_decisions.csv')
latest_decisions %>%
  select(date, fish_id, car_perc, decisions) %>%
  mutate(car_perc = round(car_perc, digits = 2)) %>%
  arrange(date, fish_id) %>%
  data.table::fwrite('selection_decisions/rep1_gen3_latest_decisions.csv')

# super basic model on differences between lines
round(coef(summary(lm(car_perc ~ selection, pics_summ))), 5)
# super basic model on heritability
m <- lmer(car_perc ~ 1 + (1 | sire) + (1 | dam), pics_summ)
#m <- lmer(car_perc ~ 1 + (1 | sire), pics_summ)
performance::icc(m, by_group = TRUE)
performance::icc(m, by_group = FALSE)

# some other plots
photo_database %>%
  group_by(date, fish_id, facing_direction) %>%
  summarise(car_perc = mean(car_perc), mel_perc = mean(mel_perc), .groups = 'drop_last') %>%
  summarise(car_perc = mean(car_perc), mel_perc = mean(mel_perc), .groups = 'drop') %>%
  {
    ggplot(., aes(car_perc, mel_perc)) + geom_point() + geom_smooth(method = 'lm') + coord_equal() +
      labs(
        x = 'Carotenoid coloration\n(% body area)',
        y = 'Melanic coloration\n(% body area)',
        subtitle = paste0('r = ', round(cor.test(.$car_perc, .$mel_perc)$estimate, 2),
                          ', p = ', round(cor.test(.$car_perc, .$mel_perc)$p.value, 2))
      )
  }

photo_database %>%
  count(fish_id) %>%
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1) +
  expand_limits(x = 0) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) + scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  labs(x = 'number of pictures per fish')

cowplot::plot_grid(
  ggplot(photo_database, aes(seq_along(fish_id), as.numeric(manual_landmarks))) +
    geom_point(shape = '|') +
    geom_smooth(
      color = 'firebrick', method = 'gam', formula = y ~ s(x, bs = "cs"),
      method.args = list(family = 'binomial')
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
    theme_minimal() +
    labs(x = 'photo', y = 'percent manual extraction', title = 'Fish extraction'),
  ggplot(photo_database, aes(
    seq_along(fish_id),
    as.numeric(file.exists(paste0('data/carotenoid_coloration_manual/', unique_id, '.png')))
  )) +
    geom_point(shape = '|') +
    geom_smooth(
      color = 'firebrick', method = 'gam', formula = y ~ s(x, bs = "cs"),
      method.args = list(family = 'binomial')
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
    theme_minimal() +
    labs(x = 'photo', y = 'percent manual extraction', title = 'Carotenoid extraction'),
  ggplot(photo_database, aes(
    seq_along(fish_id),
    as.numeric(file.exists(paste0('data/melanic_coloration_manual/', unique_id, '.png')))
  )) +
    geom_point(shape = '|') +
    geom_smooth(
      color = 'firebrick', method = 'gam', formula = y ~ s(x, bs = "cs"),
      method.args = list(family = 'binomial')
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
    theme_minimal() +
    labs(x = 'photo', y = 'percent manual extraction', title = 'Melanic extraction'),
  ggplot(photo_database, aes(
    seq_along(fish_id),
    as.numeric(file.exists(paste0('data/manual_fixed_landmarks/', unique_id, '.txt')))
  )) +
    geom_point(shape = '|') +
    geom_smooth(
      color = 'firebrick', method = 'gam', formula = y ~ s(x, bs = "cs"),
      method.args = list(family = 'binomial')
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
    theme_minimal() +
    labs(x = 'photo', y = 'percent manual extraction', title = 'Landmark detection'),
  ncol = 1
)

# proportion manual handling:
mean(photo_database$manual_landmarks)
mean(file.exists(paste0('data/carotenoid_coloration_manual/', photo_database$unique_id, '.png')))
mean(file.exists(paste0('data/melanic_coloration_manual/', photo_database$unique_id, '.png')))
mean(file.exists(paste0('data/manual_fixed_landmarks/', photo_database$unique_id, '.txt')))
# proportion terrible pics:
mean(
  photo_database$manual_landmarks &
    file.exists(paste0('data/carotenoid_coloration_manual/', photo_database$unique_id, '.png')) &
    file.exists(paste0('data/melanic_coloration_manual/', photo_database$unique_id, '.png')) &
    file.exists(paste0('data/manual_fixed_landmarks/', photo_database$unique_id, '.txt'))
)
# proportion perfect pics:
mean(
  !photo_database$manual_landmarks &
    !file.exists(paste0('data/carotenoid_coloration_manual/', photo_database$unique_id, '.png')) &
    !file.exists(paste0('data/melanic_coloration_manual/', photo_database$unique_id, '.png')) &
    !file.exists(paste0('data/manual_fixed_landmarks/', photo_database$unique_id, '.txt'))
)

