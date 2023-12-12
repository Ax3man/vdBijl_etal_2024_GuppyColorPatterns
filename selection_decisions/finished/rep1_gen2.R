# Selection decision script for Replicate 1, second selection (on the first offspring).
library(tidyverse)

photo_database <- data.table::fread('photo_database.csv') %>%
  as_tibble() %>%
  filter(replicate == 'replicate_1', generation == 'gen_2')
photo_database %>% group_by(fish_id) %>% filter(n() != 4)
photo_database %>% group_by(fish_id, facing_direction) %>% filter(n() != 2)

pedigree <- data.table::fread('data/pedigree.csv') %>%
  mutate(across(everything(), tolower))
fathers <- data.table::fread('selection_decisions/finished/rep1_parentalgen.csv') %>%
  select(fish_id, selection)
previous_decisions <- data.table::fread('selection_decisions/rep1_gen2_previous_decisions.csv')

ggplot(photo_database, aes(car_perc, fct_reorder(fish_id, car_perc), color = facing_direction)) +
  geom_point(shape = '|', size = 3) + expand_limits(x = 0) +
  scale_color_manual(values = c('firebrick',  'navy')) +
  theme_minimal() +
  labs(y = NULL, x = '% carotenoid coloration', color = 'side')

# Repeatabilities
library(lmerTest)
m <- lmer(car_perc ~ 1 + (1 | fish_id), photo_database)
performance::icc(m, by_group = TRUE)

m2 <- lmer(car_perc ~ 1 + (1 | fish_id/facing_direction), photo_database)
performance::icc(m2, by_group = TRUE)

pics_summ <- photo_database %>%
  group_by(date, fish_id, facing_direction) %>%
  summarise(car_perc = mean(car_perc), .groups = 'drop_last') %>%
  summarise(car_perc = mean(car_perc), .groups = 'drop') %>%
  left_join(pedigree, c('fish_id' = 'animal')) %>%
  left_join(fathers, c('sire' = 'fish_id'))
stopifnot(all(!is.na(pics_summ$selection)))
ggplot(pics_summ, aes(car_perc)) + geom_density() + geom_rug(sides = 't') +
  geom_boxplot(aes(y = -0.1), width = 0.1) + facet_grid(selection ~ .) + expand_limits(x = 0)

up <- filter(pics_summ, selection == 'up_selected') %>%
  # select top three per sire
  group_by(sire) %>%
  filter(rank(-car_perc) %in% 1:3) %>%
  ungroup() %>%
  # select top 30
  filter(rank(-car_perc) %in% 1:30)
down <- filter(pics_summ, selection == 'down_selected') %>%
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

data.table::fwrite(up, 'selection_decisions/rep1_gen2_up.csv')
data.table::fwrite(down, 'selection_decisions/rep1_gen2_down.csv')
data.table::fwrite(decisions, 'selection_decisions/rep1_gen2_previous_decisions.csv')
latest_decisions %>%
  select(date, fish_id, car_perc, decisions) %>%
  mutate(car_perc = round(car_perc, digits = 2)) %>%
  arrange(date, fish_id) %>%
  data.table::fwrite('selection_decisions/rep1_gen2_latest_decisions.csv')

# super basic model on differences between lines
m3 <- lmer(car_perc ~ selection + (1 | sire) + (1 | dam), pics_summ)
summary(m3)
# super basic model on heritability
m4 <- lmer(car_perc ~ 1 + (1 | sire) + (1 | dam), pics_summ)
sum(as.data.frame(VarCorr(m4))[1:2, 'vcov']) / sum(as.data.frame(VarCorr(m4))[1:3, 'vcov'])

# some other plots
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
    geom_smooth(color = 'firebrick') +
    scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
    theme_minimal() +
    labs(x = 'photo', y = 'percent manual extraction', title = 'Fish extraction'),
  ggplot(photo_database, aes(
    seq_along(fish_id),
    as.numeric(file.exists(paste0('data/carotenoid_coloration_manual/', unique_id, '.png')))
  )) +
    geom_point(shape = '|') +
    geom_smooth(color = 'firebrick') +
    scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
    theme_minimal() +
    labs(x = 'photo', y = 'percent manual extraction', title = 'Carotenoid extraction'),
  ncol = 1
)

mean(photo_database$manual_landmarks)
mean(file.exists(paste0('data/carotenoid_coloration_manual/', photo_database$unique_id, '.png')))
