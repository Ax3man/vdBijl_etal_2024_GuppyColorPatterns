# Selection decision script for Replicate 1, first selection (on parents from the stock tanks).
library(tidyverse)
library(imager)

photo_database <- data.table::fread('photo_database.csv') %>% as_tibble()

paths <- photo_database %>%
  filter(replicate == 'replicate_1', generation == 'parental_gen_1') %>%
  mutate(
    extr = file.path('data/extracted_fish', replicate, generation, paste0(unique_id, '.png')),
    car = file.path('data/carotenoid_coloration_manual', paste0(unique_id, '.png'))
  )
stopifnot(file.exists(paths$extr) & file.exists(paths$car))

count_non_transparent_pixels <- function(path) {
  if (!file.exists(path)) return(NA)
  load.image(path) %>% channel(4) %>% {sum(. > 0.5)}
}
pics <- paths %>%
  mutate(
    total_pixels = map_int(extr, count_non_transparent_pixels),
    car_pixels = map_int(car, count_non_transparent_pixels),
    car_perc = car_pixels / total_pixels * 100
  )
ggplot(pics, aes(car_perc, fct_reorder(fish_id, car_perc), color = facing_direction)) +
  geom_point(shape = '|', size = 3) + expand_limits(x = 0) +
  scale_color_manual(values = c('firebrick',  'navy')) +
  theme_minimal() +
  labs(y = NULL, x = '% carotenoid coloration', color = 'side')

# Repeatabilities
library(lme4)
m <- lmer(car_perc ~ 1 + (1 | fish_id), pics)
performance::icc(m, by_group = TRUE)
#  Group   |   ICC
# ---------------
#  fish_id | 0.952
m2 <- lmer(car_perc ~ 1 + (1 | fish_id/facing_direction), pics)
performance::icc(m2, by_group = TRUE)

pics_summ <- pics %>%
  group_by(fish_id, facing_direction) %>%
  summarise(car_perc = mean(car_perc)) %>%
  summarise(car_perc = mean(car_perc)) %>%
  mutate(
    rank = rank(desc(car_perc)),
    selection = case_when(
      rank <= 30 ~ 'up_selected',
      rank > 70 ~ 'down_selected',
      TRUE ~ 'not_selected'
    )
  )
View(pics_summ)
data.table::fwrite(pics_summ, 'selection_decisions/rep1_parentalgen.csv')


# some other plots
# photo_database %>%
#   count(fish_id) %>%
#   ggplot(aes(n)) +
#   geom_histogram(binwidth = 1) +
#   expand_limits(x = 0) +
#   scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) + scale_y_continuous(expand = c(0, 0)) +
#   theme_minimal() +
#   labs(x = 'number of pictures per fish')
# ggplot(photo_database, aes(seq_along(unique_id), as.numeric(manual_landmarks))) +
#   geom_point(shape = '|') +
#   geom_smooth(method = 'glm', method.args = list(family = 'binomial'), color = 'firebrick') +
#   scale_y_continuous(expand = c(0, 0.01), labels = scales::percent) +
#   theme_minimal() +
#   labs(x = 'photo', y = 'percent manual extraction')
