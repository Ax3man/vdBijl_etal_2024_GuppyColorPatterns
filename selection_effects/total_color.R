source('selection_decisions/compile_decisions.R')
selection <- selection %>%
  dplyr::select(selection, fish_id, generation, selection, replicate) %>%
  mutate(fish_id = tolower(fish_id))

total_color <- data.table::fread('photo_database.csv') %>%
  group_by(fish_id, facing_direction) %>%
  summarise(across(c(car_perc, mel_perc_v2), mean), .groups = 'drop_last') %>%
  summarise(across(c(car_perc, mel_perc_v2), mean), .groups = 'drop') %>%
  left_join(selection, join_by(fish_id)) %>%
  filter(generation != 'P', selection != 'not_selected')

m_car <- brm(
  car_perc ~ 1 + selection * generation + (1 | replicate:generation),
  data = total_color, family = 'normal',
  cores = 4, control = list(adapt_delta = 0.8)
)
summary(m_car)
emmeans::emmeans(m_car, pairwise ~ selection, at = list(generation = 'F3'))

ggplot(total_color, aes(generation, car_perc, fill = selection)) +
  geom_violin() +
  geom_pointrange(
    aes(y = emmean, ymin = lower.HPD, ymax = upper.HPD),
    as.data.frame(emmeans::emmeans(m_car, ~ selection | generation)),
    position = position_dodge(width = 0.9)
  ) +
  theme_bw()

m_mel <- brm(
  mel_perc_v2 ~ 1 + selection * generation + (1 | replicate:generation),
  data = total_color, family = 'normal',
  cores = 4, control = list(adapt_delta = 0.9)
)
summary(m_mel)
emmeans::emmeans(m_mel, pairwise ~ selection, at = list(generation = 'F3'))

ggplot(total_color, aes(generation, mel_perc_v2, fill = selection)) +
  geom_violin() +
  geom_pointrange(
    aes(y = emmean, ymin = lower.HPD, ymax = upper.HPD),
    as.data.frame(emmeans::emmeans(m_mel, ~ selection | generation)),
    position = position_dodge(width = 0.9)
  ) +
  theme_bw()
