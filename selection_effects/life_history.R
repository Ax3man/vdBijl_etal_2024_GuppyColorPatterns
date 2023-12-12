library(tidyverse)
library(brms)
library(emmeans)

source('selection_decisions/compile_decisions.R')
breeding <- data.table::fread('data/breeding_data.csv') %>%
  mutate(date = lubridate::as_date(date, format = '%Y/%m/%d'))

dat <- breeding %>%
  inner_join(dplyr::select(selection, fish_id, replicate, generation, selection), c('sire' = 'fish_id')) %>%
  group_by(replicate, generation) %>%
  mutate(
    days_since_first_brood = as.numeric(date - min(date, na.rm = TRUE)),
    generation = factor(generation, c('P', 'F1', 'F2'))
  )

first_broods <- dat %>%
  # drop P generation, since they weren't selected on yet, note that F3 never reproduced.
  filter(generation != 'P') %>%
  group_by(dam) %>%
  slice_min(days_since_first_brood) %>%
  ungroup()

m_time_to_first_brood <- brm(
  days_since_first_brood ~ 1 + selection * generation + (1 | replicate:generation),
  data = first_broods,
  cores = 4, control = list(adapt_delta = 0.999), iter = 4000
)
pp_check(m_time_to_first_brood)
summary(m_time_to_first_brood)
emmeans(m_time_to_first_brood, pairwise ~ selection | generation)

P_brood_timing <- ggplot(
  first_broods,
  aes(generation, days_since_first_brood, group = interaction(selection, generation))
) +
  geom_violin(aes(fill = selection), color = 1, alpha = 0.6, position = 'dodge') +
  geom_pointrange(
    aes(y = emmean, ymin = lower.HPD, ymax = upper.HPD),
    as.data.frame(emmeans(m_time_to_first_brood, ~ selection | generation)),
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = c('navy', 'firebrick'), l = c('Down-selected', 'Up-selected')) +
  scale_x_discrete(labels = c(expression(F[1]), expression(F[2]))) +
  labs(
    y = 'Time to first brood\n(days from first brood of the generation)',
    x = 'Selection',
    fill = NULL
  ) +
  theme_classic() +
  theme(legend.position = 'top')


m_fecundity <- brm(
  nr_offspring ~ selection * generation + (1 | replicate:generation),
  data = first_broods, family = negbinomial(),
  cores = 8, chains = 8, iter = 4000, control = list(adapt_delta = 0.999)
)
pp_check(m_fecundity) + xlim(0, 30)
summary(m_fecundity)
emmeans(m_fecundity, pairwise ~ selection | generation, type = 'response')

P_fecundity <- ggplot(first_broods, aes(generation, nr_offspring, group = interaction(selection, generation))) +
  geom_violin(aes(fill = selection), color = 1, alpha = 0.6, position = 'dodge') +
  geom_pointrange(
    aes(y = prob, ymin = lower.HPD, ymax = upper.HPD),
    as.data.frame(emmeans(m_fecundity, ~ selection | generation, type = 'response')),
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = c('navy', 'firebrick'), l = c('Down-selected', 'Up-selected')) +
  scale_x_discrete(labels = c(expression(F[1]), expression(F[2]))) +
  labs(
    y = 'Number of offspring of first brood',
    x = 'Generation',
    fill = NULL
  ) +
  theme_classic() +
  theme(legend.position = 'top')

((P_fecundity | P_brood_timing) &
  theme(legend.position = 'bottom')) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')

ggsave('paper_figures_supplement/life_history.png', w = 6, h = 4)
