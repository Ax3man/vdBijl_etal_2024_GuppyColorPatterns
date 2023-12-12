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
  # drop P generation, since they weren't selected on yet
  filter(generation != 'P') %>%
  group_by(dam) %>%
  slice_min(days_since_first_brood) %>%
  ungroup()

ggplot(first_broods, aes(generation, days_since_first_brood, group = interaction(selection, generation))) +
  geom_violin(aes(fill = selection), color = 1, alpha = 0.6, position = 'dodge') +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c('navy', 'firebrick'), l = c('Down-selected', 'Up-selected')) +
  labs(
    y = 'Time to first brood\n(from first brood of the generation)',
    x = 'Selection',
    fill = NULL
  ) +
  theme_classic() +
  theme(legend.position = 'top')
ggsave('visualization/basic_plots/lh_time_to_first_brood.png', w = 5, h = 5)

m_time_to_first_brood <- brm(
  days_since_first_brood ~ 1 + selection * generation + (1 | replicate:generation),
  data = first_broods,
  cores = 4, control = list(adapt_delta = 0.95), iter = 4000
)
summary(m_time_to_first_brood)
emmeans(m_time_to_first_brood, pairwise ~ selection | generation)

ggplot(first_broods, aes(generation, nr_offspring, group = interaction(selection, generation))) +
  geom_violin(aes(fill = selection), color = 1, alpha = 0.6, position = 'dodge') +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c('navy', 'firebrick'), l = c('Down-selected', 'Up-selected')) +
  labs(
    y = 'Number of offspring of first brood',
    x = 'Generation',
    fill = NULL
  ) +
  theme_classic() +
  theme(legend.position = 'top')
ggsave('visualization/basic_plots/lh_fecundity.png', w = 5, h = 5)

m_fecundity <- brm(
  nr_offspring ~ selection * generation + (1 | replicate:generation),
  data = first_broods, family = negbinomial(),
  cores = 8, chains = 8, iter = 4000, control = list(adapt_delta = 0.999)
)
pp_check(m_fecundity) + xlim(0, 30)
summary(m_fecundity)
emmeans(m_fecundity, pairwise ~ selection | generation, type = 'response')

