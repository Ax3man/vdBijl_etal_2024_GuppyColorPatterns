library(tidyverse)
source('sequencing/gwas/peak_viz/peak_viz_tools.R')

d <- data.table::fread('sequencing/yuying_parents/mapping_rates.txt') %>%
  dplyr::rename(sample_name = V1, mapping_rate = V2) %>%
  mutate(sample_name = str_remove(sample_name, '\\.bam\\:')) %>%
  left_join(yuying_parents, join_by(sample_name))

ggplot(d, aes(mapping_rate, sex, group = interaction(sex, sample_name), fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), col = 1, width = 0.5) +
  geom_text(
    aes(label = round(mapping_rate, 4)), position = position_dodge(width = 0.8),
    hjust = 1, size = 3
  ) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  theme_classic() +
  labs(x = 'Mapping rate', y = NULL)

d %>%
  group_by(sex) %>%
  summarise(mean_mapping_rate = mean(mapping_rate))

t.test(mapping_rate ~ sex, data = d)

ggsave('paper_figures_supplement/yuying_mapping_rates.png', width = 7, height = 2.5)
