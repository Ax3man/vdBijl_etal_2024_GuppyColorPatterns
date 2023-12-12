library(tidyverse)
library(brms)
library(emmeans)
library(patchwork)

d <- data.table::fread('data/behavior_morphology.csv')

m_SL <- brm(
  standard_length ~ 1 + selection + (1 | selection:replicate),
  data = d, family = 'normal',
  cores = 4, control = list(adapt_delta = 0.9)
)
pp_check(m_SL)
summary(m_SL)
emmeans(m_SL, pairwise ~ selection)

m_tail_length <- brm(
  tail_length ~ 1 + selection + (1 | selection:replicate),
  data = d, family = 'normal',
  cores = 4, control = list(adapt_delta = 0.9)
)
pp_check(m_tail_length)
summary(m_tail_length)
emmeans(m_tail_length, pairwise ~ selection)

m_tail_area <- brm(
  tail_area ~ 1 + selection + (1 | selection:replicate),
  data = d, family = 'normal',
  cores = 4, control = list(adapt_delta = 0.9)
)
pp_check(m_tail_area)
summary(m_tail_area)
emmeans(m_tail_area, pairwise ~ selection)


make_plot <- function(y, lab, model) {
  ggplot(d, aes(selection, .data[[y]])) +
    geom_violin(col = 1, linewidth = 0.1, fill = 'grey80') +
    geom_pointrange(
      aes(y = emmean, ymin = lower.HPD, ymax = upper.HPD),
      emmeans(model, ~ selection) %>% as.data.frame()
    ) +
    scale_x_discrete(limits = c('down-selected', 'up-selected')) +
    theme_classic() +
    ylab(lab)
}

(
  make_plot('standard_length', 'Standard length (cm)', m_SL) |
    make_plot('tail_length', 'Tail length (cm)', m_tail_length) |
    make_plot('tail_area', expression(Tail~area~(cm^2)), m_tail_area)
) + plot_annotation(tag_levels = 'A')

ggsave('paper_figures_supplement/gross_morphology.png', w = 10, h = 4.5)
