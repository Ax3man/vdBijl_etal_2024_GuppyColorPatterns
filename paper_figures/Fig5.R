library(patchwork)
library(ggtext)

source('paper_figures/theme.R')
source('sequencing/gwas/make_manhattan2.R')

# Define new color scale to fix the legend title
col_scale <- scale_color_manual(
  values = c('auto' = 'black', 'X' = 'firebrick', 'Y' = 'blue3', 'Not called' = 'grey40'),
  labels = c('Autosomal', 'X-linked', 'Y-linked', 'Not called'),
  name = 'Pattern of inheritance (>80% support)',
  limits = c('auto', 'X', 'Y', 'Not called'),
  guide = guide_legend(override.aes = list(size = 1.5))
)

orange_pattern <- make_manhattan_new('car_PIE', 'p_SHet') +
  theme(legend.position = 'top') +
  labs(tag = 'A') +
  col_scale

black_pattern <- make_manhattan_new('mel_PIE', 'p_SHet') +
  col_scale +
  guides(color = 'none') +
  theme(legend.position = 'top') +
  labs(y = NULL)

# It is important to load all 7 ornaments, so that the q-values are computed over the full set
orange_ornaments <- manhattan_set(
  traits = paste0('pa_car_', 1:7), name = 'orange_ornaments', pval_column = 'p_lrt',
  joined_qvalue = TRUE, fdr_level = 0.05, return_plot_object = TRUE
) + facet_wrap(vars(trait_label), ncol = 1, scales = 'free_y')
# then filter down to the selected ornaments by subsetting the plot data
orange_ornaments <- orange_ornaments %+% (
  dplyr::filter(orange_ornaments$data, trait %in% paste0('pa_car_', c(1, 5, 6))) %>%
    mutate(trait_label = factor(trait_label, levels = sort(unique(trait_label))[c(2, 1, 3)]))
) + labs(tag = 'B') +
  theme(legend.position = 'none', strip.background.y = element_blank(), strip.text.y = element_blank())

black_ornaments <- manhattan_set(
  traits = paste0('pa_mel_', 1:8), name = 'black_ornaments', pval_column = 'p_lrt',
  joined_qvalue = TRUE, fdr_level = 0.05, return_plot_object = TRUE
) + facet_wrap(vars(trait_label), ncol = 1, scales = 'free_y')
black_ornaments <- black_ornaments %+%
  dplyr::filter(black_ornaments$data, trait %in% paste0('pa_mel_', c(1, 5, 8))) +
  theme(legend.position = 'none') +
  labs(y = NULL)

# Combine panels
panelA <- ((orange_pattern | black_pattern) & theme(legend.position = 'top')) +
  plot_layout(guides = 'collect')
panelB <- (
  (orange_ornaments | black_ornaments)
)

Fig5 <- (panelA / panelB) + plot_layout(heights = c(1, 2))
