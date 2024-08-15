library(tidyverse)
source('sequencing/gwas/peak_viz/peak_viz_tools.R')

viz_range <- c(26585800, 26592800)
viz_range <- c(26570000, 26593000)
cov <- get_region_coverage(chr = 'NC_024344.1', viz_range = viz_range) |>
  mutate(sample_name = case_when(
    sample_name == '289_merged' ~ 'sample_289',
    sample_name == '323_merged' ~ 'sample_323',
    TRUE ~ sample_name
  )) |>
  inner_join(get_phenotypes('pa_car_6'), join_by(sample_name)) |>
  mutate(win_mid = (window_start + window_end) / 2)

slice_sample(cov, prop = 1) |>
  ggplot(aes(win_mid, rel_coverage, group = sample_name, color = phenotype_value == 1)) +
  geom_line(alpha = 0.3) +
  #geom_vline(xintercept = c(26587331, 26587645)) +
  geom_vline(aes(xintercept = 26587645, lty = 'Top variant')) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_color_manual(
    values = c('TRUE' = 'darkorange', 'FALSE' = 'grey10'),
    labels = c('Absent', 'Present'), limits = c(FALSE, TRUE),
    guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1))
  ) +
  scale_linetype_manual(values = 2) +
  coord_cartesian(xlim = viz_range) +
  theme_classic() +
  theme(legend.position = 'top') +
  labs(
    x = 'Position on LG14 (bp)',
    y = 'Relative coverage',
    color = 'Ornament O6',
    linetype = NULL
  )
ggsave('paper_figures_supplement/LG14_depth_by_O6.png', w = 10, h = 6)


geno <- get_genotypes(
  vcf = 'sequencing/gwas/filtered.vcf.gz', chr = 'NC_024344.1', range = c(26587645, 26587645 + 1)
) |>
  dplyr::rename(sample_name = sampleNames) |>
  inner_join(get_mean_coverage(), join_by(sample_name)) |>
  mutate(
    rel_depth = totalDepth / mean_coverage,
    sample_name = case_when(
      sample_name == '289_merged' ~ 'sample_289',
      sample_name == '323_merged' ~ 'sample_323',
      TRUE ~ sample_name
    )
  ) |>
  inner_join(get_phenotypes('pa_car_6'), join_by(sample_name))

ggplot(geno, aes(rel_depth, phenotype_value == 1, fill = phenotype_value == 1)) +
  geom_point(shape = '|', size = 3, position = position_nudge(y = -0.05)) +
  ggridges::geom_density_ridges(rel_min_height = 0.01, scale = 0.9) +
  scale_y_discrete(
    expand = expansion(c(0.1, 0.8)),
    limits = c(FALSE, TRUE),
    labels = c('Absent', 'Present')
  ) +
  scale_fill_manual(values = c('TRUE' = 'darkorange', 'FALSE' = 'grey30'), guide = 'none') +
  coord_cartesian(xlim = c(0, 3.8)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line()) +
  labs(
    x = 'Relative depth',
    y = 'Ornament O6'
  )
ggsave('paper_figures_supplement/LG14_depth_by_O6_SNP.png', w = 6, h = 4)
