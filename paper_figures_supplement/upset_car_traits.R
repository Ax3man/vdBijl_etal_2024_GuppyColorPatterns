library(tidyverse)
library(qvalue)
library(furrr)
library(ComplexUpset)

source('paper_figures/theme.R')

# first, load GWAS results, and calculate q-values:
main <- data.table::fread('sequencing/gwas/gemma_output/car_PIE.assoc.txt')
q_main <- qvalue(main$p_SHet, fdr.level = 0.05)
main$significant <- q_main$significant

main <- filter(main, significant)

orn <- c('pa_car_1', 'pa_car_2', 'pa_car_3', 'pa_car_4', 'pa_car_5', 'pa_car_6', 'pa_car_7')

plan(multisession, workers = 12)
components <- setNames(orn, orn) %>%
  future_map_dfr(
    \(trait) data.table::fread(glue::glue('sequencing/gwas/gemma_output/{trait}.assoc.txt')),
    .id = 'trait'
  ) %>%
  mutate(trait = paste0('O', parse_number(trait)))
plan(sequential)

q <- qvalue(components$p_lrt, fdr.level = 0.05)
components$significant <- q$significant
components <- filter(components, significant)

table(components$trait)
nrow(components) / nrow(main)

# now get data on the source of these SNPs
snp_inheritance <- bind_rows(
  data.table::fread('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/orange_ornaments.csv'),
  data.table::fread('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/car_PIE.csv')
) %>%
  group_by(chr, start, end, alt) %>%
  slice_max(AIC_weight, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(source = if_else(AIC_weight < 0.8, 'Not called', source)) %>%
  distinct(chr, start, end, alt, source) %>%
  dplyr::select(-end)

# pivot and combine to get the 0/1 data for upset
to_plot <- bind_rows(
  main %>% mutate(trait = 'pattern'),
  components
) %>%
  select(trait, chr, ps, allele1) %>%
  mutate(tmp = 1) %>%
  pivot_wider(names_from = trait, values_from = tmp, values_fill = 0) %>%
  left_join(snp_inheritance, join_by(chr, ps == start, allele1 == alt)) %>%
  drop_na(source)

make_upset <- function(d, min_size) {
  u <-  upset(
    as.data.frame(d),
    intersect = rev(c('pattern', paste0('O', 1:7))),
    min_size = min_size, name = NULL, sort_sets = FALSE,

    base_annotations = list(
      Inheritance = (
        ggplot(mapping = aes(fill = factor(source, rev(c('Not called', 'Y', 'X', 'auto'))))) +
          geom_bar() +
          scale_fill_manual(
            values = c(auto = 'black', 'Not called' = 'grey', 'X' = 'firebrick', 'Y' = 'blue3'),
            limits = c('auto', 'X', 'Y', 'Not called'),
            labels = c('Autosomal', 'X-linked', 'Y-linked', 'Not called'),
          ) +
          scale_y_continuous(labels = scales::label_comma(scale = 1/1000, suffix = 'k')) +
          labs(y = 'Intersection size', fill = 'Inferred pattern\nof inheritance\n(>80% support)')
      )
    )
  ) & theme_get()

  eb <- element_blank()
  u[[2]] <- u[[2]] + theme(
    panel.grid.major.y = eb, axis.text.x = eb, axis.ticks.x = eb, axis.title.x = eb, axis.line.x = eb
  )
  u[[3]]$layers[[1]] <- NULL

  u[[3]] <- u[[3]] +
    aes(fill = factor(source, rev(c('Not called', 'Y', 'X', 'auto')))) +
    scale_fill_manual(
      values = c(auto = 'black', 'Not called' = 'grey', 'X' = 'firebrick', 'Y' = 'blue3'),
      limits = c('auto', 'X', 'Y', 'Not called'),
      labels = c('Autosomal', 'X-linked', 'Y-linked', 'Not called'),
      guide = 'none'
    ) +
    theme(
      axis.text.y = eb, axis.ticks.y = eb, axis.title.y = eb, axis.line.y = eb
    ) +
    scale_y_reverse(labels = scales::label_comma(scale = 1/1000, suffix = 'k'))
  u[[4]]$layers[[1]] <- NULL
  u[[4]] <- u[[4]] + theme(
    axis.text.x = eb, axis.ticks.x = eb, axis.title.x = eb, axis.line.x = eb
  ) + ylab(NULL)

  return(u)
}

upset_all <- make_upset(to_plot, min_size = 500)
upset_auto <- make_upset(
  filter(to_plot, chr != 'NC_024342.1', str_starts(chr, 'NW', TRUE)),
  min_size = 50
)

ggsave('paper_figures_supplement/upset_car_traits.png', upset_all, width = 6, height = 4)
ggsave('paper_figures_supplement/upset_car_traits_autosomal.png', upset_auto, width = 6, height = 4)


