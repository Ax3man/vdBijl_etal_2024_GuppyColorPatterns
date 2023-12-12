library(tidyverse)
library(ggtext)

source('sequencing/genomics_helpers.R')
source('paper_figures/theme.R')

make_manhattan_new <- function(trait = NULL, pval_column, gwas_table = NULL, AIC_limit = 0.8) {
  if (is.null(gwas_table)) gwas_table <- prep_gwas_table(trait, produce_plots = FALSE, pval_column)

  if (!is.null(trait)) {
    if (any(gwas_table$significant)) {
      linkage_table <- data.table::fread(glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')) %>%
        group_by(chr, start, end, alt) %>%
        slice_max(AIC_weight) %>%
        ungroup() %>%
        mutate(source = ifelse(AIC_weight < AIC_limit, 'Not called', source))
      gwas_table <- left_join(gwas_table, linkage_table, join_by(chr, ps == start, allele1 == alt))
    } else {
      gwas_table$source <- NA
    }
  }

  scaff_labs <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
    dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
    mutate(
      scaff_mid = cum_scaff_start + 0.5 * scaff_len,
      chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                       str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                       TRUE ~ 'Un')
    ) %>% filter(!(chr2 %in% c('Un', 'MT')))

  xscale <- scale_x_continuous(expand = expansion(c(0.01, 0.01)),
                                 breaks = scaff_labs$scaff_mid, labels = scaff_labs$chr2)
  yscale <- scale_y_continuous(expand = expansion(c(0.005, 0.05)))

  ggplot(slice_sample(gwas_table, prop = 1), aes(x = cum_ps, y = -log10(.data[[pval_column]]))) +
    geom_point(
      aes(color = as.character(as.numeric(factor(chr)) %% 2)),
      data = \(d) filter(d, !significant, as.numeric(factor(chr)) %% 2 == 0),
      size = 0.6, stroke = 0, color = 'grey90'
    ) +
    geom_point(
      aes(color = as.character(as.numeric(factor(chr)) %% 2)),
      data = \(d) filter(d, !significant, as.numeric(factor(chr)) %% 2 == 1),
      size = 0.6, stroke = 0, color = 'grey80'
    ) +
    geom_point(
      aes(color = source),
      data = \(d) filter(d, significant), size = 0.6, stroke = 0
    ) +
    xscale + yscale +
    scale_color_manual(
      values = c('auto' = 'black', 'X' = 'firebrick', 'Y' = 'blue3', 'Not called' = 'grey40'),
      labels = c('Autosomal', 'X-linked', 'Y-linked', 'Not called'),
      name = 'Pattern of\ninheritance\n(>80% support)',
      limits = c('auto', 'X', 'Y', 'Not called')
    ) +
    coord_cartesian(clip = 'off') +
    labs(y = expression(-log[10](italic(p))), x = NULL) +
    theme(
      axis.ticks.x = element_blank(),
      strip.background = element_blank(), strip.placement = 'outside',
      panel.spacing = unit(2, 'points'),
      panel.background = element_rect(fill = NA) # to help avoid clipping of points
    )
}

single_manhattan <- function(trait, pval_column) {
  man <- make_manhattan_new(trait, pval_column)
  ggsave(
    glue::glue('sequencing/gwas/manhattan_plots2/{trait}.png'),
    man, units = 'cm', width = 18.4, height = 6
  )
}
manhattan_set <- function(
    traits, name, pval_column, joined_qvalue = FALSE, fdr_level = 0.05, return_plot_object = FALSE,
    AIC_limit = 0.8
) {
  if (joined_qvalue) {
    gwass <- map_dfr(
      setNames(traits, traits),
      \(tr) prep_gwas_table(tr, FALSE, pval_column, calc_q = FALSE),
      .id = 'trait'
    )
    # the joint-set q-value
    q <- qvalue(gwass[[pval_column]], fdr.level = fdr_level)
    gwass$qvalue <- q$qvalues
    gwass$significant <- q$significant

    nm <- str_remove(name, '_joinedq')
    linkage_table <- data.table::fread(
      glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{nm}.csv')
    )
  } else {
    gwass <- map_dfr(setNames(traits, traits), \(tr) prep_gwas_table(tr, FALSE, pval_column), .id = 'trait')

    linkage_table <- map_dfr(setNames(traits, traits), \(tr) {
      f <- glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{tr}.csv')
      if (file.exists(f)) data.table::fread(f) else data.frame()
    }, .id = 'trait')
  }

  linkage_table <- linkage_table %>%
    group_by(chr, start, end, alt) %>%
    slice_max(AIC_weight) %>%
    ungroup() %>%
    mutate(source = ifelse(AIC_weight < AIC_limit, 'Not called', source)) %>%
    dplyr::select(chr, start, end, alt, source) %>%
    distinct()
  gwass <- left_join(gwass, linkage_table, join_by(chr, ps == start, allele1 == alt))

  gwass$trait_label <- glue::glue(
    "<img src='ornament_analysis/ornament_figure_labels/{x}.png' width='40' />",
    x = str_remove(gwass$trait, 'pa_')
  )

  mans <- make_manhattan_new(gwas_table = gwass, pval_column = pval_column) +
    facet_wrap(vars(trait_label), dir = 'v', ncol = 2, scales = 'free_y') +
    theme(strip.text = element_markdown())

  if (return_plot_object) return(mans)

  ggsave(
    glue::glue('sequencing/gwas/manhattan_plots2/{name}.png'),
    mans, units = 'cm', width = 18.4, height = 12
  )
}

# run interactively, so if sourced you just get the functions
if (FALSE) {
  walk(
    c('car_PIE', 'mel_PIE'),
    single_manhattan, pval_column = 'p_SHet'
  )

  manhattan_set(paste0('pa_car_', 1:7), 'orange_ornaments', 'p_lrt', joined_qvalue = FALSE)
  manhattan_set(paste0('pa_mel_', 1:8), 'black_ornaments', 'p_lrt', joined_qvalue = FALSE)

  manhattan_set(paste0('pa_car_', 1:7), 'orange_ornaments_joinedq', 'p_lrt', joined_qvalue = TRUE)
  manhattan_set(paste0('pa_mel_', 1:8), 'black_ornaments_joinedq', 'p_lrt', joined_qvalue = TRUE)
}
