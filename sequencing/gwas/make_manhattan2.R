library(tidyverse)
library(ggtext)

source('sequencing/genomics_helpers.R')
source('paper_figures/theme.R')

make_manhattan_new <- function(
    trait = NULL, pval_column, gwas_table = NULL, AIC_limit = 0.8, reference = 'female'
) {
  if (is.null(gwas_table)) {
    gwas_table <- prep_gwas_table(trait, produce_plots = FALSE, pval_column, reference = reference)
  }

  if (!is.null(trait)) {
    if (reference == 'female') {
      linkage_file <- glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')
    } else {
      linkage_file <- glue::glue('sequencing/male_reference_gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')
    }
    if (any(gwas_table$significant) && file.exists(linkage_file)) {
      linkage_table <- data.table::fread(linkage_file) %>%
        group_by(chr, start, end, alt) %>%
        slice_max(AIC_weight) %>%
        ungroup() %>%
        mutate(
          source = ifelse(AIC_weight < AIC_limit, 'Not called', source),
          chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr)
        )
      if (reference == 'male') {
        # for the male reference, we update the chr12 coordinates, following Fraser et al. and relevel the chr
        source('sequencing/male_reference_gwas/update_chr12_liftover.R')
        linkage_table <- linkage_table %>% mutate(
          across(c(start, end), \(ps) ifelse(chr == '12', update_chr12(ps), ps))
        )
      }
      gwas_table <- left_join(gwas_table, linkage_table, join_by(chr, ps == start, allele1 == alt)) |>
        # relevel for the male reference
        mutate(chr = fct_relevel(factor(chr), as.character(1:23))) |>
        # drop the MT for the female reference
        filter(chr != 'NC_024238.1')
    } else {
      gwas_table$significant <- gwas_table$source <- FALSE
    }
  }

  if (reference == 'female') {
    scaff_labs <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      mutate(
        scaff_mid = cum_scaff_start + 0.5 * scaff_len,
        chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                         str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                         TRUE ~ 'Un')
      ) %>% filter(!(chr2 %in% c('Un', 'MT')))
  } else {
    scaff_labs <- data.table::fread('sequencing/male_reference_gwas/male_reference.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      dplyr::select(-V4, -V5) |>
      mutate(
        chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr),
        chr = fct_relevel(factor(chr), as.character(1:23))
      ) |>
      arrange(chr) |>
      mutate(
        cum_scaff_start = cumsum(lag(scaff_len, default = 0)),
        scaff_mid = cum_scaff_start + 0.5 * scaff_len,
        chr2 = chr
      ) |>
      filter(!str_starts(chr, '0'))
  }

  xscale <- scale_x_continuous(expand = expansion(c(0.01, 0.01)),
                                 breaks = scaff_labs$scaff_mid, labels = scaff_labs$chr2)
  yscale <- scale_y_continuous(expand = expansion(c(0.005, 0.05)))

  # add two categories for the non-signficant values (for two tone chr display)
  gwas_table <- mutate(gwas_table, category = case_when(
    significant ~ source,
    !significant & as.numeric(factor(chr2)) %% 2 == 1 ~ 'unsig1',
    TRUE ~ 'unsig2'
  )) |>
    arrange(-.data[[pval_column]])

  ggplot(gwas_table, aes(x = cum_ps, y = -log10(.data[[pval_column]]), color = category)) +
    geom_point(size = 0.6, stroke = 0) +
    xscale + yscale +
    scale_color_manual(
      values = c(
        'auto' = 'black', 'X' = 'firebrick', 'Y' = 'blue3', 'Not called' = 'grey40',
        'unsig1' = 'grey80', 'unsig2' = 'grey90'
      ),
      breaks = c('auto', 'X', 'Y', 'Not called'),
      labels = c('Autosomal', 'X-linked', 'Y-linked', 'Not called'),
      name = 'Pattern of\ninheritance\n(>80% support)',
      limits = c('auto', 'X', 'Y', 'Not called', 'unsig1', 'unsig2')
    ) +
    coord_cartesian(clip = 'off') +
    labs(y = expression(-log[10](italic(p))), x = NULL) +
    theme(
      axis.ticks.x = element_blank(),
      strip.background = element_blank(), strip.placement = 'outside',
      #strip.text.y = element_blank(),
      panel.spacing = unit(2, 'points'),
      panel.background = element_rect(fill = NA) # to help avoid clipping of points
    )
}

single_manhattan <- function(trait, pval_column, reference = 'female') {
  man <- make_manhattan_new(trait, pval_column, reference = reference)
  ggsave(
    glue::glue('sequencing/gwas/manhattan_plots2/{trait}.png'),
    man, units = 'cm', width = 18.4, height = 6
  )
}
manhattan_set <- function(
    traits, name, pval_column, joined_qvalue = FALSE, fdr_level = 0.05, return_plot_object = FALSE,
    AIC_limit = 0.8, reference = "female"
) {
  if (joined_qvalue) {
    gwass <- map_dfr(
      setNames(traits, traits),
      \(tr) prep_gwas_table(
        tr, FALSE, pval_column, calc_q = FALSE, reference = reference, adjust_chr12 = FALSE
      ),
      .id = 'trait'
    )
    # the joint-set q-value
    q <- qvalue(gwass[[pval_column]], fdr.level = fdr_level)
    gwass$qvalue <- q$qvalues
    gwass$significant <- q$significant

    nm <- str_remove(name, '_joinedq')
    if (reference == 'female') {
      linkage_table <- data.table::fread(
        glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{nm}.csv')
      )
    } else {
      linkage_table <- data.table::fread(
        glue::glue('sequencing/male_reference_gwas/snp_inheritance/gwas_snp_inheritance/{nm}.csv')
      )
    }
  } else {
    gwass <- map_dfr(
      setNames(traits, traits),
      \(tr) prep_gwas_table(tr, FALSE, pval_column, reference = reference, adjust_chr12 = FALSE),
      .id = 'trait'
    )

    linkage_table <- map_dfr(setNames(traits, traits), \(tr) {
      if (reference == 'female') {
        f <- glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{tr}.csv')
      } else {
        f <- glue::glue('sequencing/male_reference_gwas/snp_inheritance/gwas_snp_inheritance/{tr}.csv')
      }
      if (file.exists(f)) data.table::fread(f) else data.frame()
    }, .id = 'trait')
  }

  linkage_table <- linkage_table %>%
    group_by(chr, start, end, alt) %>%
    slice_max(AIC_weight) %>%
    ungroup() %>%
    mutate(
      source = ifelse(AIC_weight < AIC_limit, 'Not called', source),
      chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr),
    ) %>%
    dplyr::select(chr, start, end, alt, source) %>%
    distinct()
  gwass <- left_join(gwass, linkage_table, join_by(chr, ps == start, allele1 == alt)) |>
    mutate(chr = fct_relevel(factor(chr), as.character(1:23)))

  if (reference == 'male') {
    # for the male reference, we update the chr12 coordinates, following Fraser et al. and relevel the chr
    source('sequencing/male_reference_gwas/update_chr12_liftover.R')
    gwass <- gwass %>% mutate(ps = ifelse(chr == '12', update_chr12(ps), ps))
  }

  gwass$trait_label <- glue::glue(
    "<img src='ornament_analysis/ornament_figure_labels/{x}.png' width='40' />",
    x = str_remove(gwass$trait, 'pa_')
  )

  mans <- make_manhattan_new(gwas_table = gwass, pval_column = pval_column, reference = reference) +
    facet_wrap(vars(trait_label), dir = 'v', ncol = 2, scales = 'free_y') +
    theme(strip.text = element_markdown())

  if (return_plot_object) {
    return(mans)
  } else {
    ggsave(
      glue::glue('sequencing/gwas/manhattan_plots2/{name}.png'),
      mans, units = 'cm', width = 18.4, height = 12
    )
    return(invisible(NULL))
  }
}

# run interactively, so if sourced you just get the functions
if (FALSE) {
  walk(
    c('car_PIE', 'mel_PIE'),
    single_manhattan, pval_column = 'p_SHet'
  )

  #manhattan_set(paste0('pa_car_', 1:7), 'orange_ornaments', 'p_lrt', joined_qvalue = FALSE)
  #manhattan_set(paste0('pa_mel_', 1:8), 'black_ornaments', 'p_lrt', joined_qvalue = FALSE)

  manhattan_set(paste0('pa_car_', 1:7), 'orange_ornaments_joinedq', 'p_lrt', joined_qvalue = TRUE)
  manhattan_set(paste0('pa_mel_', 1:8), 'black_ornaments_joinedq', 'p_lrt', joined_qvalue = TRUE)
}
