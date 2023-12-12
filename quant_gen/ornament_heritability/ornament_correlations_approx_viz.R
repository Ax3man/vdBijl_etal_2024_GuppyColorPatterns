library(tidyverse)

viz_cors <- function(output) {
  # This data.frame is organized in groups of 12, 3 parameters for 4 variance components. Note that
  # these are all variances and co-variances. First, we can convert all of those to correlations:
  cors <- output %>%
    mutate(model_component = rep(seq_len(n() / 3), each = 3)) %>%
    split(., .$model_component) %>%
    map_dfr(\(df) {
      x <- df$VarComp
      V <- matrix(x[c(1, 2, 2, 3)], nc = 2)

      invalid_VCV <- any(abs(V) < .Machine$double.eps) | any(x[c(1, 3)] < 0)
      sig_variances <- all(df$p_asymp[c(1, 3)] < 0.05)

      if (sig_variances && !invalid_VCV) {
        r <- cov2cor(V)[1, 2]
      } else {
        r <- NA
      }

      out <- df[2, c('param', 'trait1', 'trait2', 'p_asymp')]
      out$r <- r

      return(out)
    })


  cors %>%
    drop_na() %>%
    mutate(
      trait1 = trait1 %>% str_replace('mel_', 'B') %>% str_replace('car_', 'O'),
      trait2 = trait2 %>% str_replace('mel_', 'B') %>% str_replace('car_', 'O'),
      param = fct_recode(
        param,
        'autosomal' = 'u:fish_id', 'X-linked' = 'u:fish_idX', 'Y-linked' = 'patriline', 'sigma' = 'units'
      )
    ) %>%
    ggplot(aes(trait2, trait1, size = -log10(p_asymp), fill = r)) +
    geom_point(shape = 21) +
    scale_fill_gradient2(limits = c(-1, 1), na.value = 'white', oob = scales::squish) +
    scale_size_area(max_size = 10, breaks = c(1, 3, 5, 7)) +
    scale_x_discrete(limits = c(paste0('O', 1:7), paste0('B', 1:8))) +
    scale_y_discrete(limits = rev(c(paste0('O', 1:7), paste0('B', 1:8)))) +
    facet_wrap(~param, scales = 'free') +
    geom_abline(slope = -1, intercept = 16, linewidth = 0.5) +
    labs(x = NULL, y = NULL, size = expression(-log[10](italic(p))), fill = expression(italic(r))) +
    theme_classic() +
    theme(strip.placement = 'outside', strip.background = element_blank())
}

size_cors <- read_rds('ornament_analysis/ornament_size_correlations.rds')
viz_cors(size_cors)
ggsave('paper_figures_supplement/ornament_size_correlations.png', width = 10, height = 8)

pa_cors <- read_rds('ornament_analysis/ornament_presence_correlations.rds')
viz_cors(pa_cors)
ggsave('paper_figures_supplement/ornament_presence_correlations.png', width = 10, height = 8)
