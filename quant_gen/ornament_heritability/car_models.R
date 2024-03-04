library(tidyverse)
library(brms)
library(furrr)

source('quant_gen/prepare_pedigrees.R')
car_out <- read_rds('ornament_analysis/car_ornaments.rds')
car_out_wide <- pivot_wider(car_out, id_cols = unique_id, names_from = ornament, values_from = pixels_present)

# Get the carotenoid values for each photo
d <- data.table::fread('photo_database.csv') %>%
  as_tibble() %>%
  select(replicate, generation, fish_id, facing_direction, unique_id) %>%
  mutate(fish_idA = fish_id, fish_idX = fish_id) %>%
  # Add the patriline information to the dataset, as well as dam and sire id's
  add_patriline(ped_df) %>%
  left_join(ped_df, by = c('fish_id' = 'animal')) %>%
  left_join(car_out_wide, 'unique_id')

## ACTUAL FITTING ----------------------------------------------------------------------------------

fit_ornament_hurdle <- function(ornament, data, A, X) {
  out_file <- paste0('quant_gen/ornament_heritability/saved_models/', ornament, '.rds')

  f <- bf(
    as.formula(paste(
      ornament,
      '~ 1 + (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) +
        (1 | patriline) + (1 | grow_tank) + (1 | fish_id)'
    )),
    hu ~ 1 + (1 |gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) +
      (1 | patriline) + (1 | grow_tank) + (1 | fish_id),
    family = 'hurdle_lognormal'
  )
  my_priors <- c(
    prior(logistic(0, 5), class = 'Intercept', dpar = 'hu'),
    prior(student_t(3, 0, 5), class = 'sd')
  )

  m <- brm(
    formula = f, data = data, data2 = list(A = A, X = X),
    prior = my_priors,
    backend = 'cmdstanr',
    chains = 4, cores = 4, warmup = 1000, iter = 5000,
    control = list(adapt_delta = 0.98, max_treedepth = 10),
    file = out_file, file_refit = 'on_change'
  )
}

plan(multisession, workers = 7)
future_walk(
  paste0('car_', 1:7),
  fit_ornament_hurdle, data = d, A = pedA, X = pedX,
  .options = furrr_options(seed = TRUE)
)
plan(sequential)

if (FALSE) { # run interactively only

  library(tidyverse)
  library(tidybayes)
  library(imager)

  finished_models <- list.files('quant_gen/ornament_heritability/saved_models', f = TRUE) %>%
    setNames(., basename(.) %>% tools::file_path_sans_ext()) %>%
    map(read_rds)

  ### Model critique -------------------------------------------------------------------------------
  map(finished_models, summary)

  map(
    finished_models,
    \(m) brms::pp_check(m, nd = 100, type = 'ecdf_overlay')
  ) %>% cowplot::plot_grid(plotlist = .)

  M <- finished_models[[2]]
  yrep <- posterior_predict(M, ndraws = 10)
  y <- get_y(M)

  yrep %>% t() %>% as.data.frame() %>% as_tibble() %>%
    mutate(
      sample = row_number(),
      y = y
    ) %>%
    pivot_longer(-sample) %>%
    ggplot(aes(x = name, fill = value == 0)) +
    geom_bar() + facet_grid(cols = vars(name == 'y'), space = 'free', scales = 'free')

  yrep %>% t() %>% as.data.frame() %>% as_tibble() %>%
    mutate(
      sample = row_number(),
      y = y
    ) %>%
    pivot_longer(-sample) %>%
    filter(value > 0) %>%
    ggplot(aes(x = value, color = name == 'y', group = name)) +
    geom_density(fill = NA) +
    scale_x_log10(breaks = c(1, 100, 10000))

  ### Plot estimates -------------------------------------------------------------------------------

  # do image prep for fancy figure labels
  bg <- (1 - (load.image('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
                channel(4))) %>% add.color()
  recolor_image <- function(im, old_cols, new_cols) {
    for (i in seq_along(old_cols)) {
      im <- colorise(
        im,
        R(im) == old_cols[[i]][1] & G(im) == old_cols[[i]][2] & B(im) == old_cols[[i]][3],
        new_cols[[i]]
      )
    }
    im
  }
  car_pixsets <- list.files('ornament_analysis/ornament_images', 'car_._new\\.png', full.names = TRUE) %>%
    map(\(.x) {
      out <- load.image(.x) %>% channel(4) %>% threshold(0.5) %>% as.cimg() %>% add.color()
      G(out) <- 0; B(out) <- 0
      as.pixset(out)
    }) %>%
    setNames(paste0('car_', seq_along(.)))
  all_car_ornaments <- reduce(car_pixsets, `+`, .init = bg)
  all_car_ornaments <- recolor_image(
    all_car_ornaments,
    list(c(1, 0, 0), c(0, 0, 0)), list(c(0.8, 0.5, 0), 'grey70')
  )
  plot(all_car_ornaments)
  save.image(all_car_ornaments, 'ornament_analysis/visualization/car_ornaments_abstract.png')
  car_ornaments <- map(car_pixsets, ~bg + .x)
  car_ornaments <- map(
    car_ornaments, recolor_image,
    list(c(1, 0, 0), c(0, 0, 0)), list(c(0.8, 0.5, 0), 'grey70')
  )
  fn <- paste0('ornament_analysis/ornament_minis/', names(car_ornaments), '.png')
  walk2(car_ornaments, fn, ~save.image(.x, .y))
  labs <- paste0("<img src='", fn, "'width='75' />") %>% setNames(paste0('car_', seq_along(fn)))
  labs_large <- paste0("<img src='", fn, "'width='125' />") %>% setNames(paste0('car_', seq_along(fn)))

  all_draws <- map_dfr(
    finished_models,
    \(model) {
      gather_draws(
        model,
        sd_fish_id__Intercept,
        sd_fish_idA__Intercept,
        sd_fish_idX__Intercept,
        sd_patriline__Intercept,
        sd_grow_tank__Intercept,

        sd_fish_id__hu_Intercept,
        sd_fish_idA__hu_Intercept,
        sd_fish_idX__hu_Intercept,
        sd_patriline__hu_Intercept,
        sd_grow_tank__hu_Intercept,

        sigma
      ) %>%
        ungroup() %>%
        mutate(
          # give clearer names
          .variable = fct_recode(
            factor(.variable),
            size_env = 'sd_fish_id__Intercept',
            size_auto = 'sd_fish_idA__Intercept',
            size_X = 'sd_fish_idX__Intercept',
            size_Y = 'sd_patriline__Intercept',
            size_tank = 'sd_grow_tank__Intercept',

            ap_env = 'sd_fish_id__hu_Intercept',
            ap_auto = 'sd_fish_idA__hu_Intercept',
            ap_X = 'sd_fish_idX__hu_Intercept',
            ap_Y = 'sd_patriline__hu_Intercept',
            ap_tank = 'sd_grow_tank__hu_Intercept',

            size_sigma = 'sigma'
          ),
          # change to variances (from standard deviations)
          .value = .value ^ 2
        ) %>%
        separate(.variable, into = c('.response', '.variable'), sep = '_') %>%
        pivot_wider(names_from = '.variable', values_from = '.value') %>%
        mutate(
          .response = fct_recode(factor(.response), 'Absence/Presence' = 'ap', 'Size' = 'size'),

          # encode sigma to be the logit link variance, but only for the presence/absence part
          sigma = coalesce(sigma, (pi ^ 2) / 3),

          Va = auto + X + Y,
          Ve = env + tank,
          Vtotal = Va + Ve + sigma,

          h2 = Va / (Va + Ve),
          h2_auto = auto / (Va + Ve),
          h2_X = X / (Va + Ve),
          h2_Y = Y / (Va + Ve),

          rel_env_effect = env / (Va + Ve),
          rel_tank_effect = tank / (Va + Ve),

          measurement_error = sigma / Vtotal
        )
    },
    .id = 'ornament'
  ) %>%
    mutate(labels = labs[ornament], large_labels = labs_large[ornament])

  draws_long <- all_draws %>%
    select(ornament, .response, labels, large_labels, h2:measurement_error) %>%
    pivot_longer(cols = h2:measurement_error)

  ggplot(draws_long, aes(value, fct_rev(name), color = .response)) +
    stat_pointinterval(position = position_dodge(width = 0.3)) +
    facet_wrap(~ large_labels) +
    labs(x = NULL, y = NULL, color = NULL) +
    theme_bw() +
    theme(
      strip.text = ggtext::element_markdown(), strip.background = element_blank(),
      legend.position = 'top'
    )

  draws_long %>%
    filter(str_starts(name, 'h2')) %>%
    ggplot(aes(value, fct_rev(name), color = fct_rev(.response))) +
    stat_pointinterval(position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 3.5, linewidth = 0.3) +
    scale_y_discrete(labels = c(
      'h2' = expression(italic(h^2)), 'h2_auto' = expression(italic(h[auto]^2)),
      'h2_X' = expression(italic(h[X]^2)), 'h2_Y' = expression(italic(h[Y]^2))
    )) +
    scale_x_continuous(lim = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(
      breaks = c('Absence/Presence', 'Size'),
      labels = c('Absence/Presence (logistic)', 'Ornament size (lognormal)'),
      values = c('black', 'grey40')
    ) +
    labs(x = NULL, y = NULL, color = "Hurdle component") +
    coord_cartesian(clip = 'off') +
    facet_wrap(~ labels, scales = 'free', ncol = 4) +
    theme_classic() +
    theme(
      strip.text = ggtext::element_markdown(), strip.background = element_blank(),
      legend.position = c(7/8, 1/4)
    )

  ggsave('ornament_analysis/visualization/car_ornament_h2.png', w = 8, h = 4, scale = 1.2)
}
