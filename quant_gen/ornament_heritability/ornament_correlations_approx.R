library(sommer)
library(tidyverse)
library(furrr)

source('quant_gen/prepare_pedigrees.R')

orn <- bind_rows(
  read_rds('ornament_analysis/car_ornaments.rds'),
  read_rds('ornament_analysis/mel_ornaments.rds')
)
pd <- data.table::fread('photo_database.csv') %>% dplyr::select(fish_id, facing_direction, unique_id)

orn_summ <- orn %>%
  left_join(pd, 'unique_id') %>%
  group_by(ornament, fish_id, facing_direction) %>%
  summarise(
    fraction_present = mean(fraction_present),
    present_10 = round(mean(present_10)),
    .groups = 'drop_last'
  ) %>%
  summarise(
    fraction_present = mean(fraction_present),
    present_10 = round(mean(present_10)),
    .groups = 'drop'
  ) %>%
  add_patriline(ped_df) %>%
  mutate(
    fish_idX = fish_id,
    fraction_present = na_if(fraction_present, 0)
  )

bivariate_correlations <- function(d, trait1, trait2) {
  #d <- orn_summ_size; trait1 <- 'car_2'; trait2 <- 'car_6'

  md <- dplyr::select(
    d,
    !!sym(trait1), !!sym(trait2),
    fish_id, fish_idX, patriline
  ) %>%
    drop_na()

  A_small <- A[md$fish_id, md$fish_id]
  X_small <- X[md$fish_id, md$fish_id]

  m <- mmer(
    as.formula(glue::glue('cbind({trait1}, {trait2}) ~ 1')),
    random = ~ vsr(fish_id, Gu = A_small) + vsr(fish_idX, Gu = X_small) + patriline,
    data = md
  )

  vc <- summary(m)$varcomp %>%
    rownames_to_column(var = 'param') %>%
    separate(param, into = c('param', 'traits'), sep = '\\.') %>%
    separate(traits, into = c('trait1', 'trait2'), sep = '-') %>%
    # add asymptotic p-values
    mutate(p_asymp = pnorm(abs(Zratio), lower.tail = FALSE) * 2)
}
safe_bivariate_correlations <- function(...) safely(bivariate_correlations)(...)$result

trait_combs <- expand_grid(trait1 = unique(orn_summ$ornament), trait2 = unique(orn_summ$ornament)) %>%
  # don't correlate traits with themselves, and only do the correlation 1 way.
  filter(
    trait1 != trait2,
    as.numeric(as.factor(trait1)) > as.numeric(as.factor(trait2))
  )

orn_summ_size <- pivot_wider(
  orn_summ, id_cols = -present_10, names_from = ornament, values_from = fraction_present
)
orn_summ_present_10 <- pivot_wider(
  orn_summ, id_cols = -fraction_present, names_from = ornament, values_from = present_10
)

plan(multisession, workers = 12)
# size_correlations <- future_map2_dfr(
#   trait_combs$trait1, trait_combs$trait2,
#   \(x, y) safe_bivariate_correlations(orn_summ_size, x, y),
#   .options = furrr_options(seed = TRUE)
# )
# write_rds(size_correlations, 'ornament_analysis/ornament_size_correlations.rds')

presence_correlations <- future_map2_dfr(
  trait_combs$trait1, trait_combs$trait2,
  \(x, y) safe_bivariate_correlations(orn_summ_present_10, x, y),
  .options = furrr_options(seed = TRUE)
)
write_rds(presence_correlations, 'ornament_analysis/ornament_presence_correlations.rds')
plan(sequential)

