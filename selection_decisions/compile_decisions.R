library(tidyverse)

gen1 <- bind_rows(
  '1' = data.table::fread('selection_decisions/finished/rep1_parentalgen.csv'),
  '2' = data.table::fread('selection_decisions/finished/rep2_parentalgen.csv'),
  '3' = data.table::fread('selection_decisions/finished/rep3_parentalgen.csv'),
  .id = 'replicate'
) %>%
  mutate(
    generation = 'P',
    replicate = paste0('replicate_', replicate),
    decisions = case_when(
      selection == 'not_selected' ~ 'sacrifice',
      selection == 'up_selected' ~ 'move to up',
      selection == 'down_selected' ~ 'move to down'
    ),
    selection = 'not_selected',
  )

gen2 <- bind_rows(
  '1' = data.table::fread('selection_decisions/finished/rep1_gen2_previous_decisions.csv'),
  '2' = data.table::fread('selection_decisions/finished/rep2_gen2_previous_decisions.csv'),
  '3' = data.table::fread('selection_decisions/finished/rep3_gen2_previous_decisions.csv'),
  .id = 'replicate'
) %>%
  mutate(generation = 'F1', replicate = paste0('replicate_', replicate))

gen3 <- bind_rows(
  '1' = data.table::fread('selection_decisions/finished/rep1_gen3_previous_decisions.csv'),
  '2' = data.table::fread('selection_decisions/finished/rep2_gen3_previous_decisions.csv'),
  '3' = data.table::fread('selection_decisions/finished/rep3_gen3_previous_decisions.csv'),
  .id = 'replicate'
) %>%
  mutate(generation = 'F2', replicate = paste0('replicate_', replicate))

gen4 <- bind_rows(
  '1' = data.table::fread('selection_decisions/finished/rep1_gen4_previous_decisions.csv'),
  '2' = data.table::fread('selection_decisions/finished/rep2_gen4_previous_decisions.csv'),
  '3' = data.table::fread('selection_decisions/finished/rep3_gen4_previous_decisions.csv'),
  .id = 'replicate'
) %>%
  mutate(generation = 'F3', replicate = paste0('replicate_', replicate))

selection <- bind_rows(as_tibble(gen1), as_tibble(gen2), as_tibble(gen3), as_tibble(gen4)) %>%
  mutate(
    across(c(fish_id, sire, dam), toupper),
    phenotype_date = lubridate::as_date(as.character(date), format = '%Y%m%d'),
    generation = factor(generation, c('P', 'F1', 'F2', 'F3')),
    selection2 = case_when(
      generation != 'P' ~ selection,
      decisions == 'move to down' ~ 'down_selected',
      decisions == 'move to up' ~ 'up_selected',
      TRUE ~ 'not_selected'
    )
  ) %>%
  dplyr::select(
    replicate, generation, fish_id, sire, dam, phenotype_date, selection, decisions,
    car_perc, selection2
  )

rm(gen1, gen2, gen3, gen4)
