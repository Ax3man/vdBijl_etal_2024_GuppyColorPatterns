library(tidyverse)

source('quant_gen/prepare_pedigrees.R')
pd <- data.table::fread('photo_database.csv', data.table = FALSE)

patrilines <- pd %>%
  dplyr::select(fish_id) %>%
  distinct() %>%
  add_patriline(ped_df)

data.table::fwrite(patrilines, 'data/patrilines.csv')
