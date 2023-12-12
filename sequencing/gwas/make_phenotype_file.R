library(tidyverse)

# Load the embeddings
em_car_ped <- read_rds(
  'dimension_reduction/triplet_loss_encoders/embeddings/car_model_ped_fullcolor_comparison/embed_dim_5.rds'
)
colnames(em_car_ped) <- paste0('car_PIE_', seq_len(ncol(em_car_ped)))

em_mel_ped <- read_rds(
  'dimension_reduction/triplet_loss_encoders/embeddings/mel_model_ped_fullcolor_comparison/embed_dim_5.rds'
)
em_mel_ped <- em_mel_ped[rownames(em_car_ped), ]
colnames(em_mel_ped) <- paste0('mel_PIE_', seq_len(ncol(em_mel_ped)))

# Load the photo database for total orange and black measures
# We'll include the selection direction (up vs down), in case we want to use it as a covariate
source('selection_decisions/compile_decisions.R')
P <- data.table::fread('photo_database.csv', data.table = FALSE) %>%
  mutate(fish_id = toupper(fish_id)) %>%
  left_join(dplyr::select(selection, fish_id, selection), 'fish_id') %>%
  mutate(selection = ifelse(selection == 'up_selected', 1L, 0L))
p <- dplyr::select(P, unique_id, selection, car_perc, mel_perc_v2)

pm <- as.matrix(p[-1])
rownames(pm) <- p$unique_id
pm <- pm[rownames(em_car_ped), ]

# Load Y-haplogroups
yhap <- data.table::fread('sequencing/Y_haplogroups.csv') %>%
  as.data.frame() %>%
  mutate(fish_id = toupper(fish_id)) %>%
  full_join(dplyr::select(P, unique_id, fish_id), join_by(fish_id)) %>%
  # recode missing values (since otherwise they are dropped)
  mutate(Y_haplogroup = factor(coalesce(Y_haplogroup, "missing"), levels = c('Y1', 'Y2', 'Y3', 'Y4', 'missing')))
rownames(yhap) <- yhap$unique_id
yhap <- yhap[rownames(em_car_ped), ]
# dummy coding the Y-haplogroups, ignore the intercept
yhap_dummy <- model.matrix(~ Y_haplogroup, data = yhap, na.rm = FALSE)[, 2:4]

# Load the ornament data, first presence/absence (pa), then ornament size.
ornaments <- bind_rows(
  read_rds('ornament_analysis/car_ornaments.rds'),
  read_rds('ornament_analysis/mel_ornaments.rds')
)
ornaments_pa <- ornaments %>%
  pivot_wider(id_cols = unique_id, names_from = ornament, values_from = present_10, names_prefix = 'pa_')
ornaments_size <- ornaments %>%
  pivot_wider(id_cols = unique_id, names_from = ornament, values_from = fraction_present, names_prefix = 'size_')
orn <- bind_cols(ornaments_pa[-1], ornaments_size[-1]) %>%
  # select the ornament measures that we actually want to use
  dplyr::select(
    pa_car_1, pa_car_2, pa_car_3, pa_car_4, pa_car_5, pa_car_6, pa_car_7,
    pa_mel_1, pa_mel_2, pa_mel_3, pa_mel_4, pa_mel_5, pa_mel_6, pa_mel_7, pa_mel_8
  ) %>%
  as.matrix()
rownames(orn) <- ornaments_pa$unique_id
orn <- orn[rownames(em_car_ped), ]

# Current structure of out (start column: phenotypes)
# 1:  selection (up vs down)
# 2:  car_perc
# 3:  mel_perc_v2
# 4:6 Y-haplogroup dummies
# 7-11: carotenoid fullcolor 5D pedigree embedding
# 12-16: melanin fullcolor 5D pedigree embedding
# 17-31: presence/absence of ornaments

out <- cbind(pm, yhap_dummy, em_car_ped, em_mel_ped, orn)
write_rds(out, 'sequencing/gwas/phenotypes_untransformed.rds')

# GWAS implementation in GEMMA is very sensitive to normalization, and so we normalize all
# phenotype columns using the inverse-rank transformation
# If not, we get p-values in the order of 10^-1000, which seems unreasonable
binary_variables <- c(TRUE, str_starts(colnames(out)[-1], 'pa_|Y_haplogroup'))
for (f in (1:ncol(out))[!binary_variables]) {
  out[,f] <- qnorm(ecdf(out[,f])(out[,f]) - 0.5/nrow(out))
}
# cowplot::plot_grid(plotlist = map(1:ncol(out), ~qplot(out[, .x])))
# pairs(out[, 4:9]); cor(out[, 5:7])
# pairs(out[, 10:14]); cor(out[, 8:10])
# pairs(out[, 15:19]); cor(out[, 11:13])
# pairs(out[, 20:24]); cor(out[, 14:16])

stopifnot(ncol(out) == (3 + 3 + 5 * 2 + 7 + 8))

write_rds(out, 'sequencing/gwas/phenotypes.rds')





# incl <- data.table::fread('photo_database.csv') %>% filter(fish_id %in% sampling$fish_id) %>%
#   dplyr::select(unique_id, fish_id)
# to_plot <- ornaments %>%
#   inner_join(incl, join_by(unique_id)) %>%
#   group_by(fish_id, ornament) %>%
#   summarise(across(pixels_present:present, mean)) %>%
#   mutate(present = round(present))
#
# cowplot::plot_grid(
#   ggplot(to_plot, aes(ornament, fill = present == 1)) + geom_bar(show.legend = FALSE),
#   ggplot(to_plot, aes(ornament, fill = fraction_present > 0.1)) + geom_bar(show.legend = FALSE),
#   ncol = 1
# )

to_plot %>%
  ggplot(aes(fraction_present, after_stat(scaled))) +
  geom_density(adjust = 1, data = . %>% filter(fraction_present > 0.1), color = 'navy') +
  #geom_density(adjust = 1, data = . %>% filter(fraction_present > 0), color = 'black') +
  facet_wrap(~ornament) + scale_x_log10()
