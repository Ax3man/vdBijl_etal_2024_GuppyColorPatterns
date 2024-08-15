library(tidyverse)
library(furrr)
library(uwot) # for UMAP
source('paper_figures/theme.R')
source('selection_decisions/compile_decisions.R')

photo_database <- data.table::fread('photo_database.csv') %>%
  select(replicate, generation, fish_id, facing_direction, unique_id) %>%
  as_tibble() %>%
  mutate(
    image_file = paste0('data/carotenoid_coloration_warped/', replicate, '/', generation, '/', unique_id, '.png')
  )
##
embedding <- read_rds(
  'dimension_reduction/triplet_loss_encoders/embeddings/car_model_ped_fullcolor_comparison/embed_dim_5.rds'
) %>%
  magrittr::set_colnames(paste0('V', 1:ncol(.))) %>%
  as.data.frame() %>%
  rownames_to_column('unique_id') %>%
  as_tibble()

umap_model <- umap(embedding[-1], ret_model = TRUE, fast_sgd = TRUE, spread = 2) #, n_neighbors = 50, spread = 5)

embeddings <- bind_cols(
  embedding,
  umap_model$embedding %>% as.data.frame() %>% setNames(paste0('UMAP', 1:2))
) %>%
  left_join(select(photo_database, fish_id, unique_id, facing_direction, image_file), 'unique_id')

emb_df <- embeddings %>%
  group_by(fish_id) %>%
  summarise(
    across(starts_with('UMAP'), mean),
    image_file = sample(image_file, 1)
  ) %>%
  left_join(selection %>% mutate(fish_id = tolower(fish_id)), 'fish_id')

embedding_selection <- ggplot(
  slice_sample(emb_df, prop = 1),
  aes(UMAP1, UMAP2, color = ifelse(generation == 'P', 'Stock', selection))
) +
  geom_point(alpha = 1, size = 0.1) +
  #stat_ellipse(linewidth = 0.2) +
  # geom_text(
  #   aes(x = 2, y = -15, label = paste('n =', n)), count(emb_df, generation),
  #   color = 'grey40', size = 3
  # ) +
  scale_color_manual(
    name = NULL,
    values = c(Stock = 'grey40', down_selected = 'navy', up_selected = 'firebrick'),
    breaks = c('Stock', 'down_selected', 'up_selected'),
    labels = c('Stock', 'Down-selected', 'Up-selected'),
    guide = guide_legend(override.aes = list(size = 0.3))
  ) +
  facet_wrap(vars(generation), nrow = 1) +
  coord_fixed() +
  #theme_classic() +
  labs(x = 'Pattern axis 1', y = 'Pattern axis 2') +
  theme(
    legend.position = 'top', strip.background = element_rect(color = NA, fill = 'grey90'),
    strip.text = element_text(size = 6), legend.box.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0), plot.margin = margin(0,0,0,0)
  )
