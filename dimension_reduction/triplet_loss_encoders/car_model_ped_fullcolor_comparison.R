library(tidyverse)
source('dimension_reduction/triplet_loss_encoders/fit_triplet_loss_model.R')

# I'd like to change the embed_dim
embed_dims <- c(1:6)

purrr::walk(
  embed_dims,
  \(embed_dim) {
    cat('Using', embed_dim, 'embedding dimensions...\n')

    output_name <- paste0('car_model_ped_fullcolor_comparison/embed_dim_', embed_dim)
    if (file.exists(paste0('dimension_reduction/triplet_loss_encoders/weights/', output_name, '.hdf5'))) {
      return(NULL)
    }

    fit_triplet_loss_model(
      img_dir = 'data/carotenoid_coloration_warped/',
      output_name = output_name,
      use_cpu = TRUE,
      embed_dim = embed_dim,
      filter_base = 8,
      margin = 1,

      method = 'pedigree',
      min_r = 1/8,

      input_dim = c(250, 70, 3),
      epochs = 50,
      batch_size = 128,
      n_validation_triplets = 500,
      n_samples = 5000,
      dropout = 0,

      lr_schedule = data.frame(epoch = c(0, 30, 40), lr = 1e-4 * c(1, 0.1, 0.01))
    )
  }
)

source('dimension_reduction/triplet_loss_encoders/triplet_loss_tools.R')
library(furrr)
library(data.table)
library(keras); reticulate::use_virtualenv('r-reticulate')

img_dir <- 'data/carotenoid_coloration_warped/'
input_dim <- c(250, 70, 3)

plan(multisession, workers = 20)
all_images <- load_images_into_array(img_dir, input_dim)
all_images <- drop(all_images)
plan(sequential)

cat('Finished loading data...\n')

source('selection_decisions/compile_decisions.R')
img_names <- dimnames(all_images)[[1]]
db <- data.table::fread('photo_database.csv') %>%
  select(replicate, generation, fish_id, unique_id, facing_direction) %>%
  as_tibble() %>%
  # drop rows for which there is (not yet) an image
  filter(unique_id %in% img_names) %>%
  # reorder to be in the same order as img
  slice(match(img_names, unique_id)) %>%
  # make a variable to denote which images are training, and which validation, also add sample id
  mutate(
    validation = fish_id %in% sample(unique(fish_id), floor(n_distinct(fish_id) * 0.2)),
    sample_id = row_number(),
    fish_id = toupper(fish_id),
    generation = factor(generation, c('parental_gen_1', 'gen_2', 'gen_3'), c('P', 'F1', 'F2'))
  ) %>%
  left_join(select(selection, fish_id, selection, car_perc), 'fish_id')

ped <- data.table::fread('data/pedigree.csv')
ped_df <- ped %>%
  select(animal, dam, sire) %>%
  nadiv::makeA() %>%
  as.matrix() %>% as.data.frame() %>%
  rownames_to_column('id1') %>%
  pivot_longer(-id1, names_to = 'id2', values_to = 'r') %>%
  as.data.table() %>%
  left_join(select(ped, animal, dam1 = dam), c('id1' = 'animal')) %>%
  left_join(select(ped, animal, dam2 = dam), c('id2' = 'animal')) %>%
  mutate(shared_dam = coalesce(dam1 == dam2, FALSE)) %>%
  select(-dam1, -dam2)

eval_fun <- function(
    weights, embed_dim,
    filter_base = 8, input_dim = c(250, 70, 3), margin = 1, dropout = 0
) {
  #weights <- 'dimension_reduction/triplet_loss_encoders/weights/car_model_ped_fullcolor_comparison/embed_dim_2.hdf5'
  #embed_dim <- 2

  model <- build_model(
    embed_dim = embed_dim, filter_base = filter_base, input_dim = input_dim, lr = 1e-3,
    margin = margin, type = 3, dropout = dropout
  )
  load_model_weights_hdf5(model, weights)
  embedder <- keras_model_sequential(layers = get_layers(model, c(1, 4:find_layer(model, 'norm'))))

  trip_gen <- \(r, n) generate_triplet_subset(
    db, n = n, all_images, method = 'pedigree', min_r = r, validation = TRUE, batch_match = FALSE,
    pedigree = ped_df, remove_maternal_effect = FALSE
  )

  map_dfr(
    setNames(1 / 2 ^ (1:6), 1 / 2 ^ (1:6)),
    \(r) trip_gen(r, n = 1e4) %>%
      get_metrics(embedder = embedder, batch_size = 512, margin = margin),
    .id = 'r'
  )
}
out <- map2_dfr(
  list.files('dimension_reduction/triplet_loss_encoders/weights/car_model_ped_fullcolor_comparison', f = TRUE),
  1:6,
  eval_fun,
  .id = 'embed_dim'
)


library(tidyverse)
ggplot(
  out,
  aes(as.numeric(embed_dim), acc, ymin = acc_lower, ymax = acc_upper, color = r)
) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_x_continuous(breaks = 1:10) +
  expand_limits(y = 1) +
  labs(x = 'Nr of embedding dimensions', y = 'Accuracy', color = 'Filter base') +
  theme_classic()
