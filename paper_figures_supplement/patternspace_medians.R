# load the embeddings etc. using the same UMAP as in figure 2
source('paper_figures/subplots/embedding_selection.R')

n_x <- 10
#n_y <- n_x * 4

x_img <- 500
y_img <- 140

x_breaks <- seq(min(embeddings$UMAP1), max(embeddings$UMAP1), l = n_x + 1)
x_labs <- ((lag(x_breaks) + x_breaks) / 2)[-1]
step_size <- mean(diff(x_labs))
y_breaks <- seq(min(embeddings$UMAP2), max(embeddings$UMAP2), step_size / x_img * y_img)
y_labs <- ((lag(y_breaks) + y_breaks) / 2)[-1]

dx <- mean(diff(x_labs)) / x_img
dy <- mean(diff(y_labs)) / y_img

combine_images <- function(files, enhance_contrast, add_background) {
  x <- lapply(files, imager::load.image)
  #a <- imager::average(x)
  a <- imager::parmed(x)
  a <- imager::mirror(a, 'y')
  tmp <- tempfile(fileext = '.png')
  imager::save.image(a, tmp)
  if (enhance_contrast) magick::image_read(tmp) %>% magick::image_contrast() %>% magick::image_write(tmp)
  if (add_background) {
    c(
      magick::image_read('data/extracted_fish_warped/replicate_3/gen_2/20210602_IMG_2097.png') %>%
        magick::image_flip() %>%
        magick::image_channel('alpha') %>%
        magick::image_transparent('#111111') %>%
        magick::image_fill('grey80', magick::geometry_point(x_img/2, y_img/2), fuzz = 90),
      magick::image_read(tmp)
    ) %>% magick::image_flatten() %>% magick::image_write(tmp)
  }

  a <- imager::load.image(tmp)
  as.data.frame(a, wide = 'c')
}

mean_images <- embeddings %>%
  mutate(
    cV1 = cut(UMAP1, breaks = x_breaks, labels = x_labs) %>% as.character() %>% as.numeric(),
    cV2 = cut(UMAP2, breaks = y_breaks, labels = y_labs) %>% as.character() %>% as.numeric()
  ) %>%
  group_by(cV1, cV2) %>%
  # only locations with at least 10 images
  filter(n() > 10) %>%
  summarise(
    mean_image = combine_images(image_file, enhance_contrast = FALSE, add_background = TRUE),
    .groups = 'drop'
  ) %>%
  unnest(mean_image) %>%
  ungroup()

eb <- element_blank()
# FOR FULL FISH
orange_emsp <- mean_images %>%
  drop_na() %>%
  mutate(x_ = cV1 + (x - x_img/2) * dx, y_ = cV2 + (y - y_img/2) * dy) %>%
  ggplot(aes(x_, y_, fill = rgb(c.1, c.2, c.3, c.4))) +
  geom_raster(show.legend = FALSE) +
  scale_fill_identity() +
  labs(x = 'Pattern axis 1', y = 'Pattern axis 2')


## black
black_embedding <- read_rds(
  'dimension_reduction/triplet_loss_encoders/embeddings/mel_model_ped_fullcolor_comparison/embed_dim_5.rds'
) %>%
  magrittr::set_colnames(paste0('V', 1:ncol(.))) %>%
  as.data.frame() %>%
  rownames_to_column('unique_id') %>%
  as_tibble()

black_umap_model <- umap(black_embedding[-1], ret_model = TRUE, fast_sgd = TRUE, spread = 2) #, n_neighbors = 50, spread = 5)

black_embeddings <- bind_cols(
  black_embedding,
  black_umap_model$embedding %>% as.data.frame() %>% setNames(paste0('UMAP', 1:2))
) %>%
  left_join(select(photo_database, fish_id, unique_id, facing_direction, image_file), 'unique_id') %>%
  mutate(image_file = str_replace(image_file, 'carotenoid_coloration_warped', 'melanic_coloration_warped_v2'))

black_mean_images <- black_embeddings %>%
  mutate(
    cV1 = cut(UMAP1, breaks = x_breaks, labels = x_labs) %>% as.character() %>% as.numeric(),
    cV2 = cut(UMAP2, breaks = y_breaks, labels = y_labs) %>% as.character() %>% as.numeric()
  ) %>%
  group_by(cV1, cV2) %>%
  # only locations with at least 10 images
  filter(n() > 10) %>%
  summarise(
    mean_image = combine_images(image_file, enhance_contrast = FALSE, add_background = TRUE),
    .groups = 'drop'
  ) %>%
  unnest(mean_image) %>%
  ungroup()

black_emsp <- black_mean_images %>%
  drop_na() %>%
  mutate(x_ = cV1 + (x - x_img/2) * dx, y_ = cV2 + (y - y_img/2) * dy) %>%
  ggplot(aes(x_, y_, fill = rgb(c.1, c.2, c.3, c.4))) +
  geom_raster(show.legend = FALSE) +
  scale_fill_identity() +
  labs(x = 'Pattern axis 1', y = 'Pattern axis 2')

combined <- ((orange_emsp + coord_fixed()) / (black_emsp + coord_fixed())) +
  plot_annotation(tag_levels = 'A')
ggsave('paper_figures_supplement/patternspace_medians.png', combined, width = 12, height = 20,
       units = 'cm', dpi = 600)
