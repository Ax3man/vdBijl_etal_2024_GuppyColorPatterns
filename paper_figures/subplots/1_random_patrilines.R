library(tidyverse)
library(magick)
source('paper_figures/theme.R')

set.seed(12345)

source('quant_gen/prepare_pedigrees.R')

pd <- data.table::fread('photo_database.csv')
F3s <- pd %>%
  filter(generation == 'gen_4') %>%
  dplyr::select(fish_id) %>%
  distinct() %>%
  pull(fish_id)

paths <- list.files('data/extracted_fish_warped', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

grab_image <- function(fish_id) {
  pd %>%
    filter(.data$fish_id == .env$fish_id, facing_direction == 'r') %>%
    slice_sample(n = 1) %>%
    pull(unique_id) %>%
    { image_read(paths[.]) } %>%
    # Our photos have been captured for analysis, avoiding clipping as much as possible. However,
    # this makes for very bland image displays. So we slightly bump the brightness by 10% and
    # the contrast by 10%.
    image_modulate(brightness = 110) %>% image_contrast(sharpen = 1.1)
}

patriline_display <- function(fish_id) {
  ims <- grab_image(fish_id)
  for (i in 2:4) {
    fish_id <- ped_df$sire[ped_df$animal == fish_id]
    ims <- c(grab_image(fish_id), ims)
  }
  image_append(ims)
}

stacked_display <- lapply(sample(F3s, 8), patriline_display) %>% do.call(c, .)
x <- image_append(stacked_display, stack = TRUE)

tmpfile <- tempfile(fileext = '.png')
image_write(x, tmpfile)

im <- imager::load.image(tmpfile) %>%
  as.data.frame(wide = 'c')
example_patrilines <- ggplot(im, aes(x, y, fill = rgb(c.1, c.2, c.3, c.4))) +
  geom_raster() +
  scale_fill_identity() +
  scale_y_reverse(name = 'Patriline') +
  scale_x_continuous(
    breaks = 250 + 0:3 * 500,
    labels = c('P', expression(F[1]), expression(F[2]), expression(F[3])),
    position = 'top'
  ) +
  coord_equal(expand = FALSE) +
  theme_void(base_size = 7, base_family = 'Arial') +
  theme(
    axis.text.x = element_text(),
    axis.title.y = element_text(angle = 90)
  )

