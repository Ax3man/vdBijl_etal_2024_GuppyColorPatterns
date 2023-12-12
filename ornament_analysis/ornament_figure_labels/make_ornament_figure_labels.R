library(magick)

make_ornament_image <- function(orn, color) {
  file <- glue::glue('ornament_analysis/ornament_images/{orn}.png')

  bg1 <- image_read('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
    image_channel('alpha')
  bg2 <- image_negate(bg1)
  bg3 <- image_colorize(bg1, 100, 'grey60') %>% image_composite(bg2, operator = 'CopyOpacity')

  im1 <- image_read(file) %>% image_channel('alpha') %>% image_threshold()
  im2 <- image_negate(im1)
  im3 <- image_colorize(im1, 100, color) %>% image_composite(im2, operator = 'CopyOpacity')
  out <- image_composite(bg3, im3)
  out_name <- str_remove(orn, '_new')
  image_write(out, glue::glue('ornament_analysis/ornament_figure_labels/{out_name}.png'))
}
walk(paste0('car_', 1:7, '_new'), make_ornament_image, color = 'orange')
walk(paste0('mel_', 1:8), make_ornament_image, color = 'black')
