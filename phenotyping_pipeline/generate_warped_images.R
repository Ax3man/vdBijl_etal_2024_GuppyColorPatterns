## Setup & prep ------------------------------------------------------------------------------------

# New resolution; starting nr of pixels is 800x300, will be downsampled to this many:
#x_bins <- 320
#y_bins <- 80

x_bins <- 500
y_bins <- 140

source('phenotyping_pipeline/process_full_landmarks.R')
library(data.table)
setDTthreads(1L)
library(tidyverse)
library(geomorph)
suppressMessages(library(imager))
library(furrr)
plan(multisession, workers = 16)

# !WARNING!: Run this to delete all previously generated results: takes a long time to re-generate
# invisible(file.remove(
#   list.files('data/extracted_fish_warped', recursive = TRUE, full.names = TRUE),
#   list.files('data/carotenoid_coloration_warped', recursive = TRUE, full.names = TRUE),
#   list.files('data/melanic_coloration_warped', recursive = TRUE, full.names = TRUE),
#   list.files('data/melanic_coloration_warped_v2', recursive = TRUE, full.names = TRUE)
# ))

# create vectors of all the file names we'll need
# Images of the extracted fish
ims <- list.files('data/extracted_fish', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

# to run small tests, randomly choose some images
#ims <- sample(ims, 200)

# Their shape files
shp <- paste0('data/full_landmarks/', names(ims), '.rds')

# remove those with missing landmarks
ims <- ims[file.exists(shp)]
shp <- shp[file.exists(shp)]

# remove those with already warped images
compl <- file.exists(str_replace(ims, 'extracted_fish', 'extracted_fish_warped')) &
  file.exists(str_replace(ims, 'extracted_fish', 'carotenoid_coloration_warped')) &
  # turned this off, since rep 3, F3 has no melanic V1
  #file.exists(str_replace(ims, 'extracted_fish', 'melanic_coloration_warped')) &
  file.exists(str_replace(ims, 'extracted_fish', 'melanic_coloration_warped_v2'))
ims <- ims[!compl]; shp <- shp[!compl]

cat(length(ims), 'images need warping.\n')

## Shape files & consensus -------------------------------------------------------------------------
combine_shapes <- function(..., shape_list = NULL) {
  if (is.null(shape_list)) shape_list <- list(...)

  fixed <- map(shape_list, 'landmarks.pixel')
  landmarks.pixel <- do.call(abind::abind, list(fixed, along = 3))
  curves.pixel <- map(shape_list, 'curves.pixel')
  out <- list(landmarks.pixel = landmarks.pixel, curves.pixel = curves.pixel)
  class(out) <- 'shapes'
  return(out)
}

# load all shape files
all_shapes <- map(shp, read_rds)
failing <- !map_lgl(all_shapes, ~all(map_dbl(.x$curves.pixel, nrow) > 3))
if (any(failing)) {
  stop(
    'These landmark files have signifant problems:\n',
    paste(shp[failing], collapse = '\n')
  )
  if (FALSE) {
    # run this to delete the affected files
    file.remove(shp[failing])
    file.remove(gsub('.rds', '.txt', gsub('full', 'fixed', shp[failing])))
  }
}

my_shape <- combine_shapes(shape_list = all_shapes)
land <- readland.shapes(my_shape, c(50, 50, 50))

# load consensus shape
ref <- read_rds('data/consensus_shape.rds')
# get the range of consensus shape coordinates, to calculate bins, give 1% extra room
ref_range <- range(ref[, 1]) * 1.01
# calculate the bins we'll need later
bins <- seq(ref_range[1], ref_range[2], length.out = x_bins)

## Helper functions --------------------------------------------------------------------------------

transform_picture <- function(pic_df, landmarks, consensus) {
  # compute transformation from landmarks to consensus
  trans <- Morpho::computeTransform(consensus, landmarks, type = "tps")
  # apply to all image coordinates
  warped <- Morpho::applyTransform(as.matrix(pic_df[1:2]), trans)
  # overwrite coordinates
  mutate(pic_df, x = warped[, 1], y = warped[, 2])
}

warp_image <- function(image_file, landmarks, mirror, consensus) {
  im <- load.image(image_file)
  #plot(im); text(landmarks, 1:nrow(landmarks))
  if (mirror) im <- mirror(im, 'x')
  img_df <- as.data.frame(im, wide = 'c')
  #img_df <- img_df[img_df$c.4 > 0, ]
  transform_picture(img_df, landmarks, consensus)
}

rasterize <- function(d, bins, method) {
  # d <- warped_im; method <- 'color_average'
  # method must be 'color_average' or 'mode'
  method <- match.arg(method, c('color_average', 'mode'))
  cut_to_middle <- function(v, bins) {
    as.numeric(as.character(cut(v, c(bins, Inf), bins + 0.5 * diff(bins[1:2]))))
  }
  cut_to_int <- function(v, bins) {
    as.numeric(as.character(cut(v, c(bins, Inf), seq_along(bins))))
  }
  require(data.table)
  dt <- setDT(d)
  dt[, X := cut_to_int(x, bins)]
  dt[, x := cut_to_middle(x, bins)]
  dt[, Y := cut_to_int(y, bins)]
  dt[, y := cut_to_middle(y, bins)]
  setkey(dt, x, y, X, Y)
  if (method == 'color_average') {
    dt <- dt[
      ,
      lapply(.SD, function(v) sqrt(mean(v ^ 2))), by = .(x, y, X, Y),
      .SDcols = paste0('c.', 1:4)
    ]
    #d[, color := rgb(c.1, c.2, c.3, c.4)]
  }
  if (method == 'mode') {
    stop('Method not implemented.')
  }
  dt %>%
    data.frame() %>%
    drop_na() %>%
    filter(c.4 > 0) %>%
    select(X, Y, c.1:c.4) %>%
    mutate(Y = Y - floor((x_bins - y_bins) / 2)) %>%
    complete(X = 1:x_bins, Y = 1:y_bins, fill = list(c.1 = 0, c.2 = 0, c.3 = 0, c.4 = 0)) %>%
    pivot_longer(c.1:c.4, names_to = 'cc', values_to = 'value') %>%
    drop_na() %>%
    mutate(z = 1, cc = as.numeric(factor(cc))) %>%
    select(x = X, y = Y, z, cc, value) %>%
    as.cimg(dims = c(x_bins, y_bins, 1, 4))
}

#warp_and_rasterize(ims[x], land$landmarks[[x]], facing_left[x], consensus = ref, bins = bins)
warp_and_rasterize <- function(image_file, landmarks, mirror, consensus, bins) {
  # xx <- which(names(ims) == '20220526_IMG_9666')
  # image_file <- ims[xx]; landmarks <- land$landmarks[[xx]]; mirror <- facing_left[xx]
  # consensus <- ref; bins <- bins

  # s <- sample(length(ims), 1)
  # image_file <- ims[s]; landmarks <- land$landmarks[[s]]; mirror <- facing_left[s];
  # consensus <- ref; bins <- bins

  # plot(load.image(image_file)); text(landmarks, labels = 1:148, col = 'red')

  # The extracted fish:
  warped_im <- warp_image(image_file, landmarks, mirror, consensus)
  # if (any(warped_im[, 1:2] < head(bins, 1)) | any(warped_im[, 1:2] > tail(bins, 1))) {
  #   stop('This file borked: ', image_file)
  # }
  rast_im <- rasterize(warped_im, bins, 'color_average')

  # ggplot(warped_im %>% filter(c.4 == 1), aes(x, y, color = rgb(c.1, c.2, c.3, c.4))) + geom_point() +
  #   scale_color_identity() + scale_y_reverse() + coord_fixed()

  # Orange coloration:
  warped_car <- warp_image(
    str_replace(image_file, 'extracted_fish', 'carotenoid_coloration'),
    landmarks, mirror, consensus
  )
  rast_car <- rasterize(warped_car, bins, 'color_average')

  # Black coloration:
  old_mel_file <- str_replace(image_file, 'extracted_fish', 'melanic_coloration')
  old_mel_exists <- file.exists(old_mel_file)
  if (old_mel_exists) {
    warped_mel <- warp_image(old_mel_file, landmarks, mirror, consensus)
    rast_mel <- rasterize(warped_mel, bins, 'color_average')
  }

  # Black coloration v2:
  warped_mel_v2 <- warp_image(
    str_replace(image_file, 'extracted_fish', 'melanic_coloration_v2'),
    landmarks, mirror, consensus
  )
  rast_mel_v2 <- rasterize(warped_mel_v2, bins, 'color_average')

  # Now that all three succeeded without error, save to disk:
  save.image(
    rast_im, quality = 1,
    str_replace(image_file, 'extracted_fish', 'extracted_fish_warped')
  )
  save.image(
    rast_car, quality = 1,
    str_replace(image_file, 'extracted_fish', 'carotenoid_coloration_warped')
  )
  if (old_mel_exists) {
    save.image(
      rast_mel, quality = 1,
      str_replace(image_file, 'extracted_fish', 'melanic_coloration_warped')
    )
  }
  save.image(
    rast_mel_v2, quality = 1,
    str_replace(image_file, 'extracted_fish', 'melanic_coloration_warped_v2')
  )

  # opa <- par(mfrow = c(3, 1))
  # plot(rast_im, re = F, axes = F); plot(rast_car, re = F, axes = F); plot(rast_mel, re = F, axes = F)
  # par(opa)
  return(invisible(NULL))
}

## Running the gauntlet ----------------------------------------------------------------------------

# load the photo database to get the facing directions, needed for mirroring the images
pd <- data.table::fread('photo_database.csv')
facing_left <- pd$facing_direction[match(names(ims), pd$unique_id)] == 'l'

if (length(ims) > 0) {
  cat('Warping images....\n')
  warped <- future_pwalk(
    list(ims, land$landmarks, facing_left),
    warp_and_rasterize, consensus = ref, bins = bins,
    .progress = TRUE, .options = furrr_options(seed = NULL)
  )

  # when having trouble, debug:
  # warp_and_rasterize(ims[[1]], land$landmarks[[1]], facing_left[[1]], ref, bins)
  # An "unable to create..." error likely means the required folder does not exist.
  # Other times, the facing direction may have been wrongly assigned.
  # On rare occasions, the full landmarks are wrong and just need to be regenerated. TSP can be off?

  cat('\nFinished.\n')
} else {
  cat('No images to warp.\n')
}

plan(sequential)

#
# x <- 1
# pwalk(
#   list(ims[x], land$landmarks[x], facing_left[x]),
#   warp_and_rasterize, consensus = ref, bins = bins
# )

