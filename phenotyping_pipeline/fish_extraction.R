suppressPackageStartupMessages(
  {
    library(imager)
    library(tensorflow)
    library(keras)
    library(tidyverse)
  }
)
# tf_version() ## avoid problems with not finding the tensorflow install on startup
# suppressMessages(gpu <- tf$config$experimental$get_visible_devices('GPU')[[1]])
# virtual_gpu <- tf$config$experimental$VirtualDeviceConfiguration(memory_limit = 5000)
# tf$config$experimental$set_virtual_device_configuration(gpu, list(virtual_gpu))

fish_extractor <- unserialize_model(
  read_rds('../color_analyses/keras_bg_removal/unet.rds'),
  custom_objects = list(dice_loss = dice_loss, dice_coef = dice_coef)
)

source('phenotyping_pipeline/tools.R')

extract_fish1 <- function(im) {
  # 1: Obtain a mask, and the correct cropped image
  mask <- predict_mask(im)
  # plot(as.cimg(mask[[2]]))
  return(mask)
}

extract_fish2 <- function(im, mask) {
  if (length(mask$best_crop) == 1) {
    crop <- crop5(im, mask$best_crop)
  } else {
    crop <- imsub(im, x %inr% mask$best_crop[c('x0', 'x1')], y %inr% mask$best_crop[c('y0', 'y1')])
  }

  segmented <- mask_and_box_image(
    crop,
    mask$mask
  )
  return(segmented)
}

extract_fish_manual <- function(im, landmarks, scale) {
  # 1: Obtain a mask from the landmarks
  mask <- produce_mask(landmarks$curves.pixel, scale = scale, w = width(im), h = height(im))
  # 2: Use the mask to extract out the fish
  segmented <- mask_and_box_image(im, mask)
  return(segmented)
}

predict_mask <- function(im, input_width = 256 * (3/2), input_height = 256) {
  # input_width and input_height refer to the dimensions the convnet accepts
  fives <- map(
    1:5,
    ~crop5(im, .x) %>% resize(input_width, input_height, interpolation_type = 5L)
  )
  get_model_prediction <- function(five) {
    arr <- aperm(five, c(3, 1, 2, 4))
    predict(fish_extractor, arr)[1, , , ]
  }
  masks <- map(fives, get_model_prediction)
  # display(map(masks, as.cimg))
  crosses_boundary <- map_lgl(masks, detect_mask_on_boundary)
  # get squared error of distance from center mass to center of the picture
  if (!all(crosses_boundary)) masks <- masks[!crosses_boundary]
  mass <- map_dbl(masks, ~sum(abs(center_of_mass(.x) - 0.5) ^ 2))
  # Very rarely, the 5 crop strategy fails. In that case, ask for a manual crop.
  if (all(crosses_boundary) | all(is.na(mass))) {
    message('Having trouble finding a cropped picture. Please drag a box over the fish:')
    l <- manual_crop(im)
    best_crop <- l$crop
    im2 <- resize(l$im, input_width, input_height)
    mask <- get_model_prediction(im2)
  } else {
    best_crop <- (1:5)[!crosses_boundary][which.min(mass)]
    mask <- masks[[which.min(mass)]]
  }

  list(
    best_crop = best_crop,
    mask = mask
  )
}

# This function takes data of a continuous curve, and uses it to produce a correct mask where the
# guppy is. Note that it is hard coded for image dimensions, currently. It allows for rescaling of
# the coordinates. This is useful, since we typically train on smaller images (e.g. 1/5 of the
# normal size).
produce_mask <- function(curve, scale, w, h) {
  curve <- do.call(rbind, curve[c('a', 'b', 'd', 'e', 'g')])
  coord_small <- round((curve - 0.5) * scale + 0.5)

  im <- imfill(dim = c(w, h, 1, 1))[, , 1, 1]
  im[coord_small] <- 1
  im <- as.cimg(im)
  im <- !bucketfill(im, 1, 1, color = 1)
}

mask_and_box_image <- function(Im, Mask, out_width = 800, out_height = 300) {
  # Im <- crop5(im, mask$best_crop)
  # Mask <- mask$mask

  # enlarge and smooth the mask
  Mask <- as.cimg(Mask) %>%
    resize(width(Im), height(Im), interpolation_type = 3) %>%
    blur_anisotropic(amplitude = 50, sharpness = 0, anisotropy = .01, gauss_prec = 1, sigma = 100)
  # ensure binary mask
  Mask[] <- as.numeric(Mask > 0.5)
  Mask <- boxblur(Mask, 3)

  #pad extra space for rotation
  Pad <- height(Im)
  Im <- pad(Im, Pad, 0, axes = 'xy', 'black')
  Mask <- pad(Mask, Pad, 0, axes = 'xy', 0)

  if (sum(Mask) == 0) return(Mask)

  bbox <- find_smallest_enclosing_rect(which(Mask == 1, arr.ind = TRUE)[, 1:2])
  bbox_angle <- angle(bbox$x[2], bbox$y[2], bbox$x[1], bbox$y[1])
  if (abs(bbox_angle) == 1/2*pi) bbox_angle <- 0
  # center point of fish to rotate around
  cp <- c(
    x = {rowMeans(Mask) > 0.001} %>% which() %>% mean(),
    y = {colMeans(Mask) > 0.001} %>% which() %>% mean()
  )
  new <- abind::abind(Im, Mask, along = 4) %>%
    as.cimg() %>%
    rotate_xy(
      angle = -bbox_angle / (2 * pi) * 360, cx = cp['x'], cy = cp['y'],
      interpolation = 2L, boundary_conditions = 0L
    )

  a <- channel(new, 4) %>% {imappend(list(., ., ., .), "c")}
  a[a < .1] <- 0
  new2 <- new * a
  # plot(new2); points(landmarks2[c('x', 'y')], col = 'blue', pch = 16)

  cr <- autocrop_coords(new2)
  w <- diff(cr)[1]; h <- diff(cr)[3]
  new3 <- imsub(new2, x %inr% cr[1:2], y %inr% cr[3:4])
  if (width(new3) < height(new3)) { # ensure image is horizontal
    new3 <- imrotate(new3, 90)
  }
  new3 <- resize(new3, -out_width / w * 100, -out_width / w * 100)
  h <- height(new3)
  new3 <- pad(new3, out_height - h, 'y', pos = 0, val = rep(0, 4))
  #clip values
  new3[new3 > 1] <- 1
  new3[new3 < 0] <- 0

  return(new3)
}

# If the 5-way crop does not correctly crop out the fish, we can do a manual crop instead. This
# function takes an image, asks the user for a box, then finds the appropriate 3/2 ratio box.
manual_crop <- function(im) {
  box <- imager::grabRect(im)
  # check if we need to expand to 3/2 ratio in the y or x direction
  if (((box['x1'] - box['x0']) / (box['y1'] - box['y0'])) > 3/2) {
    # Add this many pixels to the box in the y direction
    to_add <- (box['x1'] - box['x0']) * 2/3 - (box['y1'] - box['y0'])
    # Check and correct for boundaries
    if (box['y0'] - to_add > 0 && box['y1'] + to_add < height(im))
      box[c('y0', 'y1')] <- box[c('y0', 'y1')] + c(-to_add, to_add)
    if (box['y0'] - to_add > 0 && box['y1'] + to_add >= height(im))
      box[c('y0', 'y1')] <- box[c('y0', 'y1')] + c(-to_add, 0) + (height(im) - box['y1'])
    if (box['y0'] - to_add <= 0 && box['y1'] + to_add < height(im))
      box[c('y0', 'y1')] <- box[c('y0', 'y1')] + c(0, to_add) - (box['y0'])
  } else {
    # Add this many pixels to the box in the x direction
    to_add <- (box['y1'] - box['y0']) * 3/2 - (box['x1'] - box['x0'])
    # Check and correct for boundaries
    if (box['x0'] - to_add > 0 && box['x1'] + to_add < width(im))
      box[c('x0', 'x1')] <- box[c('x0', 'x1')] + c(-to_add, to_add)
    if (box['x0'] - to_add > 0 && box['x1'] + to_add >= width(im))
      box[c('x0', 'x1')] <- box[c('x0', 'x1')] + c(-to_add, 0) + (width(im) - box['x1'])
    if (box['y0'] - to_add <= 0 && box['y1'] + to_add < height(im))
      box[c('x0', 'x1')] <- box[c('x0', 'x1')] + c(0, to_add) - (box['x0'])
  }
  im2 <- imsub(im, x %inr% box[c('x0', 'x1')], y %inr% box[c('y0', 'y1')])
  return(list(im = im2, crop = box))
}
