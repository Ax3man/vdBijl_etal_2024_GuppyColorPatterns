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

landmark_placer <- unserialize_model(
  read_rds('../color_analyses/keras_landmark_detection/covnet.rds')
)

auto_place_landmarks <- function(ef) {
  # catch if the extracted fish image is empty, since that the predictor will fail
  if (sum(ef) == 0) {
    # simply return landmarks in the center of the image instead.
    return(data.frame(
      landmark = 1:4,
      x = width(ef) / 2,
      y = height(ef) / 2
    ))
  }
  # ef is an extracted fish image

  ## BROKEN: automatic size detection has broken on the new tf version, for some reason.
  # detect required input sizes for the landmark placer
  #input_width <- landmark_placer$get_input_shape_at(0L)[[2]]
  #input_height <- landmark_placer$get_input_shape_at(0L)[[3]]
  # Manual override for now:
  input_width <- 400
  input_height <- 144

  # resize the image to fit the
  ef_small <- resize(ef, input_width, input_height, interpolation_type = 5L)

  pr <- predict(landmark_placer, aperm(ef_small, c(3, 1, 2, 4)))
  max_coords <- apply(pr, 4, function(x) arrayInd(which.max(x[1, , ]), dim(x)[2:3]))

  df <- data.frame(
    landmark = 1:4,
    x = (max_coords[1, ] - .5) / width(ef_small) * width(ef) + .5,
    y = (max_coords[2, ] - .5) / height(ef_small) * height(ef) + .5
  )
  # plot(ef_small); points(t(max_coords), col = 'red', pch = 16)
  # plot(ef); points(df$x, df$y, col = 'red', pch = 16)
  return(df)
}
