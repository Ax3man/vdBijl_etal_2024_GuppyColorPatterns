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

source('phenotyping_pipeline/tools.R')

melanin_extractor_v2 <- unserialize_model(
  read_rds('../color_analyses/keras_melanin_detection_v2/unet.rds'),
  custom_objects = list(dice_loss = dice_loss, dice_coef = dice_coef)
)

extract_melanin_v2 <- function(extracted_fish) {
  # resize the image to the dimension needed for the neural net prediction
  x <- resize(extracted_fish, 800, 288, interpolation_type = 5) %>% channel(1:3)
  X <- aperm(as.array(x), c(3, 1, 2, 4))

  p <- predict(melanin_extractor_v2, X)

  # turn the prediction into an image, and enlarge it to the size of the original
  mask <- p[1, , , ] %>%
    as.cimg() %>%
    resize(width(extracted_fish), height(extracted_fish), interpolation_type = 5L)
  # ensure binary mask
  mask[] <- as.numeric(mask > 0.5)

  out <- extracted_fish
  R(out) <- R(out) * mask
  G(out) <- G(out) * mask
  B(out) <- B(out) * mask
  out[, , , 4] <- mask

  # plot(extracted_fish, rescale = FALSE)
  # highlight(channel(out, 4) > 0)

  return(out)
}
