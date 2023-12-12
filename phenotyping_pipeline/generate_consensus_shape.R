source('phenotyping_pipeline/process_full_landmarks.R')
library(tidyverse)
library(geomorph)

combine_shapes <- function(..., shape_list = NULL) {
  if (is.null(shape_list)) shape_list <- list(...)

  fixed <- map(shape_list, 'landmarks.pixel')
  landmarks.pixel <- do.call(abind::abind, list(fixed, along = 3))
  curves.pixel <- map(shape_list, 'curves.pixel')
  out <- list(landmarks.pixel = landmarks.pixel, curves.pixel = curves.pixel)
  class(out) <- 'shapes'
  return(out)
}

library(geomorph)
library(furrr)
plan(multisession, workers = 12)

lms <- list.files('data/full_landmarks', full.names = TRUE)
all_shapes <- future_map(lms, read_rds, .progress = TRUE)
plan(sequential)

my_shape <- combine_shapes(shape_list = all_shapes)
land <- readland.shapes(my_shape, c(50, 50, 50))

gpa <- gpagen(land, Parallel = 16)
plot(gpa, plot.param = list(pt.cex = .1, mean.cex = 1, mean.bg = 'firebrick'))
ref <- gpa$consensus

# If one of them looks really odd in the plot, use this to find which one it is:
# Morpho::find.outliers(gpa$coords, text = FALSE)

write_rds(ref, 'data/consensus_shape.rds')
