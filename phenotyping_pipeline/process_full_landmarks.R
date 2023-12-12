library(tidyverse)
suppressMessages(library(imager))
library(TSP)
library(furrr)

ims <- list.files('data/extracted_fish', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())
lms <- paste0('data/fixed_landmarks/', names(ims), '.txt')
flms <- lms %>% str_replace('fixed_landmarks', 'full_landmarks') %>% str_replace('.txt', '.rds')

# remove those with missing fixed landmarks
ims <- ims[file.exists(lms)]
flms <- flms[file.exists(lms)]
lms <- lms[file.exists(lms)]

# remove those that already have full landmarks
ims <- ims[!file.exists(flms)]
lms <- lms[!file.exists(flms)]
flms <- flms[!file.exists(flms)]

if (FALSE) {
  # x <- '20210921_IMG_0089'
  # file.remove(paste0('data/fixed_landmarks/', x, '.txt'))
  # file.remove(paste0('data/manual_fixed_landmarks/', x, '.txt'))
}

# Threshold determines what fraction of the maximum image dimension the landmarks can be away from
# the boundary. When larger than the threshold, probably something went wrong with the landmarks.
# By default, you can only be wrong by 2%.
extract_full_landmarks <- function(im, lm, file_name, threshold = 0.02) {
  # xx <- which(lms == 'data/fixed_landmarks/20210915_IMG_9454.txt')
  # xx <- which(names(ims) == '20211101_IMG_2785')
  # im <- ims[xx]; lm <- lms[xx]; file_name <- flms[xx]
  Img <- load.image(im)
  # landmarks.pixel contains the 4 fixed landmarks, curves.pixel will contain all outline pixels
  Lm <- StereoMorph::readShapes(lm)

  # If the fish is facing left, we flip everything:
  if (Lm$landmarks.pixel[1, 1] < Lm$landmarks[3, 1]) {
    Img <- mirror(Img, 'x')
    Lm$landmarks.pixel[, 1] <- width(Img) - Lm$landmarks.pixel[, 1]
  }

  # Get all the pixels around the boundary
  bound <- as.matrix(where(boundary(channel(Img, 4))))
  # Order all the pixels, so they form a large "circle" (continuous outline)
  bound <- bound[as.numeric(solve_TSP(as.TSP(dist(bound)))), ]

  if (FALSE) {  # Only run interactively for debugging
    plot(Img)
    points(bound, cex = 0.1, col = 'red')
    points(Lm$landmarks.pixel, cex = 1, col = 'green', pch = 16)
    text(Lm$landmarks.pixel, label = 1:4, cex = 1)
  }

  # Find the three boundary pixels that are closest to the fixed landmarks.
  find_closest <- function(v, m = bound) {
    which.min(apply(bound, 1, function(x) sum((v - x) ^ 2)))
  }
  snapped_lm <- Lm$landmarks.pixel
  bound_ids <- apply(snapped_lm[2:4, ], 1, find_closest)
  snapped_lm[2:4, ] <- bound[bound_ids, ]

  # Check if any of the landmarks are way off, if so, throw an error
  max_bound_distance <- max(sqrt(rowSums((snapped_lm - Lm$landmarks.pixel) ^ 2)))
  if (max_bound_distance > (max(dim(Img)) * threshold)) {
    stop(paste('Some of the landmarks are off the boundary by more than the threshold allows.\n',
               'It is this landmark file:', lm, '\n',
               'The worst landmark was',  round(max_bound_distance), 'pixels away from the boundary.'))
  }

  # We can now rearrange the bounds one more time, starting at landmark 2, then 3, then 4
  # First, rearrange so landmark 2 is the first and last pixel:
  bound <- bound[c(bound_ids['2']:nrow(bound), 1:bound_ids['2']), ]
  # All ids changed, so find bound_ids again:
  bound_ids <- apply(snapped_lm[2:4, ], 1, find_closest)
  # If landmark 4 comes earlier than landmark 3, it is going in the wrong direction, so flip it
  if (bound_ids['3'] > bound_ids['4']) {
    bound <- bound[nrow(bound):1, ]
    # All ids changed, so find bound_ids again:
    bound_ids <- apply(snapped_lm[2:4, ], 1, find_closest)
  }

  curves <- list(
    # Curve A, which is from the tip of the snout (landmark 2), along the back, to the dorsal attachment
    # of the caudal fin (landmark 3).
    A = bound[bound_ids['2']:bound_ids['3'], ],
    # Curve B, which is from dorsal attachment of the caudal fin (landmark 3) to the ventral side
    # attachment (landmark 4).
    B = bound[bound_ids['3']:bound_ids['4'], ],
    # Curve C, which is from ventral attachment of the caudal fin (landmark 4), along the belly, back
    # to the tip of the snout (landmark 2).
    C = bound[bound_ids['4']:nrow(bound), ]
  )

  Lm$curves.pixel <- curves
  Lm$landmarks.pixel <- snapped_lm
  write_rds(Lm, file_name)
  return(invisible(NULL))
}

if (length(ims) > 0) {
  cat('Processing', length(ims), 'landmark and images into full landmark sets...\n')

  plan(multisession, workers = 8)
  my_shapes <- future_pwalk(
    list(ims, lms, flms),
    extract_full_landmarks,
    .progress = TRUE,
    .options = furrr_options(seed = NULL)
  )
  plan(sequential)
  cat('\nFinished processing landmarks.\n')
} else {
  cat('No landmarks needed processing.\n')
}
