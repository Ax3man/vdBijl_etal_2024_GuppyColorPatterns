intake_photos <- function() {
  library(tidyverse)

  # loads a data.frame called `photos`, with meta-data and file names for all photos
  source('phenotyping_pipeline/index_photos.R')

  # load the saved database of already analyzed photos
  photo_database <- data.table::fread('photo_database.csv', header = TRUE) %>%
    mutate(date = as.character(date)) %>% as_tibble()
  # find the photos that have not yet been analyzed

  missing_photos <- which(!(photos$unique_id %in% photo_database$unique_id))

  # randomize order to reduce order effects
  if (length(missing_photos) > 1) {
    missing_photos <- sample(missing_photos, size = length(missing_photos))
  }

  if (length(missing_photos) == 0) {
    stop('No new photos found.')
  }

  cat(length(missing_photos), 'new photos found, out of', nrow(photos), 'total.\n')

  source('phenotyping_pipeline/tools.R')
  source('phenotyping_pipeline/fish_extraction.R')
  source('phenotyping_pipeline/auto_place_landmarks.R')
  source('phenotyping_pipeline/carotenoid_extraction.R')
  source('phenotyping_pipeline/melanin_extraction_v2.R')
  source('phenotyping_pipeline/place_manual_landmarks.R')
  library(future)
  plan(multisession, workers = 4)

  x11(width = 13, height = 13, title = 'GUPPY EXTRACTOR')
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))

  for (p in seq_along(missing_photos)) {
    ### Preparation and loading --------------------------------------------------------------------
    cat('Preparing next photo...\n')
    # current photo
    i <- missing_photos[p]
    # next photo
    j <- missing_photos[p + 1]
    # two photos from now
    k <- missing_photos[p + 2]
    cat(length(missing_photos) - p + 1, 'photos left...\n')
    # If first photo, load the current photo and do the extraction
    if (p == 1) {
      image_path <- photos$file[i]
      im <- image_path %>% load.image() %>% imresize(1 / 3)
      mask <- extract_fish1(im)
      extracted_fish <- extract_fish2(im, mask) %>% recalibrate()
      # Also start loading the next image, for iter 2.
      next_image_path2 <- photos$file[j]
      im_promise2 <- future(
        next_image_path2 %>% load.image() %>% imresize(1 / 3),
        seed = NULL
      )
    }
    # If second photo, do the extraction
    if (p == 2) {
      im <- value(im_promise2)
      mask <- extract_fish1(im)
      extracted_fish <- extract_fish2(im, mask) %>% recalibrate()
      # Prep for round 3
      next_im <- value(im_promise)
      mask <- extract_fish1(next_im)
    }
    # If not the first or second photo
    if (p > 2) {
      im <- next_im
      next_im <- value(im_promise)
      mask <- extract_fish1(next_im)

      extracted_fish <- value(next_extracted_fish)
    }
    # And if not the last or 2nd to last photo, already start loading the next photo using futures
    if (p < (length(missing_photos) - 1)) {
      next_image_path <- photos$file[k]
      im_promise <- future(
        next_image_path %>% load.image() %>% imresize(1 / 3),
        seed = NULL
      )
    }
    if ((p > 1 && p < length(missing_photos))) {
      next_extracted_fish <- future(extract_fish2(next_im, mask) %>% recalibrate(), seed = NULL)
    }
    # Prepare an entry for the photo database, we'll save it at the end if everything worked
    database_entry <- photos[i, ]

    ### Fish extraction & landmark placement -------------------------------------------------------
    layout(t(t(1:2)), heights = 2:1)
    actual_scale <- dim(im)[1] / 6000

    cat('Trying to find fixed landmarks...\n')
    fixed_landmarks <- auto_place_landmarks(extracted_fish)

    # Rotate the fish around if needed.
    if (fixed_landmarks$y[3] > fixed_landmarks$y[4]) {
      extracted_fish <- imrotate(extracted_fish, 180, interpolation = 0)
      fixed_landmarks <- auto_place_landmarks(extracted_fish)
    }

    plot(im, axes = FALSE, rescale = FALSE)
    plot(extracted_fish, axes = FALSE, rescale = FALSE)
    points(fixed_landmarks$x, fixed_landmarks$y, col = 'firebrick', pch = 16, cex = 1.5)
    text(fixed_landmarks$x, fixed_landmarks$y, labels = fixed_landmarks$landmark, cex = 0.7)

    # Ask user which photo method was used
    # while (TRUE) {
    #   photo_method <- tolower(readline('What photo method was used? (1/2) '))
    #   if (photo_method %in% c('1', '2')) break
    # }
    photo_method <- 2
    database_entry$photo_method <- as.integer(photo_method)

    # Ask user if the extraction looks okay. If not, prompt the user to do it manually
    while (TRUE) {
      correct_extraction <- tolower(readline('Does the extraction look correct? (y/n) '))
      if (correct_extraction %in% c('y', 'n')) break
    }
    if (correct_extraction == 'n') {
      # while (TRUE) {
      #   add_to_training <- tolower(readline(
      #     'Would you like to add this photo to the training set? (y/n) '
      #   ))
      #   if (add_to_training %in% c('y', 'n')) break
      # }
      add_to_training <- 'y'
      suppressWarnings( # Package has spurious warnings, I tried fixing them in a PR.
        landmarks <- place_manual_landmarks(
          database_entry, add_to_training_set = add_to_training == 'y'
        )
      )
      extracted_fish <- extract_fish_manual(im, landmarks, actual_scale) %>% recalibrate()
      fixed_landmarks <- auto_place_landmarks(extracted_fish)

      # Rotate the fish around if needed.
      if (fixed_landmarks$y[3] > fixed_landmarks$y[4]) {
        extracted_fish <- imrotate(extracted_fish, 180, interpolation = 0)
        fixed_landmarks <- mutate(
          fixed_landmarks,
          x = -x + width(extracted_fish),
          y = -y + height(extracted_fish)
        )
        print(fixed_landmarks)
      }

      plot(im, axes = FALSE, rescale = FALSE)
      plot(extracted_fish, axes = FALSE, rescale = FALSE)
      points(fixed_landmarks$x, fixed_landmarks$y, col = 'firebrick', pch = 16, cex = 1.5)
      text(fixed_landmarks$x, fixed_landmarks$y, labels = fixed_landmarks$landmark, cex = 0.7)

      database_entry <- mutate(database_entry, manual_landmarks = TRUE)
    } else {
      database_entry <- mutate(database_entry, manual_landmarks = FALSE)
    }

    # Save the extracted fish image to disk
    extracted_fish_path <- file.path(
      'data', 'extracted_fish', photos$replicate[i], photos$generation[i],
      paste0(photos$unique_id[i], '.png')
    )
    imager::save.image(extracted_fish, extracted_fish_path)

    # Ask user if the landmarks look okay. If not, prompt the user to do it manually
    while (TRUE) {
      correct_landmarks <- tolower(readline('Do the landmarks look correct? (y/n) '))
      if (correct_landmarks %in% c('y', 'n')) break
    }
    if (correct_landmarks == 'n') {
      # while (TRUE) {
      #   add_to_training <- tolower(readline(
      #     'Would you like to add these landmarks to the training set? (y/n) '
      #   ))
      #   if (add_to_training %in% c('y', 'n')) break
      # }
      add_to_training <- 'y'
      suppressWarnings(
        fixed_landmarks <- place_manual_fixed_landmarks(
          database_entry,
          add_to_training_set = add_to_training == 'y'#,
          #fixed_landmarks = fixed_landmarks
        )
      )
      # Rotate the fish around if needed.
      if (fixed_landmarks$y[3] > fixed_landmarks$y[4]) {
        extracted_fish <- imrotate(extracted_fish, 180, interpolation = 0)
        fixed_landmarks <- mutate(
          fixed_landmarks,
          x = -x + width(extracted_fish),
          y = -y + height(extracted_fish)
        )
        imager::save.image(extracted_fish, extracted_fish_path)
      }
    }
    # Save the fixed landmarks to disk
    fixed_landmark_path <- file.path(
      'data', 'fixed_landmarks',
      paste0(photos$unique_id[i], '.txt')
    )
    write_shape_file(fixed_landmarks, fixed_landmark_path)

    # Note if the picture shows the left or the right flank of the fish
    database_entry <- mutate(
      database_entry,
      facing_direction = ifelse(fixed_landmarks$x[2] > fixed_landmarks$x[3], 'r', 'l')
    )

    ### Carotenoid extraction ----------------------------------------------------------------------
    cat('Performing automatic carotenoid extraction...\n')
    layout(1)
    carotenoid <- extract_carotenoid(extracted_fish)
    plot(extracted_fish)
    highlight(channel(carotenoid, 4) == 1)

    # Ask user if the carotenoid extraction looks correct. If not, prompt to do it manually.
    while (TRUE) {
      correct_carotenoid <- tolower(readline('Does the carotenoid extraction look correct? (y/n) '))
      if (correct_carotenoid %in% c('y', 'n')) break
    }
    if (correct_carotenoid == 'n') {
      manual_car_file <- paste0(
        'data/carotenoid_coloration_manual/', database_entry$unique_id, '.png'
      )
      cat('Please perform manual carotenoid extraction on this file:\n', extracted_fish_path, '\n')
      # Wait for the user to do the manual extraction. Check every second if the manual extraction
      # file exists.
      while (!file.exists(manual_car_file)) {
        Sys.sleep(.1)
      }
      # Use the manual carotenoid extraction to build a similarly structure image as the automated
      # one.
      Sys.sleep(1) # give GIMP time to save the image properly
      manual_carotenoid <- load.image(manual_car_file)
      mask <- as.numeric(channel(manual_carotenoid, 4) > 0.5)
      carotenoid <- extracted_fish
      R(carotenoid) <- R(carotenoid) * mask
      G(carotenoid) <- G(carotenoid) * mask
      B(carotenoid) <- B(carotenoid) * mask
      carotenoid[, , , 4] <- mask
      # Show it to the user

      plot(extracted_fish)
      highlight(channel(carotenoid, 4) == 1)
      continue <- tolower(readline(
        'Showing manual carotenoid extraction. Press return to continue... '
      ))

      # Ask whether to add the manual extraction to the training data
      # while (TRUE) {
      #   add_to_training <- tolower(readline(
      #     'Would you like to add this photo to the training set? (y/n) '
      #   ))
      #   if (add_to_training %in% c('y', 'n')) break
      # }
      add_to_training <- 'y'
      if (add_to_training == 'y') {
        imager::save.image(
          extracted_fish,
          file.path(
            '~/Documents/KerasImagePool/Guppy_carotenoid_detection/raw_images',
            paste0(database_entry$unique_id, '.png')
          )
        )
        file.copy(
          from = manual_car_file,
          to = file.path(
            '~/Documents/KerasImagePool/Guppy_carotenoid_detection/manual_extractions',
            paste0(database_entry$unique_id, '.png')
          )
        )
      }
    }
    # calculate the amount of carotenoid coloration, as the % body area
    database_entry <- mutate(
      database_entry,
      car_perc = count_non_transparent_pixels(carotenoid) /
        count_non_transparent_pixels(extracted_fish) * 100
    )
    # Save the carotenoid extraction to disk
    imager::save.image(
      carotenoid,
      file.path('data', 'carotenoid_coloration', photos$replicate[i], photos$generation[i],
                paste0(photos$unique_id[i], '.png'))
    )

    ### Melanin extraction -------------------------------------------------------------------------
    cat('Performing automatic melanin extraction...\n')
    layout(1)
    melanin <- extract_melanin_v2(extracted_fish)
    plot(extracted_fish)
    highlight(channel(melanin, 4) == 1)

    # Ask user if the melanin extraction looks correct. If not, prompt to do it manually.
    while (TRUE) {
      correct_melanin <- tolower(readline('Does the melanin extraction look correct? (y/n) '))
      if (correct_melanin %in% c('y', 'n')) break
    }
    if (correct_melanin == 'n') {
      manual_mel_file <- paste0(
        'data/melanic_coloration_manual_v2/', database_entry$unique_id, '.png'
      )
      cat('Please perform manual melanin extraction on this file:\n', extracted_fish_path, '\n')
      while (!file.exists(manual_mel_file)) {
        Sys.sleep(.1)
      }
      # Use the manual carotenoid extraction to build a similarly structure image as the automated
      # one.
      Sys.sleep(1) # give GIMP time to save the image properly
      manual_melanin <- load.image(manual_mel_file)
      mask <- as.numeric(channel(manual_melanin, 4) > 0.5)
      melanin <- extracted_fish
      R(melanin) <- R(melanin) * mask
      G(melanin) <- G(melanin) * mask
      B(melanin) <- B(melanin) * mask
      melanin[, , , 4] <- mask
      # Show it to the user

      plot(extracted_fish)
      highlight(channel(melanin, 4) == 1)
      continue <- tolower(readline(
        'Showing manual melanin extraction. Press return to continue... '
      ))

      # Ask whether to add the manual extraction to the training data
      # while (TRUE) {
      #   add_to_training <- tolower(readline(
      #     'Would you like to add this photo to the training set? (y/n) '
      #   ))
      #   if (add_to_training %in% c('y', 'n')) break
      # }
      add_to_training <- 'y'
      if (add_to_training == 'y') {
        imager::save.image(
          extracted_fish,
          file.path(
            '~/Documents/KerasImagePool/Guppy_melanin_detection_v2/raw_images',
            paste0(database_entry$unique_id, '.png')
          )
        )
        file.copy(
          from = manual_mel_file,
          to = file.path(
            '~/Documents/KerasImagePool/Guppy_melanin_detection_v2/manual_extractions',
            paste0(database_entry$unique_id, '.png')
          )
        )
      }
    }
    # calculate the amount of melanin coloration, as the % body area
    database_entry <- mutate(
      database_entry,
      mel_perc_v2 = count_non_transparent_pixels(melanin) /
        count_non_transparent_pixels(extracted_fish) * 100
    )
    # Save the melanin extraction to disk
    imager::save.image(
      melanin,
      file.path('data', 'melanic_coloration_v2', photos$replicate[i], photos$generation[i],
                paste0(photos$unique_id[i], '.png'))
    )

    ### Finishing up -------------------------------------------------------------------------------
    # Add this photo to the photo database, and save it to disk
    photo_database <- bind_rows(photo_database, database_entry)
    data.table::fwrite(photo_database, 'photo_database.csv')
  }
  plan(sequential)
}
intake_photos()
