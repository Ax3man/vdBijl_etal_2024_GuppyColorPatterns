# This function both saves the landmarks to disk, but also returns them.
# Optionally, photo and landmarks can be added to the training set.
place_manual_landmarks <- function(photo, add_to_training_set) {
  require(StereoMorph)

  landmark_id <- paste0(photo$unique_id, '.txt')
  # collect all the landmarks and save them in files
  curves <- data.frame(
    id = letters[1:7],
    start = c(2, 3, 4, 5, 6, 3, 4),
    end =   c(3, 4, 5, 6, 2, 3, 5)
  ) %>%
    as.matrix()

  digitizeImages(
    photo$file,
    paste0('data/manual_landmarks/', landmark_id),
    landmarks.ref = 1:6,
    curves.ref = curves
  )

  # We don't use readShapes here, because that will only allow the reading of single files
  # if they have the .txt extension. So we use readXML4R directly
  landmarks <- readXML4R(paste0('data/manual_landmarks/', landmark_id))[[1]][2:4]
  class(landmarks) <- 'shapes'

  # If requested, add the image and landmarks to the training set for the neural net
  if (add_to_training_set) {
    file.copy(
      from = photo$file,
      to = file.path(
        '~/Documents/KerasImagePool/Guppy_bg_removal/raw_images',
        paste0(photo$unique_id, '.jpg')
      )
    )

    file.copy(
      from = paste0('data/manual_landmarks/', landmark_id),
      to = file.path(
        '~/Documents/KerasImagePool/Guppy_bg_removal/manual_landmarks/',
        paste0(photo$unique_id, '.txt')
      )
    )
  }
  return(landmarks)
}

# This function both saves the landmarks to disk, but also returns them. This version will only
# ask for the 4 fixed landmarks.
# Optionally, photo and landmarks can be added to the training set.
place_manual_fixed_landmarks <- function(photo, add_to_training_set, fixed_landmarks = NULL) {
  require(StereoMorph)

  landmark_id <- paste0(photo$unique_id, '.txt')
  lm_file_path <- paste0('data/manual_fixed_landmarks/', landmark_id)

  # Optionally, to "pre-load" predicted landmarks into the shiny app, we will save the predicted
  # landmarks to the shape file before starting the app
  if (!is.null(fixed_landmarks)) {
    write_shape_file(fixed_landmarks, lm_file_path)
  }

  # create the file path for the extracted fish image
  im_file_path <- file.path(
    'data', 'extracted_fish', photo$replicate, photo$generation,
    paste0(photo$unique_id, '.png')
  )
  # collect all the landmarks and save them in files
  digitizeImages(
    im_file_path,
    lm_file_path,
    landmarks.ref = 1:4
  )

  # We don't use readShapes here, because that will only allow the reading of single files
  # if they have the .txt extension. So we use readXML4R directly
  landmarks <- readXML4R(lm_file_path)[[1]][[2]]
  landmarks <- as.data.frame(cbind(1:4, landmarks))
  names(landmarks) <- c('landmark', 'x', 'y')

  # If requested, add the image and landmarks to the training set for the neural net
  if (add_to_training_set) {
    file.copy(
      from = im_file_path,
      to = file.path(
        '~/Documents/KerasImagePool/Guppy_landmark_detection/extracted_fish/',
        paste0(photo$unique_id, '.png')
      )
    )
    file.copy(
      from = lm_file_path,
      to = file.path(
        '~/Documents/KerasImagePool/Guppy_landmark_detection/manual_fixed_landmarks/',
        landmark_id
      )
    )
  }
  return(landmarks)
}
