# Take an image, then use a lookup table to recalibrate the colors. This is done to get the same
# color characteristics between two photography setups. Uses Lab color space.
recalibrate <- function(im) {
  if (dim(im)[4] == 1) {
    warnings('Image does not have color, returning original image.')
    return(im)
  }

  LUT <- readr::read_rds('calibration/lookuptable.rds')
  LUT <- LUT[2:(nrow(LUT) - 1), ]

  lookup <- function(x, axis) {
    breaks <- c(-Inf, LUT[[paste0(axis, 2)]], Inf)
    ref <- LUT[[paste0(axis, 1)]]
    ref <- c(ref[1], ((ref + lag(ref)) / 2)[-1], tail(ref, 1))
    ref[cut(x, breaks, labels = FALSE)]
  }
  im_df <- as.data.frame(im, wide = 'c')
  Lab <- convertColor(im_df[3:5], 'sRGB', 'Lab')
  Lab[, 'L'] <- lookup(Lab[, 'L'], 'L')
  Lab[, 'a'] <- lookup(Lab[, 'a'], 'a')
  Lab[, 'b'] <- lookup(Lab[, 'b'], 'b')
  RGB <- convertColor(Lab, 'Lab', 'sRGB')
  im_df[3:5] <- as.data.frame(RGB)
  if (!('c.4' %in% names(im_df))) {
    im_df_long <- tidyr::pivot_longer(im_df, c.1:c.3, 'cc')
  } else {
    im_df[im_df$c.4 == 0, 3:5] <- 0
    im_df_long <- tidyr::pivot_longer(im_df, c.1:c.4, 'cc')
  }
  im_df_long$cc <- as.numeric(as.factor(im_df_long$cc))
  as.cimg(im_df_long, dims = dim(im))
}

detect_mask_on_boundary <- function(m) {
  any(m[c(1, 2, nrow(m), nrow(m) - 1), ] > 0.5) |
    any(m[, c(1, 2, ncol(m), ncol(m) - 1)] > 0.5) |
    all(m == 0)
}

center_of_mass <- function(m) {
  w <- which(m > 0.5, arr.ind = TRUE)
  colMeans(w) / dim(m)
}

count_non_transparent_pixels <- function(im, threshold = 0.5) {
  if (is.character(im)) {
    im <- load.image(im)
  }
  im %>% channel(4) %>% {sum(. > threshold)}
}

# The 5 crop takes a picture, and returns one of 5 predefined crops, determined by i. It will give
# (in order): The top left corner, bottom left, the center, top right, bottom right. Provided
# factor is < 2, you will always have the full picture covered. This means you can easily provide
# crops that are similarly scaled down as those from crop_random.
crop5 <- function(X, i, w = width(X), h = height(X), factor = 1.5) {
  xs <- round(w / factor)
  ys <- round(h / factor)
  if (inherits(X, 'imager_array')) {
    out <- switch(
      i,
      imsub(X, x %inr% c(0, xs), y %inr% c(0, ys)),
      imsub(X, x %inr% c(0, xs), y %inr% c(h - ys, h)),
      imsub(X, x %inr% c((w - xs) / 2, w - (w - xs) / 2),
            y %inr% c((h - ys) / 2, h - (h - ys) / 2)),
      imsub(X, x %inr% c(w - xs, w), y %inr% c(0, ys)),
      imsub(X, x %inr% c(w - xs, w), y %inr% c(h - ys, h))
    )
  } else {
    out <- switch(
      i,
      mutate(X),
      mutate(X, y = y - (h - ys)),
      mutate(X, x = x - (w - xs) / 2, y = y - (h - ys) / 2),
      mutate(X, x = x - (w - xs)),
      mutate(X, x = x - (w - xs), y = y - (h - ys))
    )
  }
  return(out)
}

angle <- function(x1, y1, x2, y2) atan2(y2 - y1, x2 - x1)
rotate <- function(x, y, angle) {
  list(x = x * cos(angle) - y * sin(angle),
       y = x * sin(angle) + y * cos(angle))
}
find_smallest_enclosing_rect <- function(points) {
  if (!requireNamespace('geometry', quietly = TRUE))
    stop("Package \'geometry\' needs to be installed to use this function.")
  a2 <- geometry::convhulln(points, options = 'FA')
  e <- points[a2$hull[, 2], ] - points[a2$hull[, 1], ]
  norms <- apply(e, 1, function(x) sqrt(x %*% x))
  v <- diag(1 / norms) %*% as.matrix(e)
  w <- cbind(-v[, 2], v[, 1])

  vertices <- as.matrix((points) [a2$hull, 1:2])
  minmax <- function(x) c(min(x), max(x))
  x <- apply(vertices %*% t(v), 2, minmax)
  y <- apply(vertices %*% t(w), 2, minmax)
  areas <- (y[1, ] - y[2, ]) * (x[1, ] - x[2, ])
  k <- which.min(areas)

  rect <- cbind(x[c(1, 2, 2, 1), k], y[c(1, 1, 2, 2), k]) %*%
    rbind(v[k, ], w[k, ])
  rect <- as.data.frame(rect)
  names(rect) <- c('x', 'y')

  return(rect)
}
autocrop_coords <- function(x) {
  x <- x[, , 1, 1]
  c(
    xstart = nrow(x) - which.max(cumsum(rev(rowSums(x)) != 0)),
    xend = which.max(cumsum(rowSums(x) != 0)),
    ystart = ncol(x) - which.max(cumsum(rev(colSums(x)) != 0)),
    yend = which.max(cumsum(colSums(x) != 0))
  )
}

dice_coef <- function(y_true, y_pred, smooth = 1.0) {
  y_true_f <- k_flatten(y_true)
  y_pred_f <- k_flatten(y_pred)
  intersection <- k_sum(y_true_f * y_pred_f)
  (2 * intersection + smooth) / (k_sum(y_true_f) + k_sum(y_pred_f) + smooth)
}

bce_dice_loss <- function(y_true, y_pred) {
  loss_binary_crossentropy(y_true, y_pred) + (1 - dice_coef(y_true, y_pred))
}

dice_loss <- function(y_true, y_pred) {
  1 - dice_coef(y_true, y_pred)
}

write_shape_file <- function(fixed_landmarks, path) {
  # only works with a data.frame of fixed landmarks
  stopifnot(
    is.data.frame(fixed_landmarks),
    identical(names(fixed_landmarks), c('landmark', 'x', 'y')),
    nrow(fixed_landmarks) > 0
  )

  formatted_landmarks <- paste(
    apply(fixed_landmarks, 1, function(x) paste(c('\t', x), collapse = '\t')),
    collapse = '\n'
  )
  writeLines(
    c(
      '<shapes type=list>',
      '\t<scaling type=logical names=FALSE length=1 as.numeric=FALSE >NA</scaling>',
      '\t<landmarks.pixel type=matrix rownames=TRUE colnames=FALSE nrow=4 ncol=2 as.numeric=TRUE >',
      formatted_landmarks,
      '\t</landmarks.pixel>',
      '</shapes>',
      ''
    ),
    path
  )
  return(invisible(NULL))
}
