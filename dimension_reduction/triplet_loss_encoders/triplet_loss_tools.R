# Find the index of a layer by its name
find_layer <- function(model, layer_name) {
  layer_names <- map_chr(model$layers, 'name')
  which(layer_names == layer_name)
}

# Create a subsection of the model, that is only the shared layers that yields the n-dimensional
# embedding.
get_layers <- function(model, indices) {
  l <- list()
  for (i in seq_along(indices)) l[[i]] <- get_layer(model, index = indices[i])
  return(l)
}

# Identity loss: workaround Keras loss definition to use custom triplet loss.
# There is no true label: we just want to minimize the margin-based triplet loss
# to learn the embeddings
loss_identity <- function(y_true, y_pred) k_mean(y_pred - 0 * y_true)

# margin-based triplet loss. Returns a loss function with approprate margin
loss_margin_triplet <- function(margin) {
  function(x) {
    anchor <- x[[1]]
    positive <- x[[2]]
    negative <- x[[3]]
    # each is tensor with shape (batch, n_dim)

    DAP <- k_sum(k_square(anchor - positive), axis = 2, keepdims = T)
    DAN <- k_sum(k_square(anchor - negative), axis = 2, keepdims = T)
    loss <- k_maximum(x = DAP - DAN + margin, y = 0)
    # optionally ignore easy triplets:
    #loss[loss > 0]
    return(loss)
  }
}

chunk_triplets <- function(triplets, batch_size) {
  n <- nrow(triplets$anchors)
  n_batches <- ceiling(n / batch_size)

  batch_ids <- split(
    sample(n),
    rep(seq_len(n_batches), each = batch_size)[seq_len(n)]
  )

  if (length(dim(triplets[[1]])) == 3) {
    lapply(
      batch_ids,
      \(batch) lapply(triplets, \(arr) arr[batch, , , drop = FALSE])
    )
  } else {
    lapply(
      batch_ids,
      \(batch) lapply(triplets, \(arr) arr[batch, , , , drop = FALSE])
    )
  }
}

# define our own function to calculate the metrics I want. Not sure we can use the normal system for
# metrics, since I want to use the embedder
get_metrics <- function(triplets, embedder, epoch = NULL, batch_size, margin, debug = FALSE) {
  anchors <- positives <- negatives <- list()

  if (length(triplets) == 3) {
    triplets <- chunk_triplets(triplets, batch_size)
  }

  for (batch in seq_along(triplets)) {
    anchors[[batch]] <- predict(embedder, triplets[[batch]]$anchors, batch_size = batch_size)
    positives[[batch]] <- predict(embedder, triplets[[batch]]$positives, batch_size = batch_size)
    negatives[[batch]] <- predict(embedder, triplets[[batch]]$negatives, batch_size = batch_size)
  }
  anchors <- do.call(rbind, anchors)
  positives <- do.call(rbind, positives)
  negatives <- do.call(rbind, negatives)

  if (debug) {
    message('Debug enabled, making plots...')

    # # Show some triplets:
    ns <- 3
    opa <- par(mfrow = c(ns, 3), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    for (i in seq_len(dim(triplets[[1]]$anchors)[1])) {
      for (j in 1:3) {
        suppressWarnings(plot(as.cimg(triplets[[1]][[j]][i, , , ]), axes = FALSE, rescale = FALSE))
        #plot(as.raster(t(triplets[[1]][[j]][i, , , 1])))
      }
    }
    par(opa)

    bind_cols(anchors, positives, negatives) %>%
      slice(1:10) %>%
      ggplot() +
      geom_segment(aes(...1, ...2, xend = ...3, yend = ...4, color = 'positive')) +
      geom_segment(aes(...5, ...6, xend = ...3, yend = ...4, color = 'negative')) +
      geom_point(aes(...1, ...2, color = 'anchor')) +
      geom_point(aes(...3, ...4, color = 'positive')) +
      geom_point(aes(...5, ...6, color = 'negative')) +
      scale_color_manual(values = c('anchor' = 'black', positive = 'green', negative = 'firebrick')) +
      facet_wrap(vars(triplet = row_number(`...1`)))
  }

  if (length(dim(anchors)) > 2 | ncol(anchors) > 100) {
    stop('Malformed prediction. Check output layer of the embedder.')
  }
  if (embedder$layers[[length(embedder$layers)]]$name != c("norm"))  {
    stop('Last layer of embedder should be the batch normalization layer named "norm"')
  }
  if ((nrow(distinct(as.data.frame(anchors))) / nrow(anchors)) < 0.5) {
    stop('More than 50% duplicated embeddings detected. Are many images identical?')
  }

  DAP <- rowSums((anchors - positives) ^ 2)
  DAN <- rowSums((anchors - negatives) ^ 2)

  ret <- data.frame(
    loss = mean(pmax(DAP - DAN + margin, 0)),
    DAP = mean(DAP),
    DAN = mean(DAN),
    delta_dist = mean(DAP - DAN),
    acc = mean(DAP < DAN),
    acc_lower = binom.test(sum(DAP < DAN), length(DAP))$conf.int[1],
    acc_upper = binom.test(sum(DAP < DAN), length(DAP))$conf.int[2],
    acc2 = mean(DAP < (DAN - margin))
  )
  ret$epoch <- epoch
  return(ret)
}

# build and compile the model with margin-based triplet loss
build_model <- function(
    embed_dim, filter_base, lr, margin, input_dim = c(160, 40, 3), norm_momentum = 0.99, type = 1,
    dropout = 0
) {
  input_anchor <- layer_input(shape = input_dim, name = "anchors")
  input_positive <- layer_input(shape = input_dim, name = "positives")
  input_negative <- layer_input(shape = input_dim, name = "negatives")

  if (type == 1) {
    # shared layers
    #rcontr <- layer_random_contrast(factor = 0.5)
    conv1 <- layer_separable_conv_2d(
      filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu', name = 'conv1',
    )
    pool1 <- layer_max_pooling_2d(padding = 'same', name = 'pool1')
    conv2 <- layer_separable_conv_2d(
      filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu', name = 'conv2'
    )
    pool2 <- layer_max_pooling_2d(padding = 'same', name = 'pool2')
    conv3 <- layer_separable_conv_2d(
      filters = filter_base * 4, kernel_size = c(3, 3), padding = 'same', activation = 'relu', name = 'conv3'
    )
    pool3 <- layer_max_pooling_2d(padding = 'same', name = 'pool3')
    flat <- layer_flatten(name = 'flat')
    dlay <- layer_dense(units = embed_dim, name = 'embed')
    norm <- layer_batch_normalization(axis = 1L, momentum = norm_momentum, name = 'norm')

    # the triple encoders
    embed_anchor <- input_anchor %>% #rcontr() %>%
      conv1() %>% pool1() %>% conv2() %>% pool2() %>% conv3() %>% pool3() %>% flat() %>% dlay() %>% norm()
    embed_positive <- input_positive %>% #rcontr() %>%
      conv1() %>% pool1() %>% conv2() %>% pool2() %>% conv3() %>% pool3() %>% flat() %>% dlay() %>% norm()
    embed_negative <- input_negative %>% #rcontr() %>%
      conv1() %>% pool1() %>% conv2() %>% pool2() %>% conv3() %>% pool3() %>% flat() %>% dlay() %>% norm()
  }
  if (type == 2) {
    # shared layers
    #rcontr <- layer_random_contrast(factor = 0.5)
    conv1 <- layer_separable_conv_2d(
      filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu', name = 'conv1',
    )
    pool1 <- layer_max_pooling_2d(padding = 'same', name = 'pool1')
    conv2 <- layer_separable_conv_2d(
      filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu', name = 'conv2'
    )
    pool2 <- layer_max_pooling_2d(padding = 'same', name = 'pool2')

    flat <- layer_flatten(name = 'flat')
    dlay1 <- layer_dense(units = 32, name = 'dense')
    dlay2 <- layer_dense(units = embed_dim, name = 'embed')
    norm <- layer_batch_normalization(axis = 1L, momentum = norm_momentum, name = 'norm')

    # the triple encoders
    embed_anchor <- input_anchor %>% #rcontr() %>%
      conv1() %>% pool1() %>% conv2() %>% pool2() %>% flat() %>% dlay1() %>% dlay2() %>% norm()
    embed_positive <- input_positive %>% #rcontr() %>%
      conv1() %>% pool1() %>% conv2() %>% pool2() %>% flat() %>% dlay1() %>% dlay2() %>% norm()
    embed_negative <- input_negative %>% #rcontr() %>%
      conv1() %>% pool1() %>% conv2() %>% pool2() %>% flat() %>% dlay1() %>% dlay2() %>% norm()
  }
  if (type == 3) {
    # shared layers
    #rcontr <- layer_random_contrast(factor = 0.5)
    conv1 <- layer_separable_conv_2d(
      filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv1'
    )
    drop1 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop1')
    pool1 <- layer_max_pooling_2d(padding = 'same', name = 'pool1')
    conv2 <- layer_separable_conv_2d(
      filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv2'
    )
    drop2 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop2')
    pool2 <- layer_max_pooling_2d(padding = 'same', name = 'pool2')
    conv3 <- layer_separable_conv_2d(
      filters = filter_base * 4, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv3'
    )
    drop3 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop3')
    pool3 <- layer_max_pooling_2d(padding = 'same', name = 'pool3')
    conv4 <- layer_separable_conv_2d(
      filters = filter_base * 8, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv4'
    )
    drop4 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop4')
    pool4 <- layer_max_pooling_2d(padding = 'same', name = 'pool4')
    flat <- layer_flatten(name = 'flat')
    dlay1 <- layer_dense(units = pmax(8, embed_dim), name = 'dense')
    dlay2 <- layer_dense(units = embed_dim, name = 'embed')
    norm <- layer_batch_normalization(axis = 1L, momentum = norm_momentum, name = 'norm')

    # the triple encoders
    embed_anchor <- input_anchor %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      flat() %>% dlay1() %>% dlay2() %>% norm()
    embed_positive <- input_positive %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      flat() %>% dlay1() %>% dlay2() %>% norm()
    embed_negative <- input_negative %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      flat() %>% dlay1() %>% dlay2() %>% norm()
  }
  if (type == 4) {
    # shared layers
    #rcontr <- layer_random_contrast(factor = 0.5)
    conv1 <- layer_separable_conv_2d(
      filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv1'
    )
    drop1 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop1')
    pool1 <- layer_max_pooling_2d(padding = 'same', name = 'pool1')
    conv2 <- layer_separable_conv_2d(
      filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv2'
    )
    drop2 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop2')
    pool2 <- layer_max_pooling_2d(padding = 'same', name = 'pool2')
    conv3 <- layer_separable_conv_2d(
      filters = filter_base * 4, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv3'
    )
    drop3 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop3')
    pool3 <- layer_max_pooling_2d(padding = 'same', name = 'pool3')
    conv4 <- layer_separable_conv_2d(
      filters = filter_base * 8, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv4'
    )
    drop4 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop4')
    pool4 <- layer_max_pooling_2d(padding = 'same', name = 'pool4')
    flat <- layer_flatten(name = 'flat')
    dlay <- layer_dense(units = embed_dim, name = 'embed')

    norm <- layer_lambda(f = k_l2_normalize, name = 'norm')

    # the triple encoders
    embed_anchor <- input_anchor %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      flat() %>% dlay()
    embed_positive <- input_positive %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      flat() %>% dlay()
    embed_negative <- input_negative %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      flat() %>% dlay() %>% norm()
  }
  if (type == 5) {
    # This is type 3, but with an additional convolutional layer for larger images
    # shared layers
    #rcontr <- layer_random_contrast(factor = 0.5)
    conv1 <- layer_separable_conv_2d(
      filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv1'
    )
    drop1 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop1')
    pool1 <- layer_max_pooling_2d(padding = 'same', name = 'pool1')
    conv2 <- layer_separable_conv_2d(
      filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv2'
    )
    drop2 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop2')
    pool2 <- layer_max_pooling_2d(padding = 'same', name = 'pool2')
    conv3 <- layer_separable_conv_2d(
      filters = filter_base * 4, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv3'
    )
    drop3 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop3')
    pool3 <- layer_max_pooling_2d(padding = 'same', name = 'pool3')
    conv4 <- layer_separable_conv_2d(
      filters = filter_base * 8, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv4'
    )
    drop4 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop4')
    pool4 <- layer_max_pooling_2d(padding = 'same', name = 'pool4')
    conv5 <- layer_separable_conv_2d(
      filters = filter_base * 8, kernel_size = c(3, 3), padding = 'same', activation = 'relu',
      name = 'conv5'
    )
    drop5 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop5')
    pool5 <- layer_max_pooling_2d(padding = 'same', name = 'pool5')

    flat <- layer_flatten(name = 'flat')
    dlay1 <- layer_dense(units = pmax(8, embed_dim), name = 'dense')
    dlay2 <- layer_dense(units = embed_dim, name = 'embed')
    norm <- layer_batch_normalization(axis = 1L, momentum = norm_momentum, name = 'norm')

    # the triple encoders
    embed_anchor <- input_anchor %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      conv5() %>% drop5() %>% pool5() %>%
      flat() %>% dlay1() %>% dlay2() %>% norm()
    embed_positive <- input_positive %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      conv5() %>% drop5() %>% pool5() %>%
      flat() %>% dlay1() %>% dlay2() %>% norm()
    embed_negative <- input_negative %>% #rcontr() %>%
      conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
      conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
      conv5() %>% drop5() %>% pool5() %>%
      flat() %>% dlay1() %>% dlay2() %>% norm()
  }

  # margin-based triplet loss is the output
  this_loss <- loss_margin_triplet(margin)

  loss <- list(embed_anchor, embed_positive, embed_negative) %>%
    layer_lambda(this_loss, output_shape = c(1), name = 'triplet_loss')

  # define model inputs/outputs
  model <- keras_model(
    inputs = c(input_anchor, input_positive, input_negative),
    outputs = loss
  )

  # compile model
  model %>% compile(
    loss = loss_identity,
    optimizer = optimizer_adam(learning_rate = lr)
  )
  model
}

# build and compile the model with margin-based triplet loss
build_model2 <- function(
    embed_dim, filter_base, lr, margin, input_dim = c(160, 40, 3), norm_momentum = 0.99,
    dropout = 0
) {
  input_anchor <- layer_input(shape = input_dim, name = "anchors")
  input_positive <- layer_input(shape = input_dim, name = "positives")
  input_negative <- layer_input(shape = input_dim, name = "negatives")

  # This is type 3, but with an additional convolutional layer for larger images
  # shared layers
  #rcontr <- layer_random_contrast(factor = 0.5)
  conv1 <- layer_separable_conv_2d(
    filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu', depth_multiplier = 2,
    name = 'conv1'
  )
  drop1 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop1')
  pool1 <- layer_max_pooling_2d(padding = 'same', name = 'pool1')
  conv2 <- layer_separable_conv_2d(
    filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu',depth_multiplier = 2,
    name = 'conv2'
  )
  drop2 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop2')
  pool2 <- layer_max_pooling_2d(padding = 'same', name = 'pool2')
  conv3 <- layer_separable_conv_2d(
    filters = filter_base, kernel_size = c(3, 3), padding = 'same', activation = 'relu',depth_multiplier = 2,
    name = 'conv3'
  )
  drop3 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop3')
  pool3 <- layer_max_pooling_2d(padding = 'same', name = 'pool3')
  conv4 <- layer_separable_conv_2d(
    filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu',depth_multiplier = 2,
    name = 'conv4'
  )
  drop4 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop4')
  pool4 <- layer_max_pooling_2d(padding = 'same', name = 'pool4')
  # conv5 <- layer_separable_conv_2d(
  #   filters = filter_base * 2, kernel_size = c(3, 3), padding = 'same', activation = 'relu',depth_multiplier = 2,
  #   name = 'conv5'
  # )
  # drop5 <- layer_spatial_dropout_2d(rate = dropout, name = 'drop5')
  # pool5 <- layer_max_pooling_2d(padding = 'same', name = 'pool5')

  flat <- layer_flatten(name = 'flat')
  dlay1 <- layer_dense(units = pmax(16, embed_dim), name = 'dense')
  dlay2 <- layer_dense(units = embed_dim, name = 'embed')
  norm <- layer_batch_normalization(axis = 1L, momentum = norm_momentum, name = 'norm')

  # the triple encoders
  embed_anchor <- input_anchor %>% #rcontr() %>%
    conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
    conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
    #conv5() %>% drop5() %>% pool5() %>%
    flat() %>% dlay1() %>% dlay2() %>% norm()
  embed_positive <- input_positive %>% #rcontr() %>%
    conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
    conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
    #conv5() %>% drop5() %>% pool5() %>%
    flat() %>% dlay1() %>% dlay2() %>% norm()
  embed_negative <- input_negative %>% #rcontr() %>%
    conv1() %>% drop1() %>% pool1() %>% conv2() %>% drop2() %>% pool2() %>%
    conv3() %>% drop3() %>% pool3() %>% conv4() %>% drop4() %>% pool4() %>%
    #conv5() %>% drop5() %>% pool5() %>%
    flat() %>% dlay1() %>% dlay2() %>% norm()

  # margin-based triplet loss is the output
  this_loss <- loss_margin_triplet(margin)

  loss <- list(embed_anchor, embed_positive, embed_negative) %>%
    layer_lambda(this_loss, output_shape = c(1), name = 'triplet_loss')

  # define model inputs/outputs
  model <- keras_model(
    inputs = c(input_anchor, input_positive, input_negative),
    outputs = loss
  )

  # compile model
  model %>% compile(
    loss = loss_identity,
    optimizer = optimizer_adam(learning_rate = lr)
  )
  model
}

# when using method = 'pedigree', this usually returns only about 0.9 * n rows.
generate_triplet_subset <- function(
    data, n, image_array, validation = FALSE, batch_match = FALSE, method = 'self',
    pedigree, min_r = 0.25, remove_maternal_effect = FALSE
) {
  # data <- db; pedigree <- ped_df
  # n = 100; validation = TRUE; image_array = all_images; batch_match = FALSE

  if (validation) {
    data <- data[data$validation, ]
  } else {
    data <- data[!data$validation, ]
  }
  # use dtplyr for faster execution of dplyr verbs, joins in particular
  data <- dtplyr::lazy_dt(data)
  data <- select(data, replicate, generation, selection, fish_id, unique_id, sample_id)

  ### Matching images to other images of the same fish ---------------------------------------------
  if (method == 'self') {
    if (batch_match) {
      anchors <- slice_sample(data, n = n, replace = TRUE)

      pairs <- anchors %>%
        left_join(
          select(data, -replicate, -generation, -selection),
          by = c('fish_id'), suffix = c('_anchor', '_pos')
        ) %>%
        group_by(unique_id_anchor) %>%
        slice_sample(n = 1)

      triplets <- pairs %>%
        left_join(data, by = c('replicate', 'generation', 'selection')) %>%
        rename(unique_id_neg = unique_id, sample_id_neg = sample_id) %>%
        group_by(unique_id_anchor) %>%
        slice_sample(n = 1) %>%
        as_tibble()

    } else {

      data <- select(data, -replicate, -generation)

      anchors <- slice_sample(data, n = n)

      pairs <- anchors %>%
        left_join(
          data,
          by = c('fish_id'), suffix = c('_anchor', '_pos')
        ) %>%
        group_by(unique_id_anchor) %>%
        slice_sample(n = 1)

      triplets <- pairs %>%
        left_join(data, by = character()) %>% # full cross join
        rename(unique_id_neg = unique_id, sample_id_neg = sample_id) %>%
        group_by(unique_id_anchor) %>%
        slice_sample(n = 1) %>%
        as_tibble()
    }
  }
  if (method == 'pedigree') {
    if (remove_maternal_effect) stop('Currently, maternal brothers are never excluded')
    # remove females from pedigree
    pedigree <- filter(pedigree, str_starts(id1, 'M'), str_starts(id2, 'M'))

    anchors <- data %>%
      filter(generation != 'parental_gen_1') %>%
      slice_sample(n = n, replace = TRUE) %>%
      mutate(triplet_id = row_number())
    #n_distinct(as_tibble(anchors)$fish_id)

    pairs <- anchors %>%
      # we drop these, since if we batch match, we only match the positive and the negative.
      select(-replicate, -generation, -selection) %>%
      # join only pedigree values that are actually relevant (instead of filtering everything after)
      # TODO: consider whether removing full sibs is a good idea, include !shared_dam here:
      left_join(filter(pedigree, r >= min_r, r < 0.99999), by = c('fish_id' = 'id1')) %>%
      select(-r, -shared_dam) %>%
      inner_join(
        data %>% group_by(fish_id) %>% slice_sample(n = 1) %>% ungroup(),
        by = c('id2' = 'fish_id'), suffix = c('_anchor', '_pos'), relationship = 'many-to-many'
      ) %>%
      rename(fish_id_anchor = fish_id, fish_id_pos = id2) %>%
      group_by(triplet_id) %>% slice_sample(n = 1) %>% ungroup()
    #n_distinct(as_tibble(pairs)$fish_id_anchor)

    if (batch_match) {
      triplets <- pairs %>%
        left_join(
          rename(data, fish_id_neg = fish_id, unique_id_neg = unique_id, sample_id_neg = sample_id),
          by = c('replicate', 'generation', 'selection'), relationship = 'many-to-many'
        ) %>%
        inner_join(
          filter(pedigree, r == 0),
          by = c('fish_id_anchor' = 'id1', 'fish_id_neg' = 'id2')
        ) %>%
        group_by(triplet_id) %>% slice_sample(n = 1) %>% ungroup() %>%
        as_tibble()
      #n_distinct(triplets$fish_id_anchor)

    } else {
      triplets <- data %>%
        rename(fish_id_neg = fish_id, unique_id_neg = unique_id, sample_id_neg = sample_id) %>%
        select(-replicate, -generation, -selection) %>%
        right_join(pairs, by = character()) %>%
        inner_join(
          filter(pedigree, r == 0),
          by = c('fish_id_anchor' = 'id1', 'fish_id_neg' = 'id2')
        ) %>%
        group_by(triplet_id) %>% slice_sample(n = 1) %>% ungroup() %>%
        as_tibble()
    }
  }

  ### Collect images in an array to return ---------------------------------------------------------
  if (length(dim(image_array)) == 4) {
    ret <- list(
      anchors   = image_array[triplets$sample_id_anchor, , , ],
      positives = image_array[triplets$sample_id_pos   , , , ],
      negatives = image_array[triplets$sample_id_neg   , , , ]
    )
  } else {
    ret <- list(
      anchors   = image_array[triplets$sample_id_anchor, , ],
      positives = image_array[triplets$sample_id_pos   , , ],
      negatives = image_array[triplets$sample_id_neg   , , ]
    )
  }
  return(ret)
}

find_threeway_triplets <- function(
    test_mat, match_mat1, match_mat2, n_workers = parallelly::availableCores(logical = FALSE)
) {
  # source('quant_gen/prepare_pedigrees.R')
  # patrilines <- data.table::fread('photo_database.csv', data.table = FALSE) %>%
  #   select(fish_id) %>% distinct() %>% add_patriline(ped_df)
  # Y <- outer(patrilines$patriline, patrilines$patriline, \(X, Y) as.numeric(X == Y))
  # rownames(Y) <- colnames(Y) <- patrilines$fish_id
  # A <- A[rownames(Y), colnames(Y)]
  # X <- X[rownames(Y), colnames(Y)]
  # test_mat <- A; match_mat1 <- Y; match_mat2 <- X

  # i is the row number, selecting the anchor under consideration
  find_triplets_for_anchor <- function(i) {
    Test <- Rfast::Outer(test_mat[i,], test_mat[i,], '-') < 0
    match1 <- Rfast::Outer(match_mat1[i,], match_mat1[i,], '-') == 0
    match2 <- Rfast::Outer(match_mat2[i,], match_mat2[i,], '-') == 0

    passed <- which(Test & match1 & match2, arr.ind = TRUE)

    if (length(passed) == 0) return(data.frame())

    data.frame(
      anchor = rownames(test_mat)[i],
      positive = rownames(test_mat)[passed[, 'row']],
      negative = rownames(test_mat)[passed[, 'col']]
    )
  }

  future::plan(future::multisession, workers = n_workers)
  triplets <- furrr::future_map_dfr(
    1:nrow(test_mat),
    find_triplets_for_anchor,
    .progress = TRUE, .options = furrr::furrr_options(seed = NULL)
  )
  future::plan(future::sequential)

  return(triplets)
}


sample_threeway_triplets <- function(
    image_array, triplet_files, n, data, validation
) {
  # triplet_files <- list.files('dimension_reduction/triplet_loss_encoders/triplet_lists/', 'AXY', f = TRUE)
  # data <- db; validation = FALSE; n <- 100; validation_rep <- 1

  # We'll load more triplets, divided over three replicates. We need additional, since we're going
  # to drop validation (or training) data later. We need some more for training, way more for val
  if (validation) {
    n_large <- n * 1000
  } else {
    n_large <- n * 10
  }

  # if we give triplet_files as file names, load random triplets from those files
  if (is.character(triplet_files)) {
    warning('loading triplets from disk can be very slow, passing data.tables is recommended.')
    tr <- map_dfr(
      triplet_files,
      \(X) data.table::fread(cmd = paste("gshuf -n", floor(n_large / 3), X), header = FALSE)
    )
    names(tr) <- c('anchor', 'positive', 'negative')
  } else {
    # otherwise, it should be a list of data.frames from we can sample instead (faster)
    tr <- map_dfr(
      triplet_files,
      \(X) X %>% slice_sample(n = floor(n_large / 3))
    )
  }

  if (validation) {
    data <- data[data$validation, ]
  } else {
    data <- data[!data$validation, ]
  }
  data2 <- select(data, fish_id, sample_id) %>% mutate(fish_id = tolower(fish_id))

  # for each triplet, take a random combination of photos
  tr <- tr %>%
    mutate(triplet = row_number()) %>%
    inner_join(data2, by = c('anchor' = 'fish_id'), relationship = 'many-to-many') %>%
    inner_join(
      data2, by = c('positive' = 'fish_id'), suffix = c('', '_pos'), relationship = 'many-to-many'
    ) %>%
    inner_join(
      data2, by = c('negative' = 'fish_id'), suffix = c('', '_neg'), relationship = 'many-to-many'
    ) %>%
    group_by(triplet) %>%
    slice_sample(n = 1) %>%
    collect()

  tr <- tr %>%
    ungroup() %>%
    slice_sample(n = pmin(n, nrow(tr)))

  if (nrow(tr) < (0.8 * n)) {
    warning(paste(n, 'triplets requested, but only', nrow(tr), 'triplets generated.'))
  }

  if (length(dim(image_array)) == 4) {
    ret <- list(
      anchors   = image_array[tr$sample_id, , , , drop = FALSE],
      positives = image_array[tr$sample_id_pos   , , , , drop = FALSE],
      negatives = image_array[tr$sample_id_neg   , , , , drop = FALSE]
    )
  } else {
    ret <- list(
      anchors   = image_array[tr$sample_id, , , drop = FALSE],
      positives = image_array[tr$sample_id_pos   , , , drop = FALSE],
      negatives = image_array[tr$sample_id_neg   , , , drop = FALSE]
    )
  }
  return(ret)
}

generate_triplet_list <- function(
  data, validation = FALSE, batch_match = FALSE, method = 'self', pedigree, min_r = 0.25,
  remove_maternal_effect = TRUE, pairs_limit = 8
) {
  # data <- db; pedigree <- ped_df
  require(data.table)
  data <- as.data.table(data)
  if (batch_match && method != 'pedigree') stop('batch_match only implemented for pedigree triplets.')
  if (validation) {
    data <- data[data$validation, ]
  } else {
    data <- data[!data$validation, ]
  }
  if (method == 'self') {
    anchors <- data[, .(anchor_fish_id = fish_id, anchor_id = unique_id, anchor_sample_id = sample_id)]
    positives <- data[, .(pos_fish_id = fish_id, pos_id = unique_id, pos_sample_id = sample_id)]
    negatives <- data[, .(neg_fish_id = fish_id, neg_id = unique_id, neg_sample_id = sample_id)]

    pairs <- anchors[positives, on = .(anchor_fish_id = pos_fish_id), allow.cartesian = TRUE] %>%
      .[anchor_id != pos_id, ] %>%
      .[, dummy := 1]
    # limit this to 'pairs_limit' photo pairs per anchor fish
    pairs <- pairs %>%
      group_by(anchor_fish_id) %>%
      mutate(i = sample(n())) %>%
      ungroup() %>%
      filter(i <= pairs_limit) %>%
      select(-i) %>%
      as.data.table()

    triplets <- pairs[negatives[, dummy := 1], on = .(dummy), allow.cartesian = TRUE] %>%
      .[anchor_fish_id != neg_fish_id, ] %>%
      .[, c('dummy', 'anchor_fish_id', 'anchor_id', 'neg_fish_id', 'pos_id', 'neg_id') := NULL]
  }
  if (method == 'same_side') {
    anchors <- data[, .(anchor_fish_id = fish_id, anchor_id = unique_id,
                        anchor_sample_id = sample_id, anchor_side = facing_direction)
    ]
    positives <- data[, .(pos_fish_id = fish_id, pos_id = unique_id, pos_sample_id = sample_id,
                          pos_side = facing_direction)]
    negatives <- data[, .(neg_fish_id = fish_id, neg_id = unique_id, neg_sample_id = sample_id)]

    pairs <- anchors[
      positives,
      on = .(anchor_fish_id = pos_fish_id, anchor_side = pos_side),
      allow.cartesian = TRUE
    ] %>%
      .[anchor_id != pos_id, ] %>%
      .[, dummy := 1]

    triplets <- pairs[negatives[, dummy := 1], on = .(dummy), allow.cartesian = TRUE] %>%
      .[anchor_fish_id != neg_fish_id, ] %>%
      .[, c('dummy', 'anchor_fish_id', 'anchor_id', 'anchor_side', 'neg_fish_id', 'pos_id', 'neg_id') := NULL]
  }
  if (method == 'pedigree') {
    if (!remove_maternal_effect) stop('Currently, maternal brothers are always ignored.')
    if (batch_match) {
      # first, define the three tables. For each, select and rename variables, then pick only 1
      # photo per anchor fish, at random
      anchors <- data[
        ,
        .(anchor_fish_id = fish_id, anchor_sample_id = sample_id, anchor_batch = batch)
      ] %>%
        .[, .SD[sample(.N, 1)], by = anchor_fish_id]
      positives <- data[
        ,
        .(pos_fish_id = fish_id, pos_sample_id = sample_id)
      ] %>%
        .[, .SD[sample(.N, 1)], by = pos_fish_id]
      negatives <- data[
        ,
        .(neg_fish_id = fish_id, neg_sample_id = sample_id, neg_batch = batch)
      ] %>%
        .[, .SD[sample(.N, 1)], by = neg_fish_id]

      # Pair up
      pairs <- anchors[, dummy := 1] %>%
        .[positives[, dummy := 1], on = .(dummy), allow.cartesian = TRUE] %>%
        .[
          pedigree[r >= min_r & r < 0.9999 & !shared_dam, ],
          on = .(anchor_fish_id = id1, pos_fish_id = id2),
          nomatch = 0
        ] %>%
        .[, c('r', 'pos_fish_id') := NULL]

      triplets <- pairs[negatives, on = .(anchor_batch = neg_batch), allow.cartesian = TRUE] %>%
        .[pedigree[r == 0, ], on = .(anchor_fish_id = id1, neg_fish_id = id2), nomatch = 0] %>%
        .[, c('dummy', 'anchor_fish_id', 'neg_fish_id', 'anchor_batch', 'r') := NULL]
    } else {
      anchors <- data[, .(anchor_fish_id = fish_id, anchor_sample_id = sample_id)] %>%
        .[, .SD[sample(.N, 1)], by = anchor_fish_id]
      positives <- data[, .(pos_fish_id = fish_id, pos_sample_id = sample_id)] %>%
        .[, .SD[sample(.N, 1)], by = pos_fish_id]
      negatives <- data[, .(neg_fish_id = fish_id, neg_sample_id = sample_id)] %>%
        .[, .SD[sample(.N, 1)], by = neg_fish_id]

      pairs <- anchors[, dummy := 1] %>%
        .[positives[, dummy := 1], on = .(dummy), allow.cartesian = TRUE] %>%
        .[pedigree[r >= min_r & r < 0.9999, ], on = .(anchor_fish_id = id1, pos_fish_id = id2), nomatch = 0] %>%
        .[, c('r', 'pos_fish_id') := NULL]

      triplets <- pairs[negatives[, dummy := 1], on = .(dummy), allow.cartesian = TRUE] %>%
        .[pedigree[r == 0, ], on = .(anchor_fish_id = id1, neg_fish_id = id2), nomatch = 0] %>%
        .[, c('dummy', 'anchor_fish_id', 'neg_fish_id', 'r') := NULL]
    }
  }
  return(triplets)
}

# we need to generate triplets on the fly, so let's define an easy to use function
generate_triplets <- function(
  triplet_list, image_array, embeddings, n_samples = NULL, triplets = NULL, n_dim, method,
  validation = FALSE, random_gamma = FALSE, gamma_adj_level = 2
) {
  if (!(method %in% c('random', 'semi-hard'))) stop("Your method isn't a valid option.")
  # n_samples = 100; validation = FALSE; batch_match = TRUE
  if (method == 'random') {
    triplet_list <- triplet_list[sample(nrow(triplet_list), n_samples)]

    if (length(dim(image_array)) == 4) {
      ret <- list(
        anchors   = image_array[triplet_list$anchor_sample_id, , , ],
        positives = image_array[triplet_list$pos_sample_id   , , , ],
        negatives = image_array[triplet_list$neg_sample_id   , , , ]
      )
    } else {
      ret <- list(
        anchors   = image_array[triplet_list$anchor_sample_id, , ],
        positives = image_array[triplet_list$pos_sample_id   , , ],
        negatives = image_array[triplet_list$neg_sample_id   , , ]
      )
    }
  }
  if (method == 'semi-hard') {
    require(data.table)
    if (!is.null(triplets)) triplet_list <- triplet_list[triplets, ]

    em <- triplet_list[embeddings, on = c('anchor_sample_id' = 'id'), nomatch = NULL]
    setnames(em, paste0('V', seq_len(n_dim)), paste0('A', 1:n_dim))
    em <- em[embeddings, on = c('pos_sample_id' = 'id'), nomatch = NULL]
    setnames(em, paste0('V', seq_len(n_dim)), paste0('P', 1:n_dim))
    em <- em[embeddings, on = c('neg_sample_id' = 'id'), nomatch = NULL]
    setnames(em, paste0('V', seq_len(n_dim)), paste0('N', 1:n_dim))

    df <- function(dt) { rowSums((dt[, 1:n_dim] - dt[, (n_dim + 1):(2 * n_dim)]) ^ 2) }
    em[, DAP := df(.SD), .SDcols = c(paste0('A', 1:n_dim), paste0('P', 1:n_dim))]
    em[, DAN := df(.SD), .SDcols = c(paste0('A', 1:n_dim), paste0('N', 1:n_dim))]
    # drop unnesseary columns
    em[, c(paste0('A', 1:n_dim), paste0('P', 1:n_dim), paste0('N', 1:n_dim)) := NULL]

    # find only the semi-hard samples
    semi_hards <- em[(DAP < DAN) & ((DAP - DAN + margin) > 0), , ]
    if (!is.null(n_samples) & (nrow(semi_hards) > n_samples)) {
      semi_hards <- semi_hards[sample(.N, n_samples), ]
    }
    if (length(dim(image_array)) == 4) {
      ret <- list(
        anchors   = image_array[semi_hards$anchor_sample_id, , , ],
        positives = image_array[semi_hards$pos_sample_id   , , , ],
        negatives = image_array[semi_hards$neg_sample_id   , , , ]
      )
    } else {
      ret <- list(
        anchors   = image_array[semi_hards$anchor_sample_id, , ],
        positives = image_array[semi_hards$pos_sample_id   , , ],
        negatives = image_array[semi_hards$neg_sample_id   , , ]
      )
    }
  }

  if (random_gamma && !validation) {
    ns <- dim(ret$anchors)[1]
    ret$anchors   <- ret$anchors   ^ runif(ns, min = 1 / gamma_adj_level, max = gamma_adj_level)
    ret$positives <- ret$positives ^ runif(ns, min = 1 / gamma_adj_level, max = gamma_adj_level)
    ret$negatives <- ret$negatives ^ runif(ns, min = 1 / gamma_adj_level, max = gamma_adj_level)
  }
  return(ret)
}

load_images_into_array <- function(img_dirs, input_dim, n = NULL) {
  # img_dir should be a list of image-directories, concatenated as color channels
  # If img_dir has length 1, an input_dim of (r, c, 1) will extract the alpha channel, and a
  #   input_dim of (r, c, 3) will extract the rgb image
  # n is the number of images to load (first n images). If n is NULL, load all.
  #
  # img_dirs <- list('data/carotenoid_coloration_warped', 'data/melanic_coloration_warped_v2')
  require(furrr)
  require(imager)

  img <- map(img_dirs, \(i) {
    l <- list.files(i, full.names = TRUE, recursive = TRUE) %>%
      setNames(basename(.) %>% tools::file_path_sans_ext())
    # deal with folder icon files on macos (sigh)
    l <- l[!str_detect(l, 'Icon')] # using str_subset() removes names, so doing it this way
    if (!is.null(n)) {
      l <- l[seq_len(n)]
    }
    return(l)
  })
  stopifnot(n_distinct(lengths(img)) == 1)

  # check if we need full color images, or only the 1D mask (from the alpha channel)
  if (length(img_dirs) == 1 && input_dim[3] == 3) {
    # check if we are working on full sized images, or reduce it to half size
    if (input_dim[1] == 250) {
      my_load_image <- \(f) load.image(f) %>% imager::channel(1:3) %>% resize_halfXY() %>% as.array()
    } else {
      my_load_image <- \(f) load.image(f) %>% imager::channel(1:3) %>% as.array()
    }
  } else {
    # check if we are working on full sized images, or reduce it to half size
    if (input_dim[1] == 250) {
      my_load_image <- \(f) load.image(f) %>% imager::channel(4) %>% resize_halfXY() %>% as.array()
    } else {
      my_load_image <- \(f) load.image(f) %>% imager::channel(4) %>% as.array()
    }
  }

  all_images <- map(img, \(imgs) {
    # L <- future_map(
    #   imgs,
    #   \(x) my_load_image(x)[, , 1, 1],
    #   .progress = TRUE, .options = furrr_options(seed = NULL)
    # )
    L <- map(
      imgs,
      \(x) my_load_image(x),
      .progress = TRUE
    )
    abind::abind(L, along = 3) # We'll bind along 3, since that's the empty "frame" dimension
  })
  #plot(as.cimg(all_images[[1]][ , , 3, ]))
  all_images <- abind::abind(all_images, along = 5)
  #plot(as.cimg(all_images[ , , 3, , 1]))
  all_images <- aperm(all_images, c(3, 1:2, 4, 5))
  #plot(as.cimg(all_images[ 3, , , 1, 1]))

  # make images binary (but only if we asked for the mask)
  if (input_dim[3] == 1) {
    all_images[] <- as.numeric(all_images > 0.5)
  }
  return(all_images)
}
