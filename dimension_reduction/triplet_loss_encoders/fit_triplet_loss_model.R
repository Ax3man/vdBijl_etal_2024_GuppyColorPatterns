fit_triplet_loss_model <- function(
  # in and output
  img_dir,
  output_name,
  use_cpu = FALSE,
  method = 'self', # method of choosing triplets
  batch_match = TRUE, # whether to only generate triplets within replicate:generation
  min_r = 0.25,    # min relatedness (full-sib = 0.5) for positive example, if method = 'pedigree'

  # model architecture parameters
  input_dim = c(250, 70, 1), # last is color channel
  embed_dim = 2,
  filter_base = 32,
  margin = 1,
  dropout = 0,

  # training parameters
  epochs = 100,
  batch_size = 64,
  n_validation_triplets = 1000,
  n_samples = 5000,

  # use a schedule for the learning rate, data.frame of epoch and lr
  lr_schedule = data.frame(epoch = c(0, 75, 90), lr = 5e-4 * c(1, 0.1, 0.01)),

  # whether to refit or load from a previous solution
  refit = TRUE
) {
  # test arguments
  # img_dir <- list('data/carotenoid_coloration_warped/', 'data/melanic_coloration_warped_v2/')
  # output_name = paste0('double_model_comparsion/test')
  # use_cpu = TRUE; input_dim = c(250, 70, 1); embed_dim = 2; filter_base = 32; margin = 1; dropout = 0
  # epochs = 100; batch_size = 64; n_validation_triplets = 500; n_samples = 1e3; n_triplets = 1e6
  # lr_schedule = data.frame(epoch = c(0, 60, 80), lr = 1e-4 * c(1, 0.1, 0.01)); refit = TRUE

  cat('Fitting new model...\nStarting with loading data and general setup...\n')

  library(keras)
  library(tidyverse)
  library(future)
  library(furrr)
  library(imager)
  library(tensorflow)
  library(data.table)
  source('dimension_reduction/triplet_loss_encoders/triplet_loss_tools.R')

  reticulate::use_virtualenv('r-reticulate')
  tf_version()
  # gpu <- tf$config$experimental$get_visible_devices('GPU')[[1]]
  # tf$config$experimental$set_memory_growth(gpu, TRUE)
  #virtual_gpu <- tf$config$experimental$VirtualDeviceConfiguration(memory_limit = 6000)
  #tf$config$experimental$set_virtual_device_configuration(gpu, list(virtual_gpu))
  #if (tensorflow::tf$executing_eagerly()) tensorflow::tf$compat$v1$disable_eager_execution()

  if (use_cpu) {
    cpu <- tf$config$get_visible_devices()[[1]]
    tf$config$set_visible_devices(cpu)
  }

  ## Model and modelling parameters ------------------------------------------------------------------

  model <- build_model(
    embed_dim = embed_dim, filter_base = filter_base, input_dim = input_dim, lr = lr_schedule$lr[1],
    margin = margin, type = 3, dropout = dropout
  )

  ## Data loading ------------------------------------------------------------------------------------
  plan(multisession, workers = 20)
  all_images <- load_images_into_array(img_dir, input_dim)
  all_images <- drop(all_images)
  plan(sequential)

  cat('Finished loading data...\n')

  source('selection_decisions/compile_decisions.R')
  img_names <- dimnames(all_images)[[1]]
  db <- data.table::fread('photo_database.csv') %>%
    select(replicate, generation, fish_id, unique_id, facing_direction) %>%
    as_tibble() %>%
    # drop rows for which there is (not yet) an image
    filter(unique_id %in% img_names) %>%
    # reorder to be in the same order as img
    slice(match(img_names, unique_id)) %>%
    # make a variable to denote which images are training, and which validation, also add sample id
    mutate(
      validation = fish_id %in% sample(unique(fish_id), floor(n_distinct(fish_id) * 0.2)),
      sample_id = row_number(),
      fish_id = toupper(fish_id),
      generation = factor(generation, c('parental_gen_1', 'gen_2', 'gen_3'), c('P', 'F1', 'F2'))
    ) %>%
    left_join(select(selection, fish_id, selection, car_perc), 'fish_id')

  ped <- data.table::fread('data/pedigree.csv')
  ped_df <- ped %>%
    select(animal, dam, sire) %>%
    nadiv::makeA() %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column('id1') %>%
    pivot_longer(-id1, names_to = 'id2', values_to = 'r') %>%
    as.data.table() %>%
    left_join(select(ped, animal, dam1 = dam), c('id1' = 'animal')) %>%
    left_join(select(ped, animal, dam2 = dam), c('id2' = 'animal')) %>%
    mutate(shared_dam = coalesce(dam1 == dam2, FALSE)) %>%
    select(-dam1, -dam2)

  # plot(as.raster(aperm(all_images[sample(dim(all_images)[1], 1), , ], c(2:1))))

  # Show some triplets:
  # ns <- 10
  # ex <- generate_triplet_subset(
  #   db, ns, all_images, method = method, validation = FALSE, batch_match = batch_match
  # )
  # opa <- par(mfrow = c(ns, 3), oma = c(0, 0, 0, 0), mar = c(1, 1, 1, 1)/10)
  # for (i in seq_len(ns)) {
  #   for (j in 1:3) {
  #     suppressWarnings(plot(as.cimg(ex[[j]][i, , , ]), axes = FALSE, rescale = FALSE))
  #     #plot(as.raster(t(triplets[[1]][[j]][i, , , 1])))
  #   }
  # }
  # par(opa)

  if (!refit) {
    message('Loading weights from disk...')
    # stop("This function currently doesn't know where to find the correct weights. You should
    #      probably set `refit` to TRUE")
    load_model_weights_hdf5(
      model,
      paste0('dimension_reduction/triplet_loss_encoders/weights/', output_name, '.hdf5')
    )
  }

  # Create a subsection of the model, that is only the shared layers that yields the n-dimensional
  # embedding.

  embedder <- keras_model_sequential(
    layers = get_layers(model, c(1, 4:find_layer(model, 'norm')))
  )
  #print(embedder)

  # set up empty data.frames to store learning metrics
  train_metrics <- val_metrics <- data.frame()

  # train 1 epoch at a time, so we can keep generating semi-hard triplets
  for (epoch in seq_len(epochs)) {
    if (epoch %in% lr_schedule$epoch) {
      k_set_value(model$optimizer$learning_rate, lr_schedule[lr_schedule$epoch == epoch, 'lr'])
    }
    cat('Epoch:', epoch, '- Generating triplets...\n')
    training_triplets <- generate_triplet_subset(
      db, n_samples, all_images, method = method, pedigree = ped_df,
      validation = FALSE, batch_match = batch_match
    )
    # fail with memory error, are the whole triplet sets copied when sent to python??
    # Sink <- train_on_batch(model, x = training_triplets, y = rep(1, nrow(training_triplets[[1]])))

    actual_n_samples <- dim(training_triplets$anchors)[1]
    cat('Epoch:', epoch, '- Training on', actual_n_samples, 'triplets')

    training_triplets <- chunk_triplets(training_triplets, batch_size)

    for (batch in seq_along(training_triplets)) {
      y <- rep(1, times = nrow(training_triplets[[batch]][[1]]))
      Sink <- train_on_batch(model, x = training_triplets[[batch]], y = y)
      cat('.')
    }
    cat('\n')

    cat('Epoch:', epoch, '- Evaluating metrics...\n')
    train_metrics <- generate_triplet_subset(
      db, n = if (epoch == epochs) {n_samples} else {n_validation_triplets}, pedigree = ped_df,
      image_array = all_images, method = method, validation = FALSE, batch_match = FALSE
    ) %>%
      get_metrics(embedder, epoch, batch_size = batch_size, margin = margin) %>%
      bind_rows(train_metrics, .)

    val_metrics <- generate_triplet_subset(
      db, n = if (epoch == epochs) {n_samples} else {n_validation_triplets}, pedigree = ped_df,
      image_array = all_images, method = method, validation = TRUE, batch_match = FALSE
    ) %>%
      get_metrics(embedder, epoch, batch_size = batch_size, margin = margin) %>%
      bind_rows(val_metrics, .)
    bind_rows(
      training = slice(train_metrics, n()),
      validation = slice(val_metrics, n()),
      .id = 'type'
    ) %>% mutate(across(where(is.numeric), \(x) round(x, 4))) %>% print()


    # plot my progress
    metrics <- bind_rows(training = train_metrics, validation = val_metrics, .id = 'data')
    p <- ggplot(metrics, aes(epoch, color = data)) +
      geom_vline(xintercept = lr_schedule$epoch[-1], lty = 2, color = 'grey40') +
      geom_point() +
      scale_color_manual(values = c(training = 'firebrick', validation = 'navy')) +
      xlim(1, epochs) +
      theme_bw() +
      theme(plot.margin = unit(c(0, 12, 0, 12), 'points'))
    if (epoch > 5) p <- p + geom_smooth(method = 'loess', formula = 'y ~ x', se = FALSE)
    sl <- scale_y_log10(); g <- guides(color = 'none'); eb <- element_blank()
    xt <- theme(axis.text.x = eb, axis.title.x = eb, axis.ticks.x = eb)
    pc <- cowplot::plot_grid(
      p + aes(y = loss) + list(sl, g, xt) + labs(y = 'Triplet loss'),
      p + aes(y = DAP) + list(g, xt) + labs(y = 'd(A, P)'),
      p + aes(y = DAN) + xt + labs(y = 'd(A, N)'),
      p + aes(y = delta_dist) + list(g, xt) + labs(y = 'd(A, P) - d(A, N)'),
      p + aes(y = acc) + geom_errorbar(aes(ymin = acc_lower, ymax = acc_upper), width = 0.1) +
        labs(y = 'Accuracy') + list(g, xt), # expand_limits(y = 1),
      p + aes(y = acc2) + g + labs(y = 'Accuracy (incl margin)'), # expand_limits(y = 1),
      ncol = 1, align = 'hv', axis = 'lr'
    )
    print(pc)
    ggsave(
      paste0('dimension_reduction/triplet_loss_encoders/histories/', output_name, '.png'),
      pc, h = 12, w = 8
    )
    #sink <- gc() # force garbage collection, not sure it does anything.
  }

  write_rds(
    metrics,
    paste0('dimension_reduction/triplet_loss_encoders/histories/', output_name, '.rds')
  )
  save_model_weights_hdf5(
    model,
    paste0('dimension_reduction/triplet_loss_encoders/weights/', output_name, '.hdf5')
  )
  cat('Weights saved.\nCalculating embeddings...\n')
  embeddings <- predict(embedder, all_images)
  rownames(embeddings) <- rownames(all_images)
  write_rds(
    embeddings,
    paste0('dimension_reduction/triplet_loss_encoders/embeddings/', output_name, '.rds')
  )
}
