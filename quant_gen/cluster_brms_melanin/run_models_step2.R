suppressMessages(library(imager))
library(furrr)
library(tidyverse)
library(cmdstanr)

source('quant_gen/prepare_pedigrees.R')

ims <- list.files('data/melanic_coloration_warped_v2/', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

# load all images. For memory reasons, we keep just the pixels with melanin. Pixels not in df are
# by definition not melanin
plan(multisession, workers = 16)
df <- future_map_dfr(
  ims,
  ~load.image(.x) %>%
    resize_halfXY() %>% #resize_halfXY() %>%
    as.data.frame(wide = 'c') %>%
    filter(c.4 > 0.5) %>%
    dplyr::select(x, y),
  .id = 'image',
  .progress = TRUE, .options = furrr_options(seed = NULL)
) %>% as_tibble()
plan(sequential)

pedigree <- data.table::fread('data/pedigree.csv') %>% mutate(across(c(animal, sire, dam), tolower))
photo_database <- data.table::fread('photo_database.csv') %>%
  select(unique_id, fish_id, replicate, generation, facing_direction) %>%
  add_patriline(ped_df)

#################################################################################################
### IMPORTANT: We skip every 5th pixel, both x and y, since we fitted larger models for those ###
#################################################################################################

df <- filter(df, !((x %% 5) == 0 & (y %% 5) == 0))
df <- df %>% group_by(x, y) %>% group_split()
# run them in random order
df <- sample(df)

# get the coordinates of the models that were fully fit
full_models <- data.frame(file = list.files('quant_gen/cluster_brms_melanin/posteriors_1/')) %>%
  separate(file, into = c(NA, 'x', 'y', NA), convert = TRUE, remove = FALSE)

## Statistical maps
fit_second_models <- function(Data, cutoff = 0.01, verbose = FALSE) {
  # Data <- df[[1000]]
  xcoord <- Data$x[1]; ycoord <- Data$y[1]

  # find closest sampled pixels
  closest_models <- full_models %>%
    mutate(dist = sqrt((x - xcoord) ^ 2 + (y - ycoord) ^ 2)) %>%
    slice_min(dist, n = 9, with_ties = FALSE)

  stepsize <- closest_models$file %>%
    paste0('quant_gen/cluster_brms_melanin/meta_1/', .) %>%
    map(read_rds) %>%
    map(c('metadata', 'step_size_adaptation')) %>%
    map_dbl(mean) %>%
    weighted.mean(w = 1 / closest_models$dist)

  params <- closest_models$file %>%
    paste0('quant_gen/cluster_brms_melanin/posterior_summaries_1/', .) %>%
    map_dfr(~read_rds(.x)[c('variable', 'mean')], .id = 'model') %>%
    left_join(data.frame(model = as.character(1:9), w = 1 / closest_models$dist), 'model') %>%
    group_by(variable) %>%
    summarise(
      sd = sd(mean),
      mean = weighted.mean(mean, w)
    )

  post_file <- paste0(
    'quant_gen/cluster_brms_melanin/posteriors_2/pixel_',
    formatC(xcoord, width = 4, format = 'd', flag = '0'),
    '_',
    formatC(ycoord, width = 4, format = 'd', flag = '0'),
    '.rds'
  )
  post_summ_file <- gsub('posteriors', 'posterior_summaries', post_file)
  meta_file <- gsub('posteriors', 'meta', post_file)

  # reconstruct full dataset
  Data <- Data %>%
    select(image) %>%
    mutate(melanin = 1) %>%
    complete(image = setdiff(names(ims), image), fill = list(melanin = 0)) %>%
    left_join(photo_database, c('image' = 'unique_id')) %>%
    left_join(select(ped_df, animal, dam, grow_tank), by = c('fish_id' = 'animal')) %>%
    mutate(fish_idA = fish_id, fish_idX = fish_id) %>%
    # drop individuals without a dam (i.e. the parental generation)
    drop_na(dam)

  if (mean(Data$melanin) < cutoff) return(invisible(NULL))
  if (verbose) cat(
    'This data has', nrow(Data), 'obervations, with an incidence of', round(mean(Data$melanin), 3), '\n'
  )

  stan <- readLines('quant_gen/cluster_brms_melanin/model2.stan')
  stan[68] <- glue::glue_data(
    filter(params, variable == 'Intercept'),
    '  lprior += normal_lpdf(Intercept | {mean}, {sd});'
  )
  stan[69] <- glue::glue_data(
    filter(params, variable == 'sd_1[1]'),
    '  lprior += normal_lpdf(sd_1 | {mean}, {sd}) - 1 * normal_lccdf(0 | {mean}, {sd});'
  )
  stan[70] <- glue::glue_data(
    filter(params, variable == 'sd_2[1]'),
    '  lprior += normal_lpdf(sd_2 | {mean}, {sd}) - 1 * normal_lccdf(0 | {mean}, {sd});'
  )
  stan[71] <- glue::glue_data(
    filter(params, variable == 'sd_3[1]'),
    '  lprior += normal_lpdf(sd_3 | {mean}, {sd}) - 1 * normal_lccdf(0 | {mean}, {sd});'
  )
  stan[72] <- glue::glue_data(
    filter(params, variable == 'sd_4[1]'),
    '  lprior += normal_lpdf(sd_4 | {mean}, {sd}) - 1 * normal_lccdf(0 | {mean}, {sd});'
  )
  stan[73] <- glue::glue_data(
    filter(params, variable == 'sd_5[1]'),
    '  lprior += normal_lpdf(sd_5 | {mean}, {sd}) - 1 * normal_lccdf(0 | {mean}, {sd});'
  )
  stan_file <- write_stan_file(stan)
  model <- cmdstan_model(stan_file)

  data <- list(
    # total nr of observations
    N = as.integer(nrow(Data)),
    # response variable
    Y = array(Data$melanin),
    # fixed effects design matrix(?) Only an intercept, so only 1s
    X = matrix(rep_len(1, nrow(Data)), dimnames = list(seq_len(nrow(Data)), 'Intercept')),
    # Random effect design matrices, all model only an intercept, so all 1s
    Z_1_1 = array(rep_len(1, nrow(Data))), Z_2_1 = array(rep_len(1, nrow(Data))),
    Z_3_1 = array(rep_len(1, nrow(Data))), Z_4_1 = array(rep_len(1, nrow(Data))),
    Z_5_1 = array(rep_len(1, nrow(Data))),
    # The J matrices give the grouping indicator per observation. This is just a numeric coding of the factors
    J_1 = as.integer(as.factor(Data$fish_id)),
    J_2 = as.integer(as.factor(Data$fish_idA)),
    J_3 = as.integer(as.factor(Data$fish_idX)),
    J_4 = as.integer(as.factor(Data$grow_tank)),
    J_5 = as.integer(as.factor(Data$patriline)),
    # N for number of grouping levels
    N_1 = n_distinct(Data$fish_id),
    N_2 = n_distinct(Data$fish_idA),
    N_3 = n_distinct(Data$fish_idX),
    N_4 = n_distinct(Data$grow_tank),
    N_5 = n_distinct(Data$patriline),
    # M for number of coefficients per level (only intercepts for me, so all 1)
    M_1 = 1L, M_2 = 1L, M_3 = 1L, M_4 = 1L, M_5 = 1L,
    # NC for IDK, but they are all 0
    NC_1 = 0L, NC_2 = 0L, NC_3 = 0L, NC_4 = 0L, NC_5 = 0L,
    # Cholesky factor of known covariance matrices
    Lcov_2 = t(chol(A[sort(unique(Data$fish_idA)), sort(unique(Data$fish_idA))])),
    Lcov_3 = t(chol(X[sort(unique(Data$fish_idX)), sort(unique(Data$fish_idX))])),
    # K for IDK
    K = 1L,
    # whether to only sample the prior, 0 is FALSE
    prior_only = 0L
  )

  make_z <- function(z) { #z should be the string identifying the z, e.g. 'z_1'
    xx <- filter(params, str_starts(variable, z)) %>%
      select(-sd) %>%
      separate(variable, c(NA, NA, 'i', 'j', NA), convert = TRUE) %>%
      arrange(i, j)
    matrix(xx$mean, nrow = max(xx$i), ncol = max(xx$j), byrow = TRUE)
  }
  inits <- list(
    Intercept = filter(params, variable == 'Intercept')$mean,
    sd_1 = filter(params, variable == 'sd_1[1]')$mean,
    z_1 = make_z('z_1'),
    sd_2 = filter(params, variable == 'sd_2[1]')$mean,
    z_2 = make_z('z_2'),
    sd_3 = filter(params, variable == 'sd_3[1]')$mean,
    z_3 = make_z('z_3'),
    sd_4 = filter(params, variable == 'sd_4[1]')$mean,
    z_4 = make_z('z_4'),
    sd_5 = filter(params, variable == 'sd_5[1]')$mean,
    z_5 = make_z('z_5')
  )

  set.seed(NULL) # reset RNG
  fit <- model$sample(
    data = data,
    refresh = 100,
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    init_buffer = 10, # short adaptation time for step_size, use warmup for estimating the metric
    term_buffer = 10, # see init_buffer
    adapt_delta = 0.8,
    step_size = stepsize,
    init = list(inits, inits)
  )
  # full posteriors for just the interesting parameters
  draws <- fit$draws(c('Intercept', paste0('sd_', 1:5)))
  write_rds(draws, post_file)

  # posterior summaries for each parameter, that's ~20K paramaters
  summ <- fit$summary()
  write_rds(summ, post_summ_file)

  # write the meta data
  meta <- list(
    diagnostic_summary = fit$diagnostic_summary(),
    cmdstan_diagnose = fit$cmdstan_diagnose(),
    metadata = fit$metadata()
  )
  write_rds(meta, meta_file)

  if (verbose) cat('Done!\n')
  return(invisible(NULL))
}

cat('\n', length(df), 'pixels to analyze...\nStarting to fit animal models...\n\n')

fork <- parallel::makeForkCluster(parallel::detectCores() - 2)
# warning('Running on 8 cores only..\n')
# fork <- parallel::makeForkCluster(8)

sink <- parallel::clusterMap(
  fork,
  fun = fit_second_models,
  df,
  .scheduling = 'dynamic'
)

if (FALSE) { #debug only
  sink <- lapply(df, fit_second_models, verbose = TRUE)
}
