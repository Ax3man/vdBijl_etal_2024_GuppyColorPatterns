suppressMessages(library(imager))
library(furrr)
library(tidyverse)
library(cmdstanr)

# compile the model if necessary
model <- cmdstan_model('quant_gen/cluster_brms_carotenoid/model.stan')

source('quant_gen/prepare_pedigrees.R')


ims <- list.files('data/carotenoid_coloration_warped/', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

# load all images. For memory reasons, we keep just the pixels with carotenoid. Pixels not in df are
# by definition not carotenoid
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

#####################################################################################
### IMPORTANT: We only take every 5th pixel, both x and y, for full model fitting ###
#####################################################################################

df <- filter(df, (x %% 5) == 0, (y %% 5) == 0)
df <- df %>% group_by(x, y) %>% group_split()

## Statistical maps
fit_first_models <- function(Data, cutoff = 0.01, verbose = FALSE) {
  # reset seed (for forking on the cluster?)
  set.seed(NULL)

  # Data <- df[[100]]
  xcoord <- Data$x[1]; ycoord <- Data$y[1]

  post_file <- paste0(
    'quant_gen/cluster_brms_carotenoid/posteriors_1/pixel_',
    formatC(xcoord, width = 4, format = 'd', flag = '0'),
    '_',
    formatC(ycoord, width = 4, format = 'd', flag = '0'),
    '.rds'
  )
  post_summ_file <- gsub('posteriors', 'posterior_summaries', post_file)
  meta_file <- gsub('posteriors', 'meta', post_file)

  if (all(file.exists(post_file, post_summ_file, meta_file))) return(NULL)

  # reconstruct full dataset
  Data <- Data %>%
    select(image) %>%
    mutate(carotenoid = 1) %>%
    complete(image = setdiff(names(ims), image), fill = list(carotenoid = 0)) %>%
    left_join(photo_database, c('image' = 'unique_id')) %>%
    left_join(select(ped_df, animal, dam, grow_tank), by = c('fish_id' = 'animal')) %>%
    mutate(fish_idA = fish_id, fish_idX = fish_id) %>%
    # drop individuals without a dam (i.e. the parental generation)
    drop_na(dam)

  if (mean(Data$carotenoid) < cutoff) return(invisible(NULL))
  if (verbose) cat(
    'This data has', nrow(Data), 'obervations, with an incidence of', round(mean(Data$carotenoid), 3), '\n'
  )

  ## Only have to make this once, it's the same for each dataset.
  #
  # code <- brms::make_stancode(
  #   carotenoid ~ 1 +
  #     (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) + (1 | patriline) +
  #     (1 | grow_tank) + (1 | fish_id),
  #   data = Data, data2 = list(A = A, X = X),
  #   family = 'bernoulli',
  #   priors = c(
  #     prior(student_t(3, 0, 5), class = b),
  #     prior(student_t(3, 0, 5), class = sd)
  #   ),
  #   save_model = 'model.stan'
  # )

  # example, but needs brms
  # data <- make_standata(
  #   carotenoid ~ 1 +
  #     (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) + (1 | patriline) +
  #     (1 | grow_tank) + (1 | fish_id),
  #   data = Data, data2 = list(A = A, X = X),
  #   family = 'bernoulli',
  #   priors = c(
  #     prior(student_t(3, 0, 5), class = b),
  #     prior(student_t(3, 0, 5), class = sd)
  #   )
  # )

  data <- list(
    # total nr of observations
    N = as.integer(nrow(Data)),
    # response variable
    Y = array(Data$carotenoid),
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
  model <- cmdstan_model('quant_gen/cluster_brms_carotenoid/model.stan')

  set.seed(NULL) # reset RNG
  fit <- model$sample(
    data = data,
    refresh = 100,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 1500,
    adapt_delta = 0.85,
    step_size = NULL # might change this to a good starting value later (should only affect speed)
  )
  # full posteriors for just the interesting variables
  draws <- fit$draws(c('Intercept', paste0('sd_', 1:5)))
  write_rds(draws, post_file)

  # posterior summaries for each variable, that's ~20K variables
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

sink <- parallel::clusterMap(
    fork,
    fun = fit_first_models,
    df,
    .scheduling = 'dynamic'
)
