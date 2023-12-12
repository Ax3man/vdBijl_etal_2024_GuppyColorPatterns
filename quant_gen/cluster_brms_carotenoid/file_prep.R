suppressMessages(library(imager))
library(furrr)
library(tidyverse)
library(cmdstanr)
eb <- element_blank()

source('quant_gen/prepare_pedigrees.R')

#warning('Running this job on 46 cores, this better be on the HPC.')
#plan(multisession, workers = 46)
warning('Running this job on only 12, this better be local.')
plan(multisession, workers = 12)

ims <- list.files('data/carotenoid_coloration_warped/', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

# load all images. For memory reasons, we keep just the pixels with carotenoid. Pixels not in df are
# by definition not carotenoid
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

df <- df %>% group_by(x, y) %>% group_split()

pedigree <- data.table::fread('data/pedigree.csv') %>% mutate(across(c(animal, sire, dam), tolower))
photo_database <- data.table::fread('photo_database.csv') %>%
  select(unique_id, fish_id, replicate, generation, facing_direction) %>%
  add_patriline(ped_df)

## Statistical maps
generate_model_files <- function(Data, cutoff = 0.05, verbose = FALSE) {
  # Data <- df[[600]]
  xcoord <- Data$x[1]; ycoord <- Data$y[1]

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

  # code_file <- paste0(
  #   'quant_gen/cluster_brms_carotenoid/',
  #   'stan_code/',
  #   'pixel_',
  #   formatC(xcoord, width = 3, format = 'd', flag = '0'), '_',
  #   formatC(ycoord, width = 3, format = 'd', flag = '0'),
  #   '.stan'
  # )
  data_file <- paste0(
    'quant_gen/cluster_brms_carotenoid/',
    'stan_data/',
    'pixel_',
    formatC(xcoord, width = 3, format = 'd', flag = '0'), '_',
    formatC(ycoord, width = 3, format = 'd', flag = '0'),
    '.json'
  )

  ## Only have to make this once, it's the same for each dataset.
  #
#   code <- make_stancode(
#     carotenoid ~ 1 +
#       (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) + (1 | patriline) +
#       (1 | dam) + (1 | grow_tank) + (1 | fish_id) +
#       (1 | fish_id:facing_direction),
#     data = Data, data2 = list(A = A, X = X),
#     family = 'bernoulli',
#     priors = c(
#       prior(student_t(3, 0, 5), class = b),
#       prior(student_t(3, 0, 5), class = sd)
#     ),
#     save_model = code_file
#   )

  data <- make_standata(
    carotenoid ~ 1 +
      (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) + (1 | patriline) +
      (1 | dam) + (1 | grow_tank) + (1 | fish_id) +
      (1 | fish_id:facing_direction),
    data = Data, data2 = list(A = A, X = X),
    family = 'bernoulli',
    priors = c(
      prior(student_t(3, 0, 5), class = b),
      prior(student_t(3, 0, 5), class = sd)
    )
  )
  cmdstanr::write_stan_json(data, data_file)

  if (verbose) cat('Done!\n')
  return(invisible(NULL))
}

cat('\n', length(df), 'pixels to analyze...\nStarting to fit animal models...\n\n')

future_walk(
  df,
  generate_model_files,
  .progress = TRUE
)

plan(sequential)
