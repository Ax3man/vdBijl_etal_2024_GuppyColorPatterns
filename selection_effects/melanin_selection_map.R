suppressMessages(library(imager))
library(furrr)
library(lme4)
library(emmeans)
library(tidyverse)
eb <- element_blank()

# We use an example of a complete fish to throw away all pixels we don't need, as well as carry
# over the edge transparency for a more pleasing visual result.
complete_fish <- load.image('data/extracted_fish_warped/replicate_3/gen_2/20210602_IMG_2097.png') %>%
  #resize_halfXY() %>%
  as.data.frame(wide = 'c') %>% filter(c.4 > 0.1) %>% select(x, y, alpha = c.4)

photo_database <- data.table::fread('photo_database.csv', data.table = FALSE) %>%
  select(unique_id, fish_id, replicate, generation, facing_direction) %>%
  filter(generation != 'parental_gen_1') %>%
  mutate(generation = factor(generation, c('gen_2', 'gen_3', 'gen_4'), c('F1', 'F2', 'F3')))

source('selection_decisions/compile_decisions.R')
selection <- select(selection, replicate, generation, fish_id, selection) %>%
  mutate(
    fish_id = tolower(fish_id),
    generation = factor(generation, c('F1', 'F2', 'F3'))
  )

ims <- list.files('data/melanic_coloration_warped_v2', full.names = TRUE, recursive = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())
ims <- ims[names(ims) %in% photo_database$unique_id]

plan(multisession, workers = 24)
df <- future_map_dfr(
  ims,
  ~load.image(.x) %>%
    #resize_halfXY() %>%
    as.data.frame(wide = 'c') %>%
    mutate(melanin = as.numeric(c.4 > 0.5)) %>%
    select(x, y, melanin) %>%
    inner_join(select(complete_fish, x, y), by = c('x', 'y')),
  .id = 'image',
  .progress = TRUE, .options = furrr_options(seed = NULL)
) %>% as_tibble()
plan(sequential)

df <- df %>%
  group_by(x, y) %>%
  filter(mean(melanin) > 0.01) %>%
  ungroup() %>%
  inner_join(photo_database, by = c('image' = 'unique_id')) %>%
  inner_join(selection, by = c('fish_id', 'replicate', 'generation')) %>%
  select(x, y, replicate, generation, selection, image, fish_id, facing_direction, melanin)

spl <- df %>%
  group_by(x, y) %>%
  group_split()
rm(df)

selection_effect <- function(Data, verbose = FALSE) {
  # Data <- spl[[3000]]
  if (mean(Data$melanin) < 0.01) {
    return(data.frame(generation = unique(Data$generation)))
  }

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  data_small <- Data %>%
    group_by(fish_id, selection, replicate, generation) %>%
    summarize(melanin = Mode(melanin), .groups = 'drop')
  m <- suppressMessages(suppressWarnings(
    glmer(
      melanin ~ 1 + selection * generation + (1 | generation:replicate),
      data_small, 'binomial'
    )
  ))
  # get the log odds ratio of selection effect by generation, using the marginal means
  em <- suppressMessages(suppressWarnings(
    emmeans(m, pairwise ~ selection, by = 'generation')
  ))

  #cat('I fitted a model!\n')

  means <- as.data.frame(summary(em$emmeans)) %>% mutate(x = Data$x[1], y = Data$y[1])
  contrast <- as.data.frame(summary(em$contrasts)) %>% mutate(x = Data$x[1], y = Data$y[1])

  return(list(means = means, contrast = contrast))
}

cat('\n', length(spl), 'pixels to analyze...\n')

plan(multisession, workers = 24)
model_results <- future_map(
  spl,
  possibly(selection_effect, data.frame(generation = unique(spl[[1]]$generation))),
  .progress = TRUE, .options = furrr_options(seed = NULL)
)
plan(sequential)

write_rds(model_results, 'visualization/per_pixel_models/mel_selection_model_results.rds')

