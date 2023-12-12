library(tidyverse)
library(nadiv)

source('selection_decisions/compile_decisions.R')

add_patriline <- function(data, pedigree) {
  # Create new column with the patriline, empty for now
  data <- mutate(data, patriline = NA)
  # Loop through the pedigree.
  for (i in seq_len(nrow(pedigree))) {
    id <- pedigree[i, 'animal', drop = TRUE]
    # If the fish is female, we use NA
    if (pedigree[i, 'sex', drop = TRUE] == 'F') {
      current_id <- NA
    } else {
      # Initialize the pedigree chain, with the current individuals sire
      current_id <- id
      while (TRUE) {
        # Find the sire
        sire <- pedigree[pedigree$animal == current_id, 'sire', drop = TRUE]
        # if the new sire is NA, break the while loop and keep the previous id
        if (is.na(sire)) break
        # Otherwise change the id to be one level up, so from son to father and look for more sires
        current_id <- sire
      }
    }
    # assign current_id as the patriline for all photos of the current iteration
    data$patriline[data$fish_id == id] <- current_id
  }
  return(data)
}

# Load the pedigree data
ped_df <- data.table::fread('data/pedigree.csv') %>%
  as_tibble() %>%
  mutate(
    sex = str_sub(animal, 1, 1),
    grow_tank = ifelse(is.na(sire), 'source_pop', paste(sire, dam, date_of_birth, sep = '_'))
  )

# Get the carotenoid values for each photo
d <- data.table::fread('photo_database.csv') %>%
  as_tibble() %>%
  select(replicate, generation, fish_id, facing_direction, date, car_perc) %>%
  mutate(
    fish_id = toupper(fish_id),
    fish_idA = fish_id, fish_idX = fish_id,
    date = lubridate::as_date(as.character(date), format = '%Y%m%d')
  ) %>%
  # Add the patriline information to the dataset, as well as dam and sire id's
  add_patriline(ped_df) %>%
  left_join(ped_df, by = c('fish_id' = 'animal')) %>%
  # Finally, add the selection type
  left_join(selection %>% select(fish_id, selection), by = c('sire' = 'fish_id')) %>%
  left_join(selection %>% select(fish_id, selection2 = selection), by = c('fish_id')) %>%
  mutate(selection = coalesce(selection, selection2)) %>%
  select(-selection2)

# Make a new data.frame with representative females
representative_females <- d %>%
  filter(replicate == 'replicate_1', generation == 'gen_3') %>%
  group_by(dam, sire) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    fish_id = paste('representative_female', row_number(), sep = '_'),
    fish_idA = fish_id, fish_idX = fish_id,
    sex = 'F',
    car_perc = NA,
    facing_direction = NA,
    grow_tank = 'females'
  )

# add those females to our pedigree
ped_df_plus <- bind_rows(
  ped_df,
  select(representative_females, animal = fish_id, sire, dam, sex)
)

# Combine the data with our representative females
model_data <- bind_rows(d, representative_females) %>%
  # Shuffle for better threading in the model
  slice_sample(prop = 1)

## IT'S THE QUANT GEN ZONE BABY --------------------------------------------------------------------
# Make the additive relationship matrix for the autosomes
pedA <- makeA(as.data.frame(select(ped_df_plus, animal, dam, sire)))
# Make the additive relationship matrix for the X-chromosome
S <- makeS(as.data.frame(select(ped_df_plus, animal, dam, sire, sex)), 'M', 'ngdc', TRUE)
pedX <- S$S
colnames(pedX) <- rownames(pedX)

## ACTUAL FITTING ----------------------------------------------------------------------------------
library(brms)
m <- brm(
  car_perc | mi() ~ 1 + (1 | fish_id) + (1 | grow_tank) +
    (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) + (1 | patriline),
  data = model_data, data2 = list(A = pedA, X = pedX),
  backend = 'cmdstanr', chains = 6, cores = 6, threads = threading(2),
  warmup = 2000, iter = 6000, thin = 4,
  control = list(adapt_delta = 0.995, max_treedepth = 14),
  file = 'selection_decisions/models/rep1_gen3.rds'
)
#Total execution time: 7.2 days. (OUTDATED, more like 4 days?)
summary(m)
pp_check(m, nd = 50)

# only run interactively, so I can make sure the model fit well enough.
if (FALSE) {
  VC <- VarCorr(m, summary = FALSE) %>% map('sd') %>% map(`^`, 2)
  VC_df <- as.data.frame(VC) %>% setNames(names(VC)) %>%
    mutate(draw = row_number())
  VC_df_long <- pivot_longer(VC_df, -draw, names_to = 'param')
  VC_summ <- VC %>% map(posterior_summary) %>% map(as.data.frame) %>% bind_rows(.id = 'param')

  panels <- expand_grid(x = head(names(VC_df), -1), y = head(names(VC_df), -1))
  plots <- list()
  for (i in seq_len(nrow(panels))) {
    plots[[i]] <- ggplot(slice_sample(VC_df, n = 500), aes(.data[[panels$x[i]]], .data[[panels$y[i]]])) +
      geom_point()
  }
  cowplot::plot_grid(plotlist = plots, byrow = FALSE)

  library(ggdist)
  ggplot(VC_df_long, aes(param, value)) +
    ggdist::stat_ccdfinterval() +
    scale_x_discrete(
      limits = c('fish_idA' ,'fish_idX', 'patriline', 'grow_tank', 'fish_id', 'residual__'),
      labels = c('Autosome', 'X', 'Y', 'Grow\ntank', 'Environment\n(other)', 'Measurement\nerror')
    ) +
    labs(
      title = 'Where does the variance in carotenoid coloration come from?',
      subtitle = 'P, F1, F2 (33%)',
      x = 'Source of variance', y = 'Variance'
    ) +
    theme_minimal()

  cowplot::plot_grid(
    ggplot(VC_df, aes((fish_idA + fish_idX + patriline) / (fish_idA + fish_idX + patriline + fish_id + grow_tank))) +
      stat_halfeye() +
      xlim(0, 1) + labs(x = 'h2', y = 'Posterior density') + theme_minimal(),
    ggplot(VC_df, aes((patriline) / (fish_idA + fish_idX + patriline))) +
      stat_halfeye() +
      xlim(0, 1) + labs(x = 'Proportion Y-linked h2', y = 'Posterior density') + theme_minimal(),
    ggplot(VC_df, aes((fish_idX) / (fish_idA + fish_idX + patriline))) +
      stat_halfeye() +
      xlim(0, 1) + labs(x = 'Proportion X-linked h2', y = 'Posterior density') + theme_minimal(),
    ncol = 1
  )

  # predict female values
  fem <- posterior_epred(
    m,
    newdata = representative_females,
    re_formula = ~ (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X))
  )
  fem_sum <- fem %>% posterior_summary() %>% as_tibble() %>% bind_cols(representative_females)

  ggplot(fem_sum, aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = row_number(Estimate), color = selection)) +
    geom_pointrange()

  # Now let's draw females from each pair. Max 3 per pair.
  # We have to do more than 40 (even though we only need 40 in the lab), because we may not actually
  # have 3 females available from each pair. So we do 100, just in case.
  up_fem <- rep(NA, 100)
  for (i in seq_along(up_fem)) {
    draw <- sample(nrow(fem), 1)
    n <- 1
    while (is.na(up_fem[i])) {
      chosen <- tail(order(fem[draw, ]), n)[1]
      if (sum(up_fem == chosen, na.rm = TRUE) < 3 && fem_sum$selection[chosen] == 'up_selected') {
        up_fem[i] <- chosen
      } else {
        n <- n + 1
      }
    }
  }
  up_selected_females <- fem_sum %>%
    slice(up_fem) %>%
    mutate(order = row_number()) %>%
    select(order, dam, sire) %>%
    arrange(dam, order)

  # Now same for down-selected females:
  down_fem <- rep(NA, 100)
  for (i in seq_along(down_fem)) {
    draw <- sample(nrow(fem), 1)
    n <- 1
    while (is.na(down_fem[i])) {
      chosen <- tail(order(-fem[draw, ]), n)[1]
      if (sum(down_fem == chosen, na.rm = TRUE) < 3 && fem_sum$selection[chosen] == 'down_selected') {
        down_fem[i] <- chosen
      } else {
        n <- n + 1
      }
      if (n > ncol(fem)) break
    }
  }
  down_selected_females <- fem_sum %>%
    slice(down_fem) %>%
    mutate(order = row_number()) %>%
    select(order, dam, sire) %>%
    arrange(dam, order)

  data.table::fwrite(up_selected_females, 'selection_decisions/fem_rep1_gen3_up.csv')
  data.table::fwrite(down_selected_females, 'selection_decisions/fem_rep1_gen3_down.csv')
}
