library(tidyverse)
library(nadiv)

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
  # Finally, add the selection type, only for rep1 gen2 right now.
  left_join(
    data.table::fread('selection_decisions/finished/rep1_parentalgen.csv') %>%
      select(fish_id, selection) %>% mutate(fish_id = toupper(fish_id)),
    by = c('sire' = 'fish_id')
  )

# Make a new data.frame with representative females
representative_females <- d %>%
  filter(replicate == 'replicate_1', generation == 'gen_2') %>%
  group_by(dam, sire) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    fish_id = paste('representative_female', row_number(), sep = '_'),
    fish_idA = fish_id, fish_idX = fish_id,
    sex = 'F',
    car_perc = NA,
    facing_direction = NA
  )

# add those females to our pedigree
ped_df_plus <- bind_rows(
  ped_df,
  select(representative_females, animal = fish_id, sire, dam, sex)
)

# Select the subset of the data to focus on right now, and combine with representative females
model_data <- d %>%
  filter(replicate == 'replicate_1') %>%
  bind_rows(representative_females)

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
  car_perc | mi() ~ 1 + generation + (1 | fish_id) + (1 | grow_tank) +
    (1 | gr(fish_idA, cov = A)) + (1 | gr(fish_idX, cov = X)) + (1 | patriline),
  data = model_data, data2 = list(A = pedA, X = pedX),
  backend = 'cmdstanr', chains = 12, cores = 12, warmup = 2000, iter = 6000,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  file = 'quant_gen/car_heritability_model.rds'
)
# Mean chain execution time: 16047.1 seconds.
# Total execution time: 17946.8 seconds.
summary(m)
pp_check(m, ns = 50)

VC <- VarCorr(m, summary = FALSE) %>% map('sd') %>% map(`^`, 2)
VC_df <-  VC %>% map(posterior_summary) %>% map(as.data.frame) %>%  bind_rows(.id = 'param')
ggplot(VC_df, aes(param, Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_col(position = 'stack') +
  geom_errorbar(width = 0.2) +
  scale_x_discrete(
    limits = c('fish_idA' ,'fish_idX', 'patriline', 'grow_tank', 'fish_id', 'residual__'),
    labels = c('Autosome', 'X', 'Y', 'Grow\ntank', 'Environment\n(other)', 'Measurement\nerror')
  ) +
  labs(
    title = 'Where does the variance in carotenoid coloration come from?',
    subtitle = 'Very preliminary estimates, only data from replicate 1',
    x = 'Source of variance', y = 'Variance'
  ) +
  theme_minimal()

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
  select(order, dam, sire)

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
  select(order, dam, sire)

data.table::fwrite(up_selected_females, 'selection_decisions/fem_rep1_gen2_up.csv')
data.table::fwrite(down_selected_females, 'selection_decisions/fem_rep1_gen2_down.csv')
