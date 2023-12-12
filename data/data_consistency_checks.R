library(tidyverse)

# These sources of data should be consistent:
# - The photo database, derived from tank labels of individual males after phenotyping, transcribed
#   from photos of those labels, and entered as folder names.
# - The pedigree, derived from the same photos (but entered into spreadsheet).
# - The breeding data (data/breeding.csv), derived from tank labels of the breeding pairs, written
#   in the lab book, and then transcribed into a spreadsheet.
#
# Mistakes like data entry typos should show up as differences between these three data sources.
# Other data was constructed from these three sources.

# Only for add_patriline()
source('quant_gen/prepare_pedigrees.R')
source('selection_decisions/compile_decisions.R')

photo_database <- data.table::fread('photo_database.csv') %>%
  mutate(fish_id = toupper(fish_id))
pedigree <- data.table::fread('data/pedigree.csv') %>%
  as_tibble() %>%
  mutate(sex = str_sub(animal, 1, 1)) %>%
  { add_patriline(rename(., fish_id = animal), .) }
breeding <- data.table::fread('data/breeding_data.csv')

# Test 1: Every male in the pedigree should have photos
males_in_ped <- filter(pedigree, sex == 'M') %>% pull(fish_id)
setdiff(males_in_ped, photo_database$fish_id)
# M1-439 is not in the photo_database, he was photographed but a lot of tail damage and was
# therefore excluded and not analyzed.

# Test 2: Every male with photos should be in the pedigree
setdiff(photo_database$fish_id, males_in_ped)
# character(0) is expected

# Test 3: Every sire in the pedigree should have photos
setdiff(na.omit(pedigree$sire), photo_database$fish_id)
# character(0) is expected

# Test 4: The parents of every male in the pedigree should have breeding data
breeding_pairs1 <- unique(paste(breeding$sire, breeding$dam))
breeding_pairs2 <- unique(paste(pedigree$sire, pedigree$dam))
setdiff(breeding_pairs2, breeding_pairs1)
# NA NA is expected

# Note that there are breeding pairs which occur in the breeding data, but are not in the pedigree
# as dam/sire. This is because only breeding females occur in the pedigree, so a pair with only
# daughters will have no offspring in the pedigree if their daughters did not breed.
# This is expected.
setdiff(breeding_pairs1, breeding_pairs2)

# Test 5: Same as 4, but including DOB.
breeding_pairs3 <- breeding %>% with(
  unique(paste(sire, dam, str_replace_all(date, '\\/', '\\-')))
)
breeding_pairs4 <- pedigree %>%
  separate(fish_id, c('chunk1', 'chunk2'), sep = '-') %>%
  filter(!str_starts(chunk1, 'F'), parse_number(chunk2) > 100) %>%
  with(unique(paste(sire, dam, date_of_birth)))
setdiff(breeding_pairs4, breeding_pairs3)
# NA NA NA is expected
# M1-38 F1-80 2020-11-02: Breeding data on this brood is missing. Male M1-298 was likely part of
# the same brood as M1-199 till M1-202, which were born a few days earlier, but recorded separately
# by accident.

# Test 6, all males in a patriline have matching replicate numbers
pedigree %>%
  group_by(patriline) %>%
  summarise(check = n_distinct(substr(fish_id, 2, 2)) == 1) %>%
  filter(!check)
# NA FALSE is expected

# Test 7, all fish should be from the same selection line as their father (up / down)
select(pedigree, fish_id, sire) %>%
  drop_na() %>%
  left_join(select(selection, fish_id, fish_line = selection2), 'fish_id') %>%
  left_join(select(selection, sire = fish_id, sire_line = selection2), 'sire') %>%
  filter(fish_line != sire_line)
# 0 row tibble is expected

