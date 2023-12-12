library(tidyverse)

all_photos <- tibble(
  file = list.files('photos', pattern = '.JPG*', full.names = TRUE, recursive = TRUE)
)

# first we build an index of all the raw images. We can use the file structure to get fish ids.
## Note: the warnings with created NAs are the color cards and are removed
raw_photos <- all_photos %>%
  filter(str_detect(file, 'raw')) %>%
  separate(
    file,
    c('tmp', 'replicate', 'generation', 'image_type', 'date', 'fish_id', 'photo_number'),
    '/',
    fill = 'right',
    extra = 'drop'
  ) %>%
  select(-tmp, -image_type) %>%
  drop_na(photo_number)

corrected_photos <- all_photos %>%
  filter(str_detect(file, 'color_corrected')) %>%
  separate(file, c('tmp', 'replicate', 'generation', 'image_type', 'photo_number'), '/', remove = FALSE) %>%
  separate(photo_number, c('date', 'photo_number'), '-', extra = 'drop') %>%
  select(-tmp, -image_type)

photos <- inner_join(
  raw_photos,
  corrected_photos,
  by = c("replicate", "generation", "date", "photo_number")
) %>%
  mutate(unique_id = paste(date, photo_number, sep = '_') %>% str_remove('.JPG')) %>%
  relocate(file, .after = last_col())

if (nrow(photos) != nrow(raw_photos)) {
  message('Failed to find color corrected photos for all raw photos.')
}
# View(anti_join(raw_photos, photos))
# View(anti_join(photos, raw_photos))

rm(all_photos, raw_photos, corrected_photos)
