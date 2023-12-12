library(imager)
library(tidyverse)

photo_database <- data.table::fread('photo_database.csv')

# get only those photos for which I ran both methods
double_tested <- photo_database %>%
  group_by(fish_id) %>%
  filter(n_distinct(photo_method) > 1) %>%
  arrange(unique_id) %>%
  select(photo_method, fish_id, unique_id, facing_direction) %>%
  # make sure we have three photos for each fish, for each method
  group_by(fish_id, photo_method) %>%
  sample_n(3) %>%
  ungroup()

# Load the images from the calibration folder. Note that the matching images in the actual data
# have already have corrected with the LUT calculated here.
warped_ims <- with(
  double_tested,
  paste0('calibration/extracted_fish_warped/', unique_id, '.png')
) %>% map(load.image)

# plot some examples, you should be able to see that method 1 (left) is more greenish, while method
# 2 (right) is more pinkish
sp <- split(warped_ims, double_tested$fish_id)
opa <- par(mfcol = c(3, 2), oma = c(0,0,0,0), mar = c(0,0,0,0))
for (i in 1:5) walk(sp[[i]], plot, axes = FALSE, rescale = FALSE)
par(opa)

df <- warped_ims %>%
  setNames(double_tested$unique_id) %>%
  map_dfr(as.data.frame, wide = 'c', .id = 'unique_id') %>%
  rename(R = c.1, G = c.2, B = c.3, alpha = c.4) %>%
  cbind(., convertColor(as.matrix(.[4:6]), from = 'sRGB', to = 'Lab')) %>%
  as_tibble() %>%
  left_join(double_tested, 'unique_id')

# now calculate reference distributions, i.e. what the quantiles look like
LUT <- df %>%
  filter(alpha > 0.5) %>%
  group_by(photo_method) %>%
  summarise(across(c(L, a, b), quantile, seq(0, 1, l = 2048)), .groups = 'keep') %>%
  mutate(quant = row_number()) %>%
  pivot_wider(names_from = photo_method, values_from = L:b, names_sep = '')

# Save these references to disk, we'll need them lots. I'll save both .csv (for compatibility and
# archiving) and .rds (for quick access)
data.table::fwrite(LUT, 'calibration/lookuptable.csv')
write_rds(LUT, 'calibration/lookuptable.rds')

cowplot::plot_grid(
  ggplot(LUT, aes(quant/20, b1)) + geom_point() + theme_minimal() +
    labs(x = 'quantile', y = 'b (old)') + ylim(min(c(LUT$b1, LUT$b2)), max(c(LUT$b1, LUT$b2))),
  ggplot(LUT, aes(quant/20, b2)) + geom_point() + theme_minimal() +
    labs(x = 'quantile', y = 'b (new)') + ylim(min(c(LUT$b1, LUT$b2)), max(c(LUT$b1, LUT$b2))),
  ggplot(LUT, aes(b2, b1)) + geom_point() + theme_minimal() + labs(x = 'b (new)', y = 'b (old)') +
    geom_abline(slope = 1),
  nrow = 1
)
ggsave('calibration/LUT_example_plot.png', width = 12, height = 5, bg = 'white')

# These can now be used as look up tables
lookup <- function(x, axis) {
  LUT_limited <- LUT[2:(nrow(LUT) - 1), ]
  breaks <- c(-Inf, LUT_limited[[paste0(axis, 2)]], Inf)
  ref <- LUT_limited[[paste0(axis, 1)]]
  ref <- c(ref[1], ((ref + lag(ref)) / 2)[-1], tail(ref, 1))
  ref[cut(x, breaks, labels = FALSE)]
}

# apply the LUT, and the greenish/pinkish difference is no longer visible:
df_new <- df %>%
  mutate(
    L_new = if_else(photo_method == 1 | alpha < 0.5, L, lookup(L, 'L')),
    a_new = if_else(photo_method == 1 | alpha < 0.5, a, lookup(a, 'a')),
    b_new = if_else(photo_method == 1 | alpha < 0.5, b, lookup(b, 'b')),
  ) %>%
  cbind(
    .,
    setNames(as.data.frame(
      convertColor(as.matrix(.[c('L_new', 'a_new', 'b_new')]), from = 'Lab', to = 'sRGB'),
    ), c('R_new', 'G_new', 'B_new'))
  )
translated <- split(df_new, df_new$unique_id) %>%
  map(~.x[c('x', 'y', 'R_new', 'G_new', 'B_new')] %>%
        pivot_longer(c('R_new', 'G_new', 'B_new'), 'cc') %>%
        mutate(cc = setNames(1:3, c('R_new', 'G_new', 'B_new'))[cc]) %>%
        as.cimg(dims = dim(warped_ims[[1]]))
  )
sp <- split(translated, double_tested$fish_id)
pdf('calibration/calibrated_examples.pdf', width = 5, height = 2)
par(mfcol = c(3, 2), oma = c(0,0,0,0), mar = c(0,0,0,0))
for (i in seq_along(sp)) walk(sp[[i]], plot, axes = FALSE, rescale = FALSE)
dev.off()
