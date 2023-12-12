library(tidyverse)
library(imager)
source('paper_figures/theme.R')
source('selection_decisions/compile_decisions.R')
eb <- element_blank()

photo_database <- data.table::fread('photo_database.csv')

mel_out <- read_rds('ornament_analysis/mel_ornaments.rds') %>%
  mutate(
    ornament = ornament %>% parse_number() %>% as.factor()#,
    #ornament_pic_label = labs[ornament]
  ) %>%
  left_join(select(photo_database, fish_id, unique_id), 'unique_id') %>%
  mutate(fish_id = toupper(fish_id)) %>%
  left_join(select(selection, fish_id, replicate, generation, selection), 'fish_id')

Mode <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(table(match(x, ux)))]
}

pd3 <- mel_out %>%
  group_by(ornament, generation, selection, fish_id) %>%
  summarise(present = Mode(present),
            ornament_size = mean(pixels_present),
            .groups = 'drop_last'
  ) %>%
  mutate(selection = ifelse(generation == 'P', 'stock', selection))

P_mel_orn_pres <- ggplot(pd3, aes(generation, present, fill = selection)) +
  geom_bar(stat = 'summary', fun = 'mean', position = 'dodge', width = 0.7) +
  stat_summary(
    fun.data = 'mean_cl_boot', geom = 'linerange',
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_manual(
    values = c('grey60', 'navy', 'firebrick'),
    labels = c('Stock', 'Down-selected', 'Up-selected'),
    limits = c('stock', 'down_selected', 'up_selected')
  ) +
  scale_x_discrete(name = NULL, labels = expression(P, F[1], F[2], F[3])) +
  ylim(0, 1) +
  facet_wrap(ornament ~ ., nrow = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = 'Generation', y = 'Incidence', fill = NULL) +
  theme(
    strip.background = eb,
    legend.position = 'bottom',
    axis.title.x = eb, axis.ticks.x = eb, axis.text.x = eb,
    strip.text = eb, strip.background.x = eb
  )

P_mel_orn_size <- pd3 %>%
  # NOTE!!! only counting ornament sizes when it is present. NOTE!!!
  filter(present == 1) %>%
  ggplot(aes(generation, ornament_size, col = selection, group = selection)) +
  geom_line(
    stat = 'summary', fun = 'mean',
    position = position_dodge(width = 0.7), show.legend = FALSE
  ) +
  stat_summary(
    fun.data = 'mean_cl_boot', position = position_dodge(width = 0.7),
    fatten = .5, size = 1.5
  ) +
  scale_color_manual(
    values = c('navy', 'grey60', 'firebrick'),
    labels = c('Down-selected', 'Stock', 'Up-selected'),
    guide = 'none'
  ) +
  scale_x_discrete(name = NULL, labels = expression(P, F[1], F[2], F[3])) +
  facet_wrap(ornament ~ ., nrow = 1) +#, scales = 'free') +
  labs(x = 'Generation', color = 'Selection', y = 'Size (pixels)') +
  expand_limits(y = 0) +
  theme(
    strip.background = element_blank(),
    strip.text = eb, strip.background.x = eb
  )
