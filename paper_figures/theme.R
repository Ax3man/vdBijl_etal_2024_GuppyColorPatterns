library(ggplot2)

my_theme <- theme_classic(
  base_size = 7,
  base_family = 'Arial'
) +
  theme(
    legend.key.size = unit(0.4, 'lines'),
    legend.key.height = unit(0.4, 'lines'),
    strip.text = element_text(size = 7),
    plot.margin = margin(0, 3, 0, 3),
    panel.spacing = unit(5, 'points')
  )
theme_set(my_theme)
