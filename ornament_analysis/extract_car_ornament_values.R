library(tidyverse)
library(imager)
library(furrr)

template <- load.image('ornament_analysis/ornament_images/template_new.png')
car_ornaments <- list.files('ornament_analysis/ornament_images', 'car_._new', full.names = TRUE) %>%
  map(~load.image(.x) %>% channel(4) %>% threshold(thr = 0.5) %>% as.pixset()) %>%
  setNames(paste0('car_', seq_along(.)))

plot(template, rescale = FALSE, axes = FALSE)
walk(car_ornaments, highlight)

get_ornament_values <- function(image_file, ornaments) {
  # image_file <- sample(car_ims, 1); ornaments <- car_ornaments
  img <- load.image(image_file) %>% channel(4) %>% threshold(0.5) %>% as.cimg()
  # plot(img); walk(ornaments, highlight)
  counts <- map_dbl(ornaments, ~sum(img[.x]))
  fractions <- map_dbl(ornaments, ~mean(img[.x]))
  return(data.frame(
    ornament = names(ornaments),
    pixels_present = counts,
    fraction_present = fractions)
  )
}

car_ims <- list.files('data/carotenoid_coloration_warped/', full.names = TRUE, recursive = TRUE) %>%
  # remove dumb Apple icon file
  str_subset('Icon', negate = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

plan(multisession, workers = 16)
car_out <- future_map_dfr(
  car_ims, get_ornament_values, car_ornaments,
  .id = 'unique_id', .progress = TRUE, .options = furrr_options(seed = NULL)
) %>%
  mutate(
    present = as.integer(pixels_present > 0),
    present_10 = as.integer(fraction_present > 0.1)
  ) %>%
  as_tibble()
plan(sequential)

write_rds(car_out, 'ornament_analysis/car_ornaments.rds')

# for interactive use
if (FALSE) {
  car_out <- read_rds('ornament_analysis/car_ornaments.rds')

  library(imager)
  library(tidyverse)
  library(ggtext)

  # do image prep for fancy figure labels
  bg <- (1 - (load.image('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
                channel(4))) %>% add.color()
  recolor_image <- function(im, old_cols, new_cols) {
    for (i in seq_along(old_cols)) {
      im <- colorise(
        im,
        R(im) == old_cols[[i]][1] & G(im) == old_cols[[i]][2] & B(im) == old_cols[[i]][3],
        new_cols[[i]]
      )
    }
    im
  }
  car_pixsets <- list.files('ornament_analysis/ornament_images', 'car_._new\\.png', full.names = TRUE) %>%
    map(\(.x) {
      out <- load.image(.x) %>% channel(4) %>% threshold(0.5) %>% as.cimg() %>% add.color()
      G(out) <- 0; B(out) <- 0
      as.pixset(out)
    }) %>%
    setNames(paste0('car_', seq_along(.)))
  all_car_ornaments <- reduce(car_pixsets, `+`, .init = bg)
  all_car_ornaments <- recolor_image(
    all_car_ornaments,
    list(c(1, 0, 0), c(0, 0, 0)), list(c(0.8, 0.5, 0), 'grey70')
  )
  plot(all_car_ornaments)
  save.image(all_car_ornaments, 'ornament_analysis/visualization/car_ornaments_abstract.png')
  car_ornaments <- map(car_pixsets, ~bg + .x)
  car_ornaments <- map(
    car_ornaments, recolor_image,
    list(c(1, 0, 0), c(0, 0, 0)), list(c(0.8, 0.5, 0), 'grey70')
  )
  fn <- paste0('ornament_analysis/ornament_minis/', names(car_ornaments), '.png')
  walk2(car_ornaments, fn, ~save.image(.x, .y))
  labs <- paste0("<img src='", fn, "'width='75' />") %>% setNames(paste0('car_', seq_along(fn)))
  labs_large <- paste0("<img src='", fn, "'width='125' />") %>% setNames(paste0('car_', seq_along(fn)))

  # summary of the ornaments
  car_out %>% group_by(ornament) %>% summarise(across(pixels_present:present, mean))

  pd <- car_out %>%
    mutate(
      ornament = ornament %>% parse_number() %>% as.factor(),
      ornament_pic_label = labs[ornament]
    )

  ggplot(pd, aes(pixels_present)) +
    #geom_density(fill = 'grey60') +
    geom_histogram(bins = 100) +
    facet_wrap(ornament_pic_label ~ ., scales = 'free', ncol = 3) +
    #geom_vline(xintercept = 100, lty = 2) +
    scale_x_log10() +
    coord_cartesian(expand = FALSE) +
    labs(x = 'Pixels with ornament', y = 'Density') +
    theme_classic() +
    theme(strip.text = element_markdown(), strip.background = element_blank())
  #ggsave('ornament_analysis/visualization/cut_off_decision.png', w = 8, h = 5)

  panels <- expand_grid(x = names(car_out_wide)[-1], y = names(car_out_wide)[-1])
  plots <- list()
  for (i in seq_len(nrow(panels))) {
    plots[[i]] <- ggplot(
      car_out_wide, aes(x = .data[[panels$x[i]]], fill = factor(.data[[panels$y[i]]]))
    ) + geom_bar(show.legend = FALSE, position = 'fill') + theme_void()
  }
  cowplot::plot_grid(plotlist = plots, byrow = FALSE)

  # for (i in seq_len(nrow(panels))) {
  #   chisq.test(car_out_wide[[panels$x[[i]]]], car_out_wide[[panels$y[i]]]) %>% print()
  # }

  source('selection_decisions/compile_decisions.R')
  photo_database <- data.table::fread('photo_database.csv')

  pd2 <- pd %>%
    left_join(select(photo_database, fish_id, unique_id), 'unique_id') %>%
    mutate(fish_id = toupper(fish_id)) %>%
    left_join(select(selection, fish_id, replicate, generation, selection), 'fish_id')

  Mode <- function(x, na.rm = FALSE) {
    if (na.rm) x <- x[!is.na(x)]
    ux <- unique(x)
    ux[which.max(table(match(x, ux)))]
  }

  pd3 <- pd2 %>%
    group_by(ornament_pic_label, generation, selection, fish_id) %>%
    summarise(present = Mode(present),
      ornament_size = mean(pixels_present),
      .groups = 'drop_last'
    ) %>%
    mutate(selection = ifelse(generation == 'P', 'stock', selection))

  ggplot(pd3, aes(generation, present, fill = selection)) +
    geom_bar(stat = 'summary', fun = 'mean', position = 'dodge') +
    stat_summary(
      fun.data = 'mean_cl_boot', geom = 'errorbar',
      position = position_dodge(width = 0.9), width = 0.2
    ) +
    scale_fill_manual(
      values = c('navy', 'grey60', 'firebrick'),
      labels = c('Down-selected', 'Stock', 'Up-selected')) +
    ylim(0, 1) +
    facet_wrap(ornament_pic_label ~ ., ncol = 4, scales = 'free') +
    coord_cartesian(expand = FALSE) +
    labs(x = 'Generation', fill = 'Selection', y = 'Fraction of population with ornament') +
    theme_classic() +
    theme(
      strip.text = element_markdown(), strip.background = element_blank(),
      legend.position = c(7/8, 1/4)
    )
  ggsave('ornament_analysis/visualization/ornament_response_to_selection.png', w = 8, h = 4, scale = 1.2)

  pd3 %>%
    # NOTE!!! only counting ornament sizes when it is present. NOTE!!!
    filter(present == 1) %>%
    ggplot(aes(generation, ornament_size, col = selection, group = selection)) +
    geom_line(
      stat = 'summary', fun = 'mean',
      position = position_dodge(width = 0.5), show.legend = FALSE
    ) +
    stat_summary(
      fun.data = 'mean_cl_boot', position = position_dodge(width = 0.5),
      fatten = 1, size = 1.5
    ) +
    scale_color_manual(
      values = c('navy', 'grey60', 'firebrick'),
      labels = c('Down-selected', 'Stock', 'Up-selected'),
      guide = guide_legend(override.aes = list(size = 0.5))
    ) +
    facet_wrap(ornament_pic_label ~ ., ncol = 4, scales = 'free') +
    labs(x = 'Generation', color = 'Selection', y = 'Ornament size (pixels)') +
    expand_limits(y = 0) +
    theme_classic() +
    theme(
      strip.text = element_markdown(), strip.background = element_blank(),
      legend.position = c(7/8, 1/4)
    )
  ggsave('ornament_analysis/visualization/ornament_response_to_selection2.png', w = 8, h = 4, scale = 1.2)
}
