library(tidyverse)
library(uwot) # for UMAP
source('paper_figures/theme.R')
eb <- element_blank()

embedding <- read_rds(
  'dimension_reduction/triplet_loss_encoders/embeddings/mel_model_ped_fullcolor_comparison/embed_dim_5.rds'
) %>%
  magrittr::set_colnames(paste0('V', 1:ncol(.))) %>%
  as.data.frame() %>%
  rownames_to_column('unique_id') %>%
  as_tibble()

umap_model <- umap(embedding[-1], ret_model = TRUE, fast_sgd = TRUE, spread = 2) #, n_neighbors = 50, spread = 5)

#pairs(embedding[c('V1', 'V2', 'V3')])

embedding <- bind_cols(
  embedding,
  umap_model$embedding %>% as.data.frame() %>% setNames(paste0('UMAP', 1:2))
)

# plot(embedding$UMAP1, embedding$UMAP2)

pd <- data.table::fread('photo_database.csv') %>% as_tibble() %>%
  select(fish_id, unique_id) %>%
  mutate(fish_id = toupper(fish_id))

ornaments <- read_rds('ornament_analysis/mel_ornaments.rds')

P_mel_orn_emb <- inner_join(ornaments, embedding, 'unique_id') %>%
  mutate(pixels_present = na_if(pixels_present, 0)) %>%
  # shuffle for random point overlap
  slice_sample(prop = 1) %>%
  ggplot(aes(UMAP1, UMAP2, color = pixels_present / max(pixels_present, na.rm = TRUE))) +
  geom_point(size = .2, stroke = 0) +
  # sqrt trans color scale?
  scale_color_viridis_c(
    na.value = 'grey80', trans = 'sqrt', breaks = c(.01, .1, .3, .6, 1),
    guide = guide_colorbar(barheight = unit(4, 'lines'))
  ) + #limits = c(0, 1),
  #scale_color_viridis_c(na.value = 'grey80', limits = c(0, 1)) +
  coord_equal() +
  labs(color = 'Relative size') +
  facet_grid(~ornament) +
  theme(
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.title = element_text(size = rel(0.8)),
    strip.background = eb, strip.text = eb,
  )

