library(VariantAnnotation)
library(irlba)
library(ggtext)

vcf_file <- "sequencing/gwas/filtered.vcf.gz"

g <- data.table::fread('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/car_PIE.csv') %>%
  group_by(chr, start, end, alt) %>%
  slice_max(AIC_weight) %>%
  filter(AIC_weight > 0.8, source == 'Y') %>%
  ungroup()

if (FALSE) { # only run interactively, otherwise load from disk
  regions <- GRanges(
    seqnames = g$chr,
    ranges = IRanges(start = g$start, end = g$end)
  )
  geno <- readVcfAsVRanges(vcf_file, param = ScanVcfParam(which = regions)) %>%
    as.data.frame() %>%
    filter(
      !(sampleNames %in% c('NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355'))
    )

  geno_mat <- geno %>%
    semi_join(g, join_by(seqnames == chr, start, end, alt)) %>%
    mutate(
      dose = case_when(
        GT == '0/0' ~ 0,
        GT == '0/1' ~ 1,
        GT == '1/1' ~ 2
      ),
      snp_id = paste(seqnames, start, end, alt, sep = '_')
    ) %>%
    dplyr::select(sampleNames, snp_id, dose) %>%
    distinct() %>%
    pivot_wider(names_from = snp_id, values_from = dose) %>%
    as.data.frame() %>%
    column_to_rownames('sampleNames') %>%
    as.matrix()

  pca <- prcomp_irlba(geno_mat, n = 20, center = TRUE, scale. = TRUE)
  rownames(pca$x) <- rownames(geno_mat)
  write_rds(pca, 'sequencing/gwas/Y_haplotype_PCA.rds', compress = 'gz')
}
pca <- read_rds('sequencing/gwas/Y_haplotype_PCA.rds')
summary(pca)

# scree plot
scree <- summary(pca)$importance %>%
  as.data.frame() %>%
  dplyr::slice(2) %>%
  pivot_longer(everything()) %>%
  mutate(name = parse_number(name)) %>%
  ggplot(aes(name, value)) + geom_line() + geom_point() + theme_classic() +
  labs(x = 'Principal component', y = 'Proportion of variance')

# hclust
dist_matrix <- dist(pca$x) # distance matrix from the first two PCs
hc <- hclust(dist_matrix)
#plot(hc) # to visualize the dendrogram

##
source('quant_gen/prepare_pedigrees.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('selection_decisions/compile_decisions.R')
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
orn <- read_rds('ornament_analysis/car_ornaments.rds') %>%
  left_join(data.table::fread('photo_database.csv') %>% dplyr::select(unique_id, fish_id)) %>%
  group_by(fish_id, ornament) %>%
  summarise(present = Mode(present)) %>%
  pivot_wider(names_from = ornament, values_from = present)

ind_df <- as.data.frame(pca$x) %>% rownames_to_column('sample_name') %>% as_tibble() %>%
  left_join(get_sampling_structure(), join_by(sample_name)) %>%
  left_join(mutate(selection, fish_id = tolower(fish_id)), join_by(fish_id)) %>%
  left_join(orn, join_by(fish_id)) %>%
  add_patriline(ped_df)

ind_df$hclust_groups <- factor(cutree(hc, k = 4))

ind_df %>%
  mutate(Y_haplogroup = factor(hclust_groups, levels = 1:4, labels = paste0('Y', 1:4))) %>%
  dplyr::select(sample_name, fish_id, Y_haplogroup) %>%
  data.table::fwrite('sequencing/Y_haplogroups.csv')

ind_df %>%
  pivot_longer(car_1:car_7, names_to = 'ornament', values_to = 'presence') %>%
  mutate(ornament_label = glue::glue(
    "<img src='ornament_analysis/ornament_figure_labels/{ornament}.png' width='60' />"
  )) %>%
  ggplot(aes(hclust_groups, fill = factor(presence, 0:1, c('absent', 'present')))) +
  geom_bar() +
  scale_y_continuous(expand = c(0, 0), name = 'Count') +
  scale_fill_manual(values = c('grey20', 'darkorange'), name = NULL) +
  facet_grid(cols = vars(ornament_label)) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_markdown()) +
  labs(x = 'Assigned Y-haplogroup')

ggsave('paper_figures_supplement/Y_haplogroup_vs_ornaments.png', w = 8, h = 3.5)

biplot <- ggplot(ind_df, aes(PC1, PC2, color = hclust_groups)) +
  geom_point() + theme_classic() + coord_fixed() +
  scale_color_brewer(palette = "Set1", labels = paste0('Y', 1:4)) +
  labs(
    x = 'PC1 (29%)', y = 'PC2 (23%)',
    color = 'Assigned\nY-haplogroup',
    #subtitle = 'PCA on SNPs associated w/ orange pattern and Y-like inheritance'
  ) +
  theme(legend.position = 'top')

patchwork::wrap_plots(biplot, scree, ncol = 2) + patchwork::plot_annotation(tag_levels = 'A')
ggsave('paper_figures_supplement/Y_haplogroups.png', w = 8, h = 3.5)

