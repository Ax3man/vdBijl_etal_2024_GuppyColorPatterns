library(VariantAnnotation)
library(imager)
library(tidyverse)

source('sequencing/gwas/peak_viz/peak_viz_tools.R')

trait <- 'car_PIE'
pval_column <- 'p_SHet'
fdr.level <- 0.01
min_dist <- 1e5 / 2

gwas_selected <- data.frame(
  chr = c('NC_024344.1', 'NC_024333.1', 'NC_024333.1'),
  ps = c(26587645, 593580, 20103960)
)

# Load genotypes -----------------------------------------------------------------------------------
vcf_file <- 'sequencing/gwas/filtered.vcf.gz'
if (!file.exists(paste0(vcf_file, '.tbi'))) indexVcf(vcf_file)

regions <- GRanges(seqnames = gwas_selected$chr, ranges = IRanges(start = gwas_selected$ps, width = 1))
genotypes <- readVcfAsVRanges(vcf_file, param = ScanVcfParam(which = regions)) %>%
  as_tibble() %>%
  mutate(snp_id = paste(seqnames, start, end, alt, sep = '_')) %>%
  # remove dropped individual
  filter(
    !(sampleNames %in% c('NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355'))
  )

# prepare a simple genotype file and save to disk. This is necessary for GMMAT to function
genotype_flat <- genotypes %>%
  mutate(SNP = paste(seqnames, start, end, alt, sep = '/')) %>%
  dplyr::select(SNP, sampleNames, GT) %>%
  mutate(GT = case_when(GT == '0/0' ~ 0L, GT == '0/1' ~ 1L, GT == '1/1' ~ 2L, TRUE ~ NA_integer_)) %>%
  pivot_wider(names_from = sampleNames, values_from = GT)
genotype_flat_file <- tempfile(fileext = '.csv')
data.table::fwrite(genotype_flat, genotype_flat_file, row.names = FALSE, col.names = FALSE)

SAMPLE_ORDER <- colnames(genotype_flat)[-1]

# Load GRM -----------------------------------------------------------------------------------------
sample_names <- readVcfAsVRanges(
  vcf_file,
  param = ScanVcfParam(which = GRanges("NC_024342.1", IRanges(26272132)))
) %>% as.data.frame() %>% dplyr::pull(sampleNames) %>% as.character() %>%
  # remove dropped individual
  str_subset('NS.2125.002.IDT_i7_111---IDT_i5_111.280', negate = TRUE) %>%
  str_subset('NS.2145.001.IDT_i7_89---IDT_i5_89.355', negate = TRUE)

GRM <- get_GRM()
GRM <- GRM[SAMPLE_ORDER, SAMPLE_ORDER]

# Load phenotypes ----------------------------------------------------------------------------------
sampling <- get_sampling_structure()

img_files <- list.files('data/carotenoid_coloration_warped', recursive = TRUE, full.names = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

# first, match the sequenced individuals with all their color pattern images
phenotypes <- data.table::fread('photo_database.csv') %>%
  dplyr::select(fish_id, facing_direction, unique_id) %>%
  inner_join(sampling, join_by(fish_id)) %>%
  # load those images
  mutate(image = map(
    img_files[unique_id],
    \(p) load.image(p) %>% channel(4) %>% resize_halfXY() %>% threshold(thr = 0.5)
  )) %>%
  # reduce to 1 measurment per individual, by averaging the images first per side, then per fish.
  group_by(sample_name, fish_id, facing_direction) %>%
  summarise(image = list(average(image)), .groups = 'drop_last') %>%
  summarise(image = list(average(image)), .groups = 'drop') %>%
  mutate(image = map(image, as.data.frame)) %>%
  unnest(cols = image) %>%
  # since we have averaged a lot of binary values, we need to round again to binary values
  mutate(value = round(value)) %>%
  # finally, drop pixels where < 5% of sequenced males have color
  group_by(x, y) %>%
  filter(mean(value) > 0.05)
cat('Analyzing', n_groups(phenotypes), 'pixels for', nrow(genotype_flat), 'SNPs...\n')

# Run pixelwise SNP effect estimation ---------------------ht-----------------------------------------
get_snp_effects <- function(pheno, method = 'score') {
  #pheno <- group_split(phenotypes)[[1]]
  require(GMMAT)

  pheno <- pheno[match(pheno$sample_name, SAMPLE_ORDER), ]

  if (method == 'score') {
    # The score test is really fast, but can give some anti-conservative p-values
    null <- glmmkin(
      fixed = value ~ 1,
      data = pheno,
      kins = GRM,
      id = 'sample_name',
      family = gaussian()
    )
    #h2 <- null$theta['kins1'] / sum(null$theta)
    tmp_out_file <- tempfile()
    glmm.score(
      null,
      infile = genotype_flat_file,
      outfile = tmp_out_file,
      infile.sep = ',',
      select = seq_along(SAMPLE_ORDER)
    )
    out <- data.table::fread(tmp_out_file) %>%
      separate(SNP, into = c('seqnames', 'start', 'end', 'alt'), sep = '/', remove = FALSE)
  } else {
    # the Wald test is better, and gives us effect sizes (Beta + se), but is really slow:
    wald <- glmm.wald(
      fixed = value ~ 1,
      data = pheno,
      kins = GRM,
      id = 'sample_name',
      family = binomial(link = 'logit'),
      infile = genotype_flat_file,
      snps = genotype_flat$SNP,
      infile.sep = ','
    )
    out <- wald %>%
      mutate(Z = abs(BETA) / SE) %>%
      separate(SNP, into = c('seqnames', 'start', 'end'), sep = '/', remove = FALSE)
  }
  return(out)
}

library(furrr)
plan(multisession, workers = 24)
result <- phenotypes %>%
  group_nest() %>%
  mutate(out = future_map(
    data,
    \(x) safely(get_snp_effects, otherwise = data.frame())(x)$result, .progress = TRUE)) %>%
  unnest(cols = out)
plan(sequential)

# Plot the results in mini-heatmaps ----------------------------------------------------------------

# most significant cell across all variants and pixels (used for top end of color scale)
top_snp_sig <- slice_min(result, PVAL, with_ties = FALSE)$PVAL

bg <- load.image('data/extracted_fish_warped/replicate_3/gen_4/20221024_IMG_2089.png') %>%
  channel(4) %>% resize_halfXY() %>% as.data.frame() %>% filter(value > 0)

# show successful models
result %>% dplyr::select(x, y) %>% distinct() %>%
  ggplot(aes(x, y, fill = 1)) +
  geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
  geom_raster() + scale_y_reverse() + coord_fixed(expand = FALSE)

minimal_plot <- function(.x) {
  p <- ggplot(.x, aes(x, y, fill = -log10(PVAL))) +
    geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_viridis_c(
      name = expression(-log[10](italic(P))),
      limits = c(3, -log10(top_snp_sig)), na.value = 'grey80', option = 'C'
    ) +
    coord_fixed(expand = FALSE) +
    theme_void()# +
  #guides(fill = 'none')
  return(p)
}

out <- result %>%
  split(result$SNP) %>%
  map(minimal_plot)

write_rds(out, 'paper_figures/subplots/6_snp_effects.rds')
