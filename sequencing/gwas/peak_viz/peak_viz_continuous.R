peak_viz_continuous <- function(
    trait,
    chr,
    pos,
    file, # set file = NULL to return the plot object instead.
    ornament,
    fdr.level,
    include_kmers,

    pval_column = 'p_SHet',
    viz_width = 1e5,
    viz_range = NULL,
    vcf = 'sequencing/gwas/filtered.vcf.gz',
    overwrite = TRUE,

    reference = 'female'
) {
  if (!overwrite && file.exists(file)) {
    message('File exists and overwrite is FALSE, so skipping.\n')
    return(NULL)
  }

  library(VariantAnnotation)
  library(tidyverse)
  library(patchwork)
  theme_set(theme_classic())
  # This function builds a bunch of panels to further investigate GWAS peaks

  cat(trait, chr, pos, ornament, '\n')

  # source subfunctions
  source('sequencing/gwas/peak_viz/peak_viz_tools.R')
  source('sequencing/genomics_helpers.R')
  list.files('sequencing/gwas/peak_viz/panel_functions', f = TRUE) %>% map(source) %>% invisible()

  if (is.null(viz_range)) {
    viz_range <- c(pos - viz_width/2, pos + viz_width/2)
  }
  gwas <- get_gwas_results(trait, chr, pos, viz_range, fdr.level, reference = reference)
  geno <- get_genotypes(vcf, chr, viz_range)
  region_cov <- get_region_coverage(chr, pos, viz_range, reference = reference)
  linkage <- get_linkage(geno, trait, workers = 12)

  genes <- make_gene_panel(chr, viz_range, reference = reference)
  manhattan <- make_manhattan_linkage_panel(gwas, linkage, viz_range, pval_column)
  depth <- make_region_depth_panel(region_cov, viz_range)

  yuying_region_cov <- get_region_coverage_yuying(chr, pos, viz_range)
  # not currently used
  #yuying_geno <- get_genotypes('sequencing/yuying_parents/yuying_parents_filtered.vcf.gz', chr, viz_range)

  depth_yuying <- depth + aes(color = 'male') +
    geom_line(aes(color = 'female'), data = filter(yuying_region_cov, sex == 'female')) +
    scale_color_manual(
      values = c(female = 'firebrick', male = 'black'),
      name = NULL,
      guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1))
    )

  if (include_kmers) {
    kmer <- make_kmer_manhattan_panel(trait, chr, viz_range)
  } else {
    kmer <- ggplot()
  }
  #snp_effects <- get_snp_effect_heatmap(trait, chr, pos)

  allelic_depth <- make_allelic_depth_panel(geno, chr, pos, reference = reference)
  ornament_effect <- make_ornament_panel(ornament, geno, chr, pos, reference = reference)

  #pl <- genes / manhattan / kmer / depth_yuying / snp_effects / (ornament_effect | allelic_depth)
  pl <- (genes / depth_yuying / manhattan / kmer /  (allelic_depth | ornament_effect)) + plot_layout(heights = c(1, 3, 3, 3, 3))

  #partA / partB + plot_layout(heights = c(1, 1, 1, 1, 2))

  if (!is.null(file)) {
    ggsave(file, pl, width = 10, height = 10)
  } else {
    return(pl)
  }

  # Panels that are no longer used:
  #linkage_plot <- make_linkage_panel(linkage, viz_range)
  #QC1 <- make_QC_panel(geno, viz_range, type = 'site')
  #QC2 <- make_QC_panel(geno, viz_range, type = 'other')
  #excess_het <- make_excess_het_panel(geno, viz_range)
  #sig_range <- range(gwas %>% filter(significant) %>% pull(ps))
  #cov_pheno <- coverage_vs_traits(region_cov, sig_range)
  #snp_pheno <- SNP_vs_traits(geno, pos)
  #linkage_plot_cov <- coverage_linkage(region_cov, sig_range)
}
