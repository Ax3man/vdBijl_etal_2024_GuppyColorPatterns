suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(tidyverse)
})
vcf <- 'sequencing/gwas/filtered.vcf.gz'

## Define helper functions
source('sequencing/genomics_helpers.R')

source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('sequencing/gwas/peak_viz/panel_functions/make_linkage_panel.R')

run <- function(traits, name = traits, pval_column, overwrite = FALSE) {
  outfile <- glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{name}.csv')
  if (!overwrite && file.exists(outfile)) {
    message('File exists and overwrite is FALSE, so skipping')
    return(invisible(NULL))
  }

  # Load GWAS results
  if (length(traits) == 1) {
    g <- prep_gwas_table(traits, FALSE, pval_column, fdr_level = 0.05) %>%
      filter(significant)
  } else {
    g <- map_dfr(
      setNames(traits, traits),
      \(tr) prep_gwas_table(tr, FALSE, pval_column, fdr_level = 0.05, calc_q = FALSE),
      .id = 'trait'
    )
    # the joint-set q-value
    q <- qvalue(g[[pval_column]], fdr.level = 0.05)
    g$qvalue <- q$qvalues
    g$significant <- q$significant
    g <- g %>%
      dplyr::select(chr, ps, allele1, significant) %>%
      dplyr::filter(significant) %>%
      distinct()
  }

  regions <- GRanges(
    seqnames = g$chr,
    ranges = IRanges(start = g$ps, end = g$ps)
  )
  geno <- readVcfAsVRanges(vcf, param = ScanVcfParam(which = regions)) %>%
    as.data.frame() %>%
    filter(
      !(sampleNames %in% c('NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355'))
    )
  cat('\nLoaded genotypes...\n')
  link <- get_linkage(geno, workers = 4, min_MAF = 0.1)

  data.table::fwrite(link, outfile)
}

run('car_PIE', pval_column = 'p_SHet')
run('mel_PIE', pval_column = 'p_SHet')
run(glue::glue('pa_car_{x}', x = 1:7), name = 'orange_ornaments', pval_column = 'p_lrt')
run(glue::glue('pa_mel_{x}', x = 1:8), name = 'black_ornaments', pval_column = 'p_lrt')

