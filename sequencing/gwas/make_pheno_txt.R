suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
})

cat('Starting to make pheno.txt... ')

sample_list1 <- data.table::fread('data/sequencing_samples_list.csv') %>% as_tibble()
sample_list2 <- data.table::fread('data/NovaSeq_sample_list.csv') %>%
  as_tibble() %>%
  dplyr::select(Name, prefix = `Filename Prefix`)
sample_list <- full_join(sample_list1, sample_list2, c('sample_id' = 'Name')) %>%
  mutate(prefix = case_when(
    prefix == 'NS.2145.001.IDT_i7_77---IDT_i5_77.289' ~ 'sample_289',
    prefix == 'NS.2145.001.IDT_i7_144---IDT_i5_144.323' ~ 'sample_323',
    TRUE ~ prefix
  ))

em <- read_rds('sequencing/phenotypes.rds')
colnames(em) <- paste0('V', seq_len(ncol(em)))
em <- em %>%
  as.data.frame() %>%
  rownames_to_column(var = 'unique_id')
pd <- data.table::fread('photo_database.csv') %>%
  mutate(fish_id = toupper(fish_id)) %>%
  dplyr::select(fish_id, facing_direction, unique_id)
pheno <- left_join(em, pd, 'unique_id') %>%
  group_by(fish_id, facing_direction) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop_last') %>%
  summarise(across(where(is.numeric), mean))

# find which samples are in the VCF file:
sample_names <- readVcfAsVRanges(
  'sequencing/gwas/filtered.vcf.gz',
  param = ScanVcfParam(which = GRanges("NC_024342.1", IRanges(26272132)))
) %>% as.data.frame() %>% dplyr::pull(sampleNames) %>% as.character()

#	## * If you also need PLINK 1.9 to read this file, add an FID column in front,
#	##   and fill it with zeroes.
#	## * 'site' is a categorical variable.  --glm would ignore it if you loaded it
#	##   as a phenotype (multinomial logistic regression is not implemented), but
#	##   it's a valid *covariate* for --glm.
#	#IID  qt1    bmi    site
#	1110  2.3    22.22  site2
#	2202  34.12  18.23  site1
#	...

out <- left_join(sample_list, pheno, 'fish_id') %>%
  dplyr::select('#IID' = prefix, everything(), -fish_id, -sample_id) %>%
  drop_na() %>%
  arrange(`#IID`) %>%
  filter(`#IID` %in% sample_names)

write_tsv(out, 'sequencing/gwas/pheno.txt')

cat('Finished.\n')

