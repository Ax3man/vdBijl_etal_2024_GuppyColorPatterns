suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
})

cat('Starting to make pheno.txt... ')

sample_list1 <- data.table::fread('gwas2/sequencing_samples_list.csv') %>% as_tibble()
sample_list2 <- data.table::fread('gwas2/NovaSeq_sample_list.csv') %>%
  as_tibble() %>%
  select(Name, prefix = `Filename Prefix`)
sample_list <- right_join(sample_list1, sample_list2, c('sample_id' = 'Name')) %>%
  group_by(sample_id, fish_id) %>%
  summarise(prefix = ifelse(n() == 1, prefix, paste0('sample_', sample_id[[1]])), .groups = 'drop')

em <- read_rds('gwas2/phenotypes.rds')
em <- em %>%
  as.data.frame() %>%
  rownames_to_column(var = 'unique_id')
pd <- data.table::fread('gwas2/photo_database.csv') %>%
  mutate(fish_id = toupper(fish_id)) %>%
  select(fish_id, facing_direction, unique_id)
pheno <- left_join(em, pd, 'unique_id') %>%
  group_by(fish_id, facing_direction) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop_last') %>%
  summarise(across(where(is.numeric), mean))

# find which samples are in the VCF file:
# (Note that we will remove some other samples later, in 07b_run_gwas.R):
incl_samples <- system('/Linux/bcftools1.13/bin/bcftools query -l gwas2/filtered.vcf.gz', intern = TRUE)

#	## * If you also need PLINK 1.9 to read this file, add an FID column in front,
#	##   and fill it with zeroes.
#	## * 'site' is a categorical variable.  --glm would ignore it if you loaded it
#	##   as a phenotype (multinomial logistic regression is not implemented), but
#	##   it's a valid *covariate* for --glm.
#	#IID  qt1    bmi    site
#	1110  2.3    22.22  site2
#	2202  34.12  18.23  site1
#	...

joined <- left_join(sample_list, pheno, 'fish_id')
  
pheno_out <- joined %>%
  select(-selection, -starts_with('Y_haplogroup')) %>%
  select('#IID' = prefix, everything(), -fish_id, -sample_id) %>%
  arrange(`#IID`) %>%
  filter(`#IID` %in% incl_samples)

write_tsv(pheno_out, 'gwas2/pheno.txt')


dropped <- read.table('gwas2/dropped_samples.txt')$V1
covar_out <- joined %>%
  # remove individuals that are in dropped_samples.txt
  filter(!prefix %in% dropped) %>%
  mutate(intercept = 1) %>%
  select(intercept, starts_with('Y_haplogroup'))

write_tsv(covar_out, 'gwas2/covar.tsv', col_names = FALSE)

cat('Finished.\n')

