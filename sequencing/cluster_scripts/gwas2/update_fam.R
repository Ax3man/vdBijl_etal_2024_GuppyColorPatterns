#this script assumes the working directory is the folder with the .fam and supporting file,
#probably gwas2/ or gwas/

args <- commandArgs(trailingOnly = TRUE)
fam_file <- args[2]
if (length(args) < 3) {
  incl_header = FALSE 
} else {
  incl_header <- args[3]
}
if (length(args) < 4) {
  fam_file_out <- fam_file
} else {
  fam_file_out <- args[4]
}

print(fam_file)
print(incl_header)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

cat('Updating .fam file with additional phenotypes... ')

sample_list1 <- data.table::fread('sequencing_samples_list.csv') %>% as_tibble()
sample_list2 <- data.table::fread('NovaSeq_sample_list.csv') %>%
  as_tibble() %>%
  select(Name, prefix = `Filename Prefix`)
sample_list <- full_join(sample_list1, sample_list2, c('sample_id' = 'Name')) %>%
  group_by(sample_id, fish_id) %>%
  summarise(prefix = ifelse(n() == 1, prefix, paste0('sample_', sample_id[[1]])), .groups = 'drop')

em <- read_rds('phenotypes.rds')
colnames(em) <- paste0('V', seq_len(ncol(em)))
em <- em %>%
  as.data.frame() %>%
  rownames_to_column(var = 'unique_id')
pd <- data.table::fread('photo_database.csv') %>%
  mutate(fish_id = toupper(fish_id)) %>%
  select(fish_id, facing_direction, unique_id)
pheno <- left_join(pd, em, 'unique_id') %>%
  group_by(fish_id, facing_direction) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop_last') %>%
  summarise(across(where(is.numeric), mean))

to_join <- full_join(
  sample_list,
  pheno,
  'fish_id'
) %>% 
  select(-fish_id, -sample_id) %>%
  filter(!is.na(prefix)) %>%
  # drop the first phenotype column, since that's "selection", a covariate that isn't in the .fam at all
  # drop the second phenotype column, since that's already in the .fam (typically car_perc)
  # drop V4:V6, since those are also covariates, the Y-haplogroups
  select(-V1, -V2, -V4, -V5, -V6)

# Out .fam file, is like this, tab separated
# - family ID (can be all 0)
# - sample ID (like NS.2118.001.IDT_i7_100---IDT_i5_100.5, so out$prefix)
# - Sire ID (if in sample)
# - Dam ID (if in sample)
# - Sex code (1 = male, 2 = feamle, 0 = unknown)
# - Phenotypes

fam <- data.table::fread(fam_file, sep = '\t')
# set all fish to males, not actually used I think
fam$V5 <- 1

fam_out <- left_join(fam, to_join, by = c('V2' = 'prefix'))

if (incl_header) {
  colnames(fam_out) <- c('FID', 'IID', 'Sire', 'Dam', 'Sex', paste0('V', 1:(ncol(fam_out) - 5)))
  write_tsv(
    fam_out,
    fam_file_out
  )
} else {
  write_tsv(
    fam_out,
    fam_file_out,
    col_names = FALSE
  )
}

cat('Finished.\n')
