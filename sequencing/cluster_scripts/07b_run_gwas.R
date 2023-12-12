# whether the kinship matrix should be calculated this run
make_kinship <- TRUE
kinship_type <- 'GEMMA' # Either GEMMA or KING

# Note that we are dropping these samples, since its relatedness does not match the pedigree:
# NS.2125.002.IDT_i7_111---IDT_i5_111.280
# NS.2145.001.IDT_i7_89---IDT_i5_89.355
# This is listed in dropped_samples.txt, which is used in the code below

# PLINK extra options:
extra_options <- "--allow-extra-chr --max-alleles 2 --remove gwas2/dropped_samples.txt"

library(furrr)
options(datatable.showProgress = FALSE)
# load helper functions. Note that basically all these functions are run for their side effects
# (i.e. they save files on disk), and typically return the paths to their output files.
source('gwas2/GWAS_helpers.R')

# Make folders, if they do not already exist
if (!dir.exists('gwas2/results')) dir.create('gwas2/results')
if (!dir.exists('gwas2/tmp_bed')) dir.create('gwas2/tmp_bed')

# Make pheno.txt, based on phenotypes.rds and sample lists
source('gwas2/make_pheno_txt.R')

if (make_kinship) {
  if (kinship_type == "KING") {
    ## KING kinship estimation
    args <- paste(
      '--vcf gwas2/filtered_autosomes.vcf.gz',
      '--make-king square',
      '--out gwas2/king/king_kinship',
      extra_options
    )
    system2('plink2', args)
    kinship <- 'gwas2/king/king_kinship.king'
  } else {
    ## GEMMA kinship:
    cat('Making kinship bed...\n')
    kinship_bed <- plink_make_bed(
      vcf_in = 'gwas2/filtered_autosomes.vcf.gz',
      format_out = 'gwas2/tmp_bed/for_kinship',
      extra_options = extra_options
    )
    cat('Estimating kinship matrix...\n')
    kinship <- gemma_make_kinship(
      bed_in = 'gwas2/tmp_bed/for_kinship',
      name_out = 'kinship',
      dir_out = 'gwas2/results',
      gk = 1
    )['kinship']
 }
} else {
  cat('Using existing kinship matrix...\n')
  if (kinship_type == "KING") {
    kinship <- 'gwas2/king/king_kinship.king'
  } else {
    kinship <- 'gwas2/results/kinship.cXX.txt'
  }
  if (!file.exists(kinship)) stop('Kinship file does not exist. Do you need to enable kinship calculation?')
}

# Load the phenotypes file to get the column names. That way we can rely on phenotype names instead of indices
# We drop the first, because the first column has IDs
phenotypes <- colnames(data.table::fread('gwas2/pheno.txt'))[-1]

prop_missing_allowed <- 0.05
min_maf <- 0.1
hwe <- 0
#lmm <- 2			# 1: Wald, 2: LRT, 3: Score, 4: All # we set this in the my_gwas function below instead

bed_path <- plink_make_bed(
  vcf_in = 'gwas2/filtered.vcf.gz',
  format_out = 'gwas2/filtered',
  silence = FALSE
)

# The .fam file will only have phenotype column 1, so we update it using an R script that
# takes a command line argument pointing to the fam file. It assumes we are in the gwas2/ folder
currwd <- setwd('gwas2/')
args <- paste('update_fam.R --args', bed_path['fam'])
system2('/Linux/R-4.2.0/bin/Rscript', args = args, stdout = FALSE)
setwd(currwd)

my_gwas <- function(trait, name, lmm, use_covars = FALSE) {
  col <- which(phenotypes == trait)

 args <- list(
    bed_in = 'gwas2/filtered',
    kinship_in = kinship,
    name_out = name,
    dir_out = 'gwas2/results/',
    phenotype_columns = col,
    lmm = lmm, 
    miss = prop_missing_allowed, 
    maf = min_maf, 
    hwe = hwe
  )
  if (use_covars) { 
    #This enables the covariate file, that controls for Y-haplogroups    
    args <- c(args, covar_file = 'gwas2/covar.tsv')    
  }
  do.call(gemma_gwas, args)
}

trait_cols <- data.frame(
  trait = c(
    'car_perc', 
    'mel_perc_v2',
    paste0('car_PIE_', 1:5),
    paste0('mel_PIE_', 1:5),
    paste0('pa_car_', 1:7),
    paste0('pa_mel_', 1:8)
  ),
  # We need to specify the test here. For univariate, use the more reliable p-lrt, for multivariate
  # we set Wald, so we obtain Beta estimates and SE in the output, needed for Z-scores
  lmm = c(
    rep(2, 2),
    rep(1, 5),
    rep(1, 5),
    rep(2, 7),
    rep(2, 8)    
  )
)
trait_cols$name <- trait_cols$trait

plan(multisession, workers = min(nrow(trait_cols), 32))
on.exit(plan(sequential))
future_pwalk(trait_cols, my_gwas, use_covars = FALSE)

trait_cols2 <- trait_cols
trait_cols2$name <- paste0(trait_cols2$name, '_Yhap')
future_pwalk(trait_cols2, my_gwas, use_covars = TRUE)
