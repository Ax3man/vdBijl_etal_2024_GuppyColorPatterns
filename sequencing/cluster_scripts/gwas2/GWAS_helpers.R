fit_gwas_on_genome <- function(
  kinship_in,
  trait_name,
  phenotype_columns,
  lmm,
  miss,
  maf,
  covar_file = NULL,
  pheno_txt = 'gwas2/pheno.txt',
  gemma_path = '/Linux/bin/gemma',
  plink_path = '/Linux/bin/plink2',
  Rscript_path = '/Linux/R-4.2.0/bin/Rscript',
  ... # passed to plink_make_bed
) {
  require(data.table)
  
  tmp_dir <- tempdir(check = TRUE)

  # helper function to perform all steps for a GWAS on a single VCF and save the results in a tmp folder
  fit_gwas_on_vcf <- function(
    vcf_in,
    tmp_dir
  ) {
    bed_name <- paste0('gwas2/tmp_bed/', trait_name)

    # make bed files
    bed_paths <- plink_make_bed(
      format_out = bed_name, vcf_in = vcf_in, pheno_txt = pheno_txt, plink_path = plink_path, show_output = FALSE, ...
    )
    # The .fam file will only have phenotype column 1, so we update it using an R script that
    # takes a command line argument pointing to the fam file. It assumes we are in the gwas2/ folder
    currwd <- setwd('gwas2/'); on.exit(setwd(currwd))
    args <- paste('update_fam.R --args', bed_paths['fam'])
    system2(Rscript_path, args = args, stdout = FALSE)
    setwd(currwd)

    gwas_paths <- gemma_gwas(
      bed_in = gsub('\\.bed', '', bed_paths['bed']),
      name_out = trait_name,
      dir_out = tmp_dir,
      miss = miss,
      maf = maf,
      covar_file,
      kinship_in = kinship_in, phenotype_columns = phenotype_columns, lmm = lmm, gemma_path = gemma_path
    )

    setwd(currwd)

    file.remove(bed_paths)
    return(gwas_paths)
  }

  out_path <- fit_gwas_on_vcf('gwas2/filtered.vcf.gz', tmp_dir)

  assoc <- data.table::fread(out_path['assoc'])
  data.table::fwrite(assoc, paste0('gwas2/results/', trait_name, '.csv'))
}

#' Perform GWAS using GEMMA, parallelized over scaffolds/chromosomes
#'
#' @param kinship_in File path of the kinship matrix
#' @param trait_name Name used when saving the output in the results folder.
#' @param phenotype_columns Have to be numbers, 1 is refers to the first phenotype column of the
#'  .fam file, which I think is actually the 6th column in that file.
#' @param lmm As in GEMMA
#' @param pheno_txt
#' @param gemma_path
#' @param plink_path
#' @param Rscript_path
#'
#' @return
fit_gwas_on_scaffolds <- function(
    kinship_in,
    trait_name,
    phenotype_columns,
    lmm,
    miss,
    maf,
    pheno_txt = 'gwas2/pheno.txt',
    gemma_path = '/Linux/bin/gemma',
    plink_path = '/Linux/bin/plink2',
    Rscript_path = '/Linux/R-4.2.0/bin/Rscript',
    ... # passed to plink_make_bed
) {
  require(furrr)
  require(purrr)
  require(data.table)

  tmp_dir <- tempdir(check = TRUE)

  # helper function to perform all steps for a GWAS on a single VCF and save the results in a tmp folder
  fit_gwas_on_vcf <- function(
    vcf_in,
    tmp_dir,
    ...
  ) {
    scaff <- tools::file_path_sans_ext(basename(vcf_in))
    bed_name <- paste0('gwas2/tmp_bed/', trait_name, scaff)

    # make bed files
    bed_paths <- plink_make_bed(
      format_out = bed_name, vcf_in = vcf_in, pheno_txt = pheno_txt, plink_path = plink_path, show_output = FALSE, ...
    )
    # The .fam file will only have phenotype column 1, so we update it using an R script that
    # takes a command line argument pointing to the fam file. It assumes we are in the gwas2/ folder
    currwd <- setwd('gwas2/'); on.exit(setwd(currwd))
    args <- paste('update_fam.R --args', bed_paths['fam'])
    system2(Rscript_path, args = args, stdout = FALSE)
    setwd(currwd)    

    gwas_paths <- gemma_gwas(
      bed_in = gsub('\\.bed', '', bed_paths['bed']),
      name_out = scaff,
      dir_out = tmp_dir,
      miss = miss,
      maf = maf,
      kinship_in = kinship_in, phenotype_columns = phenotype_columns, lmm = lmm, gemma_path = gemma_path
    )

    setwd(currwd)

    file.remove(bed_paths)
    return(gwas_paths)
  }

  # get lists of the autosomes
  auto <- data.table::fread('reference/autosomes.list')[[1]]
  vcfs <- c(
    # All the autosome vcfs
    paste0('gwas2/scaffold_vcfs_filtered/', auto, '.vcf.gz'),
    # The VCF for the sex-chromosomes (LG12)
    'gwas2/scaffold_vcfs_filtered/NC_024342.1.vcf.gz',
    # And the combined VCF of all the unplaced scaffolds
    'gwas2/Un.vcf.gz'
  )
  stopifnot((length(vcfs) == 24) && (n_distinct(vcfs) == 24))

  out_paths <- future_map(
    vcfs, \(x) fit_gwas_on_vcf(x, tmp_dir = tmp_dir, ...),
    .progress = FALSE, .options = furrr_options(chunk_size = 1)
  )
  all_assoc <- purrr::map_dfr(map(out_paths, 'assoc'), data.table::fread)
  data.table::fwrite(all_assoc, paste0('gwas2/results/', trait_name, '.csv'))
}

#' Prepare BED file structure using PLINK2
#'
#' @param vcf_in
#' @param format_out
#' @param pheno_txt
#' @param extra_options
#' @param plink_path
#' @param show_output Print stdout to console?
#'
#' @return The paths of the created files, i.e. bed, bim, fam & log
plink_make_bed <- function(
    vcf_in,
    format_out,
    pheno_txt = 'gwas2/pheno.txt',
    extra_options = '--allow-extra-chr --max-alleles 2 --remove gwas2/dropped_samples.txt',
    plink_path = '/Linux/bin/plink2',
    show_output = '', # '' means print to R console
    silence = TRUE
) {
  in_comm <- paste('--vcf', vcf_in)
  out_comm <- paste('--out', format_out)
  pheno_comm <- paste("--pheno 'iid-only'", pheno_txt)

  args <- paste(
    '--vcf', vcf_in,
    '--make-bed', 
    '--out', format_out,
    extra_options,
    "--pheno 'iid-only'", pheno_txt
  )
  if (silence) args <- paste(args, '--silent') 

  #system2(plink_path, args = args, stdout = show_output)
  command <- paste(plink_path, args)
  system(command, wait = TRUE)

  outfiles <- paste0(format_out, c('.bed', '.bim', '.fam', '.log'))
  stopifnot(all(file.exists(outfiles)))
  return(setNames(file.path(getwd(), outfiles), c('bed', 'bim', 'fam', 'log')))
}

#' Estimate kinship matrix using GEMMA
#'
#' @param bed_in
#' @param name_out
#' @param gk 1 is centered kinship matrix, 2 is standardized kinship matrix
#' @param gemma_path
#'
#' @return The paths of the created files
gemma_make_kinship <- function(
  bed_in,
  name_out,
  dir_out,
  gk,
  gemma_path = '/Linux/bin/gemma'
) {
  args <- paste(
    '-bfile', bed_in, 
    '-gk', gk, 
    '-o', name_out,
    '-outdir', dir_out
  )
  system2(gemma_path, args = args, stdout = TRUE)

  if (gk == 1) outfiles <- file.path(dir_out, paste0(name_out, c('.cXX', '.log'), '.txt'))
  if (gk == 2) outfiles <- file.path(dir_out, paste0(name_out, c('.sXX', '.log'), '.txt'))
  stopifnot(all(file.exists(outfiles)))

  return(setNames(outfiles, c('kinship', 'log')))
}

#' Perform GWAS using GEMMA
#'
#' @param bed_in
#' @param kinship_in
#' @param name_out
#' @param phenotype_columns Have to be numbers, 1 is refers to the first phenotype column of the
#'  .fam file, which I think is actually the 6th column in that file.
#' @param lmm
#' @param miss
#' @param maf
#' @param silence If TRUE, don't output status to std out
#' @param gemma_path
#'
#' @return
gemma_gwas <- function(
  bed_in,
  kinship_in,
  name_out,
  dir_out,
  phenotype_columns,
  lmm,
  miss = 0.05,
  maf = 0.01,
  hwe = 0,
  covar_file = NULL,
  silence = FALSE,
  gemma_path = '/Linux/bin/gemma'
) {
  args <- paste(
    '-bfile', bed_in,
    '-k', kinship_in,
    '-n', paste(phenotype_columns, collapse = ' '),
    '-lmm', lmm,
    '-miss', miss,
    '-maf', maf,
    '-hwe', hwe,
    '-o', name_out,
    '-outdir', dir_out
  )
  if (silence) args <- paste(args, '-silence')
  if (!is.null(covar_file)) args <- paste(args, '-c', covar_file)
  
  command <- paste(gemma_path, args)
  system(command, wait = TRUE)

  outfiles <- file.path(dir_out, paste0(name_out, c('.assoc', '.log'), '.txt'))
  stopifnot(all(file.exists(outfiles)))

  return(setNames(outfiles, c('assoc', 'log')))
}

#' Concatenate multiple VCF files into one
#'
#' @param vcfs_in A vector or list of file paths, to .vcf(.gz) files
#' @param vcf_out A file path for the vcf.gz file to be saved
#' @param bcftools_path
#'
#' @return
bcftools_concat <- function(
    vcfs_in,
    vcf_out,
    bcftools_path = '/Linux/bcftools1.13/bin/bcftools'
) {
  args <- paste(
    paste(vcfs_in, collapse = ' '),
    'Oz >', 'vcf_out'
  )
  system2(bcftools_path, args = args)

  stopifnot(file.exists(vcf_out))

  return(setNames(file.path(getwd(), vcf_out), 'vcf'))
}
