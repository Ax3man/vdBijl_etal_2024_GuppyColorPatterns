# get the .vcf.gz files
f <- list.files('gvcf', pattern = '\\.vcf\\.gz$', full.names = TRUE)
# get the sample names by removing the folders and extension
n <- gsub('\\.g\\.vcf\\.gz', '', basename(f))

# get the .tbi files (we only want to include .vcf files with a .tbi, since those are finished)
i <- list.files('gvcf', pattern = '\\.tbi')
# get the name (to match to f)
i_n <- gsub('\\.g\\.vcf\\.gz\\.tbi', '', basename(i))

f <- f[n %in% i_n]
n <- n[n %in% i_n]

write.table(
	x = cbind(n, f),
	file = 'genomicsDB/guppies.sample_map',
	sep = '\t',
	quote = FALSE,
	row.names = FALSE,
	col.names = FALSE
)
