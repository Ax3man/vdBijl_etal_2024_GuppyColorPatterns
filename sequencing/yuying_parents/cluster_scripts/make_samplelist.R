f <- list.files('untrimmed_fastq', pattern = '_R.')
f <- gsub('_R.', '', f)
f <- gsub('\\.fastq\\.gz', '', f)
f <- unique(f)

# exlude files that already have a dedupped bam with index (the index means it is completed)
#e <- list.files('bam_dedup', pattern = '\\.bam\\.bai')
#e <- gsub('\\.bam\\.bai', '', e)

if (exists('e')) f <- f[!(f %in% e)]

cat(f, file = 'samplelist.txt', sep = '\n')
