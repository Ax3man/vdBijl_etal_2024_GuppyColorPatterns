library(tidyverse)
library(RcppRoll)
library(furrr)

get_coverage <- function(region, window_size = 1000, window_dist = 100) {
  command <- c("depth", "-a", "-q 30", "-r", region, bam_files)
  output <- system2("samtools", command, stdout = TRUE)
  
  # parse the output very quickly by abusing fread
  dt <- data.table::fread(paste(output, collapse = '\n'))
  colnames(dt) <- c('chr', 'pos', names(bam_files))

  df <- bind_cols(
    data.frame(
      chr = dt$chr[1],
      window_start = roll_min(dt$pos, n = window_size, by = window_dist),
      window_end = roll_max(dt$pos, n = window_size, by = window_dist)
    ),
    map_dfc(dt[, -(1:2)], roll_mean, n = window_size, by = window_dist)
  )
  
  data.table::fwrite(df, paste0('coverage_windows/', region, '.csv'))
}

bam_files <- list.files('bam_dedup', pattern = '\\.bam$', full.names = TRUE) %>%
  setNames(., basename(.) %>% tools::file_path_sans_ext())

intervals <- read.table('../gwas2/gwas_intervals.list')[[1]]

plan(multisession, workers = 32)
future_walk(intervals, get_coverage)
plan(sequential)
