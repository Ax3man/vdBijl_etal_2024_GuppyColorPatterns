library(tidyverse)

f <- list.files('coverage', full.names=TRUE)
f <- f[!(f %in% c("coverage/average_coverages.txt", "coverage/coverage_summary.csv"))]
names(f) <- basename(f) %>% str_remove('_average_coverage.txt')

out <- map_dfr(f, data.table::fread, .id = 'sample_name') %>%
  rename(mean_coverage = V1)

write.csv(out, 'coverage/coverage_summary.csv', row.names = FALSE)
