Rs4=/Linux/R-4.2.0/bin/Rscript

#bash 01_trim.sh
bash 02_map_to_bam.sh
bash 03_mark_dupes.sh
bash 04_coverage.sh

bash 05b_bcftools_call.sh
bash 06b_combine_filter_vcf.sh
#$Rs4 07b_run_gwas.R
