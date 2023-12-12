#!/bin/bash

BAM_DIR="bam_dedup"
OUTPUT_DIR="coverage"
mkdir -p ${OUTPUT_DIR}

# Run the calculate_coverage script in parallel
find ${BAM_DIR} -name "*.bam" | parallel -j 42 --eta ./calc_coverage.sh {}

echo "Done!"


## Old method using gatk
#gatk=/Linux/GATK4/./gatk

## -R is the reference
## -L is the intervals, I'm just using the whole chromosomes from the reference
## -I is input
## -O is output
##
## --omit-depth-output-at-each-base disables per base statistics, and speeds up running time

#bam_list=$(find bam_dedup/ -name "*.bam" -exec basename {} \;)

#echo "$bam_list" | parallel -j 36 \
#        "$gatk DepthOfCoverage \
#                -R reference/GCF_000633615.1/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna \
#                -L reference/intervals.list \
#                -I bam_dedup/{} \
#                -O coverage/{/.}.cov"
