#!/bin/bash

BAM_DIR="bam_dedup"
OUTPUT_DIR="coverage"
mkdir -p ${OUTPUT_DIR}

# Run the calculate_coverage script in parallel
find ${BAM_DIR} -name "*.bam" | parallel -j 12 --eta ../calc_coverage.sh {}

echo "Done!"
