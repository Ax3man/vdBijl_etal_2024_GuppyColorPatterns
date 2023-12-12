#!/bin/bash

SAM="/Linux/samtools-1.9/bin/samtools"
BAM=$1
OUTPUT_DIR="coverage"
SAMPLE_NAME=$(basename ${BAM} .bam)
OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_average_coverage.txt"

# Check if output exists and decide what to do
if [ -f "${OUTPUT}" ]; then
    echo "Output for ${SAMPLE_NAME} already exists. Overwriting."
fi

# Calculate average coverage
$SAM depth -q 30 ${BAM} | awk '{sum+=$3} END {print sum/NR}' > ${OUTPUT}

