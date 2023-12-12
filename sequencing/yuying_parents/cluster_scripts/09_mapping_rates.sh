#!/bin/bash

# Path to samtools
samtools_path="/Linux/samtools-1.9/bin/samtools"

# Directory containing BAM files
bam_dir="bam_dedup"
# Output file for mapping rates
output_file="mapping_rates.txt"

# Check if the output file already exists and remove it to start fresh
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# Iterate over each BAM file in the directory
for bam_file in "$bam_dir"/*.bam; do
    # Extract the base name of the file for reporting
    file_name=$(basename "$bam_file")
    
    # Get alignment statistics using samtools flagstat
    stats=$("$samtools_path" flagstat "$bam_file")
    
    # Extract total number of reads and number of mapped reads
    total_reads=$(echo "$stats" | grep "in total" | cut -d ' ' -f 1)
    mapped_reads=$(echo "$stats" | grep "mapped (" | cut -d ' ' -f 1)
    
    # Calculate mapping rate
    if [ $total_reads -gt 0 ]; then
        mapping_rate=$(echo "scale=8; $mapped_reads / $total_reads" | bc)
    else
        mapping_rate="NA"
    fi

    # Output the result to the file
    echo "$file_name: $mapping_rate" >> "$output_file"
done

echo "Mapping rates calculated. Results are saved in $output_file"

