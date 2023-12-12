#!/bin/bash

bcf=/Linux/bcftools1.13/bin/bcftools

# Generate a comma-separated list of BAM files with folder path
bam_list=$(find bam_dedup/ -name "*.bam")

# Perform mpileup and call variants for all BAM files
# mpileup is poor at parallelization, so we parallelize over scaffolds ourselves

while read name; do
	## mpileup options:
	## NOTUSED	--count-orphans / -A 	: Do not skip anomalous read pairs in variant calling
	## 	--adjust-MQ 60 		: Downgrade mapping quality for reads containing excessive mismatches
	##	--min-MQ 30		: Min mapping quality score for an alignment to be used
	## call options:
	## 	--consensus-caller	: Use the concesus caller, calling everything as biallelic
	##	--variants-only		: Only keep variant sites
	$bcf mpileup -Ou \
		--fasta-ref ../reference/GCF_000633615.1/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna \
		--adjust-MQ 60 \
		--min-MQ 30 \
		--regions ${name} \
		--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
		$bam_list | \
	$bcf call --consensus-caller --variants-only --annotate GQ -Oz -o vcfs/scaffold_vcfs/${name}.vcf.gz &

        # Check how many background jobs there are, and if it is equal to the number of cores, wait for anyone 
	# to finish before continuing.
	while :; do
	    background=( $(jobs -p))
	    if (( ${#background[@]} < 46 )); then
	        break
	    fi
	    sleep 1
	done
done < ../reference/scaffolds.list

# Wait for all jobs to finish
wait
