#!/bin/bash

## use gatk4 to mark duplicate reads (e.g. due to PCR)

gatk=/Linux/GATK4/./gatk
sam=/Linux/samtools-1.9/bin/samtools

while read name; do
	$gatk --java-options "-Xmx100g -Xms100g" MarkDuplicates \
		-I bam/${name}.bam \
		-O bam_dedup/${name}.bam \
		-M dedup_logs/${name}.txt
	$sam index bam_dedup/${name}.bam
done < samplelist.txt
