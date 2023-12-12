#!/bin/bash

bwa=/Linux/bin/bwa-0.7.15
gatk=/Linux/GATK4/./gatk

## bwa alignment to the reference, reference needs to be indexed first
## -M : Mark shorter split hits as secondary (for Picard compatibility).
## -t : nr of threads
## -@ : nr of threads
##
## Note that bwa mem has many default arugments that control the performance of the algorithm.
##
## I directly pipe into samtools to make a bam and sort
##
## I am assigning a single read group to each sample

while read name; do
	$bwa mem -t 16 -M \
		../reference/GCF_000633615.1/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna \
		trimmed_fastq/${name}_1P.fastq.gz \
		trimmed_fastq/${name}_2P.fastq.gz | \
	samtools sort -@16 -o tmp/${name}.bam
        $gatk AddOrReplaceReadGroups \
                -I tmp/${name}.bam \
                -O bam/${name}.bam \
                -LB library1 \
                -PL ILLUMINA \
                -PU barcode1 \
                -SM ${name}
	rm tmp/${name}.bam
done < samplelist.txt
	
