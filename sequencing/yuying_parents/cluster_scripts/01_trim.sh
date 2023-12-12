#!/bin/bash

## trimming illumina adapters
## Note the use of a custom fasta file with the adapters
## Run 8 of these at the same time, each using 4 cores

trimm=/Linux/Trimmomatic/trimmomatic.jar

cat samplelist.txt | parallel -j 6 \
	"java -jar $trimm PE -threads 2 \
	untrimmed_fastq/{}_R1.fastq.gz \
	untrimmed_fastq/{}_R2.fastq.gz \
	-baseout trimmed_fastq/{}.fastq.gz \
	ILLUMINACLIP:../my_adapters.fa:2:30:10:2:True" 
