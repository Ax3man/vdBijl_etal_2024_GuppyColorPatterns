#!/bin/bash

## trimming illumina adapters
## Note the use of a custom fasta file with the adapters
## Run 8 of these at the same time, each using 4 cores

trimm=/Linux/Trimmomatic/trimmomatic.jar

cat samplelist.txt | parallel -j 20 \
	"java -jar $trimm PE -threads 2 \
	untrimmed_fastq/{}_R1.fastq.gz \
	untrimmed_fastq/{}_R2.fastq.gz \
	-baseout trimmed_fastq/{}.fastq.gz \
	-trimlog trimlogs/{}.txt \
	ILLUMINACLIP:my_adapters.fa:2:30:10:2:True" 

##while read name; do
##	java -jar $trimm PE -threads 4 \
##	untrimmed_fastq/${name}_R1.fastq.gz \
##     	untrimmed_fastq/${name}_R2.fastq.gz \
##  	-baseout trimmed_fastq/${name}.fastq.gz \
##	-trimlog trimlogs/${name}.txt \
##	ILLUMINACLIP:my_adapters.fa:2:30:10:2:True &
##
##	## At most as number of CPU cores
##	## [ $( jobs | wc -l ) -ge 8 ] && wait
##done < samplelist.txt

