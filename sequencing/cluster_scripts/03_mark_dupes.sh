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

## At this time, we merge the bam files of the two samples that got sequenced twice.
## (The doubles are because the sequencing center messed up)
# First we merge:
$sam merge \
	bam_dedup/289_merged.bam \
	bam_dedup/NS.2145.001.IDT_i7_77---IDT_i5_77.289.bam \
	bam_dedup/NS.2118.002.IDT_i7_19---IDT_i5_19.289.bam
# Since we used the file name as the sample name in our read groups, we need to overwrite those. So we keep
# the read groups, but make sure both of them have the same sample name. This is necessary for mpileup
# and call to generate a vcf with the correct number of samples.
$sam view -H bam_dedup/289_merged.bam | \
    awk -v OFS='\t' -v sample_name="sample_289" '/^@RG/ { sub(/SM:[^\t]+/, "SM:" sample_name) } { print }' | \
    samtools reheader - bam_dedup/289_merged.bam > tmp.bam
mv tmp.bam bam_dedup/289_merged.bam
# Now index the bam file
$sam index bam_dedup/289_merged.bam
# And now we remove the other files, since we don't want to include the separate files in the downstream scripts
rm bam_dedup/NS.2145.001.IDT_i7_77---IDT_i5_77.289.bam bam_dedup/NS.2118.002.IDT_i7_19---IDT_i5_19.289.bam

# And finally, we take all the exact same steps for sample 323:
$sam merge \
	bam_dedup/323_merged.bam \
	bam_dedup/NS.2145.001.IDT_i7_144---IDT_i5_144.323.bam \
	bam_dedup/NS.2118.002.IDT_i7_155---IDT_i5_155.323.bam
$sam view -H bam_dedup/323_merged.bam | \
    awk -v OFS='\t' -v sample_name="sample_323" '/^@RG/ { sub(/SM:[^\t]+/, "SM:" sample_name) } { print }' | \
    samtools reheader - bam_dedup/323_merged.bam > tmp.bam
mv tmp.bam bam_dedup/323_merged.bam
$sam index bam_dedup/323_merged.bam
rm bam_dedup/NS.2145.001.IDT_i7_144---IDT_i5_144.323.bam bam_dedup/NS.2118.002.IDT_i7_155---IDT_i5_155.323.bam
