#!/bin/bash

bcf=/Linux/bcftools1.13/bin/bcftools
export BCFTOOLS_PLUGINS=/Linux/bcftools1.13/plugins

# We do variant level filtering using bcftools filter.
#	QUAL < 30 			: Phred score of the variant call, exclude (-e) sites with a score < 30, error rate >0.1%
#       F_PASS(FORMAT/DP < 5) > 0.2	: Exclude sites where more than 20% of individuals have a read depth less than 5
#	F_PASS(FORMAT/GQ < 10) > 0.2	: Exclude sites where more than 20% of individuals have a quality score of less than 10 (10% error)
cat ../reference/scaffolds.list | parallel -j 44 \
    "$bcf filter vcfs/scaffold_vcfs/{}.vcf.gz -Oz -o vcfs/scaffold_vcfs_filtered/{}.vcf.gz \
	-e 'QUAL<30|F_PASS(FORMAT/DP < 5) > 0.2|F_PASS(FORMAT/GQ < 10) > 0.2' "

# Note that we are not filtering on:
#	HWE : We applied strong selection, and expect (potentially mismapped) sex-linked variants, so we except violations of HWE
#	missingness : We will apply this filter in gemma (default < 5%)
#	minor allele frequency : We will apply this filter in gemma (default > 0.01)

# Combine the scaffolds into a single vcf (useful later, e.g. to load variants into R)
# We'll make a filtered and an unfiltered one
echo "### bcf concat 1"
$bcf concat -Oz --threads 12 vcfs/scaffold_vcfs_filtered/*.vcf.gz > vcfs/filtered.vcf.gz

# Make a vcf of just the unplaced scaffolds, we can use that in our GWAS instead of running GEMMA on each scaffold
echo '### bcf concat 2'
$bcf concat -Oz --threads 12 vcfs/scaffold_vcfs_filtered/NW*.vcf.gz > vcfs/Un.vcf.gz

# Make a version of the VCF without the unplaced scaffolds, MT, and sex chromosome (LG12).
# this is intented to be used for the estimation of kinship matrix.
echo '### bcf view'
$bcf view --targets-file ../reference/autosomes.list --threads 12 vcfs/filtered.vcf.gz -Oz -o vcfs/filtered_autosomes.vcf.gz
