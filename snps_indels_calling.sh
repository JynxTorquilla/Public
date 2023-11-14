#!/bin/bash

# Input VCF file
input_vcf="*.vcf.gz"
# Output files for indels and SNPs
output_indels="indels.vcf.gz"
output_snps="snps.vcf.gz"
# Range in base pairs
range=50
# Create tmp directoty
mkdir tmp
# Create output directory
mkdir output
            
# Create .txt file for chromosome renaming            
echo "chr1 1
chr2 2
chr3 3
chr4 4
chr5 5
chr6 6
chr7 7
chr8 8
chr9 9
chr10 10
chr11 11
chr12 12
chr13 13
chr14 14
chr15 15
chr16 16
chr17 17
chr18 18
chr19 19
chr20 20
chr21 21
chr22 22
chrY Y
chrX X
chrMT MT" > mapping.txt

# indels and snps calling
bcftools index $input_vcf
bcftools view -v snps $input_vcf -Oz -o ./tmp/$output_snps
bcftools view -v indels $input_vcf -Oz -o ./tmp/$output_indels

# rename output sns and indels chromosomes before normalization
bcftools annotate --rename-chr mapping.txt ./tmp/$output_snps -Oz -o ./tmp/rnmd-chr.snps.vcf.gz 
bcftools annotate --rename-chr mapping.txt ./tmp/$output_indels -Oz -o ./tmp/rnmd-chr.indels.vcf.gz

# normalization
bcftools norm -f *.fa ./tmp/rnmd-chr.snps.vcf.gz -Oz -o ./tmp/norm_snps.vcf.gz -c s
bcftools norm -f *.fa ./tmp/rnmd-chr.indels.vcf.gz -Oz -o ./tmp/norm_indels.vcf.gz -c s

# Create VCF files for indels and SNPs within the range
bcftools filter --IndelGap ${range} ./tmp/norm_indels.vcf.gz -Oz -o ./output/norm.filt-indels.vcf.gz
bcftools filter --SnpGap ${range} ./tmp/norm_snps.vcf.gz -Oz -o ./output/norm.filt-snps.vcf.gz

echo "Indels and SNPs extracted, normalized, and saved in separate VCF files."
