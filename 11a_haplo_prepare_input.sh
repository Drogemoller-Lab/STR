#!/bin/bash
  
# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
PLINK="/opt/plink-1.09/plink"


##subset to get the specific region (more efficient)
#bcftools view -r chr5:72000000-75000000 $base_dir/chr5/input/filtered_imputed_STR_SNPs_chr5.vcf.gz \
#       	-Oz -o $base_dir/haplotype/filtered_imputed_STR_SNPs_chr5_subset.vcf.gz
#
##get SNP + STR file (imputed from beagle) and convert to bi-allelic
#bcftools norm \
#       -m -both $base_dir/haplotype/filtered_imputed_STR_SNPs_chr5_subset.vcf.gz \
#       -Oz -o $base_dir/haplotype/filtered_imputed_STR_SNPs_chr5_subset_biallelic.vcf.gz
#
##index the output file
#tabix -p vcf $base_dir/haplotype/filtered_imputed_STR_SNPs_chr5_subset_biallelic.vcf.gz


# Extract sample IDs from VCF
bcftools query -l $base_dir/haplotype/filtered_imputed_STR_SNPs_chr5_subset_biallelic.vcf.gz > $base_dir/haplotype/temp_samples.txt
