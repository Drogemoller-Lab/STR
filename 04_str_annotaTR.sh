#!/bin/bash

# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"

#activate TRTools environment

# Loop over chromosomes 1 to 22
for chr in {1..22}; do

   echo "Starting STR annotation for chr$chr..."
   annotaTR \
	   --vcf $base_dir/chr${chr}/input/filtered_imputed_STR_SNPs_chr${chr}.vcf.gz \
	   --out $base_dir/chr${chr}/input/annotated_filtered_imputed_STR_SNPs_chr${chr} \
	   --ref-panel $base_dir/str_ref/chr$chr/ensembletr_refpanel_v4_chr${chr}.vcf.gz \
	   --outtype vcf \
           --vcftype hipstr \
	   --vcf-outtype z

   echo "STR annotation for chr$chr is completed..."
   echo "-----------------------------------"
done

echo "All chromosomes annotated."
