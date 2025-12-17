#!/bin/bash
# Description: Full pipeline to prepare genotype files for chromosomes 1 to 22

# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"


##qc
# Loop over chromosomes 1 to 22
for chr in {4..9}; do

   echo "Starting STR QC for chr$chr..."
   mkdir -p $base_dir/chr${chr}/qc_plots
   
   
   qcSTR --vcf $base_dir/chr${chr}/input/annotated_filtered_imputed_STR_SNPs_chr${chr}.vcf.gz \
	--out $base_dir/chr${chr}/qc_plots/qc-report --vcftype hipstr

done
