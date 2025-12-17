#!/bin/bash

# Description: Full pipeline to prepare genotype files for chromosomes 1 to 22

# Paths
clsa_data="/home/projects/archive/previous_projects/plinkFiles"
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
PLINK="/opt/plink-1.09/plink"

#keep only samples with covariates
##create vcf file for each chromosome
##Loop over chromosomes 1 to 22
for chr in {1..22}; do
  echo "exctracting chr$chr..."

  # Step 1: Extract genotype data using PLINK
  $PLINK --bfile $clsa_data/clsa_imp_v3_clean \
   --keep $base_dir/samples/sample_ids_to_include_2ids.txt \
   --chr $chr --recode vcf bgz \
         --out $base_dir/chr${chr}/input/genotype_chr${chr}

  echo "Completed chr$chr extraction."
  echo "---------------------------------------------"

done

#fix the chr column to include chr string
for chr in 1; do
	echo "fixing chr ${chr} ..."
	bcftools annotate --rename-chrs <(echo -e "${chr}\tchr${chr}") \
  	-o $base_dir/chr${chr}/input/fixed_genotype_chr${chr}.vcf.gz \
  	-O z $base_dir/chr${chr}/input/genotype_chr${chr}.vcf.gz
done
