#!/bin/bash

# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
beagle="/home/projects/temp/beagle.17Dec24.224.jar"


# Loop over chromosomes 1 to 22
#for chr in 9 ; do
#
#   echo "Starting STR imputation for chr$chr..."
#
##perform imputation
#java -Xmx100g -jar $beagle \
#	gt=$base_dir/chr${chr}/input/fixed_genotype_chr${chr}.vcf.gz \
#       ref=$base_dir/str_ref/chr${chr}_new/ensembletr_refpanel_v3_chr${chr}.bref3 \
#       out=$base_dir/chr${chr}/input/imputed_STR_SNPs_chr${chr}
#
#   echo "STR imputation for chr$chr is completed..."
#   echo "-----------------------------------"
#done
#
#echo "All chromosomes imputed."

#remove low quality STR
#QC DR2 >= 0.7
for chr in {7..9}; do
	echo "Starting STR filtering for chr$chr..."
	bcftools filter -i 'DR2>=0.7' $base_dir/chr${chr}/input/imputed_STR_SNPs_chr${chr}.vcf.gz \
		-Oz -o $base_dir/chr${chr}/input/filtered_imputed_STR_SNPs_chr${chr}.vcf.gz
	echo "Done STR filtering for chr$chr..."
	
	#generate index file for the imputed vcf
   	
	echo "Indexing chr$chr..."
	tabix -f -p vcf $base_dir/chr${chr}/input/filtered_imputed_STR_SNPs_chr${chr}.vcf.gz
	echo "Done indexing chr$chr..."

done
