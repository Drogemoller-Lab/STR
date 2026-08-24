#!/bin/bash
# Description: Full pipeline to prepare genotype files for chromosomes 1 to 22

# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"


##qc
# Loop over chromosomes 1 to 22
for chr in {1..22}; do

   echo "Starting STR QC for chr$chr..."
   mkdir -p $base_dir/chr${chr}/qc_plots
  
  
   qcSTR --vcf $base_dir/chr${chr}/input/annotated_filtered_imputed_STR_SNPs_chr${chr}.vcf.gz \
	--out $base_dir/chr${chr}/qc_plots/qc-report --vcftype hipstr

done

#QC for all chromosomes
echo "Starting STR QC for all CHRs"
qcSTR --vcf $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str.vcf.gz \
	--out $base_dir/qc-report_comb --vcftype hipstr

