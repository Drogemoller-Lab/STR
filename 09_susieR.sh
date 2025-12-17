#!/bin/bash
  
# Define the base directory
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
PLINK2="/opt/plink2/plink2"
clsa_data="/home/projects/archive/previous_projects/plinkFiles"

#select STRs within 1Mbp of the top STR chr5_73782133

#bcftools view \
#  -r chr5:73782133-75782133 \
#  $base_dir/chr5/input/annotated_filtered_imputed_STR_SNPs_chr5.vcf.gz \
#  -Oz -o  $base_dir/str_finemapping_susieR/chr5_window_1mb.vcf.gz



#bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' $base_dir/str_finemapping_susieR/chr5_window_1mb.vcf.gz > $base_dir/str_finemapping_susieR/chr5_window_1mb_dosage.tsv


#select the SNPs within the same region and get the dosages

$PLINK2 --bfile $clsa_data/clsa_imp_v3_clean \
	--keep $base_dir/samples/sample_ids_to_include_2ids.txt \
       --chr 5 --from-bp 72782133 --to-bp 74782133 \
       --export A --maf 0.05 \
       --out $base_dir/str_finemapping_susieR/chr5_snps_1mb

