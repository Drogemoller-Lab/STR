#!/bin/bash
  
# Script: download_ref_files.sh
# Description: Downloads STR reference files for each chromosome and organizes them.

# Absolute path for the base directory
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str/str_ref"

# Loop through chromosomes 1 to 22
for chr in 9; do
    echo "Downloading reference files for chr$chr..."

    # Create the chromosome-specific directory if it doesn't exist
    mkdir -p $base_dir/chr$chr

    # Define URL templates inside the loop (after chr is defined)
    vcf_url="https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v4/ensembletr_refpanel_v4_chr${chr}.vcf.gz"
    tbi_url="https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v4/ensembletr_refpanel_v4_chr${chr}.vcf.gz.tbi"
    bref_url="https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v4/ensembletr_refpanel_v4_chr${chr}.bref3"

    # Download the .vcf.gz file
    wget -O $base_dir/chr$chr/ensembletr_refpanel_v4_chr${chr}.vcf.gz $vcf_url

    # Download the .vcf.gz.tbi file
    wget -O $base_dir/chr$chr/ensembletr_refpanel_v4_chr${chr}.vcf.gz.tbi $tbi_url

    # Download the .bref3 file
    wget -O $base_dir/chr$chr/ensembletr_refpanel_v4_chr${chr}.bref3 $bref_url

done

echo "All reference files downloaded and organized successfully!"

