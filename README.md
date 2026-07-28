# STR Analysis Pipeline

This repository contains the scripts used for the genome-wide STR analyses.
The pipeline includes STR imputation, quality control, association testing, heritability estimation, fine-mapping, haplotype analysis, burden testing, and functional annotation.

## Pipeline Overview

The workflow consists of the following steps:

1. Download STR reference panel
2. Prepare CLSA genotype data
3. Impute STRs using Beagle
4. Annotate STRs
5. Perform quality control
6. Association testing
7. Estimate SNP and STR heritability
8. Fine-mapping
9. Haplotype analysis
10. STR burden analysis
11. Functional prediction using AlphaGenome

## Scripts

00_create_str_dirs.sh - Create project directories
01_download_str_refs.sh - Download STR reference files
02_prepare_clsa_genotypes.sh Prepare genotype data
03_str_imputation_beagle.sh - STR imputation
04_str_annotaTR.sh - STR annotation
05_str_qc.sh - STR quality control
06a_convert_to_npy.py - Convert phenotype files to NumPy arrays
06b_str_assoc_associaTR.sh - STR association analysis
06c_str_assoc_associaTR_rs6453022.sh - Conditional association analysis
07_prepare_heritability_files.sh - Prepare heritability files
08a_heritability_testing_new.sh - Heritability estimation
08b_calculate_SNPbased_LDscore.R - Calculate LD scores
09_susieR.sh - Fine mapping
10_susieR.R - Fine mapping
11a_haplo_prepare_input.sh - Prepare haplotype input
11b_haplotyping.R - Haplotype construction
11c_haplo_assoc_plink.sh - Haplotype association analysis
12_burden_analysis_loop_count.R - STR burden analysis
13a_predict_alphagenome_score_batch_STRs.py - AlphaGenome prediction
13b_prepare_alphagenome_files.R - Prepare AlphaGenome input
