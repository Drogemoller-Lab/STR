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
11. Functional annotation and prediction

## Scripts

### 00_create_str_dirs.sh
Creates the chromosome-specific directory structure used for the STR analyses.

### 01_download_str_refs.sh
Downloads the EnsembleTR reference panel used for STR imputation.

### 02_prepare_clsa_genotypes.sh
Prepares chromosome-specific CLSA genotype data for STR imputation.

### 03_str_imputation_beagle.sh
Performs STR imputation using Beagle.

### 04_str_annotaTR.sh
Annotates the imputed STRs using TRTools AnnotaTR.

### 05_str_qc.sh
Performs post-imputation STR quality control.

### 06a_convert_to_npy.py
Prepares phenotype and covariate files for AssociaTR.

### 06b_str_assoc_associaTR.sh
Performs genome-wide STR association testing for the metabolic and sensory phenotypes.

### 07_prepare_heritability_files.sh
Prepares SNP and STR genotype files for heritability analyses.

### 08a_heritability_testing_new.sh
Performs SNP-only and SNP+STR heritability analyses using GCTA.

### 08b_calculate_SNPbased_LDscore.R
Prepares SNP-based LD information for the heritability analyses.

### 09_susieR.sh
Prepares the chr5 locus for SuSiE fine-mapping.

### 10_susieR.R
Performs SuSiE fine-mapping of the chr5 metabolic locus.

### 11a_haplo_prepare_input.sh
Prepares chr5 genotype data for haplotype analysis.

### 11b_haplotyping.R
Constructs haplotypes at the chr5 locus.

### 11c_haplo_assoc_plink.sh
Performs haplotype association testing.

### 12_burden_analysis_loop_count.R
Performs STR burden analyses.

### 13a_predict_alphagenome_score_batch_STRs.py
Generates AlphaGenome predictions for selected STRs.

### 13b_prepare_alphagenome_files.R
Processes AlphaGenome prediction results.

### str_annotation.R
Performs functional annotation using TR-xQTL resources.

### gencode.sh
Performs GENCODE-based genomic annotation.

### 16_cpg_annotation.R
Performs CpG annotation of associated STR loci.
