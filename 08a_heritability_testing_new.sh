#!/bin/bash
  
# Paths
clsa_data="/home/projects/archive/previous_projects/plinkFiles"
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
PLINK="/opt/plink-1.09/plink"
GCTA="/opt/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static"

############################

#for the SNP only analysis

#$PLINK --bfile $clsa_data/clsa_imp_v3_clean \
#          --keep $base_dir/samples/sample_ids_to_include_2ids.txt \
#          --make-bed \
#          --out $base_dir/heritability/snp_heritability/input/chrAll_snp

##Step 1: segment based LD score
#$GCTA --bfile $base_dir/heritability/snp_heritability/input/chrAll_snp \
#       --ld-score-region 200 \
#       --out $base_dir/heritability/snp_heritability/input/chrAll_snp \
#       --thread-num 36
#
##Step 2 (option #1): stratify the SNPs by LD scores of individual SNPs in R
#Rscript calculate_SNPbased_LDscore.R snp_heritability chrAll_snp.score.ld
#
##Step 3: making GRMs using SNPs stratified into different groups
#for gr in {1..4}; do
#       $GCTA --bfile $base_dir/heritability/snp_heritability/input/chrAll_snp \
#               --extract $base_dir/heritability/snp_heritability/input/snp_group$gr.txt \
#               --make-grm --out $base_dir/heritability/snp_heritability/input/chrAll_snp_group$gr \
#	       --thread-num 36
#done
#
##Step 4: REML analysis with multiple GRMs

#for phen in "met" "sen"
#do
#        $GCTA --reml --mgrm $base_dir/heritability/snp_heritability/input/multi_GRMs_snp.txt \
#                --pheno $base_dir/phenotype/${phen}_phenotype_better.txt \
#                --out $base_dir/heritability/snp_heritability/results/${phen}_snp_heritability \
#		--thread-num 36
#done


#for SNPs plus STRs
#prepapre STR file and combine with SNP file
#for data with STRs try maf 0.01 and geno 0.05
#$PLINK --vcf $base_dir/heritability/str_heritability/input/chrAll_str_biallelic.vcf.gz \
#       --maf 0.01 --geno 0.05 --make-bed --recode vcf bgz \
#       --out $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str_biallelic_filtered
#
#$PLINK --bfile $clsa_data/clsa_imp_v3_clean \
#          --keep $base_dir/samples/sample_ids_to_include_2ids.txt \
#          --recode vcf bgz \
#          --out $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp
#
###create index for the files
#tabix -p vcf $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str_biallelic_filtered.vcf.gz

#fix the samples in the STR
#bcftools concat -a \
#       -o $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb.vcf.gz \
#       -O z $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp.vcf.gz \
#       $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str_biallelic_filtered.vcf.gz
#echo "concat is completed"
#str=bcftools view $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str_biallelic_filtered.vcf.gz | grep -v -c "^#"
#snp=bcftools view $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp.vcf.gz | grep -v -c "^#" 
#snp_str=bcftools view $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb.vcf.gz | grep -v -c "^#" 
#
#echo "we have $str_var STRs, $snp SNPs and $snp_str SNPs + STRs"
#
###create index for the files
#tabix -p vcf $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb.vcf.gz
#
##convert to plink format 
##for data with STRs try maf 0.01 and geno 0.05
#$PLINK --vcf $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb.vcf.gz \
#       --make-bed --threads 36 \
#       --out $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb_plink


###Step 1: segment based LD score
$GCTA --bfile $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb_plink \
       --ld-score-region 200 \
       --out $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb_plink \
       --thread-num 36

###Step 2 (option #1): stratify the SNPs by LD scores of individual SNPs in R
Rscript 08b_calculate_SNPbased_LDscore.R snp_str_comb_heritability_new chrAll_snp_str_comb_plink.score.ld

##Step 3: making GRMs using SNPs stratified into different groups
for gr in {1..4}; do
       $GCTA --bfile $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_str_comb_plink \
               --extract $base_dir/heritability/snp_str_comb_heritability_new/input/snp_group$gr.txt \
               --make-grm --out $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_snp_group$gr \
              --thread-num 36
done

#Step 4: REML analysis with multiple GRMs

for phen in "met" "sen"
do
        $GCTA --reml --mgrm $base_dir/heritability/snp_str_comb_heritability_new/input/multi_GRMs_snp.txt \
                --pheno $base_dir/phenotype/${phen}_phenotype_better.txt \
                --out $base_dir/heritability/snp_str_comb_heritability_new/results/${phen}_snp__str_heritability \
               --thread-num 36
done


