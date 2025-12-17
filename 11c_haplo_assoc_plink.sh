#!/bin/bash
  
# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
PLINK="/opt/plink-1.09/plink"
clsa_data="/home/projects/hearing_loss/clsaARHL_SA/clsaFiles"

##add metabolic phenotypes
#$PLINK --bfile $base_dir/haplotype/chr5_window \
#	--update-sex $clsa_data/updated_sex_v2.txt \
#	--pheno $clsa_data/phen_met_better.txt \
#	--make-bed --out $base_dir/haplotype/chr5_window_met
#
#$PLINK --bfile $base_dir/haplotype/chr5_window_met \
#       --linear --ci 0.95 --covar $clsa_data/covar_all_v2.txt \
#       --covar-name sen.better.rank, age, diabetes-PC10\
#       --out $base_dir/haplotype/chr5_window_met_gwas


#sort and filter
#for phenotype in "chr5_window_met_gwas"
#do
#        sort -g -k12 $base_dir/haplotype/${phenotype}.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $base_dir/haplotype/${phenotype}_sorted_noNA_add
#done
#

#try the 20kbp window with 5kbp slide
##add metabolic phenotypes
#
#$PLINK --bfile $base_dir/haplotype/chr5_window.20kbp \
#       --update-sex $clsa_data/updated_sex_v2.txt \
#       --pheno $clsa_data/phen_met_better.txt \
#       --make-bed --out $base_dir/haplotype/chr5_window.20kbp_met
#
#$PLINK --bfile $base_dir/haplotype/chr5_window.20kbp_met \
#       --linear --ci 0.95 --covar $clsa_data/covar_all_v2.txt \
#       --covar-name sen.better.rank, age, diabetes-PC10\
#       --out $base_dir/haplotype/chr5_window_met.20kbp_gwas
#
#
##sort and filter
#for phenotype in "chr5_window_met.20kbp_gwas"
#do
#        sort -g -k12 $base_dir/haplotype/${phenotype}.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $base_dir/haplotype/${phenotype}_sorted_noNA_add
#done
#


#################
##try blocks with double the marker size (5)
##get the formated .fam file
#mv chr5_window.5mrk.fam old_chr5_window/
#cp chr5_window.fam chr5_window.5mrk.fam
#
###add metabolic phenotypes
#$PLINK --bfile $base_dir/haplotype/chr5_window.5mrk \
#       --update-sex $clsa_data/updated_sex_v2.txt \
#       --pheno $clsa_data/phen_met_better.txt \
#       --make-bed --out $base_dir/haplotype/chr5_window.5mrk_met
#
#$PLINK --bfile $base_dir/haplotype/chr5_window.5mrk_met \
#       --linear --ci 0.95 --covar $clsa_data/covar_all_v2.txt \
#       --covar-name sen.better.rank, age, diabetes-PC10\
#       --out $base_dir/haplotype/chr5_window.5mrk_met_gwas
#
#
###sort and filter
#for phenotype in "chr5_window.5mrk_met_gwas"
#do
#        sort -g -k12 $base_dir/haplotype/${phenotype}.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $base_dir/haplotype/${phenotype}_sorted_noNA_add
#done
#


##try one block with only top variants
##get the formated .fam file
#mv top_var.fam old_chr5_window/
#cp chr5_window.fam top_var.fam
#
###add metabolic phenotypes
#$PLINK --bfile $base_dir/haplotype/top_var \
#       --update-sex $clsa_data/updated_sex_v2.txt \
#       --pheno $clsa_data/phen_met_better.txt \
#       --make-bed --out $base_dir/haplotype/top_var_met
#
#$PLINK --bfile $base_dir/haplotype/top_var_met \
#       --linear --ci 0.95 --covar $clsa_data/covar_all_v2.txt \
#       --covar-name sen.better.rank, age, diabetes-PC10\
#       --out $base_dir/haplotype/top_var_met_gwas
#
#
###sort and filter
#for phenotype in "top_var_met_gwas"
#do
#        sort -g -k12 $base_dir/haplotype/${phenotype}.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $base_dir/haplotype/${phenotype}_sorted_noNA_add
#done
#



#try 1 STR alt allele at a time

#get the formated .fam file
for i in {0..5}
do 
	mv top_var_maf01.${i}.fam old_chr5_window/
	cp chr5_window.fam top_var_maf01.${i}.fam
	
	##add metabolic phenotypes
	$PLINK --bfile $base_dir/haplotype/top_var_maf01.${i} \
		--update-sex $clsa_data/updated_sex_v2.txt \
       		--pheno $clsa_data/phen_met_better.txt \
       		--make-bed --out $base_dir/haplotype/top_var_maf01.${i}_met

	$PLINK --bfile $base_dir/haplotype/top_var_maf01.${i}_met \
       		--linear --ci 0.95 --covar $clsa_data/covar_all_v2.txt \
       		--covar-name sen.better.rank, age, diabetes-PC10\
       		--out $base_dir/haplotype/top_var_maf01.${i}_met_gwas

	##sort and filter
	sort -g -k12 $base_dir/haplotype/top_var_maf01.${i}_met_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $base_dir/haplotype/top_var_maf01.${i}_met_gwas_sorted_noNA_add

done

