#/bin/bash/
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
PLINK2="/opt/plink2/plink2"
PLINK="/opt/plink-1.09/plink"

##merge filtered imputed vcf files for heritability testing

#bcftools concat \
#       -o $base_dir/heritability/snp_str_comb_heritability/input/chrAll_snp.vcf.gz \
#       -O z $base_dir/chr*/input/filtered_imputed_STR_SNPs_chr*.vcf.gz 

#for str files
vcf_list="$base_dir/heritability/vcfs_to_concat.txt"
#> $vcf_list
#for chr in {1..22}; do
#  echo "$base_dir/chr${chr}/input/annotated_filtered_imputed_STR_SNPs_chr${chr}.vcf.gz" >> $vcf_list
#done

bcftools concat --file-list $vcf_list -O z -o $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str.vcf.gz

#index the file 
tabix -p vcf $base_dir/heritability/snp_str_comb_heritability_new/input/chrAll_str.vcf.gz

#variants=$(bcftools view -c1 "$base_dir/heritability/snp_str_comb_heritability/input/chrAll_snp_str.vcf.gz" | grep -v -c '^#')
#echo "Total variant are $variants"


###create index for the merged file
#tabix -p vcf $base_dir/heritability/str_heritability/input/chrAll_str.vcf.gz
#
##convert multi-allelic to bi-allelic 
#bcftools norm \
#	-m -both $base_dir/heritability/str_heritability/input/chrAll_str.vcf.gz \
#	-Oz -o $base_dir/heritability/str_heritability/input/chrAll_str_biallelic.vcf.gz 

#convert to plink
#for data with STRs try maf 0.01 and geno 0.05
#$PLINK --vcf $base_dir/heritability/snp_str_comb_heritability/input/chrAll_snp_str_biallelic.vcf.gz \
#       --maf 0.01 --geno 0.05 --make-bed \
#       --out $base_dir/heritability/snp_str_comb_heritability/input/chrAll_snp_str_biallelic

