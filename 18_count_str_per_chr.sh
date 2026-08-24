#/bin/bash

base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"
output_file="$base_dir/str_variant_counts.txt"
> "$output_file"  # Clear the file if it exists

for chr in {1..22}; do
    count=$(bcftools view \
        "$base_dir/chr${chr}/input/annotated_filtered_imputed_STR_SNPs_chr${chr}.vcf.gz" | \
        grep -v -c '^#')

    echo "chr${chr}: $count" >> "$output_file"
done

