#!/bin/bash

# Paths
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"

#activate TRTools environment
# Initialize combined files
combined_met_file="$base_dir/results/combined_met_association_results_rs6453022.tsv"

 # Loop over chromosomes 1 to 22
 for chr in {1..22}; do
 
 	for phen in "met"
 	do
 		echo "Starting STR association for $phen phenotype for chr$chr..."
 		
 		# Directories for this chromosome
 		intermediate_dir="$base_dir/chr$chr/intermediate"
 
    		associaTR \
 	   		$base_dir/chr$chr/results/chr${chr}_${phen}_association_results_rs6453022.tsv \
 	   		$base_dir/chr$chr/input/annotated_filtered_imputed_STR_SNPs_chr${chr}.vcf.gz \
 	   		${phen}_phenotype_rs6453022 \
 	   		$base_dir/phenotype/${phen}_merge_filt_condition_rs6453022.npy \
 	   		--same-samples 
 
       		echo "STR association for $phen phenotype for chr$chr is completed..."
 	
           echo "-----------------------------------"
         
     done
 done
 
# Clear the files if they exist (to avoid appending to old files)

> "$combined_met_file"

for chr in {1..22}; do
    for phen in "met"; do
        # Directories for this chromosome
        results_dir="$base_dir/chr$chr/results"
        file="$results_dir/chr${chr}_${phen}_association_results_rs6453022.tsv"

        # Append results to the combined file for each phenotype
        if [ "$phen" = "met" ]; then
            if [ "$chr" -eq 1 ]; then
                cat "$file" >> "$combined_met_file"  # Include header for chr1
            else
                tail -n +2 "$file" >> "$combined_met_file"  # Skip header for other chromosomes
            fi
        fi
    done
done

echo "All chromosomes tested."
echo "All met phenotype results combined into $combined_met_file"


