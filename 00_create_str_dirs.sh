#!/bin/bash

# Define the base directory
base_dir="/home/projects/hearing_loss/clsaARHL_SA/str"

# Create the main directories
#mkdir -p $base_dir/str_ref
mkdir -p $base_dir/phenotype

# Create chromosome-specific directories (chr1 to chr22)
for chr in {1..22}; do
    # Create directories for each chromosome
    mkdir -p $base_dir/chr$chr/input
    mkdir -p $base_dir/chr$chr/intermediate
    mkdir -p $base_dir/chr$chr/results

    # Create reference directories for each chromosome
    mkdir -p $base_dir/str_ref/chr$chr
done

# Print confirmation
echo "Directory structure created successfully!"
echo "Base directory: $base_dir"
