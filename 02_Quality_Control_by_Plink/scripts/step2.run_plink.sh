#!/bin/bash

# Check if the correct number of arguments is provided
if [[ $# -gt 1 ]]; then
    echo "Usage: $0 [output_directory]"
    exit 1
fi

# Get the output directory from the command line arguments, default to current directory if not provided
output_directory=${1:-$(pwd)}/plink_results

# Define input and output file prefixes
input_prefix="input_data"
output_prefix="filtered_data"

module load roslin/plink/1.90p

# Check if input files exist
if [[ ! -f "${input_prefix}.ped" || ! -f "${input_prefix}.map" ]]; then
    echo "Error: input_data.ped or input_data.map not found."
    exit 1
fi

# Initial filtering: remove monomorphic SNPs and SNPs with MAF less than 0.01, also filter for genotype and sample missingness
plink --file $input_prefix --maf 0.01 --geno 0.05 --mind 0.05 --make-bed --out ${output_prefix}_maf --allow-extra-chr
if [[ $? -ne 0 ]]; then
    echo "Error: PLINK command failed at MAF filtering step."
    exit 1
fi

# HWE filtering: remove SNPs with Hardy-Weinberg equilibrium p-value less than 0.01
plink --bfile ${output_prefix}_maf --hwe 0.001 --make-bed --out ${output_prefix}_hwe --allow-extra-chr
if [[ $? -ne 0 ]]; then
    echo "Error: PLINK command failed at HWE filtering step."
    exit 1
fi

# Final filtering: check for missing genotype data again
plink --bfile ${output_prefix}_hwe --geno 0.05 --mind 0.05 --make-bed --out ${output_prefix}_final --allow-extra-chr
if [[ $? -ne 0 ]]; then
    echo "Error: PLINK command failed at final filtering step."
    exit 1
fi

# Output the results
initial_individuals=$(wc -l < ${input_prefix}.ped)
initial_snps=$(wc -l < ${input_prefix}.map)
final_individuals=$(wc -l < ${output_prefix}_final.fam)
final_snps=$(wc -l < ${output_prefix}_final.bim)

echo "Filtering complete, result file prefix is ${output_prefix}_final"
echo "Initial number of individuals: $initial_individuals"
echo "Initial number of SNPs: $initial_snps"
echo "Final number of individuals: $final_individuals"
echo "Final number of SNPs: $final_snps"

# Create the output directory if it doesn't exist
mkdir -p $output_directory

# Move all files starting with filtered_data to the specified output directory
mv ${output_prefix}* $output_directory/

# Define the final output file prefix
final_output_prefix="${output_directory}/filtered_data_final"

# Extract individual IDs from .fam file
cut -d ' ' -f 1,2 ${output_directory}/${output_prefix}_final.fam > ${final_output_prefix}_individuals.txt

# Extract SNP IDs from .bim file
cut -f 2 ${output_directory}/${output_prefix}_final.bim > ${final_output_prefix}_snps.txt
