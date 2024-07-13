#!/bin/bash

module load roslin/plink/1.90p

# Perform SNP association analysis
plink --bfile filtered_data_final \
      --linear \
      --out assoc_results \
      --allow-no-sex \
      --adjust

# Check if the association analysis command was successful
if [ $? -ne 0 ]; then
    echo "PLINK association analysis failed!"
    exit 1
fi

echo "GWAS analysis completed successfully!"
