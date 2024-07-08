#!/usr/bin/env python3

import pandas as pd
import sys

def filter_cnv(cnv_file_path, phenotypes_file_path, output_file_path):
    # Read the phenotypes file
    phenotypes_df = pd.read_csv(phenotypes_file_path, sep='\t')
    # Extract all unique Sample_IDs
    sample_ids = phenotypes_df['sample id'].unique()

    # Read the formatted_cnv file
    cnv_df = pd.read_csv(cnv_file_path, sep='\t')

    # Filter rows where Sample_ID is present in phenotypes file
    filtered_cnv_df = cnv_df[cnv_df['sample id'].isin(sample_ids)]

    # Save the filtered data to a new file
    filtered_cnv_df.to_csv(output_file_path, sep='\t', index=False)

    print(f"Filtered CNV data has been saved to {output_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <cnv_file_path> <phenotypes_file_path> <output_file_path>")
        sys.exit(1)
    
    cnv_file_path = sys.argv[1]
    phenotypes_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    filter_cnv(cnv_file_path, phenotypes_file_path, output_file_path)
