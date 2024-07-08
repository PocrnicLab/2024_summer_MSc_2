#!/usr/bin/env python3

import pandas as pd
import sys

def filter_cnv_in_cnvr(cnv_file, cnvr_file, output_file):
    # Read the input files
    cnv_df = pd.read_csv(cnv_file, sep="\t")
    cnvr_df = pd.read_csv(cnvr_file, sep="\t")

    # Add a column to store the filter results
    cnv_df['in_cnvr'] = False

    # Iterate over each row in the CNVR file
    for index, row in cnvr_df.iterrows():
        chr_cnvr = row['chr']
        posStart_cnvr = row['posStart']
        posEnd_cnvr = row['posEnd']
        
        # Update the in_cnvr column for CNVs that match the current CNVR range
        cnv_df.loc[
            (cnv_df['chr'] == f'chr{chr_cnvr}') &
            (cnv_df['start'] >= posStart_cnvr) &
            (cnv_df['end'] <= posEnd_cnvr),
            'in_cnvr'
        ] = True

    # Filter the CNV dataframe to keep only rows that are within any CNVR
    filtered_cnv_df = cnv_df[cnv_df['in_cnvr']].copy()
    filtered_cnv_df.drop(columns=['in_cnvr'], inplace=True)

    # Remove the 'chr' prefix in the 'chr' column for the filtered DataFrame
    filtered_cnv_df['chr'] = filtered_cnv_df['chr'].str.replace('chr', '')

    # Save the filtered results to a new file
    filtered_cnv_df.to_csv(output_file, sep="\t", index=False)

    print(f"Filtered results have been saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <cnv_file_path> <cnvr_file_path> <output_file_path>")
        sys.exit(1)
    
    cnv_file = sys.argv[1]
    cnvr_file = sys.argv[2]
    output_file = sys.argv[3]

    filter_cnv_in_cnvr(cnv_file, cnvr_file, output_file)
