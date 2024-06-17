#!/usr/bin/env python3

import pandas as pd
import argparse

def main(final_cnvr_types_path, cnv_clean_path, final_cnv_path, final_cnv_filtered_out_path, phenotype_path, filtered_phenotype_path):
    # Read files
    final_cnvr_types = pd.read_csv(final_cnvr_types_path, sep='\t')
    cnv_clean = pd.read_csv(cnv_clean_path, sep='\t')
    phenotype_data = pd.read_csv(phenotype_path, delim_whitespace=True)

    # Print column names for debugging
    #print("Columns in phenotype data file:", phenotype_data.columns)

    # Extract idanim from phenotype data
    idanim_list = phenotype_data['idanim'].unique()
    
    # Filter cnv_clean based on idanim
    filtered_cnv_clean = cnv_clean[cnv_clean['Sample_ID'].isin(idanim_list)]

    # Convert final_cnvr_types to a dictionary with keys as tuples (chr, posStart, posEnd) and values as CNVR_ID
    cnvr_dict = {}
    for _, row in final_cnvr_types.iterrows():
        key = (row['chr'], row['posStart'], row['posEnd'])
        cnvr_dict[key] = row['CNVR_ID']

    # Initialize lists for kept and filtered out CNVs
    kept_cnv = []
    filtered_out_cnv = []

    # Optimize filtering
    for _, cnv_row in filtered_cnv_clean.iterrows():
        matched = False
        for (chr_, pos_start, pos_end), cnvr_id in cnvr_dict.items():
            if (cnv_row['chr'] == chr_ and
                cnv_row['posStart'] >= pos_start and
                cnv_row['posEnd'] <= pos_end):
                kept_cnv.append(cnv_row)
                matched = True
                break
        if not matched:
            filtered_out_cnv.append(cnv_row)

    # Convert lists to DataFrames
    kept_cnv_df = pd.DataFrame(kept_cnv, columns=filtered_cnv_clean.columns)
    filtered_out_cnv_df = pd.DataFrame(filtered_out_cnv, columns=filtered_cnv_clean.columns)

    # Save results
    kept_cnv_df.to_csv(final_cnv_path, sep='\t', index=False)
    filtered_out_cnv_df.to_csv(final_cnv_filtered_out_path, sep='\t', index=False)
    phenotype_data.to_csv(filtered_phenotype_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter CNVs based on phenotype data and CNVR types, then save the results.')
    parser.add_argument('final_cnvr_types_path', type=str, help='Path to the final CNVR types file')
    parser.add_argument('cnv_clean_path', type=str, help='Path to the cleaned CNV file')
    parser.add_argument('final_cnv_path', type=str, help='Path to save the kept CNVs')
    parser.add_argument('final_cnv_filtered_out_path', type=str, help='Path to save the filtered out CNVs')
    parser.add_argument('phenotype_path', type=str, help='Path to the phenotype data file')
    parser.add_argument('filtered_phenotype_path', type=str, help='Path to save the filtered phenotype data')

    args = parser.parse_args()

    main(args.final_cnvr_types_path, args.cnv_clean_path, args.final_cnv_path, args.final_cnv_filtered_out_path, args.phenotype_path, args.filtered_phenotype_path)
