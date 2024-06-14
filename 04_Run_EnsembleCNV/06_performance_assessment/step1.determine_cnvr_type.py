#!/usr/bin/env python3

import pandas as pd
import argparse

def determine_cnvr_type(cn_values):
    cn_values_set = set(cn_values)
    if cn_values_set.issubset({0, 1}):
        return 'Loss'
    elif cn_values_set.issubset({3}):
        return 'Gain'
    elif cn_values_set.intersection({0, 1}) and cn_values_set.intersection({3}):
        return 'Mixed'
    else:
        return 'Undefined'

def main(cnvr_file_path, cnv_file_path, output_file_path):
    # Load files with error handling
    try:
        cnvr_df = pd.read_csv(cnvr_file_path, sep='\t', encoding='utf-8')
        cnv_df = pd.read_csv(cnv_file_path, sep='\t', encoding='utf-8')
    except Exception as e:
        print(f"Error reading files: {e}")
        raise

    # Extract necessary columns
    cnv_df = cnv_df[['chr', 'posStart', 'posEnd', 'CN']]

    # Initialize results list
    result = []

    # Process each CNVR
    for index, cnvr_row in cnvr_df.iterrows():
        cnvr_chr = cnvr_row['chr']
        cnvr_start = cnvr_row['posStart']
        cnvr_end = cnvr_row['posEnd']

        # Find all CNVs within the range of CNVR
        cnv_in_cnvr = cnv_df[(cnv_df['chr'] == cnvr_chr) & 
                             (cnv_df['posStart'] >= cnvr_start) & 
                             (cnv_df['posEnd'] <= cnvr_end)]

        # If CNVs are present, determine the type
        if not cnv_in_cnvr.empty:
            cn_values = cnv_in_cnvr['CN'].tolist()
            cnvr_type = determine_cnvr_type(cn_values)
        else:
            cnvr_type = 'Undefined'  # If no CNVs, type is undefined

        # Add results to the list
        result.append([cnvr_row['CNVR_ID'], cnvr_row['chr'], cnvr_row['arm'], cnvr_row['posStart'], cnvr_row['posEnd'], cnvr_type])

    # Convert results to DataFrame and save
    result_df = pd.DataFrame(result, columns=['CNVR_ID', 'chr', 'arm', 'posStart', 'posEnd', 'Type'])
    try:
        result_df.to_csv(output_file_path, index=False, sep='\t', encoding='utf-8')
        print("CNVR types have been successfully determined and saved to:", output_file_path)
    except Exception as e:
        print(f"Error saving results: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process CNVR and CNV files to determine CNVR types.')
    parser.add_argument('cnvr_file_path', type=str, help='Path to the CNVR file')
    parser.add_argument('cnv_file_path', type=str, help='Path to the CNV file')
    parser.add_argument('output_file_path', type=str, help='Path to the output file')

    args = parser.parse_args()

    main(args.cnvr_file_path, args.cnv_file_path, args.output_file_path)
