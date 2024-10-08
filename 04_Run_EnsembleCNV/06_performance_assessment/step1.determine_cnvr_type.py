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

def merge_overlapping_cnvr(cnvr_df):
    merged_cnvr = []
    current_cnvr = None

    for index, row in cnvr_df.iterrows():
        if current_cnvr is None:
            current_cnvr = row
        else:
            if row['chr'] == current_cnvr['chr'] and row['posStart'] <= current_cnvr['posEnd']:
                current_cnvr['posEnd'] = max(current_cnvr['posEnd'], row['posEnd'])
                current_cnvr['end_snp'] = row['end_snp']
            else:
                merged_cnvr.append(current_cnvr)
                current_cnvr = row

    if current_cnvr is not None:
        merged_cnvr.append(current_cnvr)

    return pd.DataFrame(merged_cnvr)

def main(cnvr_file_path, cnv_file_path, output_file_path):
    try:
        cnvr_df = pd.read_csv(cnvr_file_path, sep='\t', encoding='utf-8')
        cnv_df = pd.read_csv(cnv_file_path, sep='\t', encoding='utf-8')
    except Exception as e:
        print(f"Error reading files: {e}")
        raise

    cnv_df = cnv_df[['chr', 'posStart', 'posEnd', 'CN']]
    
    # Sort CNVRs by chromosome and start position
    cnvr_df = cnvr_df.sort_values(by=['chr', 'posStart'])

    # Merge overlapping CNVRs
    cnvr_df = merge_overlapping_cnvr(cnvr_df)

    result = []

    for index, cnvr_row in cnvr_df.iterrows():
        cnvr_chr = cnvr_row['chr']
        cnvr_start = cnvr_row['posStart']
        cnvr_end = cnvr_row['posEnd']
        start_snp = cnvr_row['start_snp']
        end_snp = cnvr_row['end_snp']

        cnv_in_cnvr = cnv_df[(cnv_df['chr'] == cnvr_chr) & 
                             (cnv_df['posStart'] >= cnvr_start) & 
                             (cnv_df['posEnd'] <= cnvr_end)]

        if not cnv_in_cnvr.empty:
            cn_values = cnv_in_cnvr['CN'].tolist()
            cnvr_type = determine_cnvr_type(cn_values)
        else:
            cnvr_type = 'Undefined'

        result.append([cnvr_row['CNVR_ID'], cnvr_row['chr'], cnvr_row['arm'], cnvr_row['posStart'], cnvr_row['posEnd'], start_snp, end_snp, cnvr_type])

    result_df = pd.DataFrame(result, columns=['CNVR_ID', 'chr', 'arm', 'posStart', 'posEnd', 'start_snp', 'end_snp', 'Type'])
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
