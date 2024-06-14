#!/usr/bin/env python3

import pandas as pd
import argparse

def main(final_cnvr_types_path, cnv_clean_path, final_cnv_path, final_cnv_filted_out_path):
    # Read files
    final_cnvr_types = pd.read_csv(final_cnvr_types_path, sep='\t')
    cnv_clean = pd.read_csv(cnv_clean_path, sep='\t')

    # Convert final_cnvr_types to a dictionary with keys as tuples (chr, posStart, posEnd) and values as CNVR_ID
    cnvr_dict = {}
    for _, row in final_cnvr_types.iterrows():
        key = (row['chr'], row['posStart'], row['posEnd'])
        cnvr_dict[key] = row['CNVR_ID']

    # Initialize DataFrames for kept and filtered out CNVs
    kept_cnv = []
    filtered_out_cnv = []

    # Optimize filtering
    for _, cnv_row in cnv_clean.iterrows():
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
    kept_cnv_df = pd.DataFrame(kept_cnv, columns=cnv_clean.columns)
    filtered_out_cnv_df = pd.DataFrame(filtered_out_cnv, columns=cnv_clean.columns)

    # Save results
    kept_cnv_df.to_csv(final_cnv_path, sep='\t', index=False)
    filtered_out_cnv_df.to_csv(final_cnv_filted_out_path, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter CNVs based on CNVR types and save the results.')
    parser.add_argument('final_cnvr_types_path', type=str, help='Path to the final CNVR types file')
    parser.add_argument('cnv_clean_path', type=str, help='Path to the cleaned CNV file')
    parser.add_argument('final_cnv_path', type=str, help='Path to save the kept CNVs')
    parser.add_argument('final_cnv_filted_out_path', type=str, help='Path to save the filtered out CNVs')

    args = parser.parse_args()

    main(args.final_cnvr_types_path, args.cnv_clean_path, args.final_cnv_path, args.final_cnv_filted_out_path)
