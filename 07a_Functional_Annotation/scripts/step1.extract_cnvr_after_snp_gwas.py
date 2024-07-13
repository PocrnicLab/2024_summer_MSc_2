#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def main(snp_path, cnvr_final_path, output_path, phenotype):
    try:
        # Read input files
        snp_df = pd.read_csv(snp_path, sep='\t')
        cnvr_final_df = pd.read_csv(cnvr_final_path, sep='\t')

        # Initialize result list
        results = []

        # Match SNP data with cnvr_final data
        for _, snp_row in snp_df.iterrows():
            chr = snp_row['CHR']
            bp = snp_row['BP']

            matched_cnvr_rows = cnvr_final_df[
                (cnvr_final_df['chr'] == chr) & 
                (cnvr_final_df['posStart'] <= bp) & 
                (cnvr_final_df['posEnd'] >= bp)
            ]

            for _, cnvr_row in matched_cnvr_rows.iterrows():
                result = {
                    'CNVR_ID': cnvr_row['CNVR_ID'],
                    'Trait': phenotype,  # Use the phenotype provided by the user
                    'CHR': cnvr_row['chr'],
                    'BP1': cnvr_row['posStart'],
                    'BP2': cnvr_row['posEnd'],
                    'First_marker_in_the_window': cnvr_row['start_snp'],
                    'Last_marker_in_the_window': cnvr_row['end_snp']
                }
                results.append(result)

        # Convert results to DataFrame and output, keeping unique CNVR_IDs
        result_df = pd.DataFrame(results).drop_duplicates(subset=['CNVR_ID'])
        result_df.to_csv(output_path, sep='\t', index=False)
        print(f"Results saved to {output_path}")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process SNP and cnvr_final files.')
    parser.add_argument('snp_path', type=str, help='Path to the SNP file')
    parser.add_argument('cnvr_final_path', type=str, help='Path to the cnvr_final file')
    parser.add_argument('output_path', type=str, help='Path to save the output file')
    parser.add_argument('phenotype', type=str, help='Phenotype to include in the output')

    args = parser.parse_args()
    main(args.snp_path, args.cnvr_final_path, args.output_path, args.phenotype)
