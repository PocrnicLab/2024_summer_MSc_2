#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def main(segs_pvalue_path, cnvr_final_path, output_path, pvalue_threshold):
    try:
        # Read input files
        segs_pvalue_df = pd.read_csv(segs_pvalue_path, sep='\t')
        cnvr_final_df = pd.read_csv(cnvr_final_path, sep='\t')

        # Filter segs_pvalue where MinPvalueAdjusted < threshold
        filtered_segs_pvalue_df = segs_pvalue_df[segs_pvalue_df['MinPvalueAdjusted'] < pvalue_threshold]

        # Initialize result list
        results = []

        # Match filtered data with cnvr_final data
        for _, seg_row in filtered_segs_pvalue_df.iterrows():
            chr = seg_row['seqnames']
            start = seg_row['start']
            end = seg_row['end']
            phenotype = seg_row['Phenotype']

            matched_cnvr_rows = cnvr_final_df[
                (cnvr_final_df['chr'] == chr) & 
                (cnvr_final_df['posStart'] <= start) & 
                (cnvr_final_df['posEnd'] >= end)
            ]

            for _, cnvr_row in matched_cnvr_rows.iterrows():
                result = {
                    'CNVR_ID': cnvr_row['CNVR_ID'],
                    'First.marker.in.the.window': cnvr_row['start_snp'],
                    'Last.marker.in.the.window': cnvr_row['end_snp'],
                    'Trait': phenotype,
                    'CHR': cnvr_row['chr'],
                    'BP1': cnvr_row['posStart'],
                    'BP2': cnvr_row['posEnd']
                }
                results.append(result)

        # Convert results to DataFrame and output, keeping unique CNVR_IDs
        result_df = pd.DataFrame(results).drop_duplicates(subset=['CNVR_ID'])
        result_df.to_csv(output_path, sep='\t', index=False)
        print(f"Results saved to {output_path}")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process segs_pvalue and cnvr_final files.')
    parser.add_argument('segs_pvalue_path', type=str, help='Path to the segs_pvalue file')
    parser.add_argument('cnvr_final_path', type=str, help='Path to the cnvr_final file')
    parser.add_argument('output_path', type=str, help='Path to save the output file')
    parser.add_argument('--pvalue_threshold', type=float, default=0.05, help='Threshold for MinPvalueAdjusted')

    args = parser.parse_args()
    main(args.segs_pvalue_path, args.cnvr_final_path, args.output_path, args.pvalue_threshold)
