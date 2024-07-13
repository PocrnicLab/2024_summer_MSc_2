#!/usr/bin/env python3

import pandas as pd
import argparse

def read_and_merge_files(file1, file2):
    # Read the files
    df1 = pd.read_csv(file1, delim_whitespace=True)
    df2 = pd.read_csv(file2, delim_whitespace=True)
    
    # Merge the files on CHR and SNP columns
    merged_df = pd.merge(df1, df2, on=["CHR", "SNP"])
    return merged_df

def filter_and_save(merged_df, threshold, output_file):
    # Filter rows where FDR_BH is less than the threshold
    filtered_df = merged_df[merged_df['FDR_BH'] < threshold]
    
    # Select required columns and rename them
    output_df = filtered_df[['SNP', 'CHR', 'BP']]
    output_df.columns = ['SNP reference', 'CHR', 'BP']
    
    # Save to a file
    output_df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Filter GWAS results based on FDR_BH threshold.")
    parser.add_argument('file1', type=str, help="Path to assoc_results.assoc.linear")
    parser.add_argument('file2', type=str, help="Path to assoc_results.assoc.linear.adjusted")
    parser.add_argument('threshold', type=float, help="FDR_BH threshold value")
    parser.add_argument('output', type=str, help="Output file path")
    args = parser.parse_args()

    merged_df = read_and_merge_files(args.file1, args.file2)
    filter_and_save(merged_df, args.threshold, args.output)

if __name__ == "__main__":
    main()
