#!/usr/bin/env python3

import pandas as pd
import argparse

def main(cnvr_file_path, snp_file_path, final_samples_file_path, output_file_path):
    # Read CNVR file
    cnvr_df = pd.read_csv(cnvr_file_path, sep='\t')

    # Read SNP file, skipping the header
    snp_df = pd.read_csv(snp_file_path, sep='\t', comment='[', skiprows=9)

    # Filter out SNPs not on chromosomes 1-26
    snp_df = snp_df[snp_df['Chr'].astype(str).isin([str(i) for i in range(1, 27)])]

    # Read final_samples.txt to get the list of sample IDs
    with open(final_samples_file_path, 'r') as file:
        final_samples = file.read().splitlines()
    final_samples = final_samples[1:]  # Skip the header

    # Filter SNPs to keep only those in the final sample list
    snp_df = snp_df[snp_df['Sample ID'].isin(final_samples)]

    # Initialize an empty DataFrame to store the filtered SNPs
    filtered_snps_df = pd.DataFrame()

    # Filter SNPs using vectorized operations
    for _, cnvr in cnvr_df.iterrows():
        chr_match = snp_df['Chr'].astype(str) == str(cnvr['chr'])
        pos_in_range = (snp_df['Position'] >= cnvr['posStart']) & (snp_df['Position'] <= cnvr['posEnd'])
        filtered_snps = snp_df[chr_match & pos_in_range]
        filtered_snps_df = pd.concat([filtered_snps_df, filtered_snps])

    # Drop duplicates if any SNPs fall within multiple CNVRs
    filtered_snps_df = filtered_snps_df.drop_duplicates()

    # Save the filtered SNPs to a new file
    filtered_snps_df.to_csv(output_file_path, sep='\t', index=False)

    print(f"Filtered SNPs have been saved to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter SNPs based on CNVRs and sample IDs, and save the results.')
    parser.add_argument('cnvr_file_path', type=str, help='Path to the CNVR file')
    parser.add_argument('snp_file_path', type=str, help='Path to the SNP file')
    parser.add_argument('final_samples_file_path', type=str, help='Path to the final samples file')
    parser.add_argument('output_file_path', type=str, help='Path to the output file')

    args = parser.parse_args()

    main(args.cnvr_file_path, args.snp_file_path, args.final_samples_file_path, args.output_file_path)
