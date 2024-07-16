#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse

def calculate_summary_stats(df, snp_df):
    summary = {}
    summary['Total length (Mb)'] = df['posEnd'].sub(df['posStart']).sum() / 1e6
    summary['Total number of CNVR'] = len(df)
    summary['Number of CNVR (< 10 kb)'] = len(df[(df['posEnd'] - df['posStart']) < 10 * 1e3])
    summary['Number of CNVR (10-50 kb)'] = len(df[(df['posEnd'] - df['posStart']).between(10 * 1e3, 50 * 1e3)])
    summary['Number of CNVR (50-100 kb)'] = len(df[(df['posEnd'] - df['posStart']).between(50 * 1e3, 100 * 1e3)])
    summary['Number of CNVR (100-500 kb)'] = len(df[(df['posEnd'] - df['posStart']).between(100 * 1e3, 500 * 1e3)])
    summary['Number of CNVR (500 kb-1 Mb)'] = len(df[(df['posEnd'] - df['posStart']).between(500 * 1e3, 1 * 1e6)])
    summary['Number of CNVR (>= 1 Mb)'] = len(df[(df['posEnd'] - df['posStart']) >= 1 * 1e6])
    
    # Calculate the average number of SNPs per CNVR
    snp_counts = []
    for _, cnvr in df.iterrows():
        chr_match = snp_df['Chr'].astype(str) == str(cnvr['chr'])
        pos_in_range = (snp_df['Position'] >= cnvr['posStart']) & (snp_df['Position'] <= cnvr['posEnd'])
        snp_count = snp_df[chr_match & pos_in_range].shape[0]
        snp_counts.append(snp_count)
    summary['Average number of SNPs per CNVR'] = np.mean(snp_counts)
    
    summary['Minimum size of CNVR (kb)'] = (df['posEnd'] - df['posStart']).min() / 1e3
    summary['Maximum size of CNVR (kb)'] = (df['posEnd'] - df['posStart']).max() / 1e3
    summary['Average CNVR size (kb)'] = (df['posEnd'] - df['posStart']).mean() / 1e3
    summary['Standard deviation of CNVR size (kb)'] = (df['posEnd'] - df['posStart']).std() / 1e3
    return summary

def main(cnvr_file_path, snp_file_path, output_file_path):
    # Read the CNVR data
    cnvr_data = pd.read_csv(cnvr_file_path, delimiter="\t")

    # Read the SNP data
    snp_data = pd.read_csv(snp_file_path, delimiter="\t")

    # Select SNPs from the first individual
    first_individual = snp_data['Sample ID'].unique()[0]
    individual_snp_data = snp_data[snp_data['Sample ID'] == first_individual]

    # Calculate for total, gain, loss, and gain/loss
    total_stats = calculate_summary_stats(cnvr_data, individual_snp_data)
    gain_stats = calculate_summary_stats(cnvr_data[cnvr_data['Type'] == 'Gain'], individual_snp_data)
    loss_stats = calculate_summary_stats(cnvr_data[cnvr_data['Type'] == 'Loss'], individual_snp_data)
    gain_loss_stats = calculate_summary_stats(cnvr_data[cnvr_data['Type'] == 'Mixed'], individual_snp_data)

    # Combine the results into a DataFrame
    summary_df = pd.DataFrame([total_stats, gain_stats, loss_stats, gain_loss_stats], index=['Total', 'Gain', 'Loss', 'Mixed'])

    # Transpose the DataFrame
    transposed_summary_df = summary_df.transpose()

    # Save the transposed summary statistics to a CSV file
    transposed_summary_df.to_csv(output_file_path)

    print("Transposed summary statistics table generated and saved successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate summary statistics for CNVR and SNP data.')
    parser.add_argument('cnvr_file_path', type=str, help='Path to the CNVR file')
    parser.add_argument('snp_file_path', type=str, help='Path to the SNP file')
    parser.add_argument('output_file_path', type=str, help='Path to the output file')

    args = parser.parse_args()

    main(args.cnvr_file_path, args.snp_file_path, args.output_file_path)