#!/usr/bin/env python3

import pandas as pd
import argparse

def main(input_file_path, output_file_path):
    # Load data
    data = pd.read_csv(input_file_path, sep='\t')

    # Select data only from PennCNV
    penncnv_data = data[data['method'] == 'PennCNV']

    # Compute statistics
    result = []

    # Calculate statistics by CNV_type
    for cnv_type, group in penncnv_data.groupby('CNV_type'):
        num_cnvs = group.shape[0]
        mean_length = group['length'].mean()
        min_length = group['length'].min()
        max_length = group['length'].max()
        result.append([cnv_type, num_cnvs, mean_length, min_length, max_length])

    # Compute overall statistics
    total_num_cnvs = penncnv_data.shape[0]
    total_mean_length = penncnv_data['length'].mean()
    total_min_length = penncnv_data['length'].min()
    total_max_length = penncnv_data['length'].max()
    result.append(['Total', total_num_cnvs, total_mean_length, total_min_length, total_max_length])

    # Create DataFrame with the results
    columns = ['CNV_type', 'Number_of_CNVs', 'Mean_Length', 'Min_Length', 'Max_Length']
    result_df = pd.DataFrame(result, columns=columns)

    # Save the results to a CSV file
    result_df.to_csv(output_file_path, index=False)

    print(f"Statistics have been saved to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process CNV data.')
    parser.add_argument('input_file', type=str, help='Path to the input file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    args = parser.parse_args()

    main(args.input_file, args.output_file)
