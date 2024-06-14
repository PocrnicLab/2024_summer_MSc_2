#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np

def read_and_process_data(report_path, map_path, output_path, missing_data_output_path):
    try:
        # Read the first 9 lines of headers
        with open(report_path, 'r') as file:
            header_lines = [next(file) for _ in range(9)]
        
        # Reading the report file starting from the 10th line
        report_df = pd.read_csv(report_path, sep='\t', skiprows=9)
        report_df = report_df[report_df.columns.drop(list(report_df.filter(regex='Allele1 - AB|Allele2 - AB')))]
    except Exception as e:
        print(f"Error reading report file: {e}")
        return

    try:
        # Reading the map file
        map_df = pd.read_csv(map_path, sep='\t', on_bad_lines='skip')
    except Exception as e:
        print(f"Error reading map file: {e}")
        return

    # Merging data
    merged_df = pd.merge(left=report_df, right=map_df, left_on="SNP Name", right_on="Name", how="inner")

    # Define key columns
    key_columns = ["Sample ID", "Chromosome", "Position", "SNP Name", "Log R Ratio", "B Allele Freq"]

    # Identify rows with missing values in key columns
    missing_data_df = merged_df[key_columns].copy()
    missing_data_df = missing_data_df[missing_data_df.isnull().any(axis=1)]

    # Sort the data
    final_df = merged_df[key_columns]
    final_df.rename(columns={"Chromosome": "Chr"}, inplace=True)
    final_df.sort_values(by=["Sample ID", "Chr", "SNP Name"], inplace=True)

    # Remove spaces in Sample ID
    final_df["Sample ID"] = final_df["Sample ID"].str.replace(" ", "")

    # Ensure missing values are represented as 'NaN'
    final_df.replace(['', 'NA', 'N/A', 'nan', 'NaN'], np.nan, inplace=True)
    missing_data_df.replace(['', 'NA', 'N/A', 'nan', 'NaN'], np.nan, inplace=True)

    # Convert NaN values to string 'NaN' for the output
    final_df = final_df.fillna('NaN')
    missing_data_df = missing_data_df.fillna('NaN')

    # Save processed data and missing data report
    with open(output_path, 'w') as f:
        f.writelines(header_lines)  # Write the header lines first
        final_df.to_csv(f, sep='\t', index=False)  # Then append the final DataFrame

    missing_data_df.to_csv(missing_data_output_path, sep='\t', index=False)

    print(f"Final report generated and saved to {output_path}")
    print(f"Missing data report saved to {missing_data_output_path}")

def main():
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('--report_path', required=True, help='Path to the report file')
    parser.add_argument('--map_path', required=True, help='Path to the map file')
    parser.add_argument('--output_path', required=True, help='Path to save the final report')
    parser.add_argument('--missing_data_output_path', required=True, help='Path to save the missing data report')

    args = parser.parse_args()

    read_and_process_data(args.report_path, args.map_path, args.output_path, args.missing_data_output_path)

if __name__ == "__main__":
    main()
