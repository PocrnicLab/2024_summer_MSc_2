#!/usr/bin/env python3

import pandas as pd
import argparse

def update_positions(file1_path, file2_path, output_path):
    # Read file1 and file2
    file1 = pd.read_csv(file1_path, sep='\t')  # Assume file1 is tab-separated
    file2 = pd.read_csv(file2_path, sep='\t')  # Assume file2 is tab-separated
    
    # Merge based on SNP_name from file1 and Name from file2
    merged_df = pd.merge(file1, file2, left_on='SNP_name', right_on='Name')
    
    # Update Position and Chromosome in file2 with values from file1
    merged_df['Position'] = merged_df['position']
    merged_df['Chromosome'] = merged_df['chromosome']
    
    # Select necessary columns, keeping only rows present in file1
    result_df = merged_df[file2.columns]
    
    # Sort the result by Index
    result_df = result_df.sort_values(by='Index')
    
    # Filter out rows where Position is 0
    result_df = result_df[result_df['Position'] != 0]
    
    # Write the result to the output file
    result_df.to_csv(output_path, sep='\t', index=False)

def main():
    # Define command-line argument parser
    parser = argparse.ArgumentParser(description='Update positions and chromosomes in file2 based on file1')
    parser.add_argument('file1', type=str, help='Path to the first input file')
    parser.add_argument('file2', type=str, help='Path to the second input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    
    args = parser.parse_args()
    
    # Update positions and chromosomes in file2 based on file1
    update_positions(args.file1, args.file2, args.output)

if __name__ == '__main__':
    main()
