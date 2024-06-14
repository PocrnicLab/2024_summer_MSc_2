#!/usr/bin/env python3

import pandas as pd
import argparse

def main(input_file, output_file):
    # Read TXT file
    df = pd.read_csv(input_file, sep='\t')

    # Filter rows where Type is not 'Undefined'
    filtered_df = df[df['Type'] != 'Undefined']

    # Save the results to a new TXT file
    filtered_df.to_csv(output_file, index=False, sep='\t')

    print(f"Filtered data has been saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter CNVR types and save the results.')
    parser.add_argument('input_file', type=str, help='Path to the input file')
    parser.add_argument('output_file', type=str, help='Path to the output file')

    args = parser.parse_args()

    main(args.input_file, args.output_file)
