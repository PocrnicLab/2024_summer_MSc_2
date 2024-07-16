#!/usr/bin/env python3

import pandas as pd
import argparse

def compute_statistics(data):
    N = len(data.dropna())
    mean = data.mean()
    minimum = data.min()
    maximum = data.max()
    sd = data.std()
    return N, mean, minimum, maximum, sd

def select_preferred_record(df):
    # Sort by idanim and parityC so that parityC==2 comes first
    df = df.sort_values(by=['idanim', 'parityC'], ascending=[True, False])
    # Drop duplicates, keeping the first occurrence (i.e., the one with parityC==2 if available)
    df = df.drop_duplicates(subset='idanim', keep='first')
    return df

def main(input_file, output_file):
    # Read the file into a pandas DataFrame, assuming the first row contains column names
    df = pd.read_csv(input_file, sep=' ', header=0)
    
    # Select preferred records for each individual
    df = select_preferred_record(df)
    
    # Get the list of columns to compute statistics for, excluding non-numeric columns
    numeric_columns = df.select_dtypes(include='number').columns.tolist()
    
    # Open the output file for writing
    with open(output_file, 'w') as f:
        f.write("Trait\tN\tMean\tMinimum\tMaximum\tSD\n")
        for trait in numeric_columns:
            N, mean, minimum, maximum, sd = compute_statistics(df[trait])
            f.write(f"{trait}\t{N}\t{mean}\t{minimum}\t{maximum}\t{sd}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute statistics for sheep phenotypes data.')
    parser.add_argument('input_file', type=str, help='Path to the input file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    
    args = parser.parse_args()
    
    main(args.input_file, args.output_file)
