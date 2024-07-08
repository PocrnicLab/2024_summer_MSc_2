#!/usr/bin/env python3

import pandas as pd
import sys

def main():
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python script.py <phenotype_file> <output_file> <phenotype_column> [<sample_file>]")
        sys.exit(1)

    phenotype_file = sys.argv[1]
    output_file = sys.argv[2]
    phenotype_column = sys.argv[3]
    sample_file = sys.argv[4] if len(sys.argv) == 5 else None

    # Read input files
    phenotype_df = pd.read_csv(phenotype_file, delim_whitespace=True)

    if sample_file:
        sample_df = pd.read_csv(sample_file, header=None, names=['Sample_ID'])
        # Filter phenotype data to include only samples present in the sample file
        filtered_df = phenotype_df[phenotype_df['idanim'].isin(sample_df['Sample_ID'])]
    else:
        filtered_df = phenotype_df

    # Check if the phenotype column exists
    if phenotype_column not in filtered_df.columns:
        print(f"Phenotype column '{phenotype_column}' not found in phenotype file.")
        sys.exit(1)

    # Function to get the earliest record based on parityC
    def get_earliest_record(group):
        parity2_records = group[group['parityC'] == 2]
        if not parity2_records.empty:
            return parity2_records.iloc[0]
        parity1_records = group[group['parityC'] == 1]
        if not parity1_records.empty:
            return parity1_records.iloc[0]
        return group.iloc[0]

    # Apply the function to each idanim group
    result_df = filtered_df.groupby('idanim').apply(get_earliest_record).reset_index(drop=True)

    # Select only the idanim and phenotype columns
    result_df = result_df[['idanim', phenotype_column]]

    # Rename columns to match the output format
    result_df.columns = ['sample id', phenotype_column]

    # Insert new columns 'fam' and 'sex' with values set to 'NA'
    result_df.insert(1, 'fam', 'NA')
    result_df.insert(2, 'sex', 'NA')

    # Format the second column to two decimal places
    result_df[phenotype_column] = result_df[phenotype_column].map("{:.2f}".format)

    # Save the result to the output file with header
    result_df.to_csv(output_file, sep='\t', index=False, header=True)

if __name__ == "__main__":
    main()
