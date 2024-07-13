#!/usr/bin/env python3

import pandas as pd
import argparse

def update_fam_file(fam_file_path, phenotype_file_path):
    # Read the .fam file
    fam_df = pd.read_csv(fam_file_path, sep=" ", header=None)
    fam_df.columns = ["FID", "IID", "PID", "MID", "Sex", "Pheno"]

    # Read the phenotype file
    phenotype_df = pd.read_csv(phenotype_file_path, sep="\t")
    phenotype_df.columns = ["SampleID", "Fam", "Sex", "Phenotype"]

    # Create a dictionary for phenotypes
    phenotype_dict = dict(zip(phenotype_df["SampleID"], phenotype_df["Phenotype"]))

    # Update the phenotype column in the .fam file only for matching samples
    fam_df["Pheno"] = fam_df.apply(
        lambda row: phenotype_dict[row["IID"]] if row["IID"] in phenotype_dict else row["Pheno"], axis=1
    )

    # Ensure that missing phenotypes remain as -9 (integer)
    fam_df["Pheno"] = fam_df["Pheno"].apply(lambda x: -9 if pd.isnull(x) else x)

    # Convert all values to string and specifically ensure -9 is an integer
    fam_df["Pheno"] = fam_df["Pheno"].astype(str).replace("-9.0", "-9")

    # Save the updated .fam file
    fam_df.to_csv(fam_file_path, sep=" ", header=False, index=False, float_format='%.0f')

    print("Phenotype update completed!")

def main():
    parser = argparse.ArgumentParser(description="Update phenotype in .fam file based on a phenotype file.")
    parser.add_argument("fam_file_path", type=str, help="Path to the .fam file")
    parser.add_argument("phenotype_file_path", type=str, help="Path to the phenotype file")
    args = parser.parse_args()

    update_fam_file(args.fam_file_path, args.phenotype_file_path)

if __name__ == "__main__":
    main()
