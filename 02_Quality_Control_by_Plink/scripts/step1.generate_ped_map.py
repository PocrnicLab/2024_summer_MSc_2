#!/usr/bin/env python3

import pandas as pd
import argparse

def main(genotype_file, snp_map_file, ped_file, map_file):
    # Load SNP mapping data
    snp_map_df = pd.read_csv(snp_map_file, sep='\t')
    
    # Filter SNPs to include only those on chromosomes 1-26
    valid_chromosomes = [str(i) for i in range(1, 27)]
    snp_map_df = snp_map_df[snp_map_df['Chromosome'].isin(valid_chromosomes)]
    
    snp_map = snp_map_df[['Name', 'Chromosome', 'Position']]

    # Load genotype data
    genotype_df = pd.read_csv(genotype_file, sep='\t', skiprows=9)
    genotype_data = genotype_df[['Sample ID', 'SNP Name', 'Allele1 - Top', 'Allele2 - Top']]

    # Remove spaces from Sample ID and underscores
    genotype_data['Sample ID'] = genotype_data['Sample ID'].str.replace(' ', '').str.replace('_', '')

    # Merge genotype data with SNP map to avoid repetitive lookups
    merged_data = genotype_data.merge(snp_map, left_on='SNP Name', right_on='Name', how='left')

    # Generate .map file
    snp_map.to_csv(map_file, sep=' ', header=False, index=False, columns=['Chromosome', 'Name', 'Position'])

    # Generate .ped file
    with open(ped_file, 'w') as ped_f:
        sample_ids = genotype_data['Sample ID'].unique()
        for sample_id in sample_ids:
            # Write initial columns: family ID, individual ID, father ID, mother ID, sex (0 for unknown), phenotype (-9 for unknown)
            ped_f.write(f"{sample_id} {sample_id} 0 0 0 -9 ")

            # Filter sample genotypes once
            sample_genotypes = merged_data[merged_data['Sample ID'] == sample_id]

            # Use a dict for fast lookup
            snp_dict = sample_genotypes.set_index('SNP Name')[['Allele1 - Top', 'Allele2 - Top']].to_dict('index')

            for snp_name in snp_map['Name']:
                if snp_name in snp_dict:
                    alleles = snp_dict[snp_name]
                    ped_f.write(f"{alleles['Allele1 - Top']} {alleles['Allele2 - Top']} ")
                else:
                    ped_f.write("0 0 ")

            ped_f.write('\n')

    print("PLINK files successfully generated.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate PLINK files from genotype and SNP map data.")
    parser.add_argument('--genotype_file', type=str, required=True, help="Path to the genotype file.")
    parser.add_argument('--snp_map_file', type=str, required=True, help="Path to the SNP map file.")
    parser.add_argument('--ped_file', type=str, required=True, help="Path to the output .ped file.")
    parser.add_argument('--map_file', type=str, required=True, help="Path to the output .map file.")

    args = parser.parse_args()
    main(args.genotype_file, args.snp_map_file, args.ped_file, args.map_file)