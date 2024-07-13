#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import argparse

def read_file(file_path, sep, names=None, comment=None):
    try:
        return pd.read_csv(file_path, sep=sep, names=names, comment=comment)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def extract_name(attributes):
    for attr in attributes.split(';'):
        if attr.startswith('Name='):
            return attr.split('=')[1]
    return None

def main(qtls_file, gff_file, output_file):
    qtls_df = read_file(qtls_file, sep='\t')
    gff_cols = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes']
    gff_df = read_file(gff_file, sep='\t', names=gff_cols, comment='#')

    if qtls_df is not None and gff_df is not None:
        
        gff_df['Name'] = gff_df['Attributes'].apply(extract_name)
        db_counts = gff_df.groupby(['Name']).size().reset_index(name='Count')
        total_gff = db_counts['Count'].sum()
        db_counts['Proportion'] = db_counts['Count'] / total_gff

        if 'Name' in qtls_df.columns:

            bubble_df = qtls_df.copy()
            qtls_df = qtls_df.drop_duplicates(subset=['QTL_ID'])

            qtls_df['Name'] = qtls_df['Name'].astype(str)
            bubble_df['Name'] = bubble_df['Name'].astype(str)

            name_counts = qtls_df.groupby(['Name']).size().reset_index(name='Observed Number of QTLs')
            bubble_counts = bubble_df.groupby(['Name']).size().reset_index(name='Bubble Size')

            merged_df = pd.merge(db_counts, name_counts, left_on='Name', right_on='Name', how='inner').fillna(0)
            bubble_merged_df = pd.merge(merged_df, bubble_counts, left_on='Name', right_on='Name', how='inner').fillna(0)

            total_qtls = merged_df['Observed Number of QTLs'].sum()

            p_values = []
            richness_factors = []

            for index, row in merged_df.iterrows():
                observed_qtls = row['Observed Number of QTLs']
                k = observed_qtls
                N = total_qtls
                n = row['Proportion'] * gff_df.shape[0]
                M = gff_df.shape[0]
                p_value = hypergeom.sf(k-1, M, n, N)
                p_values.append(p_value)
                
                richness_factor = k / n if n > 0 else np.nan
                richness_factors.append(richness_factor)

            if len(p_values) == len(merged_df) and len(richness_factors) == len(merged_df):
                merged_df['P_value'] = p_values
                merged_df['Richness Factor'] = richness_factors
                merged_df['P_value'] = merged_df['P_value'].replace(0, 1e-300)
                _, fdr_corrected_p_values, _, _ = multipletests(merged_df['P_value'], method='fdr_bh')
                merged_df['FDR_P_value'] = fdr_corrected_p_values
                merged_df['-log10(FDR_P_value)'] = -np.log10(merged_df['FDR_P_value'])

                final_df = pd.merge(merged_df, bubble_merged_df[['Name', 'Bubble Size']], left_on='Name', right_on='Name', how='inner')
                
                plt.figure(figsize=(15, 10))
                scatter = plt.scatter(
                    final_df['Richness Factor'],
                    final_df['Name'],
                    s=final_df['Bubble Size'] * 100, # Size based on Bubble Size
                    c=final_df['-log10(FDR_P_value)'],
                    cmap='coolwarm',
                    alpha=0.6,
                    edgecolors='w',
                    linewidth=0.5
                )
                plt.xlabel('Richness Factor')
                plt.title('Enrichment Analysis Bubble Plot')

                for i in range(final_df.shape[0]):
                    plt.text(
                        final_df['Richness Factor'].iloc[i],
                        final_df['Name'].iloc[i],
                        str(final_df['Bubble Size'].iloc[i]), # Display number of QTLs without deduplication
                        fontsize=8,
                        ha='center',
                        va='center',
                        color='white'
                    )

                cbar = plt.colorbar(scatter)
                cbar.set_label('-log10(FDR_P_value)')

                plt.grid(True, axis='x')
                plt.savefig(output_file)  # Save the plot as a PNG file
            else:
                print("Error: Length of p_values or richness_factors does not match length of merged_df.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enrichment Analysis Bubble Plot")
    parser.add_argument("qtls_file", type=str, help="Path to the QTLs file")
    parser.add_argument("gff_file", type=str, help="Path to the GFF file")
    parser.add_argument("output_file", type=str, help="Path to the output PNG file")

    args = parser.parse_args()
    main(args.qtls_file, args.gff_file, args.output_file)
