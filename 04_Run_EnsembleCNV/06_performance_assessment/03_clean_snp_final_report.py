import pandas as pd

# File paths
cnvr_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnvr_types.txt"
snp_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/data/final_report.txt"
output_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/cleaned_snp.txt"

# Read CNVR file
cnvr_df = pd.read_csv(cnvr_file_path, sep='\t')

# Read SNP file, skipping the header
snp_df = pd.read_csv(snp_file_path, sep='\t', comment='[', skiprows=9)

# Filter out SNPs not on chromosomes 1-26
snp_df = snp_df[snp_df['Chr'].astype(str).isin([str(i) for i in range(1, 27)])]

# Initialize an empty DataFrame to store the filtered SNPs
filtered_snps_df = pd.DataFrame()

# Filter SNPs using vectorized operations
for _, cnvr in cnvr_df.iterrows():
    chr_match = snp_df['Chr'].astype(str) == str(cnvr['chr'])
    pos_in_range = (snp_df['Position'] >= cnvr['posStart']) & (snp_df['Position'] <= cnvr['posEnd'])
    filtered_snps = snp_df[chr_match & pos_in_range]
    filtered_snps_df = pd.concat([filtered_snps_df, filtered_snps])

# Drop duplicates if any SNPs fall within multiple CNVRs
filtered_snps_df = filtered_snps_df.drop_duplicates()

# Save the filtered SNPs to a new file
filtered_snps_df.to_csv(output_file_path, sep='\t', index=False)
