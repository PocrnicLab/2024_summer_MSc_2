import pandas as pd

# File path for the filtered SNPs
filtered_snp_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/cleaned_snp.txt"

# Read the filtered SNP file
filtered_snp_df = pd.read_csv(filtered_snp_file_path, sep='\t')

# Group by Sample ID and count the number of SNPs for each individual
snp_counts = filtered_snp_df.groupby('Sample ID').size()

# Print the SNP counts for each individual
print("SNP counts for each individual:")
print(snp_counts)

# Check if all individuals have the same number of SNPs
if snp_counts.nunique() == 1:
    print("\nAll individuals have the same number of SNPs.")
else:
    print("\nIndividuals have different numbers of SNPs.")

# Optionally, print individuals with their SNP counts
#print("\nDetailed SNP counts for each individual:")
#for sample_id, count in snp_counts.items():
#    print(f"Sample ID: {sample_id}, SNP Count: {count}")
