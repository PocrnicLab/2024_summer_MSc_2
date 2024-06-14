import pandas as pd

# Define the input and output file paths
input_file = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnv.txt"
output_file = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnv_bed.txt"

# Read the input file
df = pd.read_csv(input_file, sep='\t')

# Extract the necessary columns and rename them for BED format
bed_df = df[['chr', 'posStart', 'posEnd', 'CNV_ID']]
bed_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name']

# Add the "chr" prefix to the chromosome column
bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)

# Write to the output file in BED format
bed_df.to_csv(output_file, sep='\t', header=False, index=False)

print(f"Conversion complete, output file saved as {output_file}")
