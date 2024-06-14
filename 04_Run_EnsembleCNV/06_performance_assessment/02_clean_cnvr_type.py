import pandas as pd

# File paths
input_file = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/cnvr_types.txt"
output_file = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnvr_types.txt"

# Read TXT file
df = pd.read_csv(input_file, sep='\t')

# Filter rows where Type is not 'Undefined'
filtered_df = df[df['Type'] != 'Undefined']

# Save the results to a new TXT file
filtered_df.to_csv(output_file, index=False, sep='\t')

print(f"Filtered data has been saved to {output_file}")
