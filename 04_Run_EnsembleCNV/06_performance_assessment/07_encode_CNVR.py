import pandas as pd

# Read the data files
cnvr_types_path = '/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnvr_types.txt'
cnv_path = '/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnv.txt'

cnvr_df = pd.read_csv(cnvr_types_path, sep='\t')
cnv_df = pd.read_csv(cnv_path, sep='\t')

# Filter out CNVRs of type 'Mixed'
cnvr_df = cnvr_df[cnvr_df['Type'] != 'Mixed']

# Create CNVR position dictionary
cnvr_dict = {(row['chr'], row['posStart'], row['posEnd']): row['CNVR_ID'] for _, row in cnvr_df.iterrows()}

# Find all unique CNVR_IDs and Sample_IDs
all_cnvr_ids = cnvr_df['CNVR_ID'].unique()
all_sample_ids = cnv_df['Sample_ID'].unique()

# Initialize result dictionary with all samples and CNVRs set to 0
encoded_results_dict = {sample: {cnvr: 0 for cnvr in all_cnvr_ids} for sample in all_sample_ids}

# Function to encode CN values
def encode_cn_value(cn_values):
    if set(cn_values).issubset({0, 1}):  # Both 0 and 1 should be treated as -1
        return -1
    if len(cn_values) == 1:
        cn_value = cn_values[0]
        if cn_value in [0, 1]:
            return -1
        elif cn_value == 2:
            return 0
        else:
            return 1
    return None  # Indicate inconsistency

# Vectorized operation to map CNVRs
cnv_df['encoded_value'] = cnv_df.groupby(['chr', 'posStart', 'posEnd'])['CN'].transform(lambda x: encode_cn_value(x.tolist()))

inconsistent_cnvr_samples = []

# Encode CNVRs
for _, cnvr_row in cnvr_df.iterrows():
    mask = (cnv_df['chr'] == cnvr_row['chr']) & (cnv_df['posStart'] >= cnvr_row['posStart']) & (cnv_df['posEnd'] <= cnvr_row['posEnd'])
    filtered_cnv_df = cnv_df[mask]

    if not filtered_cnv_df.empty:
        grouped = filtered_cnv_df.groupby('Sample_ID')
        for sample_id, group in grouped:
            encoded_value = encode_cn_value(group['CN'].unique())
            if encoded_value is not None:  # Only set if consistent
                encoded_results_dict[sample_id][cnvr_row['CNVR_ID']] = encoded_value

# Check consistency
for _, cnvr_row in cnvr_df.iterrows():
    mask = (cnv_df['chr'] == cnvr_row['chr']) & (cnv_df['posStart'] >= cnvr_row['posStart']) & (cnv_df['posEnd'] <= cnvr_row['posEnd'])
    filtered_cnv_df = cnv_df[mask]

    if not filtered_cnv_df.empty:
        grouped = filtered_cnv_df.groupby('Sample_ID')
        for sample_id, group in grouped:
            cn_values = group['CN'].unique()
            if len(cn_values) > 1 and set(cn_values).issubset({0, 1}):
                # Both 0 and 1 present, treat as -1
                encoded_results_dict[sample_id][cnvr_row['CNVR_ID']] = -1
            elif len(cn_values) > 1:
                inconsistent_cnvr_samples.append((cnvr_row['CNVR_ID'], sample_id))

# Save inconsistent CNVRs and samples to file
inconsistent_df = pd.DataFrame(inconsistent_cnvr_samples, columns=['CNVR_ID', 'Sample_ID'])
inconsistent_df.to_csv('inconsistent_cnvr_samples.txt', sep='\t', index=False)

# Transform the encoded results into a DataFrame
encoded_results = [(sample, cnvr, encoded_value) for sample, cnvrs in encoded_results_dict.items() for cnvr, encoded_value in cnvrs.items()]
encoded_df = pd.DataFrame(encoded_results, columns=['Sample_ID', 'CNVR_ID', 'Encoded_Value'])
encoded_df.to_csv('encoded_results.txt', sep='\t', index=False)

# Save the sample list to final_samples.txt
sample_list_df = pd.DataFrame(all_sample_ids, columns=['Sample_ID'])
sample_list_df.to_csv('final_samples.txt', sep='\t', index=False)

print("Processing completed. Check 'inconsistent_cnvr_samples.txt', 'encoded_results.txt', and 'final_samples.txt' for results.")