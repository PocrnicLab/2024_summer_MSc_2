import pandas as pd

# File paths
final_cnvr_types_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnvr_types.txt"
cnv_clean_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/03_create_CNVR/cnv_create.txt"
final_cnv_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnv.txt"
final_cnv_filted_out_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnv_filted_out.txt"

# Read files
final_cnvr_types = pd.read_csv(final_cnvr_types_path, sep='\t')
cnv_clean = pd.read_csv(cnv_clean_path, sep='\t')

# Convert final_cnvr_types to a dictionary with keys as tuples (chr, posStart, posEnd) and values as CNVR_ID
cnvr_dict = {}
for _, row in final_cnvr_types.iterrows():
    key = (row['chr'], row['posStart'], row['posEnd'])
    cnvr_dict[key] = row['CNVR_ID']

# Initialize DataFrames for kept and filtered out CNVs
kept_cnv = []
filtered_out_cnv = []

# Optimize filtering
for _, cnv_row in cnv_clean.iterrows():
    matched = False
    for (chr_, pos_start, pos_end), cnvr_id in cnvr_dict.items():
        if (cnv_row['chr'] == chr_ and
            cnv_row['posStart'] >= pos_start and
            cnv_row['posEnd'] <= pos_end):
            kept_cnv.append(cnv_row)
            matched = True
            break
    if not matched:
        filtered_out_cnv.append(cnv_row)

# Convert lists to DataFrames
kept_cnv_df = pd.DataFrame(kept_cnv, columns=cnv_clean.columns)
filtered_out_cnv_df = pd.DataFrame(filtered_out_cnv, columns=cnv_clean.columns)

# Save results
kept_cnv_df.to_csv(final_cnv_path, sep='\t', index=False)
filtered_out_cnv_df.to_csv(final_cnv_filted_out_path, sep='\t', index=False)