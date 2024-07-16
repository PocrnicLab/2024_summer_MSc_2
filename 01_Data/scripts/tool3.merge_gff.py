import os

def read_gff(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def extract_qtl_ids(lines):
    qtl_ids = set()
    for line in lines:
        if not line.startswith("#") and line.strip():  # Skip headers and empty lines
            parts = line.split("\t")
            attributes = parts[-1]
            try:
                qtl_id = next(attr.split("=")[1] for attr in attributes.split(";") if attr.startswith("QTL_ID"))
                qtl_ids.add(qtl_id)
            except StopIteration:
                print(f"QTL_ID not found in line: {line.strip()}")
    return qtl_ids

def filter_new_entries(lines, existing_qtl_ids):
    new_entries = []
    for line in lines:
        if not line.startswith("#") and line.strip():  # Skip headers and empty lines
            parts = line.split("\t")
            attributes = parts[-1]
            try:
                qtl_id = next(attr.split("=")[1] for attr in attributes.split(";") if attr.startswith("QTL_ID"))
                if qtl_id not in existing_qtl_ids:
                    new_entries.append(line)
            except StopIteration:
                print(f"QTL_ID not found in line: {line.strip()}")
    return new_entries

def sort_key(line):
    parts = line.split("\t")
    chr_part = parts[0]
    if chr_part == "Chr.X":
        return (0, chr_part)
    elif chr_part.startswith("Chr."):
        try:
            return (1, int(chr_part.split(".")[1]))
        except ValueError:
            return (2, chr_part)
    else:
        return (2, chr_part)

def merge_and_write_files(file1_lines, new_entries, output_path):
    merged_lines = file1_lines + new_entries
    headers = [line for line in merged_lines if line.startswith("#")]
    data_lines = [line for line in merged_lines if not line.startswith("#")]

    sorted_data_lines = sorted(data_lines, key=sort_key)

    with open(output_path, 'w') as output_file:
        for line in headers:
            output_file.write(line)
        for line in sorted_data_lines:
            output_file.write(line)

# File paths
file1_path = '/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/Animal_QTLdb_release53_sheepOAR3.gff'
file2_path = '/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/updated_QTLdb_sheepOAR4.gff'
output_path = '/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/merged_QTLdb_sheepOAR3.gff'

# Read files
file1_lines = read_gff(file1_path)
file2_lines = read_gff(file2_path)

# Extract QTL_IDs from the first file
existing_qtl_ids = extract_qtl_ids(file1_lines)

# Filter new entries from the second file
new_entries = filter_new_entries(file2_lines, existing_qtl_ids)

# Merge and write to the new GFF file
merge_and_write_files(file1_lines, new_entries, output_path)

print(f'Merged file has been saved to: {output_path}')
