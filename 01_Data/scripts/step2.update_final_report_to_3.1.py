#!/usr/bin/env python3

import pandas as pd
import sys
import io

def update_final_report(new_snp_map_file, final_report_file, output_file):
    # Read the new_snp_map file
    new_snp_map = pd.read_csv(new_snp_map_file, sep='\t')
    
    # Read the final_report file
    with open(final_report_file, 'r') as file:
        lines = file.readlines()
    
    # Find the starting index of the [Data] section
    data_start_index = lines.index('[Data]\n') + 1
    header_lines = lines[:data_start_index]
    data_lines = lines[data_start_index:]
    
    # Parse the [Data] section into a DataFrame
    final_report_data = pd.read_csv(io.StringIO(''.join(data_lines)), sep='\t')
    
    # Update Num SNPs value
    num_snps = new_snp_map['Name'].nunique()
    
    # Update the header with the new Num SNPs value
    for i, line in enumerate(header_lines):
        if line.startswith('Num SNPs'):
            header_lines[i] = f'Num SNPs\t{num_snps}\n'
    
    # Merge the two DataFrames
    updated_final_report_data = final_report_data[final_report_data['SNP Name'].isin(new_snp_map['Name'])]
    
    # Write the output file
    with open(output_file, 'w') as file:
        file.writelines(header_lines)
        updated_final_report_data.to_csv(file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python update_final_report.py <new_snp_map_file> <final_report_file> <output_file>")
        sys.exit(1)
    
    new_snp_map_file = sys.argv[1]
    final_report_file = sys.argv[2]
    output_file = sys.argv[3]
    
    update_final_report(new_snp_map_file, final_report_file, output_file)
