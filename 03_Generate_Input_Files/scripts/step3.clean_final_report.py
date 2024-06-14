#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def filter_data(final_report_path, qc_individuals_path, qc_snps_path, output_path):
    # Read file content
    with open(final_report_path, 'r') as file:
        lines = file.readlines()

    # Separate Header and Data sections
    header = []
    data_start_index = 0
    for i, line in enumerate(lines):
        header.append(line)
        if line.startswith('[Data]'):
            data_start_index = i + 1  # The next line after [Data] is the column names
            break

    # Column names
    columns = lines[data_start_index].strip().split("\t")
    data_lines = lines[data_start_index + 1:]

    # Read QC passed individual IDs
    with open(qc_individuals_path, 'r') as file:
        qc_individuals = set(line.strip().split()[0] for line in file)
    num_samples = len(qc_individuals)

    # Read QC passed SNPs
    with open(qc_snps_path, 'r') as file:
        qc_snps = set(line.strip() for line in file)
    num_snps = len(qc_snps)

    # Update Num SNPs and Num Samples values in the header
    new_header = []
    for line in header:
        if line.startswith("Num SNPs"):
            new_header.append(f"Num SNPs        {num_snps}\n")
        elif line.startswith("Num Samples"):
            new_header.append(f"Num Samples     {num_samples}\n")
        else:
            new_header.append(line)

    # Filter data
    filtered_data_lines = []
    for line in data_lines:
        parts = line.split("\t")
        sample_id = parts[0]
        snp_name = parts[3]
        if sample_id in qc_individuals and snp_name in qc_snps:
            filtered_data_lines.append(line)

    # Write filtered results
    with open(output_path, 'w') as file:
        file.writelines(new_header)
        file.write("\t".join(columns) + "\n")
        file.writelines(filtered_data_lines)

def main():
    parser = argparse.ArgumentParser(description='Process and filter final report data.')
    parser.add_argument('--final_report_path', required=True, help='Path to the final report file')
    parser.add_argument('--qc_individuals_path', required=True, help='Path to the QC passed individuals file')
    parser.add_argument('--qc_snps_path', required=True, help='Path to the QC passed SNPs file')
    parser.add_argument('--output_path', required=True, help='Path to save the filtered report')

    args = parser.parse_args()

    filter_data(args.final_report_path, args.qc_individuals_path, args.qc_snps_path, args.output_path)

if __name__ == "__main__":
    main()
