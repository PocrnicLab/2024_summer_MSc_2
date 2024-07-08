#!/usr/bin/env python3

import sys

def process_snp_map(input_file_path, output_file_path):
    # Open and read file content
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    # Parse file content and extract required information
    output_lines = ['Name\tChr\tPosition']
    for line in lines:
        if line.strip() == "":
            continue
        parts = line.strip().split()
        if len(parts) < 4 or not parts[2].isdigit() or not parts[3].isdigit():
            continue
        name = parts[1]
        chromosome = parts[2]
        position = parts[3]
        output_lines.append(f"{name}\t{chromosome}\t{position}")

    # Write the result to the output file
    with open(output_file_path, 'w') as output_file:
        for output_line in output_lines:
            output_file.write(output_line + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file_path> <output_file_path>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    process_snp_map(input_file_path, output_file_path)
