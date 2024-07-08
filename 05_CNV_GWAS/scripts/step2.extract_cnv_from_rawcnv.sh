#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

# Define input and output file paths from command line arguments
input_file="$1"
output_file="$2"

# Write the header to the output file
echo -e "chr\tstart\tend\tsample id\tstate\tnum snps\tstart probe\tend probe" > "$output_file"

# Read the input file line by line
while IFS= read -r line; do
    # Extract fields
    chr=$(echo "$line" | awk '{print $1}' | cut -d':' -f1)
    range=$(echo "$line" | awk '{print $1}' | cut -d':' -f2)
    start=$(echo "$range" | cut -d'-' -f1)
    end=$(echo "$range" | cut -d'-' -f2)
    state=$(echo "$line" | grep -oP 'cn=\K\d+')
    ne_id=$(echo "$line" | grep -oP 'HR[0-9]+')
    
    # Ensure state and num_snps are treated as integers
    state=$(printf "%d" "$state")
    num_snps=$(echo "$line" | grep -oP 'numsnp=\K\d+')
    num_snps=$(printf "%d" "$num_snps")

    # Extract the start probe, and end probe
    start_probe=$(echo "$line" | grep -oP 'startsnp=\K[^\s]+')
    end_probe=$(echo "$line" | grep -oP 'endsnp=\K[^\s]+')

    # Write the formatted output to the file
    echo -e "$chr\t$start\t$end\t$ne_id\t$state\t$num_snps\t$start_probe\t$end_probe" >> "$output_file"
done < "$input_file"
