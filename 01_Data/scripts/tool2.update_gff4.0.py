import re

# Define file paths
bed_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/hglft_genome_2a394_539f30.bed"
gff_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/QTLdb_sheepOAR4.gff"
output_gff_file_path = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/updated_QTLdb_sheepOAR4.gff"

# Read BED file and parse data
bed_data = {}
with open(bed_file_path, 'r') as bed_file:
    for line in bed_file:
        fields = line.strip().split()
        chrom, start, end, qtl_id = fields[0], int(fields[1]), int(fields[2]), re.search(r'QTL_ID=(\d+)', fields[3]).group(1)
        bed_data[qtl_id] = (chrom, start, end)

# Read GFF file, update and write to a new GFF file
with open(gff_file_path, 'r') as gff_file, open(output_gff_file_path, 'w') as output_gff_file:
    for line in gff_file:
        if line.startswith('#'):
            output_gff_file.write(line)
            continue

        fields = line.strip().split('\t')
        
        # Get the last field which contains the attributes
        attributes = fields[-1]
        match = re.search(r'QTL_ID=(\d+)', attributes)
        
        if match:
            qtl_id = match.group(1)
            if qtl_id in bed_data:
                chrom, start, end = bed_data[qtl_id]
                fields[3] = str(start)
                fields[4] = str(end)
                output_gff_file.write('\t'.join(fields) + '\n')

# Confirm the script has completed
print(f"Updated GFF file has been saved to {output_gff_file_path}")
