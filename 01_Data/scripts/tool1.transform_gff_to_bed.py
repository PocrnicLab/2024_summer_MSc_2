import re

def gff_to_bed(gff_file, bed_file):
    with open(gff_file, 'r') as gff, open(bed_file, 'w') as bed:
        for line in gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0].replace('Chr.', '')
            start = fields[3]
            end = fields[4]
            attributes = fields[8]
            
            match = re.search(r'QTL_ID=(\d+)', attributes)
            if match:
                qtl_id = match.group(1)
                formatted_chrom = f"chr{chrom}"
                bed.write(f"{formatted_chrom}\t{start}\t{end}\tQTL_ID={qtl_id}\n")

gff_file = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/QTLdb_sheepOAR4.gff"
bed_file = "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/2024_summer_MSc_2/01_Data/data/QTLdb_sheepOAR4.bed"
gff_to_bed(gff_file, bed_file)
