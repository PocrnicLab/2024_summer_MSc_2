#!/usr/bin/env Rscript

# Load required libraries
library(CNVRanger)
library(readxl)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 5) {
  stop("Usage: script.R <cnv_loc> <phen_loc> <map_loc> <chrom_size_xlsx> <output_file>")
}

# Assign command line arguments to variables
cnv.loc <- args[1]
phen.loc <- args[2]
map.loc <- args[3]
chrom_size_xlsx <- args[4]
output_file <- args[5]

# Read the CNV calls file with tab delimiter
cnv.calls <- read.delim(cnv.loc, as.is=TRUE)
cnv.calls <- makeGRangesListFromDataFrame(cnv.calls, 
    split.field="sample.id", keep.extra.columns=TRUE)
sort(cnv.calls)

# Read the phenotype data
phen.df <- read.delim(phen.loc, as.is=TRUE)

# Create RaggedExperiment object
re <- RaggedExperiment(cnv.calls, colData=phen.df)

# Read the mapping data
map.df <- read.delim(map.loc, as.is=TRUE)

# Extract the directory from the output file path
output_dir <- dirname(output_file)

# Set up CNV GWAS
phen.info <- setupCnvGWAS("CNVRanger", cnv.out.loc=re, map.loc=map.loc, folder=output_dir)
all.paths <- phen.info$all.paths

# Perform CNV GWAS
segs.pvalue.gr <- cnvGWAS(phen.info, method.m.test = "fdr")

# Convert segs.pvalue.gr to a data frame
segs.pvalue.df <- as.data.frame(segs.pvalue.gr)

# Save segs.pvalue.df to a txt file
write.table(segs.pvalue.df, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)

# Read the Excel file with chromosome sizes
chromosome_data <- read_excel(chrom_size_xlsx)

# Filter data for chromosomes 1-26
chromosome_data <- chromosome_data %>%
  filter(as.numeric(Chromosome) %in% 1:26)

# Convert chromosome number to character and size to integer
simplified_data <- chromosome_data %>%
  mutate(Chromosome = as.character(Chromosome),
         Size_bp = as.integer(gsub(",", "", `Size (bp)`))) %>%
  select(Chromosome, Size_bp)

# Print the simplified data frame
print(simplified_data)

# Save as a CSV file
output_csv_path <- file.path(output_dir, "simplified_chromosome_data.csv")
write.csv(simplified_data, output_csv_path, row.names = FALSE)

# Define the chromosome order in the plot
order.chrs <- as.character(1:26)

# Read chromosome sizes from CSV file
chr.size.file <- output_csv_path
chr.sizes <- read.csv(chr.size.file, stringsAsFactors = FALSE)

# Create a data frame for chromosome order and sizes
chr.size.order <- data.frame(chr = order.chrs, sizes = chr.sizes$Size_bp, stringsAsFactors = FALSE)

# Plot a PDF file with a Manhattan plot within the 'Results' workfolder
plotManhattan(all.paths, segs.pvalue.gr, chr.size.order, plot.pdf = FALSE)
