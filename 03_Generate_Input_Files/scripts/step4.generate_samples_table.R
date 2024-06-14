#!/usr/bin/env Rscript

# Load necessary library
library(data.table)

# Function to parse command line arguments
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2) {
        stop("Usage: script_name.R <input_file> <output_file>")
    }
    list(input_file = args[1], output_file = args[2])
}

# Parse command line arguments
args <- parse_args()
input_file <- args$input_file
output_file <- args$output_file

# Read the content of the file
lines <- readLines(input_file)

# Find the starting position of the [Data] section
data_start <- which(lines == "[Data]") + 1

# Check if [Data] section exists
if (length(data_start) == 0) {
    stop("[Data] section not found in the input file.")
}

# Read the data section
data_lines <- lines[data_start:length(lines)]

# Convert data lines to data.table
data <- fread(text = paste(data_lines, collapse = "\n"), sep = "\t", header = TRUE)

# Check if necessary columns exist
if (!all(c("Sample ID", "Chr") %in% colnames(data))) {
    stop("Input file does not contain required columns: Sample ID, Chr.")
}

# Rename columns for consistency
setnames(data, c("Sample ID", "Chr"), c("Sample.ID", "Chr"))

# Determine gender based on presence of Y chromosome for each sample
data$Gender <- ifelse(data$Chr == "Y", "Male", "Female")

# Summarize gender for each sample
samples_gender <- data[, .(Gender = ifelse(any(Chr == "Y"), "Male", "Female")), by = Sample.ID]

# Sort by Sample ID
samples_sorted <- samples_gender[order(Sample.ID)]

# Rename columns
setnames(samples_sorted, c("Sample.ID", "Gender"))

# Write the output to a file
fwrite(samples_sorted, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")

cat("Sample information has been successfully extracted and saved to", output_file, "\n")
