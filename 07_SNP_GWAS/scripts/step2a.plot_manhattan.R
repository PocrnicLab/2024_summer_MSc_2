#!/usr/bin/env Rscript

# Load necessary packages
library(data.table)
library(qqman)
library(readxl)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 4) {
  stop("Usage: script.R <assoc_linear_file> <assoc_linear_adjusted_file> <chromosome_size_file> <output_file>")
}

# Assign arguments to variables
assoc_linear_file <- args[1]
assoc_linear_adjusted_file <- args[2]
chromosome_size_file <- args[3]
output_file <- args[4]

# Read data
assoc_linear <- fread(assoc_linear_file)
assoc_linear_adjusted <- fread(assoc_linear_adjusted_file)
chromosome_size <- read_excel(chromosome_size_file)

# Ensure that CHR columns are of the same type
assoc_linear$CHR <- as.character(assoc_linear$CHR)
assoc_linear_adjusted$CHR <- as.character(assoc_linear_adjusted$CHR)

# Merge data based on SNP column
merged_data <- merge(assoc_linear, assoc_linear_adjusted, by = "SNP")

# Check the structure of merged data
# str(merged_data)

# Prepare data for Manhattan plot
manhattan_data <- data.frame(
  CHR = as.numeric(merged_data$CHR.x),  # Use CHR from the first dataset
  BP = as.numeric(merged_data$BP),
  SNP = merged_data$SNP,
  P = as.numeric(merged_data$FDR_BH)
)

# Remove rows with NA values in CHR or P
manhattan_data <- na.omit(manhattan_data)

# Ensure chromosome labels are character vector
chromosome_labels <- as.character(chromosome_size$Chromosome)

# Plot Manhattan plot
png(output_file, width = 12, height = 8, units = "in", res = 300)
manhattan(manhattan_data, main = "Manhattan Plot", col = c("blue4", "orange3"), chrlabs = chromosome_labels)

# Add horizontal lines for p-value thresholds 0.05 and 0.1
abline(h = -log10(0.05), col = "red", lty = 2)  # p-value threshold 0.05
abline(h = -log10(0.1), col = "green", lty = 2)  # p-value threshold 0.1

dev.off()

