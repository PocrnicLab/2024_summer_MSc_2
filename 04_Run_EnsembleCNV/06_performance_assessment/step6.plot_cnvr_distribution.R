#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(optparse)

# Create argument parser
option_list <- list(
  make_option(c("-c", "--chromosome_data_path"), type = "character", default = NULL, 
              help = "Path to the chromosome size Excel file", metavar = "character"),
  make_option(c("-n", "--cnvr_data_path"), type = "character", default = NULL, 
              help = "Path to the final CNVR types file", metavar = "character"),
  make_option(c("-o", "--output_plot_path"), type = "character", default = NULL, 
              help = "Path to save the output plot", metavar = "character"),
  make_option(c("-s", "--output_stats_path"), type = "character", default = NULL, 
              help = "Path to save the output statistics file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$chromosome_data_path) || is.null(opt$cnvr_data_path) || is.null(opt$output_plot_path) || is.null(opt$output_stats_path)) {
  print_help(opt_parser)
  stop("All arguments must be supplied (chromosome_data_path, cnvr_data_path, output_plot_path, output_stats_path).", call. = FALSE)
}

# Function to read data safely
read_data <- function(path, type) {
  if (type == "excel") {
    tryCatch({
      data <- read_excel(path)
      return(data)
    }, error = function(e) {
      stop("Error reading Excel file: ", path)
    })
  } else if (type == "csv") {
    tryCatch({
      data <- read.csv(path, sep = "\t")
      return(data)
    }, error = function(e) {
      stop("Error reading CSV file: ", path)
    })
  }
}

# Read chromosome data from Excel file
chromosome_data <- read_data(opt$chromosome_data_path, "excel")

# Read CNVR data from TXT file
cnvr_data <- read_data(opt$cnvr_data_path, "csv")

# Ensure chromosome numbers in CNVR data are numeric
cnvr_data <- cnvr_data %>%
  mutate(chr = as.numeric(gsub("chr", "", chr)))

# Set colors for CNVR types
cnvr_colors <- c("Gain" = "red", "Loss" = "green", "Mixed" = "purple")

# Create chromosome length dataframe
chromosome_data <- chromosome_data %>%
  mutate(Chromosome = factor(Chromosome, levels = Chromosome))

# Create CNVR dataframe and add color column
cnvr_data <- cnvr_data %>%
  mutate(color = cnvr_colors[Type])

# Calculate CNVR statistics
cnvr_stats <- cnvr_data %>%
  group_by(chr) %>%
  summarise(
    Gain = sum(Type == "Gain"),
    Loss = sum(Type == "Loss"),
    Mixed = sum(Type == "Mixed"),
    Total = n()
  )

# Save CNVR statistics to a file
write.csv(cnvr_stats, opt$output_stats_path, row.names = FALSE)

# Plot the CNVR distribution
p <- ggplot() +
  geom_segment(data = chromosome_data, aes(x = 0, xend = `Size (bp)`, y = Chromosome, yend = Chromosome), color = "gray", linewidth = 5) +
  geom_segment(data = cnvr_data, aes(x = posStart, xend = posEnd, y = factor(chr, levels = chromosome_data$Chromosome), yend = factor(chr, levels = chromosome_data$Chromosome), color = Type), linewidth = 5) +
  scale_color_manual(values = cnvr_colors) +
  labs(x = "Position (bp)", y = "Chromosome", title = "CNVR Distribution Across Chromosomes", color = "CNVR Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "lightblue"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "lightblue"))

# Save the plot to a file
tryCatch({
  ggsave(opt$output_plot_path, plot = p, width = 12, height = 10, units = "in", dpi = 300)
}, error = function(e) {
  stop("Error saving the plot: ", opt$output_plot_path)
})
