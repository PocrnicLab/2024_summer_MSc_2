#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(readxl)
library(dplyr)

# Read chromosome data from Excel file
chromosome_data <- read_excel("/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/chromosome_size.xlsx")
# Read CNVR data from TXT file
cnvr_data <- read.csv("/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnvr_types.txt", sep="\t")

# Ensure chromosome numbers in CNVR data are numeric
cnvr_data <- cnvr_data %>%
  mutate(chr = as.numeric(gsub("chr", "", chr)))

# Set colors for CNVR types
cnvr_colors <- c("Gain" = "purple", "Loss" = "green", "Mixed" = "darkblue")

# Create chromosome length dataframe
chromosome_data <- chromosome_data %>%
  mutate(Chromosome = factor(Chromosome, levels = Chromosome))

# Create CNVR dataframe and add color column
cnvr_data <- cnvr_data %>%
  mutate(color = cnvr_colors[Type])

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
ggsave("cnvr_distribution.png", plot = p, width = 12, height = 10, units = "in", dpi = 300)

# Display the plot
print(p)
