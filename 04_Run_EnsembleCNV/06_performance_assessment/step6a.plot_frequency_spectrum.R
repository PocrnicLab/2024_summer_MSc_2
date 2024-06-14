#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales) # For breaks_pretty
library(tidyr) # For replace_na and complete
library(optparse)

# Create argument parser
option_list <- list(
  make_option(c("-c", "--cnvr_types_path"), type = "character", default = NULL, 
              help = "Path to the final CNVR types file", metavar = "character"),
  make_option(c("-e", "--encoded_results_path"), type = "character", default = NULL, 
              help = "Path to the encoded results file", metavar = "character"),
  make_option(c("-d", "--output_data_path"), type = "character", default = NULL, 
              help = "Path to save the merged CNVR and MAC data file", metavar = "character"),
  make_option(c("-o", "--output_plot_path"), type = "character", default = NULL, 
              help = "Path to save the plot", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$cnvr_types_path) || is.null(opt$encoded_results_path) || is.null(opt$output_data_path) || is.null(opt$output_plot_path)) {
  print_help(opt_parser)
  stop("All arguments must be supplied (cnvr_types_path, encoded_results_path, output_data_path, output_plot_path).", call. = FALSE)
}

# Read data
read_data <- function(path) {
  tryCatch({
    data <- read.table(path, header = TRUE, sep = "\t")
    return(data)
  }, error = function(e) {
    stop("Error reading file: ", path)
  })
}

cnvr_types <- read_data(opt$cnvr_types_path)
encoded_results <- read_data(opt$encoded_results_path)

# Function to calculate Minor Allele Count (MAC) for each CNVR
calculate_mac <- function(data) {
  allele_counts <- table(data$Encoded_Value)
  if (length(allele_counts) > 1) {
    sorted_counts <- sort(allele_counts, decreasing = TRUE)
    return(sorted_counts[2]) # Return the second most frequent allele count (Minor Allele)
  } else {
    return(0) # If there's only one allele, return 0
  }
}

# Calculate MAC for each CNVR
cnvr_mac <- encoded_results %>%
  group_by(CNVR_ID) %>%
  summarize(MAC = calculate_mac(pick(everything())))

# Merge CNVR types and MAC data
cnvr_data <- merge(cnvr_types, cnvr_mac, by = "CNVR_ID")

# Calculate proportions of each CNVR type considering the MAC
cnvr_proportions <- cnvr_data %>%
  group_by(Type) %>%
  mutate(Type_Count = n()) %>% # Count the number of CNVRs for each type
  group_by(Type, MAC) %>%
  summarize(Proportion = n() / first(Type_Count), .groups = 'drop') %>% # Calculate the proportion within each type
  complete(Type, MAC, fill = list(Proportion = 0)) # Ensure all combinations are present

# Save the merged data including CNVR_ID
write.table(cnvr_data, opt$output_data_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Create a plot with improved readability and dynamic x-axis breaks
plot <- ggplot(cnvr_proportions, aes(x = MAC, y = Proportion, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = breaks_pretty(n = 15)) + # Dynamic breaks
  labs(x = "Minor Allele Count", y = "Proportions of CNVRs", fill = "Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  coord_cartesian(xlim = c(1, max(cnvr_proportions$MAC))) # Adjust x-axis limit

# Adjust the size and resolution of the image
tryCatch({
  ggsave(opt$output_plot_path, plot, width = 10, height = 6, dpi = 300)
}, error = function(e) {
  stop("Error saving the plot: ", opt$output_plot_path)
})

# Display the plot
print(plot)
