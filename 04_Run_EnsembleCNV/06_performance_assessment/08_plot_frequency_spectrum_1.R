# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales) # For breaks_pretty
library(tidyr) # For replace_na and complete

# Read data
cnvr_types <- read.table("/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/final_cnvr_types.txt", header = TRUE, sep = "\t")
encoded_results <- read.table("/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/encoded_results.txt", header = TRUE, sep = "\t")

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
write.table(cnvr_data, "/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/cnvr_mac_data_with_id_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

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
ggsave("/exports/cmvm/eddie/eb/groups/PocrnicLab/2024_ms_cnv/Wkdir/06_performance_assessment/cnvr_mac_plot_1.png", plot, width = 10, height = 6, dpi = 300)
