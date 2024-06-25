#!/usr/bin/env Rscript

# Load necessary libraries
if (!requireNamespace("tools", quietly = TRUE)) {
  install.packages("tools")
}
library(tools)

# Function to find all R scripts in the directory
find_R_scripts <- function(directory) {
  list.files(directory, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
}

# Function to extract package names from R scripts
extract_packages <- function(file) {
  lines <- readLines(file, warn = FALSE)
  packages <- unique(grep("^\\s*(library|require)\\((\"|')?\\w+(\"|')?\\)", lines, value = TRUE))
  packages <- gsub("^\\s*(library|require)\\((\"|')?(\\w+)(\"|')?\\).*", "\\3", packages)
  return(packages)
}

# Function to get the version of a package
get_package_version <- function(package) {
  info <- installed.packages()
  if (package %in% rownames(info)) {
    return(info[package, "Version"])
  } else {
    return(NA)
  }
}

# Main function to gather package information
gather_package_info <- function(directory) {
  scripts <- find_R_scripts(directory)
  package_info <- data.frame(Package = character(), Version = character(), stringsAsFactors = FALSE)
  
  for (script in scripts) {
    packages <- extract_packages(script)
    for (package in packages) {
      version <- get_package_version(package)
      package_info <- rbind(package_info, data.frame(Package = package, Version = version, stringsAsFactors = FALSE))
    }
  }
  
  package_info <- unique(package_info)
  return(package_info)
}

# Save package information to a file
save_package_info <- function(package_info, file) {
  write.csv(package_info, file, row.names = FALSE)
}

# Run the functions and save the output
directory <- "."  # Current directory
output_file <- "R_package_info.csv"
package_info <- gather_package_info(directory)
save_package_info(package_info, output_file)

cat("Package information has been saved to", output_file, "\n")
