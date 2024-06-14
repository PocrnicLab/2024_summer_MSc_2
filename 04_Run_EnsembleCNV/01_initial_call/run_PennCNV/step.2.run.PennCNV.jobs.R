#!/usr/bin/env Rscript

suppressMessages({
  require(optparse, quietly = TRUE)
})

# Define command-line options
option_list <- list(
  make_option(c("-p", "--penncnv"), action = "store", default = NA, type = "character",
              help = "Path to PennCNV installation folder."),  
  make_option(c("-a", "--data"), action = "store", default = NA, type = "character",
              help = "Path to tab-delimit text data files for each sample."),
  make_option(c("-d", "--wkdir"), action = "store", default = NA, type = "character",
              help = "Working directory."),
  make_option(c("-f", "--pfb"), action = "store", default = NA, type = "character",
              help = "PFB file."),
  make_option(c("-m", "--hmm"), action = "store", default = NA, type = "character",
              help = "HMM model file.")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

path_penncnv <- opt$penncnv
path_data    <- opt$data
path_wkdir   <- opt$wkdir
file_pfb     <- opt$pfb
file_hmm     <- opt$hmm

# Check if all parameters are supplied
if (any(is.na(c(path_data, path_wkdir, file_pfb, file_hmm)))) {
  stop("All parameters must be supplied. (--help for details)")
}

# Create necessary directories if they do not exist
path_list  <- file.path(path_wkdir, "list")
path_res   <- file.path(path_wkdir, "res")

if (!dir.exists(path_list)) {
  dir.create(path = path_list, showWarnings = FALSE, recursive = TRUE)
}
if (!dir.exists(path_res)) {
  dir.create(path = path_res, showWarnings = FALSE, recursive = TRUE)
}

# Generate list.txt for each sample
sample_files <- list.files(path = path_data)

cat("Number of samples:", length(sample_files), "\n")

for (i in 1:length(sample_files)) {
  sample_file <- sample_files[i]
  sample_list <- sub("\\.txt$", ".list", sample_file)
  
  dat1 <- data.frame(file_name = file.path(path_data, sample_file), ## Add whole path information
                     stringsAsFactors = FALSE)
  write.table(dat1, file = file.path(path_list, sample_list), sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Function to create PennCNV command
cmd_PennCNV <- function(file_hmm, file_pfb,
                        filename_sample, path_list, path_res_sample) {

  file_list <- file.path(path_list, sub("\\.txt$", ".list", filename_sample))

  samplename <- gsub(pattern = "\\.txt$", replacement = "", filename_sample)
  
  file_log   <- file.path(path_res_sample, paste0(samplename, ".log"))
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))

  cmd <- paste(file.path(path_penncnv, "detect_cnv.pl"),
               "-test --confidence",
               "-hmm", file_hmm,
               "-pfb", file_pfb,
               "-list", file_list,
               "-log", file_log,
               "-out", file_rawcnv,
               "-lastchr 26")  # Specify the last chromosome as 26
  cmd
}

# Function to check if SNP IDs in the sample file match those in the PFB file
check_snp_ids <- function(sample_file, pfb_file) {
  sample_data <- read.table(sample_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  pfb_data <- read.table(pfb_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  sample_ids <- sample_data$Name
  pfb_ids <- pfb_data$Name
  
  missing_ids <- setdiff(sample_ids, pfb_ids)
  
  if (length(missing_ids) > 0) {
    cat("Warning: The following SNP IDs are in the sample file but not in the PFB file:\n")
    print(missing_ids)
  } else {
    cat("All SNP IDs in the sample file are present in the PFB file.\n")
  }
}

# Function to print the first few rows of the sample file for inspection
print_sample_file_head <- function(sample_file) {
  sample_data <- read.table(sample_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("First few rows of the sample file:\n")
  print(head(sample_data))
}

# Main loop to process each sample
for (i in 1:length(sample_files)) {
  sample_file <- sample_files[i]
  samplename <- gsub(pattern = "\\.txt$", replacement = "", sample_file)

  path_res_sample <- file.path(path_res, samplename)
  dir.create(path = path_res_sample, showWarnings = FALSE, recursive = TRUE)

  cat("Sample_ID:", samplename, "\n")

  # Check SNP IDs before running PennCNV
  check_snp_ids(file.path(path_data, sample_file), file_pfb)
  
  # Print the first few rows of the sample file for inspection
  print_sample_file_head(file.path(path_data, sample_file))

  cmd.sample <- cmd_PennCNV(file_hmm = file_hmm,
                            file_pfb = file_pfb,
                            filename_sample = sample_file,
                            path_list = path_list,
                            path_res_sample = path_res_sample)

  cat("Running command:", cmd.sample, "\n")
  
  # Execute the command and capture the output
  result <- try(system(cmd.sample, intern = TRUE), silent = TRUE)
  
  if (inherits(result, "try-error")) {
    cat("Error in processing sample:", samplename, "\n")
    cat(result, "\n")
  } else {
    cat("Completed processing sample:", samplename, "\n")
    cat(result, "\n")  # Print the result for debugging purposes
  }
  
  Sys.sleep(0.1)
}