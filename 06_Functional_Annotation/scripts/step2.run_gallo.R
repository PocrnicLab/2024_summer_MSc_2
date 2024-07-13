#!/usr/bin/env Rscript

# Load required packages
library(GALLO)
library(data.table)
library(optparse)
library(dplyr)

# Define command line options
option_list <- list(
  make_option(c("-c", "--cnvr"), type = "character", default = NULL,
              help = "Path to the CNVR file", metavar = "character"),
  make_option(c("-q", "--qtl"), type = "character", default = NULL,
              help = "Path to the QTL file", metavar = "character"),
  make_option(c("-g", "--gene"), type = "character", default = NULL,
              help = "Path to the gene file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Directory to save output files", metavar = "character")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Function to load and validate a file
load_and_validate_file <- function(file_path, file_type) {
  if (file_type == "txt") {
    data <- fread(file_path, header = TRUE, sep = "\t")
    if (nrow(data) == 0) {
      stop(paste("File", file_path, "is empty or not loaded correctly"))
    }
  } else if (file_type %in% c("gff", "gtf")) {
    data <- import_gff_gtf(db_file = file_path, file_type = file_type)
    if (is.null(data)) {
      stop(paste("File", file_path, "is not loaded correctly"))
    }
  } else {
    stop("Unsupported file type")
  }
  return(data)
}

# Check if all file paths are provided
if (is.null(opt$cnvr) | is.null(opt$qtl) | is.null(opt$gene)) {
  stop("All file paths must be provided. Use --help for more information.")
}

# Load files
QTLwindows <- load_and_validate_file(opt$cnvr, "txt")

qtl.inp <- load_and_validate_file(opt$qtl, "gff")
qtf.inp <- load_and_validate_file(opt$gene, "gtf")

# Fix start_pos and end_pos if needed
qtl.inp <- as.data.table(qtl.inp)
qtl.inp[start_pos > end_pos, `:=`(start_pos = end_pos, end_pos = start_pos)]

# Running QTL annotation
out.qtls <- find_genes_qtls_around_markers(db_file = qtl.inp,
                                           marker_file = QTLwindows,
                                           method = "qtl",
                                           marker = "haplotype",
                                           interval = 500000)

# Filter to retain only rows where CHR and chr column values are equal
out.qtls <- out.qtls %>%
  filter(CHR == chr)

# Running gene annotation
out.genes <- find_genes_qtls_around_markers(db_file = qtf.inp,
                                            marker_file = QTLwindows,
                                            method = "gene",
                                            marker = "haplotype",
                                            interval = 500000)

# Define output file paths
out_qtls_path <- file.path(opt$output_dir, "output_qtls.txt")
out_genes_path <- file.path(opt$output_dir, "output_genes.txt")

# Save the outputs
write.table(out.qtls, out_qtls_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(out.genes, out_genes_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Plotting the percentage of each QTL class annotated and saving as JPEG
jpeg(file.path(opt$output_dir, "qtl_type_plot.jpg"), width = 1200, height = 800, res = 100)
par(mar=c(1,9,7,0))
plot_qtl_info(out.qtls, qtl_plot = "qtl_type", cex=2)
dev.off()

# Setting margin parameter to better fit the axis labels
QTL_classes <- unique(out.qtls$QTL_type)

for (c in QTL_classes) {
  tmp_file_name <- file.path(opt$output_dir, paste(c, ".png", sep = ""))
  jpeg(tmp_file_name, width = 1200, height = 800, res = 150)  # Use dynamically generated file name
  par(mar = c(6, 12, 3, 3))  # Adjust margins for better fit
  plot_qtl_info(out.qtls, qtl_plot = "qtl_name", qtl_class = c)  # Use current QTL class
  dev.off()
}

# Uncomment and modify below lines if enrichment analysis is needed
# out.enrich <- qtl_enrich(qtl_db = qtl.inp, 
#                          qtl_file = out.qtls, qtl_type = "Name",
#                          enrich_type = "genome", 
#                          padj = "fdr")
#
# Creating a new ID composed by the trait and the chromosome
# out.enrich$ID <- paste(out.enrich$QTL, " - ", "CHR", out.enrich$CHR, sep = "")
#
# Match the QTL classes and filtering the Reproduction related QTLs
# out.enrich.filtered <- out.enrich[which(out.enrich$adj.pval < 0.05),]
#
# Save the enrichment results
# out_enrich_path <- file.path(opt$output_dir, "output_enrich.txt")
# write.table(out.enrich, out_enrich_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#
# Plotting the enrichment results for the QTL enrichment analysis
# jpeg(file.path(opt$output_dir, "enrichment_plot.jpg"), width = 1600, height = 1000, res = 150)
# par(mar = c(10, 15, 5, 5))  # Adjust margins for better fit
# QTLenrich_plot(out.enrich.filtered, x = "ID", pval = "adj.pval")
# dev.off()
