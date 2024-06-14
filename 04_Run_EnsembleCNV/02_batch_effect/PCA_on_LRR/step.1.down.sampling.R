#!/usr/bin/env Rscript

# generate the list of all SNPs
# Input is SNP_pos.txt
# including 3 columns: Name, Chr, Position

suppressMessages({
  require(data.table, quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
file_snps <- args[1]   ## SNP_pos.txt
path_output <- args[2] ## path to save SNPs

## Reading the SNP data
dat_snps <- fread(input = file_snps)
dat_snps <- as.data.frame(dat_snps, stringsAsFactors = FALSE)
dat_snps <- subset(dat_snps, Chr %in% 1:26)

## Select all SNPs
snps.selected <- dat_snps$Name

write.table(snps.selected, file = file.path(path_output, "snps.down.sample.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)