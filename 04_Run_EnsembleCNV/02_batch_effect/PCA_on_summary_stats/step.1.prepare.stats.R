#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

path_penncnv <- args[1]
path_output  <- args[2]

suppressMessages({
  require(data.table)
})

# penncnv sample-level -----------------------------------------------------
cat("Processing PennCNV results ...\n")
dat_penncnv <- read.table(file = file.path(path_penncnv, "CNV.PennCNV_qc_new.txt"),
                          sep = "\t",
                          header = TRUE,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)
dat_penncnv$File <- gsub("\\.txt$", "", dat_penncnv$File, perl = TRUE) 
dat_penncnv$WF <- abs(dat_penncnv$WF)

fp <- c("LRR_SD", "BAF_SD", "BAF_drift", "WF", "NumCNV")
dat_penncnv <- dat_penncnv[, c("File", fp)]
names(dat_penncnv) <- c("Sample_ID", paste("PennCNV", fp, sep = "."))

dat_stats_penncnv <- dat_penncnv

write.table(dat_stats_penncnv, 
            file = file.path(path_output, "IPQ.stats.txt"),
            quote = F, row.names = F, sep = "\t")
cat("Done.\n")