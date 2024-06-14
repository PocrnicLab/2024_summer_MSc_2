#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(plyr))

option_list <- list(
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving results.")
)

opt <- parse_args(OptionParser(option_list = option_list))

pars <- c(opt$resultpath)
if (any(is.na(pars))) {
  stop("All parameters must be supplied. (--help for detail)")
}

path_result <- opt$resultpath

# Combine refinement results
path_refine <- file.path(path_result, "res_refine")
folders_chr <- list.files(path = path_refine, pattern = "^chr[0-9]+$", full.names = TRUE)

res_refinement <- data.frame()
for (folder.chr in folders_chr) {
  
  chr <- gsub("^chr", "", basename(folder.chr), perl = TRUE)
  file.chr <- paste0("CNVR_refine_chr_", chr, "_detail.rds")
  path.chr.data <- file.path(folder.chr, "data")
  
  file_path <- file.path(path.chr.data, file.chr)
  if (file.exists(file_path)) {
    res.chr <- readRDS(file = file_path)
    res_refinement <- rbind(res_refinement, res.chr)
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# Merge CNVRs with identical boundaries after refinement
res_refinement$identicalID <- paste(res_refinement$Chr,
                                    res_refinement$snp.start.refine,
                                    res_refinement$snp.end.refine, sep = "___")

res_refinement_same <- subset(res_refinement, type.overlap.based.on.raw == "same")
res_refinement_refine <- subset(res_refinement, type.overlap.based.on.raw != "same")

res_refinement_refine <- subset(res_refinement_refine, 
                                !identicalID %in% res_refinement_same$identicalID)

# De-duplicate CNVR
res_refinement_refine <- res_refinement_refine[!duplicated(res_refinement_refine$identicalID), ]
cat("Number of CNVRs with refined boundaries:", nrow(res_refinement_refine), "\n")

cnvrID_refine_same <- res_refinement_same$CNVR_ID

# Clean
dat_cnvr_keep <- read.delim(file = file.path(path_result, "cnvr_keep.txt"), as.is = TRUE)
dat_cnvr_keep$identicalID <- paste(dat_cnvr_keep$chr,
                                   dat_cnvr_keep$start_snp,
                                   dat_cnvr_keep$end_snp, sep = "___")

dat_cnvr_refine <- read.delim(file = file.path(path_result, "cnvr_refine.txt"), as.is = TRUE)
dat_cnvr_refine$identicalID <- paste(dat_cnvr_refine$chr,
                                     dat_cnvr_refine$start_snp,
                                     dat_cnvr_refine$end_snp, sep = "___")

cnvrID_keep <- dat_cnvr_keep$CNVR_ID
cnvrID_keep_final <- union(cnvrID_keep, cnvrID_refine_same)

dat_cnvr_keep_after_refine <- rbind(dat_cnvr_keep, dat_cnvr_refine)
dat_cnvr_keep_after_refine <- subset(dat_cnvr_keep_after_refine, CNVR_ID %in% cnvrID_keep_final)

# CNVRs with refined boundaries
res_refinement_refine_clean <- subset(res_refinement_refine, 
                                      !identicalID %in% dat_cnvr_keep$identicalID)

# CNVRs to be regenotyped after updating boundary information
dat_cnvr_regt <- subset(dat_cnvr_refine, CNVR_ID %in% res_refinement_refine_clean$CNVR_ID)
dat_cnvr_regt <- rename(dat_cnvr_regt, 
                        c("posStart" = "posStart.round1",
                          "posEnd" = "posEnd.round1",
                          "start_snp" = "start_snp.round1",
                          "end_snp" = "end_snp.round1",
                          "batch" = "batch.round1",
                          "genotype" = "genotype.round1",
                          "Freq" = "Freq.round1",
                          "identicalID" = "identicalID.round1"))
dat_cnvr_regt <- merge(dat_cnvr_regt, 
                       res_refinement_refine_clean[, 
                         c("CNVR_ID", "identicalID", "snp.posStart.refine", "snp.posEnd.refine", 
                           "snp.start.refine", "snp.end.refine")],
                       by = "CNVR_ID")
stopifnot(nrow(dat_cnvr_regt) == nrow(res_refinement_refine_clean))

dat_cnvr_regt <- rename(dat_cnvr_regt,
                        c("snp.posStart.refine" = "posStart", 
                          "snp.posEnd.refine" = "posEnd", 
                          "snp.start.refine" = "start_snp", 
                          "snp.end.refine" = "end_snp"))

write.table(dat_cnvr_keep_after_refine, 
            file = file.path(path_result, "cnvr_kept_after_refine.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

write.table(res_refinement_refine_clean, 
            file = file.path(path_result, "cnvr_refined_after_refine.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

write.table(dat_cnvr_regt, 
            file = file.path(path_result, "cnvr_regenotype_after_refine.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
