#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))

## functions
##================================================================================
## overlap.idx
overlap.idx <- function(start, end, start1, end1)
{  
  idx1 <- start >= start1 & start <= end1
  idx2 <- end   >= start1 & end   <= end1
  idx3 <- start <= start1 & end   >= end1
  idx4 <- start >= start1 & end   <= end1
    
  return(which(idx1 | idx2 | idx3 | idx4))
}


## overlap.props
overlap.props <- function(start, end, start1, end1)
{  
  overlap.intersert <- pmax(pmin(end, end1) - pmax(start, start1) + 1, 0)
  overlap.union     <- pmax(end, end1) - pmin(start, start1) + 1
   
  return(overlap.intersert / overlap.union)
}


## initialize.cnvr
initialize.cnvr <- function(cnv, prop.cutoff = 0.3)
{
  cnv <- cnv[order(cnv$start, cnv$end), ]
  
  cnv$CNVR_ID <- NA
  cnv$CNVR_ID[1] <- paste0("CNVR_", 1)
  
  cnvr <- data.frame(CNVR_ID     = paste0("CNVR_", 1), 
                     outer.start = cnv$start[1], 
                     outer.end   = cnv$end[1],
                     nCNV        = 1,
                     stringsAsFactors = FALSE)
  
  idx <- 1
 
  for (i in 2:nrow(cnv)) {
    cnv.id    <- cnv$CNV_ID[i]
    cnv.start <- cnv$start[i]
    cnv.end   <- cnv$end[i]
    
    idxs.cnvr <- overlap.idx(start = cnvr$outer.start, end = cnvr$outer.end,
                             start1 = cnv.start, end1 = cnv.end)
    
    if (length(idxs.cnvr) == 0) {
      idx <- idx + 1
      cnvr1 <- data.frame(CNVR_ID     = paste0("CNVR_", idx), 
                          outer.start = cnv.start, 
                          outer.end   = cnv.end,
                          nCNV        = 1,
                          stringsAsFactors = FALSE)
      cnvr <- rbind(cnvr, cnvr1)
      cnv$CNVR_ID[i] <- paste0("CNVR_", idx)
    } else {
      dt.cnvr <- data.frame(CNVR_ID = cnvr$CNVR_ID[idxs.cnvr], 
                            prop    = 0,
                            stringsAsFactors = FALSE)
      for (k in 1:nrow(dt.cnvr)) {
        cnv1 <- cnv[which(cnv$CNVR_ID == dt.cnvr$CNVR_ID[k]), ]
        props <- overlap.props(start = cnv1$start, end = cnv1$end,
                               start1 = cnv.start, end1 = cnv.end)
        dt.cnvr$prop[k] <- sum(props) / nrow(cnv1)
      }
      if (any(dt.cnvr$prop >= prop.cutoff)) {
        idx.max <- which.max(dt.cnvr$prop)[1]
        cnv$CNVR_ID[i] <- dt.cnvr$CNVR_ID[idx.max]
        idx.mid <- which(cnvr$CNVR_ID == dt.cnvr$CNVR_ID[idx.max])
        cnvr$outer.start[idx.mid] <- min(cnvr$outer.start[idx.mid], cnv.start)
        cnvr$outer.end[idx.mid]   <- max(cnvr$outer.end[idx.mid],   cnv.end)
        cnvr$nCNV[idx.mid]        <- cnvr$nCNV[idx.mid] + 1
      } else {
        idx <- idx + 1
        cnv$CNVR_ID[i] <- paste0("CNVR_", idx)
        cnvr1 <- data.frame(CNVR_ID     = paste0("CNVR_", idx), 
                            outer.start = cnv.start,
                            outer.end   = cnv.end,
                            nCNV        = 1,
                            stringsAsFactors = FALSE)
        cnvr <- rbind(cnvr, cnvr1)
      }
    }
  }
  return(list(cnv = cnv, cnvr = cnvr))
}


## reassign.cnvr
reassign.cnvr <- function(cnv, cnvr, prop.cutoff = 0.3)
{  
  cnv  <- cnv[ order(cnv$start, cnv$end), ]
  cnvr <- cnvr[ order(cnvr$outer.start, cnvr$outer.end), ]
  
  cnvrs.sub <- cnvr$CNVR_ID[ which(cnvr$nCNV > 1) ]
  
  if (length(cnvrs.sub) == 0) {
    return(list(cnv = cnv, cnvr = cnvr))
  } else {
    cnvr.matrix <- list() 
    for (i in 1:length(cnvrs.sub)) {
      cnvr.id <- cnvrs.sub[i]
      cnv1 <- cnv[ which(cnv$CNVR_ID == cnvr.id), ]
      
      n1 <- nrow(cnv1)
      m1 <- matrix(0, nrow = n1, ncol = n1)
      rownames(m1) <- cnv1$CNV_ID
      colnames(m1) <- cnv1$CNV_ID
      for (j in 1:(n1 - 1)) {
        m1[j, (j+1):n1] <- overlap.props(start = cnv1$start[(j+1):n1], 
                                         end   = cnv1$end[(j+1):n1],
                                         start1 = cnv1$start[j], 
                                         end1   = cnv1$end[j])
      }
      m <- m1 + t(m1)
      cnvr.matrix[[ cnvr.id ]] <- m
      avg.props <- colSums(m) / (nrow(m) - 1)
      
      if (all(avg.props >= prop.cutoff)) {
        next
      } else {
        while (nrow(m) >= 2 & any(avg.props < prop.cutoff)) {
          idx.min <- which.min(avg.props)
          cnv$CNVR_ID[which(cnv$CNV_ID %in% cnv1$CNV_ID[idx.min])] <- NA
          idx.cnvr.update <- which(cnvr$CNVR_ID == cnvr.id)
          m <- m[-idx.min, -idx.min, drop = FALSE]
          if (nrow(m) > 0) {
            cnvr.matrix[[ cnvr.id ]] <- m
            avg.props <- colSums(m) / (nrow(m) - 1)
            cnv1 <- cnv1[-idx.min, ]  
            cnvr$outer.start[idx.cnvr.update] <- min(cnv1$start)
            cnvr$outer.end[idx.cnvr.update]   <- max(cnv1$end)
            cnvr$nCNV[idx.cnvr.update]        <- cnvr$nCNV[idx.cnvr.update] - length(idx.min)
          } else {
            cnvr <- cnvr[-idx.cnvr.update, ]
            cnvr.matrix[[ cnvr.id ]] <- NULL
          }
        }
      }
    }
    cnv.out <- cnv[ which(is.na(cnv$CNVR_ID)), ]   
    if (nrow(cnv.out) == 0) {
      return(list(cnv = cnv, cnvr = cnvr))
    } else {  
      for (k in 1:nrow(cnv.out)) {
        cnv.id    <- cnv.out$CNV_ID[k]
        cnv.start <- cnv.out$start[k]
        cnv.end   <- cnv.out$end[k]
        idxs.cnvr <- overlap.idx(start = cnvr$outer.start, end = cnvr$outer.end,
                                 start1 = cnv.start, end1 = cnv.end)
        if (length(idxs.cnvr) == 0) {
          next
        } else {  
          dt.cnvr <- data.frame(CNVR_ID = cnvr$CNVR_ID[idxs.cnvr], 
                                prop    = 0,
                                stringsAsFactors = FALSE)
          for (l in 1:nrow(dt.cnvr)) {
            dt.cnv <- cnv[ which(cnv$CNVR_ID == dt.cnvr$CNVR_ID[l]), ]
            props <- overlap.props(start = dt.cnv$start, end = dt.cnv$end,
                                   start1 = cnv.start, end1 = cnv.end)
            dt.cnvr$prop[l] <- sum(props) / nrow(dt.cnv)
          }
          if (any(dt.cnvr$prop >= prop.cutoff)) {
            idx.max <- which.max(dt.cnvr$prop)[1]
            cnvr.id.max <- dt.cnvr$CNVR_ID[idx.max]
            dt.cnv <- cnv[ which(cnv$CNVR_ID == cnvr.id.max), ]
            rownames(dt.cnv) <- dt.cnv$CNV_ID
            if (nrow(dt.cnv) == 1) {
              m <- matrix(0, 1, 1)
              rownames(m) <- dt.cnv$CNV_ID
              colnames(m) <- dt.cnv$CNV_ID
            } else {
              m <- cnvr.matrix[[ cnvr.id.max ]]
              dt.cnv <- dt.cnv[ rownames(m), ]
            }
            props1 <- overlap.props(start = dt.cnv$start, end = dt.cnv$end,
                                    start1 = cnv.start, end1 = cnv.end)
            n1 <- nrow(m)
            m <- cbind(rbind(m, props1), c(props1, 0))
            rownames(m)[n1 + 1] <- cnv.id
            colnames(m)[n1 + 1] <- cnv.id
            avg.props <- colSums(m) / (nrow(m) - 1)
            if (all(avg.props >= prop.cutoff)) {
              cnv$CNVR_ID[which(cnv$CNV_ID == cnv.id)] <- cnvr.id.max
              cnvr.matrix[[ cnvr.id.max ]] <- m
              idx.mid <- which(cnvr$CNVR_ID == cnvr.id.max)
              cnvr$outer.start[idx.mid] <- min(cnvr$outer.start[idx.mid], cnv.start)
              cnvr$outer.end[idx.mid]   <- max(cnvr$outer.end[idx.mid],   cnv.end)
              cnvr$nCNV[idx.mid]        <- cnvr$nCNV[idx.mid] + 1
            }
          }
        }
      }
      return(list(cnv = cnv, cnvr = cnvr))
    }
  }
}


## recheck.singleton.cnvr
recheck.singleton.cnvr <- function(cnv, cnvr, prop.cutoff = 0.3)
{
  cnv  <- cnv[ order(cnv$start, cnv$end), ]
  cnvr <- cnvr[ order(cnvr$outer.start, cnvr$outer.end), ]
  cnvrs.sgl <- cnvr$CNVR_ID[ which(cnvr$nCNV == 1) ]
  if (length(cnvrs.sgl) == 0) {
    return(list(cnv = cnv, cnvr = cnvr))
  } else {
    cnv.sgl <- cnv[cnv$CNVR_ID %in% cnvrs.sgl, ]
  }
  for (i in 1:nrow(cnv.sgl)) {
    cnv.id    <- cnv.sgl$CNV_ID[i]
    cnv.start <- cnv.sgl$start[i]
    cnv.end   <- cnv.sgl$end[i] 	
    cnvr.id   <- cnv.sgl$CNVR_ID[i]
    cnvr1 <- cnvr[-which(cnvr$CNVR_ID == cnvr.id), ]
    idxs.cnvr <- overlap.idx(start = cnvr1$outer.start, end = cnvr1$outer.end,
                             start1 = cnv.start, end1 = cnv.end)
    if (length(idxs.cnvr) == 0) {
      next
    } else {  
      dt.cnvr <- data.frame(CNVR_ID = cnvr1$CNVR_ID[idxs.cnvr], 
                            prop    = 0,
                            stringsAsFactors = FALSE)
      for (k in 1:nrow(dt.cnvr)) {
        cnv1 <- cnv[ which(cnv$CNVR_ID == dt.cnvr$CNVR_ID[k]), ]
        props <- overlap.props(start = cnv1$start, end = cnv1$end,
                               start1 = cnv.start, end1 = cnv.end)
        dt.cnvr$prop[k] <- sum(props) / nrow(cnv1)
      }
      if (any(dt.cnvr$prop >= prop.cutoff)) {
        idx.max <- which.max(dt.cnvr$prop)[1]
        cnvr.id.max <- dt.cnvr$CNVR_ID[idx.max]
        cnv1 <- cnv[ which(cnv$CNVR_ID %in% c(cnvr.id.max, cnvr.id)), ]
        n1 <- nrow(cnv1)
        m1 <- matrix(0, nrow = n1, ncol = n1)
        for (j in 1:(n1 - 1)) {
          m1[j, (j+1):n1] <- overlap.props(start = cnv1$start[(j+1):n1], 
                                           end   = cnv1$end[(j+1):n1],
                                           start1 = cnv1$start[j], 
                                           end1   = cnv1$end[j])
        }
        m <- m1 + t(m1)
        avg.props <- colSums(m) / (nrow(m) - 1)
        if (all(avg.props >= prop.cutoff)) {
          idx.mid <- which(cnvr$CNVR_ID == cnvr.id.max)
          cnvr$outer.start[idx.mid] <- min(cnvr$outer.start[idx.mid], cnv.start) 
          cnvr$outer.end[idx.mid]   <- max(cnvr$outer.end[idx.mid],   cnv.end) 
          cnvr$nCNV[idx.mid]        <- cnvr$nCNV[idx.mid] + 1
          cnvr <- cnvr[-which(cnvr$CNVR_ID == cnvr.id), ]
          cnv$CNVR_ID[cnv$CNV_ID == cnv.id] <- cnvr.id.max
        }
      }
    }
  }
  return(list(cnv = cnv, cnvr = cnvr))
}


## create.cnvr.chr.arm
create.cnvr.chr.arm <- function(cnv, prop.cutoff = 0.3)
{  
  cnv.current <- cnv
  flag <- 0
  round <- 1
  cnv.new <- NULL
  cnvr.new <- NULL
  
  while (flag == 0) {
    create.out <- initialize.cnvr(cnv = cnv.current, prop.cutoff = prop.cutoff)
    reassign.out <- reassign.cnvr(cnv = create.out$cnv, cnvr = create.out$cnvr, prop.cutoff = prop.cutoff)
    cnv.round <- reassign.out$cnv
    idxs.part2 <- which(is.na(cnv.round$CNVR_ID))
    if (length(idxs.part2) == 0) {
        cnv.part1 <- cnv.round
        cnv.part1$CNVR_ID <- paste0(cnv.part1$CNVR_ID, "_r", round)
        cnv.part2 <- NULL        
    } else {
        cnv.part1 <- cnv.round[-idxs.part2, ]
        cnv.part1$CNVR_ID <- paste0(cnv.part1$CNVR_ID, "_r", round)
        cnv.part2 <- cnv.round[idxs.part2, ]
    }
    cnvr.round <- reassign.out$cnvr
    cnvr.round$CNVR_ID <- paste0(cnvr.round$CNVR_ID, "_r", round)
    cnv.new <- rbind(cnv.new, cnv.part1)
    cnvr.new <- rbind(cnvr.new, cnvr.round)
    if (is.null(cnv.part2)) {
      cnv.out <- cnv.new
      cnvr.out <- cnvr.new
      flag <- 1
    } else if (nrow(cnv.part2) == 1) {
      cnv.part2$CNVR_ID[1] <- paste0("CNVR_", 1, "_r", round + 1)
      cnv.out <- rbind(cnv.new, cnv.part2)
      cnvr1 <- data.frame(CNVR_ID     = paste0("CNVR_", 1, "_r", round + 1), 
                          outer.start = cnv.part2$start[1],
                          outer.end   = cnv.part2$end[1], 
                          nCNV        = 1,
                          stringsAsFactors = FALSE)
      cnvr.out <- rbind(cnvr.new, cnvr1)
      flag <- 1
    } else {  
      cnv.current <- cnv.part2
      round <- round + 1
    }
  }
  output <- recheck.singleton.cnvr(cnv = cnv.out, cnvr = cnvr.out, prop.cutoff = prop.cutoff)
  return(list(cnv = output$cnv, cnvr = output$cnvr))
}


## create.cnvr
create.cnvr <- function(cnv, centromere, prop.cutoff = 0.3)
{
  chrs <- sort(unique(cnv$chr))
  cnv.out  <- NULL
  cnvr.out <- NULL
	
  for (chr1 in chrs) {
    cat("chromsome:", chr1, "\n")
    chr.cen <- centromere$position[which(centromere$chr == chr1)]
    cnv.chr <- cnv[which(cnv$chr == chr1), ]

    if (nrow(cnv.chr) > 0) {
      cnv.chr.p <- cnv.chr[which(cnv.chr$posEnd < chr.cen), ]
      cat("p-arm: creating CNVR for", nrow(cnv.chr.p), "CNVs ...\n")
      if (nrow(cnv.chr.p) >= 2) {
        p.out <- create.cnvr.chr.arm(cnv = cnv.chr.p, prop.cutoff = prop.cutoff)
        p.out.cnv  <- p.out$cnv
        p.out.cnv$CNVR_ID  <- paste0(p.out.cnv$CNVR_ID, "_chr", chr1, "_p")
        p.out.cnvr <- p.out$cnvr
        p.out.cnvr$CNVR_ID <- paste0(p.out.cnvr$CNVR_ID, "_chr", chr1, "_p")
        p.out.cnvr$chr     <- chr1
        p.out.cnvr$arm     <- "p"
      } else if (nrow(cnv.chr.p) == 0) {
        p.out.cnv  <- NULL
        p.out.cnvr <- NULL
      } else if (nrow(cnv.chr.p) == 1) {
        p.out.cnv <- cnv.chr.p
        p.out.cnv$CNVR_ID <- paste0("CNVR_1_r1_chr", chr1, "_p")
        p.out.cnvr <- data.frame(CNVR_ID     = paste0("CNVR_1_r1_chr", chr1, "_p"), 
                                 outer.start = cnv.chr.p$start,
                                 outer.end   = cnv.chr.p$end,
                                 nCNV        = 1,
                                 chr = chr1,
                                 arm = "p",
                                 stringsAsFactors = FALSE)
      }
      cat("p-arm: done.\n")
  
      cnv.chr.q <- cnv.chr[which(cnv.chr$posStart > chr.cen), ]
      cat("q-arm: creating CNVR for", nrow(cnv.chr.q), "CNVs ...\n")
      if (nrow(cnv.chr.q) >= 2) {
        q.out <- create.cnvr.chr.arm(cnv = cnv.chr.q, prop.cutoff = prop.cutoff)
        q.out.cnv  <- q.out$cnv
        q.out.cnv$CNVR_ID  <- paste0(q.out.cnv$CNVR_ID, "_chr", chr1, "_q")
        q.out.cnvr <- q.out$cnvr
        q.out.cnvr$CNVR_ID <- paste0(q.out.cnvr$CNVR_ID, "_chr", chr1, "_q")
        q.out.cnvr$chr     <- chr1
        q.out.cnvr$arm     <- "q"
      } else if (nrow(cnv.chr.q) == 0) {
        q.out.cnv  <- NULL
        q.out.cnvr <- NULL
      } else if (nrow(cnv.chr.q) == 1) {
        q.out.cnv <- cnv.chr.q
        q.out.cnv$CNVR_ID <- paste0("CNVR_1_r1_chr", chr1, "_q")
        q.out.cnvr <- data.frame(CNVR_ID     = paste0("CNVR_1_r1_chr", chr1, "_q"), 
                                 outer.start = cnv.chr.q$start,
                                 outer.end   = cnv.chr.q$end,
                                 nCNV        = 1,
                                 chr = chr1,
                                 arm = "q",
                                 stringsAsFactors = FALSE)
      }
      cat("q-arm: done.\n")  
      cnv.out  <- rbind(cnv.out, p.out.cnv, q.out.cnv)
      cnvr.out <- rbind(cnvr.out, p.out.cnvr, q.out.cnvr)
    }
  }
  cnvr.n <- table(cnv.out$CNVR_ID)
  if (all(cnvr.n[cnvr.out$CNVR_ID] == cnvr.out$nCNV)) {
    cat("CNVR ID in cnv and cnvr tables are consistent!\n")
  } else {
    cat("CNVR ID in cnv and cnvr tables are NOT consistent!\n")
  }
  return(list(cnv = cnv.out, cnvr = cnvr.out))
}


## preprocess.CNV
preprocess.CNV <- function(cnv, snp, centromere, min.numSNP = 3)
{
  cnv$CNV_ID <- paste0(cnv$Sample_ID, "_chr", cnv$chr, "_", cnv$posStart, "_", cnv$posEnd,
                       "_CNV_type_", cnv$CNV_type, "_method_", cnv$method)
  if (nrow(cnv) != length(unique(cnv$CNV_ID))) {
    stop("The CNV_ID column value is not unqiue.")
  }
  
  cnv <- cnv[order(cnv$chr, cnv$posStart, cnv$posEnd), ]
  snp <- snp[order(snp$Chr, snp$Position, snp$Name), ]
  
  chrs <- sort(unique(cnv$chr))
  cnv.tmp <- NULL
  for (chr1 in chrs) {
  	cnv.chr <- cnv[which(cnv$chr == chr1), ]
  	snp.chr <- snp[which(snp$Chr == chr1), ]
  	cnv.chr$start <- match(cnv.chr$posStart, snp.chr$Position)
  	cnv.chr$end   <- match(cnv.chr$posEnd,   snp.chr$Position)
  	cnv.chr$snpStart = snp.chr$Name[cnv.chr$start]
  	cnv.chr$snpEnd   = snp.chr$Name[cnv.chr$end]
  	if (any(cnv.chr$start == 0) |  any(cnv.chr$end == 0)) {
  	  stop("Positions in SNP data and in CNV data are not consistency.")
  	}
    cnv.tmp <- rbind(cnv.tmp, cnv.chr)
  }
  cnv <- cnv.tmp
  idx.snp <- which(cnv$numSNP < min.numSNP)
  if (length(idx.snp) > 0) {
    cat(length(idx.snp), "CNVs with <", min.numSNP, "probes", "are excluded.\n")
    cnv <- cnv[-idx.snp, ]
  }
  for (chr1 in chrs) {
    chr.cen <- centromere$position[which(centromere$chr == chr1)]
    idx.cen <- which(cnv$chr==chr1 & cnv$posStart <= chr.cen & cnv$posEnd >= chr.cen)
    if (length(idx.cen) > 0) {
      cat("chr", chr1, ":", length(idx.cen), "CNVs crossing centromere are excluded.\n")
      cnv <- cnv[-idx.cen, ] 
    }  
  }
  cnv <- cnv[order(cnv$chr, cnv$start, cnv$end), ]
  return(cnv)
}


## compute.boundary
compute.boundary <- function(cnv, freq.cutoff = 0.5)
{
  knots.all <- sort(unique(c(cnv$start, cnv$end)))
  nk <- length(knots.all)
  idxs.count <- NULL  
  idxs.start <- match(cnv$start, knots.all)
  idxs.end   <- match(cnv$end,   knots.all) 
  for (i in 1:nrow(cnv)) {
    idxs.count <- c(idxs.count, idxs.start[i]:idxs.end[i])
  }
  idxs.freq <- table(idxs.count)[as.character(1:nk)] / nrow(cnv)
  idxs.candidate <- as.integer(names(idxs.freq)[idxs.freq >= freq.cutoff])
  start <- knots.all[min(idxs.candidate)]
  end   <- knots.all[max(idxs.candidate)]
  nc <- length(idxs.candidate)
  idxs.diff <- idxs.candidate[-1] - idxs.candidate[-nc]
  npeak <- sum(idxs.diff > 1) + 1      
  return(list(start = start, end = end, npeak = npeak))
}


## create.CNVR.boundary
create.CNVR.boundary <- function(cnv, cnvr, snp, freq.cutoff = 0.5)
{
  snp <- snp[order(snp$Chr, snp$Position, snp$Name), ]
  
  cnvr$posStart <- NA
  cnvr$posEnd   <- NA
  for (i in 1:nrow(cnvr)) {
    cnvr.id <- cnvr$CNVR_ID[i]
    chr     <- cnvr$chr[i]
    cnv1 <- subset(cnv, CNVR_ID == cnvr.id)
    cnv1 <- cnv1[order(cnv1$start, cnv1$end), ]
    bd <- compute.boundary(cnv = cnv1, freq.cutoff = 0.5)
    snp.chr <- subset(snp, Chr == chr)
    cnvr$posStart[i] <- snp.chr$Position[ bd$start ]
    cnvr$posEnd[i]   <- snp.chr$Position[ bd$end ]
    cnvr$start_snp[i]  <- snp.chr$Name[ bd$start ]
    cnvr$end_snp[i]    <- snp.chr$Name[ bd$end ]
    cnvr$nPeak[i]    <- bd$npeak
    if (i %% 100 == 0) cat("Idx:",i, "in total:", nrow(cnvr), "CNVR_ID:", cnvr.id, "\n")
  }
  return(cnvr)
}


## merge_IPQ_in_one_cnvr
merge_IPQ_in_one_cnvr <- function(cnv1)
{
  cnvr_id <- unique(cnv1$CNVR_ID)
  cnv1$alg <- "P"
  tbl1 <- table(cnv1$Sample_ID)
  idxs.single   <- which(tbl1 == 1)
  idxs.multiple <- which(tbl1 >= 2)
  
  if (length(idxs.multiple) == 0) {
    #cat("CNVR_ID:", cnvr_id, "has", length(idxs.single), 
    #    "samples with CNV called by single method.\n")
    res <- cnv1[, c("CNVR_ID", "chr", "Sample_ID", "CN", "alg")]
    res$inner_start <- cnv1$posStart
    res$inner_end   <- cnv1$posEnd
    res$outer_start <- cnv1$posStart
    res$outer_end   <- cnv1$posEnd
    
  } else {
    #cat("CNVR_ID:", cnvr_id, "has", length(idxs.single), 
    #    "samples with CNV called by single method.\n")
    if (length(idxs.single) > 0) {
      samples.single <- names(tbl1)[idxs.single]
      idx.tmp <- which(cnv1$Sample_ID %in% samples.single)
      res1 <- cnv1[idx.tmp, c("CNVR_ID", "chr", "Sample_ID", "CN", "alg")]
      res1$inner_start <- cnv1$posStart[idx.tmp]
      res1$inner_end   <- cnv1$posEnd[idx.tmp]
      res1$outer_start <- cnv1$posStart[idx.tmp]
      res1$outer_end   <- cnv1$posEnd[idx.tmp]
    } else { res1 <- NULL }
    #cat("CNVR_ID:", cnvr_id, "has", length(idxs.multiple), 
    #    "samples with CNV called by multiple methods.\n")
    samples.multiple <- names(tbl1)[idxs.multiple]
    res2 <- NULL
    for (i in 1:length(samples.multiple)) {
      sample.id <- samples.multiple[i]
      idxs1 <- which(cnv1$Sample_ID == sample.id)
      cns1  <- unique(cnv1$CN[idxs1])
      if (all(cns1 < 2)) {
        res2.1 <- data.frame(CNVR_ID = cnvr_id, 
                             chr     = unique(cnv1$chr), 
                             Sample_ID = sample.id, 
                             CN = min(cns1), 
                             alg = "P",
                             inner_start = max(cnv1$posStart[idxs1]),
                             inner_end   = min(cnv1$posEnd[idxs1]),
                             outer_start = min(cnv1$posStart[idxs1]),
                             outer_end   = max(cnv1$posEnd[idxs1]),
                             stringsAsFactors = FALSE)
        res2 <- rbind(res2, res2.1)
      } else if (all(cns1 > 2)) {
        res2.1 <- data.frame(CNVR_ID   = cnvr_id, 
                             chr       = unique(cnv1$chr), 
                             Sample_ID = sample.id, 
                             CN        = max(cns1), 
                             alg       = "P",
                             inner_start = max(cnv1$posStart[idxs1]),
                             inner_end   = min(cnv1$posEnd[idxs1]),
                             outer_start = min(cnv1$posStart[idxs1]),
                             outer_end   = max(cnv1$posEnd[idxs1]),
                             stringsAsFactors = FALSE)
        res2 <- rbind(res2, res2.1)
      } else { next }
    }
    res <- rbind(res1, res2)
  }
  if (all(table(res$Sample_ID) == 1) | is.null(res)) {
    return(res)
  } else {
    stop(paste0("ERROR in merging IPQ CNVs in CNVR: ", cnvr_id, ".\n"))
  }
}

merge_IPQ <- function(cnv, cnvr)
{
  cnv.new <- NULL
  cnvr.new <- cnvr
  
  for (i in 1:nrow(cnvr)) {
    cnvr_id <- cnvr$CNVR_ID[i]
    cnv.new1 <- merge_IPQ_in_one_cnvr(cnv1 = subset(cnv, CNVR_ID == cnvr_id))
    cnv.new <- rbind(cnv.new, cnv.new1)
    if (is.null(cnv.new1)) {
      idx <- which(cnvr.new$CNVR_ID == cnvr_id)
      if (length(idx) > 0) { cnvr.new <- cnvr.new[-idx, ] }
    }
    if (i %% 100 == 0) cat("i =", i, "\n")
  }
  
  tbl.freq <- table(cnv.new$CNVR_ID)
  cnvr.freq <- as.data.frame(tbl.freq)
  names(cnvr.freq)[1] <- "CNVR_ID"
  cnvr.freq$CNVR_ID <- as.character(cnvr.freq$CNVR_ID)
  cnvr.new <- merge(cnvr.new, cnvr.freq)
  
  return(list(cnv = cnv.new, cnvr = cnvr.new))
}


# optionparser ------------------------------------------------------------

option_list <- list(
  make_option(c("-p", "--pcnv"), action = "store", default = NA, type = "character",
              help = "cnv file from PennCNV"),
  make_option(c("-s", "--snp"), action = "store", default = NA, type = "character",
              help = "snp position file"),
  make_option(c("-c", "--centromere"), action = "store", default = NA, type = "character",
              help = "centromere information file"),
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = "result output path.")
)


opt = parse_args(OptionParser(option_list = option_list))

if (is.na(opt$pcnv) | is.na(opt$snp) | is.na(opt$centromere) | is.na(opt$output)) {
  stop("All arguments must be supplied. Type --help for details.\n")
}

path_output <- opt$output

# read in centromere (.rds format)
centromere <- read.delim(file = opt$centromere, as.is = TRUE)

# read in snp (use the raw name from PennCNV)
snp <- read.table(file = opt$snp, sep = "\t", comment.char = "",
                  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# read from .rds ----------------------------------------------------------
pcnv <- read.delim(file = opt$pcnv, as.is = TRUE)

cat("Total number of PennCNV CNV calls:", nrow(pcnv), "\n")

ecnv <- pcnv

cat("CNV preprocessing ...\n")
cnv <- preprocess.CNV(cnv = ecnv,
                      snp = snp,
                      centromere = centromere,
                      min.numSNP = 3)

## create CNVR
cat("Create CNVR ...\n")
output <- create.cnvr(cnv = cnv,
                      centromere = centromere,
                      prop.cutoff = 0.3)

cnv.create <- output$cnv
cnvr.create <- output$cnvr

write.table(cnv.create, 
            file = file.path(path_output, "cnv_create.txt"),
            quote = F, row.names = F, sep = "\t")
write.table(cnvr.create, 
            file = file.path(path_output, "cnvr_create.txt"),
            quote = F, row.names = F, sep = "\t")

## CNVR boundary 
cat("Calculate CNVR boundary ...\n")
cnvr.boundary <- create.CNVR.boundary(cnv = cnv.create, cnvr = cnvr.create, snp = snp, freq.cutoff = 0.5)

peak <- table(cnvr.boundary$nPeak)
print(peak)

write.table(cnvr.boundary, 
            file = file.path(path_output, "cnvr_boundary.txt"),
            quote = F, row.names = F, sep = "\t")

## Assign CNV calls from individual methods for each CNVR
cat("Assign CNV calls from individual methods to each CNVR ...\n")
cleanCNV <- merge_IPQ(cnv = cnv.create, cnvr = cnvr.boundary)

cnv.cleanCNV <- cleanCNV$cnv
cnvr.cleanCNV <- cleanCNV$cnvr

cat("CNV number:", nrow(cnv.create), nrow(cnv.cleanCNV), "\n")
cat("CNVR number:", nrow(cnvr.boundary), nrow(cnvr.cleanCNV), "\n")

write.table(cnv.cleanCNV, 
            file = file.path(path_output, "cnv_clean.txt"),
            quote = F, row.names = F, sep = "\t")
write.table(cnvr.cleanCNV, 
            file = file.path(path_output, "cnvr_clean.txt"),
            quote = F, row.names = F, sep = "\t")