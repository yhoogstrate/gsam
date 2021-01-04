#!/usr/bin/env R

# ---- load libraries ----

library(tidyverse)

# ---- just counts / sample ----

gsam.viii.rnaseq <- 'output/tables/v3_extract_readcounts.txt' %>%
  read.table(header=T,stringsAsFactor=F) %>%
  dplyr::mutate(sample = gsub("/.+$","",gsub("^.+/alignments/","", sample) ) ) %>%
  dplyr::mutate(sid = gsub("-replicate", "", sample, fixed=T)) %>%
  dplyr::mutate(n.reads = vIII.reads.v3 + wt.reads.v3) %>%
  dplyr::mutate(egfrviii.pct = vIII.reads.v3 / n.reads * 100) %>%
  dplyr::mutate(egfrviii.pct = ifelse(n.reads < 10, NA, egfrviii.pct)) %>%
  dplyr::rename(vIII.reads =  vIII.reads.v3) %>%
  dplyr::rename(wt.reads = wt.reads.v3)


# ---- rotated table; per patient ----

# not updated recently? x-check for old stuff??
# commented out to find possible issues

# gsam.viii.rnaseq.rotated <- read.delim("data/gsam/output/tables/v3_extract_readcounts_rotated.txt")
# rownames(gsam.viii.rnaseq.rotated) <- gsam.viii.rnaseq.rotated$sid
# gsam.viii.rnaseq.rotated$resection.1.sum <- gsam.viii.rnaseq.rotated$resection.1.wt + gsam.viii.rnaseq.rotated$resection.1.v3
# gsam.viii.rnaseq.rotated$resection.2.sum <- gsam.viii.rnaseq.rotated$resection.2.wt + gsam.viii.rnaseq.rotated$resection.2.v3
# gsam.viii.rnaseq.rotated <- gsam.viii.rnaseq.rotated[,!colnames(gsam.viii.rnaseq.rotated) %in% c("X","X.1","sid")]
# gsam.viii.rnaseq.rotated <- gsam.viii.rnaseq.rotated[,order(colnames(gsam.viii.rnaseq.rotated))]
# 
# gsam.viii.rnaseq.rotated$resection.1.pos  <- gsam.viii.rnaseq.rotated$resection.1.v3 > 0
# gsam.viii.rnaseq.rotated$resection.2.pos  <- gsam.viii.rnaseq.rotated$resection.2.v3 > 0
# 
# gsam.viii.rnaseq.rotated$v3.stat <- NA
# gsam.viii.rnaseq.rotated[!is.na(gsam.viii.rnaseq.rotated$resection.1.pos) & !is.na(gsam.viii.rnaseq.rotated$resection.2.pos) & gsam.viii.rnaseq.rotated$resection.1.pos == FALSE & gsam.viii.rnaseq.rotated$resection.2.pos == FALSE,]$v3.stat <- "none"
# gsam.viii.rnaseq.rotated[!is.na(gsam.viii.rnaseq.rotated$resection.1.pos) & !is.na(gsam.viii.rnaseq.rotated$resection.2.pos) & gsam.viii.rnaseq.rotated$resection.1.pos == TRUE & gsam.viii.rnaseq.rotated$resection.2.pos == TRUE,]$v3.stat <- "both"
# gsam.viii.rnaseq.rotated[!is.na(gsam.viii.rnaseq.rotated$resection.1.pos) & !is.na(gsam.viii.rnaseq.rotated$resection.2.pos) & gsam.viii.rnaseq.rotated$resection.1.pos == TRUE & gsam.viii.rnaseq.rotated$resection.2.pos == FALSE,]$v3.stat <- "resection1"
# gsam.viii.rnaseq.rotated[!is.na(gsam.viii.rnaseq.rotated$resection.1.pos) & !is.na(gsam.viii.rnaseq.rotated$resection.2.pos) & gsam.viii.rnaseq.rotated$resection.1.pos == FALSE & gsam.viii.rnaseq.rotated$resection.2.pos == TRUE,]$v3.stat <- "resection2"
# gsam.viii.rnaseq.rotated$v3.stat <- as.factor(gsam.viii.rnaseq.rotated$v3.stat)
# 
# gsam.viii.rnaseq.rotated$delta.percentage <- (gsam.viii.rnaseq.rotated$resection.2.v3 / gsam.viii.rnaseq.rotated$resection.2.sum * 100.0) - (gsam.viii.rnaseq.rotated$resection.1.v3 / gsam.viii.rnaseq.rotated$resection.1.sum * 100.0)
# gsam.viii.rnaseq.rotated$resection.1.p <- gsam.viii.rnaseq.rotated$resection.1.v3 / gsam.viii.rnaseq.rotated$resection.1.sum * 100
# gsam.viii.rnaseq.rotated$resection.2.p <- (gsam.viii.rnaseq.rotated$resection.2.v3 / gsam.viii.rnaseq.rotated$resection.2.sum) * 100
# 


