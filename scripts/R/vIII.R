#!/usr/bin/env R

# ---- just counts / sample ----

# these info should taken from gsam.rna.metadata ( scripts/R/gsam_metadata.R )
#vIII <- read.table('output/tables/v3_extract_readcounts.txt',header=T,stringsAsFactor=F)
#vIII$sample <- gsub("_.+$","",gsub("^[^_]+_","",vIII$sample))
#rownames(vIII) <- vIII$sample
#vIII$vIII.percentage <- vIII$vIII.reads.v3 / (vIII$vIII.reads.v3 + vIII$wt.reads.v3) * 100


# ---- rotated table; per patient ----

vIII.rot <- read.delim("output/tables/v3_extract_readcounts_rotated.txt")
rownames(vIII.rot) <- vIII.rot$sid
vIII.rot$resection.1.sum <- vIII.rot$resection.1.wt + vIII.rot$resection.1.v3
vIII.rot$resection.2.sum <- vIII.rot$resection.2.wt + vIII.rot$resection.2.v3
vIII.rot <- vIII.rot[,!colnames(vIII.rot) %in% c("X","X.1","sid")]
vIII.rot <- vIII.rot[,order(colnames(vIII.rot))]

vIII.rot$resection.1.pos  <- vIII.rot$resection.1.v3 > 0
vIII.rot$resection.2.pos  <- vIII.rot$resection.2.v3 > 0

vIII.rot$v3.stat <- NA
vIII.rot[!is.na(vIII.rot$resection.1.pos) & !is.na(vIII.rot$resection.2.pos) & vIII.rot$resection.1.pos == FALSE & vIII.rot$resection.2.pos == FALSE,]$v3.stat <- "none"
vIII.rot[!is.na(vIII.rot$resection.1.pos) & !is.na(vIII.rot$resection.2.pos) & vIII.rot$resection.1.pos == TRUE & vIII.rot$resection.2.pos == TRUE,]$v3.stat <- "both"
vIII.rot[!is.na(vIII.rot$resection.1.pos) & !is.na(vIII.rot$resection.2.pos) & vIII.rot$resection.1.pos == TRUE & vIII.rot$resection.2.pos == FALSE,]$v3.stat <- "resection1"
vIII.rot[!is.na(vIII.rot$resection.1.pos) & !is.na(vIII.rot$resection.2.pos) & vIII.rot$resection.1.pos == FALSE & vIII.rot$resection.2.pos == TRUE,]$v3.stat <- "resection2"
vIII.rot$v3.stat <- as.factor(vIII.rot$v3.stat)

vIII.rot$delta.percentage <- (vIII.rot$resection.2.v3 / vIII.rot$resection.2.sum * 100.0) - (vIII.rot$resection.1.v3 / vIII.rot$resection.1.sum * 100.0)
vIII.rot$resection.1.p <- vIII.rot$resection.1.v3 / vIII.rot$resection.1.sum * 100
vIII.rot$resection.2.p <- (vIII.rot$resection.2.v3 / vIII.rot$resection.2.sum) * 100

