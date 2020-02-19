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
vIII.rot$pid <- rownames(vIII.rot)

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


# ---- add vIII RT-qPCR ----

# GSAM qPCRs were not done in duplo

tmp <- read.csv('data/RNA/Final_qPCR_EGFR_GSAM.csv',stringsAsFactors = F)
tmp <- tmp[,colnames(tmp) %in% c("patientID", "EGFRCt002", "vIIICt002", "recurrent_EGFRCt002", "recurrent_vIIICt002") ]

colnames(tmp) <- paste0("qPCR.",colnames(tmp))
tmp$qPCR.percentageEGFRvIII <- 100 - (1/(1 + 2 ^ (tmp$qPCR.EGFRCt002 - tmp$qPCR.vIIICt002)) * 100)
tmp$qPCR.percentageEGFRvIII[tmp$qPCR.vIIICt002 == 40] <- 0
tmp$qPCR.recurrent_percentageEGFRvIII <- 100 - (1/(1 + 2 ^ (tmp$qPCR.recurrent_EGFRCt002 - tmp$qPCR.recurrent_vIIICt002)) * 100)
tmp$qPCR.recurrent_percentageEGFRvIII[tmp$qPCR.recurrent_vIIICt002 == 40] <- 0
tmp$qPCR.delta_percentage <- tmp$qPCR.recurrent_percentageEGFRvIII - tmp$qPCR.percentageEGFRvIII

print(tmp)
vIII.rot <- merge(x=vIII.rot, y = tmp, by.x="pid" , by.y = "qPCR.patientID",
                  all.x = T, all.y=T)
#rm(tmp)
#dim(vIII.rot)



