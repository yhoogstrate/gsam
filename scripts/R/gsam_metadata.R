#!/usr/bin/env R

library(dplyr) # for distinct() function


# ---- patient metadata ----
# three CNV samples samples not in metadata: "AMA" "AMA" "BAO" "BAO" "FAF" "FAF"
#gsam.patient.metadata <- read.csv('data/administratie/dbGSAM_PUBLIC_VERSION.csv',stringsAsFactors=F)
gsam.patient.metadata <- read.csv('data/administratie/GSAM_combined_clinical_molecular.csv',stringsAsFactors=F)
gsam.patient.metadata <- gsam.patient.metadata[order(gsam.patient.metadata$studyID),] # reorder

# ---- exome-seq CNV metadata ----
gsam.cnv.metadata <- read.delim('data/output/tables/cnv_copynumber-ratio.cnr_log2_all.txt',stringsAsFactors=F)
gsam.cnv.metadata <- gsam.cnv.metadata[,-c(1,2,3,4)] # remove columns with geneid, start, end etc.
gsam.cnv.metadata <- data.frame(
    cnv.table.id = colnames(head(gsam.cnv.metadata)),
    sid = gsub("\\.b[12]$","",colnames(head(gsam.cnv.metadata))),
    batch = as.factor(gsub("^[^\\.]+\\.","",colnames(head(gsam.cnv.metadata))))
    )

# Add CNV -> regular patient identifiers and some DNA metrics
# gender from this table is incomplete, add it from gsam.patient.metadata later on
tmp <- read.delim("data/DNA/sample codes sanger gsam.txt",stringsAsFactors=FALSE)
tmp$pid <- gsub("[1-2]$","",tmp$donor_ID)
tmp <- tmp[,match(c("donor_ID", "pid", "PD_ID", "donor_sex", "donor_age_at_diagnosis","Concentration.at.QC..ng.ul.","Volume.at.QC..ul.","Amount.at.QC..ug."), colnames(tmp))]
gsam.cnv.metadata <- merge(gsam.cnv.metadata,tmp,by.x="sid",by.y="PD_ID")
rm(tmp)

# previous genders are incomplete
gsam.cnv.metadata <- merge(gsam.cnv.metadata, gsam.patient.metadata, by.x = "pid",by.y = "studyID",all.x=TRUE)
gsam.cnv.metadata$donor_sex <- NULL # incomplete, the one from patient metadata = complete


# --- RNA-seq metadata ----

# duplicate reads
gsam.metadata <- read.delim("data/RNA/output/tables/duplicate_reads_sambamba_stats.txt",header=T,sep="\t",stringsAsFactors=F)
colnames(gsam.metadata)[2:ncol(gsam.metadata)] <- paste0("sambamba.",colnames(gsam.metadata)[2:ncol(gsam.metadata)])
gsam.metadata$pid <- gsub("^([A-Z]+)[0-9]+$","\\1",gsam.metadata$sample)
gsam.metadata$resection <- as.factor(gsub("^[A-Z]+([0-9]+)$","\\1",gsam.metadata$sample))


# low complexity reads
tmp <- read.delim("data/RNA/output/tables/low_complexity_reads.txt",header=T,sep="\t",stringsAsFactors = F)
tmp$percentage.low.complexity.reads <- as.numeric(gsub("%","",tmp$percentage.low.complexity.reads,fixed=T))
colnames(tmp)[2:ncol(tmp)] <- paste0("fastp.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp,by.x="sample.id", by.y = "sample")
rm(tmp)


# mapping stats STAR
tmp <- read.delim("data/RNA/output/tables/star_percentage_mapped.txt",header=T,sep="\t",stringsAsFactors = F)
tmp$percentage.uniquely.mapped <- as.numeric(gsub("%","",tmp$percentage.uniquely.mapped,fixed=T))
tmp$percentage.multimap <- as.numeric(gsub("%","", tmp$percentage.multimap,fixed=T))
colnames(tmp)[2:ncol(tmp)] <- paste0("star.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp,by.x="sample.id", by.y = "sample")
rm(tmp)


# featureCounts stats
tmp <- read.delim("data/RNA/output/tables/featureCounts_gsam_1st96.exon-level.txt.summary",header=T,sep="\t",stringsAsFactors = F,row.names=1)
colnames(tmp) <- gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(tmp))
tmp <- t(tmp) # transpose
tmp <- data.frame(tmp)
tmp$read.count <- as.numeric(rowSums(tmp))
tmp$sample <-  rownames(tmp)
tmp$percentage.exon.counts <- tmp$Assigned/ tmp$read.count * 100
sel <- colnames(tmp) != "sample"
colnames(tmp)[sel] <- paste0("featurecounts.",colnames(tmp)[sel])
gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample")
rm(tmp, sel)


# gc 
tmp <- read.delim("data/RNA/output/tables/gc_content_rmse.txt")
colnames(tmp)[2:ncol(tmp)] <- paste0("gc.control.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample.id")
rm(tmp)


# egfr-v3 (RNA-seq)
tmp <- read.delim("data/RNA/output/tables/v3_extract_readcounts.txt")
tmp$sample <- gsub("^[^_]+_([^_]+)_.+$","\\1",tmp$sample)
tmp$v3.percentage <- tmp$vIII.reads.v3 / (tmp$vIII.reads.v3 + tmp$wt.reads.v3) * 100.0
tmp[(tmp$vIII.reads.v3 + tmp$wt.reads.v3) <= 6,]$v3.percentage <- NA
colnames(tmp)[2:ncol(tmp)] <- paste0("v3.rnaseq.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample")
rm(tmp)



# all metadata, but per patient and not per resection (collapsed)
tmp <- read.delim("data/administratie/GSAM_combined_clinical_molecular.csv", sep="," , stringsAsFactors=F)
gsam.metadata <- merge(gsam.metadata, tmp, by.x="pid", by.y="studyID", all.x = T, all.y=F, no.dups=F)


# --- RNA-seq metadata [full] ----
# STAR alignment statistics + patient / sample identifiers
gsam.rna.metadata <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt.summary",stringsAsFactors = F,comment="#",row.names=1)
colnames(gsam.rna.metadata) <- gsub("^[^\\.]+\\.([^\\.]+)\\..+$","\\1",colnames(gsam.rna.metadata),fixed=F)
gsam.rna.metadata <- t(gsam.rna.metadata)
gsam.rna.metadata <- gsam.rna.metadata[,colSums(gsam.rna.metadata) != 0] # removal of columns with only 0
colnames(gsam.rna.metadata) <- paste0("STAR.",colnames(gsam.rna.metadata))
gsam.rna.metadata <- data.frame(gsam.rna.metadata)
gsam.rna.metadata$sid <- as.character(rownames(gsam.rna.metadata))
gsam.rna.metadata$pid <- as.factor(gsub("[0-9]$","",rownames(gsam.rna.metadata)))
gsam.rna.metadata$resection <- as.factor(gsub("^.+([0-9])$","r\\1",rownames(gsam.rna.metadata)))

# @ todo add batch [1/2]

## add IDH
#source("scripts/R/dna_idh.R")
#tmp <- dna_idh[dna_idh$duplicate == FALSE,]
#tmp$Sample <- NULL
#tmp$duplicate <- NULL
#print(dim(gsam.rna.metadata))
#gsam.rna.metadata <- merge(gsam.rna.metadata, tmp, by.x = "sid", by.y = "donor_ID", all.x = TRUE)
#print(dim(gsam.rna.metadata))








