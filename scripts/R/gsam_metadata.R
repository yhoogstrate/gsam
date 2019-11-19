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

## duplicate reads
#gsam.metadata <- read.delim("data/RNA/output/tables/duplicate_reads_sambamba_stats.txt",header=T,sep="\t",stringsAsFactors=F)
#colnames(gsam.metadata)[2:ncol(gsam.metadata)] <- paste0("sambamba.",colnames(gsam.metadata)[2:ncol(gsam.metadata)])
#gsam.metadata$pid <- gsub("^([A-Z]+)[0-9]+$","\\1",gsam.metadata$sample)
#gsam.metadata$resection <- as.factor(gsub("^[A-Z]+([0-9]+)$","\\1",gsam.metadata$sample))


## low complexity reads
#tmp <- read.delim("data/RNA/output/tables/low_complexity_reads.txt",header=T,sep="\t",stringsAsFactors = F)
#tmp$percentage.low.complexity.reads <- as.numeric(gsub("%","",tmp$percentage.low.complexity.reads,fixed=T))
#colnames(tmp)[2:ncol(tmp)] <- paste0("fastp.",colnames(tmp)[2:ncol(tmp)])
#gsam.metadata <- merge(gsam.metadata, tmp,by.x="sample.id", by.y = "sample")
#rm(tmp)


## mapping stats STAR
#tmp <- read.delim("data/RNA/output/tables/star_percentage_mapped.txt",header=T,sep="\t",stringsAsFactors = F)
#tmp$percentage.uniquely.mapped <- as.numeric(gsub("%","",tmp$percentage.uniquely.mapped,fixed=T))
#tmp$percentage.multimap <- as.numeric(gsub("%","", tmp$percentage.multimap,fixed=T))
#colnames(tmp)[2:ncol(tmp)] <- paste0("star.",colnames(tmp)[2:ncol(tmp)])
#gsam.metadata <- merge(gsam.metadata, tmp,by.x="sample.id", by.y = "sample")
#rm(tmp)


## featureCounts stats
#tmp <- read.delim("data/RNA/output/tables/featureCounts_gsam_1st96.exon-level.txt.summary",header=T,sep="\t",stringsAsFactors = F,row.names=1)
#colnames(tmp) <- gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(tmp))
#tmp <- t(tmp) # transpose
#tmp <- data.frame(tmp)
#tmp$read.count <- as.numeric(rowSums(tmp))
#tmp$sample <-  rownames(tmp)
#tmp$percentage.exon.counts <- tmp$Assigned/ tmp$read.count * 100
#sel <- colnames(tmp) != "sample"
#colnames(tmp)[sel] <- paste0("featurecounts.",colnames(tmp)[sel])
#gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample")
#rm(tmp, sel)


## gc 
#tmp <- read.delim("data/RNA/output/tables/gc_content_rmse.txt")
#colnames(tmp)[2:ncol(tmp)] <- paste0("gc.control.",colnames(tmp)[2:ncol(tmp)])
#gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample.id")
#rm(tmp)


## egfr-v3 (RNA-seq)
#tmp <- read.delim("data/RNA/output/tables/v3_extract_readcounts.txt")
#tmp$sample <- gsub("^[^_]+_([^_]+)_.+$","\\1",tmp$sample)
#tmp$v3.percentage <- tmp$vIII.reads.v3 / (tmp$vIII.reads.v3 + tmp$wt.reads.v3) * 100.0
#tmp[(tmp$vIII.reads.v3 + tmp$wt.reads.v3) <= 6,]$v3.percentage <- NA
#colnames(tmp)[2:ncol(tmp)] <- paste0("v3.rnaseq.",colnames(tmp)[2:ncol(tmp)])
#gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample")
#rm(tmp)



## all metadata, but per patient and not per resection (collapsed)
#tmp <- read.delim("data/administratie/GSAM_combined_clinical_molecular.csv", sep="," , stringsAsFactors=F)
#gsam.metadata <- merge(gsam.metadata, tmp, by.x="pid", by.y="studyID", all.x = T, all.y=F, no.dups=F)


# --- RNA-seq metadata [full] ----
# STAR alignment statistics + patient / sample identifiers
gsam.rna.metadata <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt.summary",stringsAsFactors = F,comment="#",row.names=1)
colnames(gsam.rna.metadata) <- gsub("^[^\\.]+\\.([^\\]+)\\.Aligned.sorted.+$","\\1",colnames(gsam.rna.metadata),fixed=F)
gsam.rna.metadata <- t(gsam.rna.metadata)
gsam.rna.metadata <- gsam.rna.metadata[,colSums(gsam.rna.metadata) != 0] # removal of columns with only 0
colnames(gsam.rna.metadata) <- paste0("STAR.",colnames(gsam.rna.metadata))
gsam.rna.metadata <- data.frame(gsam.rna.metadata)
gsam.rna.metadata$sid <- gsub(".","-",as.character(rownames(gsam.rna.metadata)),fixed=T)
gsam.rna.metadata$pid <- as.factor(  gsub("[0-9]$","",gsub(".replicate","",rownames(gsam.rna.metadata)))  )
gsam.rna.metadata$resection <- as.factor(gsub("^.+([0-9])$","r\\1",gsub(".replicate$","",rownames(gsam.rna.metadata))))
gsam.rna.metadata$pct.duplicate.STAR <- gsam.rna.metadata$STAR.Unassigned_Duplicate / rowSums(gsam.rna.metadata[,gsub("^(.....).+$","\\1",colnames(gsam.rna.metadata)) == "STAR."]) * 100
gsam.rna.metadata$pct.multimap.STAR <- gsam.rna.metadata$STAR.Unassigned_MultiMapping / rowSums(gsam.rna.metadata[,gsub("^(.....).+$","\\1",colnames(gsam.rna.metadata)) == "STAR."]) * 100
gsam.rna.metadata$pct.nofeature.STAR <- gsam.rna.metadata$STAR.Unassigned_NoFeatures / rowSums(gsam.rna.metadata[,gsub("^(.....).+$","\\1",colnames(gsam.rna.metadata)) == "STAR."]) * 100
gsam.rna.metadata$duplicate.fold.STAR <- 1 / (1 - (gsam.rna.metadata$pct.duplicate.STAR / 100)) # 75% duplicate means 4 fold duplication

# add chromosomal distribution of rRNA containing alternate loci
tmp <- read.delim("output/tables/qc/idxstats/samtools.indexstats.matrix.txt",stringsAsFactors=F,row.names=1)
tmp$ref.len <- NULL
tmp <- t(tmp)
# rowSums(tmp / rowSums(tmp) * 100) == 100
tmp <- tmp / rowSums(tmp) * 100
rownames(tmp) <- gsub(".samtools.idxstats.txt","",rownames(tmp),fixed=T)
tmp <- tmp[,colnames(tmp) == "chrUn_gl000220"]
tmp <- data.frame(pct.rRNA.by.chrUn.gl000220 = tmp)
tmp$sid <- gsub(".","-",rownames(tmp),fixed=T)
gsam.rna.metadata <- merge(gsam.rna.metadata , tmp, by.x = "sid" , by.y = 'sid' , all.x = T)
rm(tmp)
#plot(gsam.rna.metadata$pct.rRNA.by.chrUn.gl000220 , gsam.rna.metadata$pct.multimap.STAR)
#plot(gsam.rna.metadata$pct.rRNA.by.chrUn.gl000220 , gsam.rna.metadata$pct.duplicate.STAR)
#plot(gsam.rna.metadata$pct.rRNA.by.chrUn.gl000220 , gsam.rna.metadata$pct.nofeature.STAR)
#plot(gsam.rna.metadata$pct.nofeature.STAR , gsam.rna.metadata$pct.duplicate.STAR)




# tin.py qc metrics
tmp <- read.table('output/tables/qc/tin.py/tin.py.matrix.txt',stringsAsFactors=F,header=T)
tmp$sid <- gsub(".bam","",tmp$Bam_file,fixed=T)
rownames(tmp) <- tmp$sid
tmp$Bam_file <- NULL
gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x = "sid", by.y = "sid")




# vIII rna-seq counts
tmp <- read.table('output/tables/v3_extract_readcounts.txt',header=T,stringsAsFactor=F)
tmp$sample <- gsub("^.+/alignments/([^/]+)/.+$","\\1",tmp$sample)
rownames(tmp) <- tmp$sample


sel <- tmp$wt.reads.v3 + tmp$vIII.reads.v3 > 15
tmp$vIII.percentage <- NA
tmp$vIII.percentage[sel] <- tmp$vIII.reads.v3[sel] / (tmp$wt.reads.v3[sel] + tmp$vIII.reads.v3[sel]) * 100
gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x = "sid", by.y = "sample")




# vIII qPCR percentage
tmp <- read.csv('data/RNA/Final_qPCR_EGFR_GSAM.csv',stringsAsFactors = F)
tmp <- tmp[,colnames(tmp) %in% c('EGFRCt002', 'vIIICt002', 'recurrent_EGFRCt002', 'recurrent_vIIICt002','recurrent_patientID')]

tmp.1 <- tmp[,match( c( 'EGFRCt002', 'vIIICt002',  'recurrent_patientID' ) , colnames(tmp) ) ]
tmp.1$resection <- 'r1'
colnames(tmp.1)[colnames(tmp.1) == "EGFRCt002"] <- 'qPCR.ct.EGFR.wt'
colnames(tmp.1)[colnames(tmp.1) == "vIIICt002"] <- 'qPCR.ct.EGFR.vIII'

tmp.2 <- tmp[,match( c( 'recurrent_EGFRCt002', 'recurrent_vIIICt002','recurrent_patientID' ) , colnames(tmp) ) ]
tmp.2$resection <- 'r2'
colnames(tmp.2)[colnames(tmp.2) == "recurrent_EGFRCt002"] <- 'qPCR.ct.EGFR.wt'
colnames(tmp.2)[colnames(tmp.2) == "recurrent_vIIICt002"] <- 'qPCR.ct.EGFR.vIII'

tmp <- rbind(tmp.1, tmp.2)
tmp$resection <- as.factor(tmp$resection)
rm(tmp.1, tmp.2)

tmp$qPCR.percent.EGFR.vIII <- 100 - (1/(1 + 2 ^ (tmp$qPCR.ct.EGFR.wt - tmp$qPCR.ct.EGFR.vIII)) * 100)
tmp[tmp$qPCR.ct.EGFR.vIII >= 40,]$qPCR.percent.EGFR.vIII <- 0

gsam.rna.metadata$tmp.id <- paste0(gsam.rna.metadata$pid, ".", gsam.rna.metadata$resection)
tmp$tmp.id <- paste0(tmp$recurrent_patientID, '.', tmp$resection)
tmp$resection <- NULL
tmp$recurrent_patientID <- NULL

gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x="tmp.id" , by.y = "tmp.id" , all.x = T) 
gsam.rna.metadata$tmp.id <- NULL
# x-check replictes CAO1, FAB2, GAS2 => works




# blacklist by heavy DNA contamination
blacklist.too.low.assigned <- c("AKA1", "CAC1", "AAB2", "GAS1", "KAE1", "CCZ1", "GAO2", "JAN1", "BAU2", "EAV2", "AAP1", "AZE1", "HAF1", "GAM1", "HAG1", "BAX2", "EAN1", "CBQ1", "AAD2", "HAK1", "CBG2", "BAI2", "HAE1", "CDH1", "HAI1", "KAB2", "GAE1", "BAN1", "KAC1", "KAA1", "ABA1")
blacklist.heavy.dna.contamination  <- c("CAV1", "BAT2", "EBW1", "HAE1", "BAU1", "EBO1", "GAE1", "CDH1", "KAC2", "ABA1", "KAA1", "KAC1", "BAN1", "KAB2", "GAJ2", "HAI1")

gsam.rna.metadata$blacklist.too.low.assigned <- gsam.rna.metadata$sid %in% blacklist.too.low.assigned
gsam.rna.metadata$blacklist.heavy.dna.contamination <- gsam.rna.metadata$sid %in% blacklist.heavy.dna.contamination

rm(blacklist.too.low.assigned)
rm(blacklist.heavy.dna.contamination)


# add batches
tmp <- read.table('data/administratie/plate.layout.table.txt',stringsAsFactors=F,header=T)
gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x="sid" , by.y = "sid")
rm(tmp)
gsam.rna.metadata$plate <- as.factor(gsam.rna.metadata$plate)
gsam.rna.metadata$storage.box <- as.factor(gsam.rna.metadata$storage.box )



