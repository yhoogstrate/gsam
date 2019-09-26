#!/usr/bin/env R

# duplicate reads
gsam.metadata <- read.delim("output/tables/duplicate_reads_sambamba_stats.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(gsam.metadata)[2:ncol(gsam.metadata)] <- paste0("sambamba.",colnames(gsam.metadata)[2:ncol(gsam.metadata)])
gsam.metadata$pid <- gsub("^([A-Z]+)[0-9]+$","\\1",gsam.metadata$sample)
gsam.metadata$resection <- as.factor(gsub("^[A-Z]+([0-9]+)$","\\1",gsam.metadata$sample))

# low complexity reads
tmp <- read.delim("output/tables/low_complexity_reads.txt",header=T,sep="\t",stringsAsFactors = F)
tmp$percentage.low.complexity.reads <- as.numeric(gsub("%","",tmp$percentage.low.complexity.reads,fixed=T))
colnames(tmp)[2:ncol(tmp)] <- paste0("fastp.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp,by.x="sample.id", by.y = "sample")
rm(tmp)

# mapping stats STAR
tmp <- read.delim("output/tables/star_percentage_mapped.txt",header=T,sep="\t",stringsAsFactors = F)
tmp$percentage.uniquely.mapped <- as.numeric(gsub("%","",tmp$percentage.uniquely.mapped,fixed=T))
tmp$percentage.multimap <- as.numeric(gsub("%","", tmp$percentage.multimap,fixed=T))
colnames(tmp)[2:ncol(tmp)] <- paste0("star.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp,by.x="sample.id", by.y = "sample")
rm(tmp)

# featureCounts stats
tmp <- read.delim("output/tables/featureCounts_gsam_1st96.exon-level.txt.summary",header=T,sep="\t",stringsAsFactors = F,row.names=1)
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
tmp <- read.delim("output/tables/gc_content_rmse.txt")
colnames(tmp)[2:ncol(tmp)] <- paste0("gc.control.",colnames(tmp)[2:ncol(tmp)])
gsam.metadata <- merge(gsam.metadata, tmp, by.x="sample.id", by.y = "sample.id")
rm(tmp)

# egfr-v3 qPCR

