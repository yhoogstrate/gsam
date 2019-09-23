#!/usr/bin/env R

gsam.metadata <- read.delim("output/tables/v3_dedup.txt",header=T,stringsAsFactors=F)
gsam.metadata$sample <- gsub('^[^_]+_','',gsam.metadata$sample)
gsam.metadata$sample <- gsub('_.+$','',gsam.metadata$sample)

gsam.metadata$pid <- gsub("^([A-Z]+)[0-9]+$","\\1",gsam.metadata$sample)
gsam.metadata$resection <- gsub("^[A-Z]+([0-9]+)$","\\1",gsam.metadata$sample)


gsam.metadata <- gsam.metadata[order(gsam.metadata$resection, gsam.metadata$pid),]


#print(gsam.metadata)
