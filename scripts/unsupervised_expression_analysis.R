#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)

# ---- load data ----
source("scripts/R/ligands.R")
#source("scripts/R/gsam_metadata.R")

e <- read.delim("output/tables/featureCounts_gsam_1st96.exon-level.txt",stringsAsFactors = F,comment="#")

colnames(d) <- gsub("X.data.users.youri.mnt.neurogen.ro.gsam.RNA.alignments.","",colnames(d),fixed=T)
colnames(d) <- gsub(".sambamba.dedup.bam","",colnames(d),fixed=T)

rownames(d) <- d$Geneid
d$Geneid <- NULL
d$Chr <- NULL
d$Start <- NULL
d$End <- NULL
d$Strand <- NULL
d$Length <- NULL

plot(sort(colSums(d)), log="y")
#abline(h=log(2000000))
# min 2mljn reads to take the low-Q samples out
#d <- d[,colSums(d) > 2000000]
#plot(sort(log(colSums(d))))
#abline(h=log(2000000))

