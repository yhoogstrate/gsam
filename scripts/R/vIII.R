#!/usr/bin/env R


vIII <- read.table('output/tables/v3_extract_readcounts.txt',header=T,stringsAsFactor=F)
vIII$sample <- gsub("_.+$","",gsub("^[^_]+_","",vIII$sample))
rownames(vIII) <- vIII$sample
vIII$vIII.percentage <- vIII$vIII.reads.v3 / (vIII$vIII.reads.v3 + vIII$wt.reads.v3) * 100
