#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)

# ---- load data ----
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")

# ---- unsupervised DESeq2 PCA [ all genes, all samples ] ----
e <- expression_matrix

colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id

egfr_high <- colnames(e.vst) [ e.vst[match("ENSG00000146648" ,gsub("\\..+$","",rownames(e.vst))),] > 12.5 ]
cond <- as.factor(gsub("FALSE","EGFR.low",gsub("TRUE","EGFR.high",as.character(colnames(e) %in% egfr_high),fixed=T),fixed=T))

cdk4_high <- colnames(e.vst) [ e.vst[match("ENSG00000135446" ,gsub("\\..+$","",rownames(e.vst))),] > 8.64 ]
cond <- as.factor(colnames(e) %in% cdk4_high)

# v3
cond <- as.character(gsam.metadata$v3.rnaseq.v3.percentage > 1.0)
cond[is.na(cond)] <- "NA"
cond <- as.factor(cond)

# exon read count
#cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$featurecounts.Assigned < 2000000))
#cond <- as.factor(gsub("TRUE", "low.Q", gsub("FALSE", "high.Q", cond, fixed=T), fixed=T))

# resection
#cond <- factor(gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$resection)


dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))

pc <- plotPCA(vst(dds), intgroup=c("cond"))



outliers = rownames(pc$data)[ pc$data[,1] > 12.5 ] 
cond <- as.factor(colnames(e) %in% outliers)


dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)

# Keep only on average 3 reads per sample
print(dim(dds))
keep <- (rowSums(counts(dds)) >= ncol(e) * 3)
dds <- dds[keep,]
print(dim(dds))

dds <- DESeq(dds)
res <- results(dds)
dim(res)
res <- res[!is.na(res$padj),] # remove those that have failed, probably 0,0,0,0, high, high, high, high?
res <- res[order(res$padj),]
dim(res)

res <- res[order(res$padj),]
sig <- res[res$padj < 0.01,]




