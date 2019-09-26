#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)
library(ggplot2)

# ---- load data ----
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)

source("scripts/R/dna_idh.R")

# ---- unsupervised DESeq2 PCA [ all genes, all samples ] ----
e <- expression_matrix

colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id

# EGFR
#egfr_high <- colnames(e.vst) [ e.vst[match("ENSG00000146648" ,gsub("\\..+$","",rownames(e.vst))),] > 12.5 ]
#cond <- as.factor(gsub("FALSE","EGFR.low",gsub("TRUE","EGFR.high",as.character(colnames(e) %in% egfr_high),fixed=T),fixed=T))

# cdk4
#cdk4_high <- colnames(e.vst) [ e.vst[match("ENSG00000135446" ,gsub("\\..+$","",rownames(e.vst))),] > 8.64 ]
#cond <- as.factor(colnames(e) %in% cdk4_high)


# v3
#cond <- as.character(gsam.metadata$v3.rnaseq.v3.percentage > 1.0)
#cond[is.na(cond)] <- "NA"
#cond <- as.factor(cond)

# exon read count
#cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$featurecounts.Assigned < 2000000))
#cond <- as.factor(gsub("TRUE", "low.Q", gsub("FALSE", "high.Q", cond, fixed=T), fixed=T))

# resection
cond <- factor(gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$resection)

# gender
cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$donor_sex))
cond[is.na(cond)] <- "NA"
cond <- as.factor(cond)


cond <- as.factor(gsub("FALSE","DNA.wt", gsub("TRUE","DNA.IDH1",as.character(colnames(e) %in% is_idh))))

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))

#pc <- plotPCA(vst(dds), intgroup=c("cond"))
p <- plotPCA(vst(dds), intgroup=c("cond"))
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
print(p)




# do DE

outliers = rownames(p$data)[ p$data[,1] > 12.5 ] 
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

o <- match(rownames(sig) , gene_matrix$Geneid)
sig$chr <- gene_matrix[o,]$Chr
sig$start <- paste0(round(as.numeric(gene_matrix[o,]$Start) / 1000000), 'M')

write.table(sig, "/tmp/sig-right-area.txt")








outliers2 = rownames(p$data)[ p$data[,1] < 0 & p$data[,2] > 12.5] 
cond <- as.factor(colnames(e) %in% outliers2)


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




o <- match(rownames(sig) , gene_matrix$Geneid)
sig$chr <- gene_matrix[o,]$Chr
sig$start <- paste0(round(as.numeric(gene_matrix[o,]$Start) / 1000000), 'M')

write.table(sig, "/tmp/sig-upper-corner-cluster.txt")



