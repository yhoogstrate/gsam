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
source("scripts/R/ensembl_to_geneid.R")
source("scripts/R/chrom_sizes.R")

source("scripts/R/dna_idh.R")

# ---- obtain expression matrix and sizeFactors ----
e <- expression_matrix

colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id
cond <- as.factor(round(runif(ncol(e))) + 1)

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
keep <- (rowSums(counts(dds)) >= ncol(e) * 3) # post dds removal
dds <- dds[keep,]

sample_size_factors <- colData(estimateSizeFactors(dds))$sizeFactor


rm(dds, keep, cond)

# ---- get ligand table ----


sel <- match( gsub("\\..+$","",tt1) , gsub("\\..+$","",rownames(e))) 
f <- e[sel,]
rownames(f) <- ens_ids[match( gsub("\\..+$","",rownames(f)), gsub("\\..+", "" ,ens_ids$ens.id) ),]$gene.id

f_normalised <- t(t(f) * sample_size_factors)

# @todo: ADD vIII ~ does not need to be normalized as it is relative to wt EGFR
order <- match(vIII$sample , colnames(f_normalised))
#colnames(f_normalised) == vIII[order,]$sample # matches
#vIII[order,]$vIII.percentage


f_normalised <- t(t(f) * sample_size_factors)
f_normalised <- f_normalised + 0.0001


v3 <- vIII[order,]$vIII.percentage
v3[is.na(v3)] <- 0.0
v3 <- v3 + 0.0001 
g <- rbind(f_normalised , vIII.percentage = v3)

# ----- find those samples that have 2 resections -----
two_resections_rna <- gsam.metadata[ duplicated(gsam.metadata$pid),]$pid
#two_resections_rna <- gsam.metadata[  gsam.metadata$pid %in% two_resections_rna ,]$sample.id


# ---- make new table with the fold changes between res1 and res2 ----
g <- g[,gsub("[0-9]+$","",colnames(g)) %in% two_resections_rna]
g1 <- g[,gsub("^[A-Z]+","",colnames(g)) == "1"]
g2 <- g[,gsub("^[A-Z]+","",colnames(g)) == "2"]

d <- log2(g2/g1)


corrplot(cor(t(d),method="spearman"))

