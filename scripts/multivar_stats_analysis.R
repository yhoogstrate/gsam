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

# @todo: ADD vIII

f_normalised <- t(t(f) * sample_size_factors)


# ----- find those samples that have 2 resections -----
two_resections_rna <- gsam.metadata[ duplicated(gsam.metadata$pid),]$pid
#two_resections_rna <- gsam.metadata[  gsam.metadata$pid %in% two_resections_rna ,]$sample.id






