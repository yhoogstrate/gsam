# Set wd

# ---- initialization setup ----

setwd("~/projects/gsam")


# ---- load: libs ----

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(fgsea)
library(limma)
library(tidyverse)

# ---- load: functions ----

#source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
#ensembl_genes <- get_ensembl_hsapiens_gene_ids()


# ---- load: data ----
# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")

source("scripts/R/dna_idh.R")
source("scripts/R/chrom_sizes.R")

source("scripts/R/job_gg_theme.R")

# ---- unsupervised DESeq2 PCA / MDS / clustering / t-SNE ----

e <- expression_matrix_full_new
#%>%   dplyr::filter(rowSums(.) > ncol(.) * 10)
#dim(e[,colSums(e) > 5000000])

e <- e[order(-rowSums(e)),]
e <- e[! grepl("chrM", rownames(e) ),]
e <- e[,colSums(e) > 1000000]

e <- e[,1:196]

cond <- as.factor(rep(c("r1","r2"), 500)[1:ncol(e)])
dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
dds <- DESeq(dds)
plotMA(dds)





res <- results(dds)
plotMA(res)



