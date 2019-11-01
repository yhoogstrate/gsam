#!/usr/bin/env R

expression_matrix <- read.delim("data/RNA/output/tables/featureCounts_gsam_1st96.exon-level.txt",stringsAsFactors = F,comment="#")
colnames(expression_matrix) <- gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(expression_matrix),fixed=F)

gene_matrix <- expression_matrix[,1:6]

rownames(expression_matrix) <- expression_matrix$Geneid
expression_matrix$Geneid <- NULL
expression_matrix$Chr <- NULL
expression_matrix$Start <- NULL
expression_matrix$End <- NULL
expression_matrix$Strand <- NULL
expression_matrix$Length <- NULL



# ---- full dataset ----

expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt",stringsAsFactors = F,comment="#")
colnames(expression_matrix_full) <- gsub("^[^\\.]+\\.([^\\.]+)\\..+$","\\1",colnames(expression_matrix_full),fixed=F)

rownames(expression_matrix_full) <- expression_matrix_full$Geneid
expression_matrix_full$Geneid <- NULL
expression_matrix_full$Chr <- NULL
expression_matrix_full$Start <- NULL
expression_matrix_full$End <- NULL
expression_matrix_full$Strand <- NULL
expression_matrix_full$Length <- NULL


