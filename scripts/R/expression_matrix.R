#!/usr/bin/env R

expression_matrix <- read.delim("output/tables/featureCounts_gsam_1st96.exon-level.txt",stringsAsFactors = F,comment="#")
colnames(expression_matrix) <- gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(expression_matrix),fixed=F)

gene_matrix <- expression_matrix[,1:6]

rownames(expression_matrix) <- expression_matrix$Geneid
expression_matrix$Geneid <- NULL
expression_matrix$Chr <- NULL
expression_matrix$Start <- NULL
expression_matrix$End <- NULL
expression_matrix$Strand <- NULL
expression_matrix$Length <- NULL

