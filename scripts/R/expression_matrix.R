#!/usr/bin/env R

#expression_matrix <- read.delim("data/RNA/output/tables/featureCounts_gsam_1st96.exon-level.txt",stringsAsFactors = F,comment="#")
#colnames(expression_matrix) <- gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(expression_matrix),fixed=F)
#
#gene_matrix <- expression_matrix[,1:6]
#
#rownames(expression_matrix) <- expression_matrix$Geneid
#expression_matrix$Geneid <- NULL
#expression_matrix$Chr <- NULL
#expression_matrix$Start <- NULL
#expression_matrix$End <- NULL
#expression_matrix$Strand <- NULL
#expression_matrix$Length <- NULL

# ---- load libraries ----

library(tidyverse)

# ---- full dataset ----


gencode.31 <- read.delim("ref/star-hg19/gencode.v31lift37.annotation.gtf", comment.char="#",stringsAsFactors = F,header=F) %>%
  dplyr::filter(V3 == "gene") %>%
  dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1", V9)) %>%
  dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1", V9)) %>%
  dplyr::mutate(V9 = NULL)



# file is currently not at server and @ wrong permissions @ ccbc
#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
expression_matrix_full_new <- read.delim("output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#") %>%
  `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1", colnames(.) ,fixed=F)) %>%
  dplyr::left_join(gencode.31, by=c('Geneid' = 'ENSG') ) %>%
  dplyr::mutate(rn = paste0 ( Geneid,"|", GENE , "|", unlist(lapply(str_split(Chr,";") , unique)), ':', unlist(lapply(str_split(Start,";") , min)), '-', unlist(lapply(str_split(End,";") , max)) , '(', unlist(lapply(str_split(Strand,";") , unique)), ')' ) ) %>%
  dplyr::mutate(Chr=NULL, Start = NULL, End = NULL, Strand = NULL, Length=NULL, Geneid=NULL) %>%
  tibble::column_to_rownames("rn")








