#!/usr/bin/env R

# ---- load libs ----

library(tidyverse)


# ---- Load data ----

source('scripts/R/ligands.R')

if(!exists('expression_matrix_full_new')) {
  source('scripts/R/gsam_expression_matrix.R')
}

if(!exists('gsam.patient.metadata') | !exists('gsam.rna.metadata') ) {
  source('scripts/R/gsam_metadata.R')
}

if(!exists('gsam.viii.rnaseq')) {
  source('scripts/R/gsam_rnaseq_egfrviii_expression.R')
}


# ---- subsample - take out low-Q ----


if(file.exists("tmp/expression_matrix_full_new.vst.Rds")) {

  expression_matrix_full_new.vst <- readRDS("tmp/expression_matrix_full_new.vst.Rds")

} else {

  exclude <- gsam.rna.metadata %>%
      dplyr::filter(blacklist.pca == T | batch == "old")  %>%
      dplyr::pull('sid') %>%
      # "EBP1", "FAH2", hebben geen paar
      c("KAC1-new","KAE1-new") %>%
      # ,"AAC1","AAC2" laag tumour percentage?
      unique()


  expression_matrix_full_new.vst <- expression_matrix_full_new %>%
    dplyr::select( -exclude )
  
  viii.tmp <- gsam.viii.rnaseq %>%
    tibble::column_to_rownames('sample') %>%
    dplyr::select('vIII.reads') %>%
    t() %>%
    data.frame() %>%
    `colnames<-`(gsub ('.', '-', colnames(.) , fixed=T)) %>% # bla.new -> bla-new
    dplyr::select(colnames(expression_matrix_full_new.vst)) %>%
    `rownames<-`(c("EGFRvIII"))
  
  stopifnot(ncol(viii.tmp) == ncol(expression_matrix_full_new.vst))
  
  
  expression_matrix_full_new.vst <- expression_matrix_full_new.vst %>%
    dplyr::bind_rows(viii.tmp) %>%
    dplyr::filter(rowSums(.) >= 3 * ncol(.) | rownames(.) %in% c(tt1, "EGFRvIII") ) %>%
    DESeq2::DESeqDataSetFromMatrix(S4Vectors::DataFrame(tmp=as.factor(1:ncol(.))), ~tmp) %>%
    DESeq2::varianceStabilizingTransformation(blind=T) %>%
    SummarizedExperiment::assay() %>%
    data.frame(stringsAsFactors = F) %>%
    `colnames<-`(gsub ('.', '-', colnames(.) , fixed=T)) # bla.new -> bla-new
  
  # check duplicates?
  stopifnot(sum(duplicated(gsub("-new","", colnames(expression_matrix_full_new.vst) , fixed = T))) == 0)
  
  
  rm(exclude, viii.tmp)
  
  
  saveRDS(expression_matrix_full_new.vst, file="tmp/expression_matrix_full_new.vst.Rds")
}



