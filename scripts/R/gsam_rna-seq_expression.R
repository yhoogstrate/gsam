#!/usr/bin/env R

# obtain metadata to exclude poor samples (blacklist.pca == T) ----


if(!exists("gsam.rna.metadata")) {
  source('scripts/R/gsam_metadata.R')
}

if(!exists("gsam.viii.rnaseq")) {
  source('scripts/R/gsam_rnaseq_egfrviii_expression.R')
}



#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
#Load gene annotation data
if(!exists("gencode.31")) {
  gencode.31 <- read.delim("data/gsam/ref/star-hg19/gencode.v31lift37.annotation.gtf", comment.char="#",stringsAsFactors = F,header=F) %>%
    dplyr::filter(V3 == "gene") %>%
    dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
    dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
    dplyr::mutate(V9 = NULL)
}


# full dataset ----

blacklist <- gsam.rna.metadata %>%
  dplyr::filter(blacklist.pca == T) %>%
  dplyr::pull('sid') %>%
  c("CAO1-replicate", "GAS2-replicate","FAB2") %>% # replicates
  c("KAC2-new") %>% # dna contaminiation
  unique()


gsam.rnaseq.expression <- "data/gsam/output/tables/gsam_featureCounts_readcounts_new.txt" %>%
  read.delim(stringsAsFactors = F,comment="#") %>%
  `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1", colnames(.) ,fixed=F)) %>%
  `colnames<-`(gsub(".","-", colnames(.) ,fixed=T)) %>%
  dplyr::left_join(gencode.31, by=c('Geneid' = 'ENSG') ) %>%
  dplyr::mutate(rn = paste0 ( 
    Geneid,
    "|",
    GENE ,
    "|",
    unlist(lapply(stringr::str_split(Chr,";") , unique)),
    ':',
    unlist(lapply(stringr::str_split(Start,";") , min)),
    '-',
    unlist(lapply(stringr::str_split(End,";") , max)) ,
    '(',
    unlist(lapply(stringr::str_split(Strand,";") , unique)),
    ')' ) ) %>%
  dplyr::mutate(Chr=NULL, Start = NULL, End = NULL, Strand = NULL, Length=NULL, Geneid=NULL, GENE=NULL) %>%
  dplyr::select(-paste0("V",1:8)) %>%
  tibble::column_to_rownames("rn") %>%
  dplyr::select( -c(blacklist) )



colnames(gsam.rnaseq.expression) %in% gsam.viii.rnaseq$sample

rm(blacklist)


# Santoesha:
# - CAO1.replicate, GAS2.replicate, FAB2 want dat waren replicated
# - blacklistp.pca
# - KAC2-new << van nieuwe samples


## vst ----


if(!exists('gsam.viii.rnaseq')) {
  source('scripts/R/gsam_rnaseq_egfrviii_expression.R')
}

if(!exists('gsam.rnaseq.expression.vst')) {
  tmp.1 <- gsam.rnaseq.expression %>%
      dplyr::filter(rowSums(.) >= ncol(.))

  tmp.2 <- gsam.viii.rnaseq %>%
      dplyr::mutate(wt.reads = NULL) %>%
      dplyr::mutate(n.reads = NULL) %>%
      dplyr::mutate(egfrviii.pct = NULL) %>%
      dplyr::mutate(sid = NULL) %>%
      dplyr::rename(EGFRvIII = vIII.reads) %>%
      tibble::column_to_rownames('sample') %>%
      t() %>%
      as.data.frame(stringsAsFactors = F) %>%
      dplyr::select(colnames(tmp.1))

    tmp <- rbind(tmp.1, tmp.2)
  
  rm(tmp.1, tmp.2)
  
  cond <- as.factor(paste0('r',sample(c(1,2),ncol(tmp), T)))
  tmp <- DESeq2::DESeqDataSetFromMatrix(tmp,  DataFrame(cond), ~cond)
  gsam.rnaseq.expression.vst <- assay(vst(tmp,blind=T))
  rm(cond, tmp)
}



# ## egfr high/low ----
# 
# tmp <- sort(gsam.rnaseq.expression.vst[gsub("\\..+$","",rownames(gsam.rnaseq.expression.vst)) == "ENSG00000146648",])
# 
# #plot(tmp)
# 
# splt <- 189
# 
# #abline(v=splt)
# #abline(v=splt + 1)
# 
# cutoff <- (tmp[splt] + tmp[splt + 1]) / 2
# #abline(h=cutoff,col="red")
# 
# gsam.egfr.high.expressed <- names(tmp[tmp >= cutoff])
# gsam.egfr.not.high.expressed <- names(tmp[tmp < cutoff])
# 
# rm(cutoff, tmp, splt)
# 
# 
# ## MDM2 high/low ----
# # ENSG00000135679
# 
# tmp <- sort(gsam.rnaseq.expression.vst[gsub("\\..+$","",rownames(gsam.rnaseq.expression.vst)) == "ENSG00000135679",])
# 
# #plot(tmp)
# 
# splt <- 315
# 
# #abline(v=splt)
# #abline(v=splt + 1)
# 
# cutoff <- (tmp[splt] + tmp[splt + 1]) / 2
# 
# #abline(h=cutoff,col="red")
# 
# gsam.mdm2.high.expressed <- names(tmp[tmp >= cutoff])
# gsam.mdm2.not.high.expressed <- names(tmp[tmp < cutoff])
# 
# rm(cutoff, tmp, splt)
# 
# 
# 
# 
# ## ---- CDK4 high/low ----
# 
# ## ENSG00000135446
# 
# tmp <- sort(gsam.rnaseq.expression.vst[gsub("\\..+$","",rownames(gsam.rnaseq.expression.vst)) == "ENSG00000135446",])
# 
# #plot(tmp)
# 
# splt <- 300
# 
# #abline(v=splt)
# #abline(v=splt + 1)
# 
# cutoff <- (tmp[splt] + tmp[splt + 1]) / 2
# 
# #abline(h=cutoff,col="red")
# 
# gsam.cdk4.high.expressed <- names(tmp[tmp >= cutoff])
# gsam.cdk4.not.high.expressed <- names(tmp[tmp < cutoff])
# 
# rm(cutoff, tmp, splt)
# 



