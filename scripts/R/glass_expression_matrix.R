#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(DESeq2)


## read counts from GLASS samples ----

# possibly some match w/ Wang dataset?


# The samples are taken from a data portal:
# www.synapse.org

# expression values are deterimined by Kallisto, per transcript


# hacking this back from GLASS git code refers to gencode v19...
# gencode.v19.chr_patch_hapl_scaff.annotation.gtf << almost no mismatches [1193]



tmp <- 'data/gsam/data/GLASS_GBM_R1-R2/gencode.v19.chr_patch_hapl_scaff.annotation.translate-table.txt' %>%
  read.table(header=F, stringsAsFactors = F) %>%
  dplyr::rename(gene_symbol=V1) %>%
  dplyr::rename(gene_id=V2) %>%
  dplyr::rename(transcript_id=V3) 


glass.gbm.rnaseq.expression <- 'data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt' %>%
  read.delim(stringsAsFactors = F) %>%
  dplyr::mutate(length = NULL) %>% # not needed
  dplyr::left_join(tmp, by=c('target_id' = 'transcript_id')) %>%
  dplyr::filter(!is.na(gene_symbol) ) %>% # 1193 transcript id's not matching gtf Ensembl 64
  dplyr::mutate(target_id = NULL) %>% # aggregate @ gene_id level
  dplyr::mutate(gene_symbol = NULL) %>%
  dplyr::group_by(gene_id) %>%
  summarise(across(everything(), list(sum))) %>%
  tibble::rownames_to_column('tmp') %>% 
  dplyr::mutate(tmp=NULL) %>%
  tibble::column_to_rownames('gene_id') %>%
  round()

rm(tmp)



# Load metadata ----


glass.gbm.rnaseq.metadata <- data.frame(sid = colnames(glass.gbm.rnaseq.expression),  stringsAsFactors = F) %>%
  dplyr::mutate(sid = gsub(".","-",sid,fixed=T)) %>% # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sid = gsub("_1$","",sid,fixed=F)) %>% # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(pid = gsub("^(............).+$","\\1",sid)) %>%
  dplyr::arrange(pid) %>%
  dplyr::mutate(resection = as.factor(gsub("^.............(..).+$","\\1",sid))) %>% # TP is primary tumour? https://github.com/fpbarthel/GLASS
  dplyr::mutate(dataset =  as.factor(gsub("^(....).+$","\\1",sid)) ) %>%
  dplyr::left_join(
    read.table('data/gsam/data/GLASS_GBM_R1-R2/GLASS.GBM.subtypes.from.Synapse.portal.tsv', header=T,stringsAsFactors = F) %>% # https://www.synapse.org/#!Synapse:syn21441635/tables/
      dplyr::select(c('aliquot_barcode','subtype')) %>%
      dplyr::mutate(subtype = as.factor(subtype)) %>%
      dplyr::rename(GBM.transcriptional.subtype.Synapse = subtype)
    , by     = c('sid' = 'aliquot_barcode' )  )

stopifnot(!is.na(glass.gbm.rnaseq.metadata $GBM.transcriptional.subtype.Synapse))




# VST transform [n >= 3] ----


if(!exists('glass.gbm.rnaseq.expression.vst')) {
  tmp <- glass.gbm.rnaseq.expression %>%
    dplyr::filter(rowSums(.) >= (ncol(.) * 1) )
  
  cond <- as.factor(paste0('r',sample(c(1,2),ncol(tmp), T)))
  tmp <- DESeq2::DESeqDataSetFromMatrix(tmp, S4Vectors::DataFrame(cond), ~cond)
  glass.gbm.rnaseq.expression.vst <- SummarizedExperiment::assay(DESeq2::vst(tmp, blind=T))
  rm(cond, tmp)
}





