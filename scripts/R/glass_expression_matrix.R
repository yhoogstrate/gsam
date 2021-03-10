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
  dplyr::filter(gene_id %in% c("ENSG00000198804", "ENSG00000198886", "ENSG00000198938","ENSG00000198712",
                               "ENSG00000198727") == F ) %>% # extreme high counts from odd genes
  tibble::column_to_rownames('gene_id') %>%
  round() %>%
  `colnames<-`( gsub(".","-",colnames(.), fixed=T) ) %>%
  `colnames<-`( gsub("_1$","",colnames(.), fixed=F) )



  # Using the following will pick the max() transcript per gene
  # I tried this as well in the hope that it would remove the batch
  # effect but it didn't:'
  # 
  # dplyr::arrange(gene_id) %>%
  # as.data.frame() %>%
  # mutate(rs = rowSums(dplyr::select( ., !dplyr::contains("gene_id") ))) %>%
  # tibble::as_tibble() %>%
  # dplyr::group_by(gene_id) %>%
  # top_n(n=1, wt = rs) %>%
  # dplyr::mutate(rs=NULL) %>%
  # dplyr::arrange(gene_id) %>%
  # dplyr::filter(duplicated(gene_id) == F) %>% # exclude ties
  # tibble::rownames_to_column('tmp') %>% 
  # dplyr::mutate(tmp=NULL) %>%
  # tibble::column_to_rownames('gene_id') %>%
  # round() %>%
  # `colnames<-`( gsub(".","-",colnames(.), fixed=T) ) %>%
  # `colnames<-`( gsub("_1$","",colnames(.), fixed=F) )


  
rm(tmp)


tmp = glass.gbm.rnaseq.expression.vst[rownames(glass.gbm.rnaseq.expression.vst) == "ENSG00000198727",]
type = as.factor(gsub("^(.).*$","\\1",names(tmp)) )
o = order(tmp)
plot(tmp[o], pch=19,cex=0.95,col=as.numeric(type[o]) + 1)
 

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
    , by     = c('sid' = 'aliquot_barcode' )  )  %>%
  dplyr::mutate(sid.label = gsub("^(.)...(-..-....-).(.).*$","\\1\\2\\3",sid) ) %>%
  dplyr::mutate(dataset = gsub("^(....).*$","\\1",sid) )

stopifnot(!is.na(glass.gbm.rnaseq.metadata $GBM.transcriptional.subtype.Synapse))




# VST transform ----


cond <- as.factor(gsub("^(.).*$","\\1",colnames(glass.gbm.rnaseq.expression)))
glass.gbm.rnaseq.expression.vst <- DESeq2::DESeqDataSetFromMatrix(glass.gbm.rnaseq.expression, S4Vectors::DataFrame(cond), ~cond)
glass.gbm.rnaseq.expression.vst <- SummarizedExperiment::assay(DESeq2::vst(glass.gbm.rnaseq.expression.vst, blind=F))
rm(cond)








