#!/usr/bin/env R

# obtain metadata to exclude poor samples (blacklist.pca == T) ----


if(!exists("gsam.rna.metadata")) {
  source('scripts/load_G-SAM_metadata.R')
}

if(!exists("gsam.viii.rnaseq")) {
  source('scripts/R/gsam_rnaseq_egfrviii_expression.R')
}



#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
#Load gene annotation data
if(!file.exists('tmp/gencode.31.Rds')) {
  gencode.31 <- read.delim("data/gsam/ref/star-hg19/gencode.v31lift37.annotation.gtf", comment.char="#",stringsAsFactors = F,header=F) %>%
    dplyr::filter(V3 == "gene") %>%
    dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
    dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
    dplyr::mutate(gene_type = gsub("^.+gene_type ([^;]+);.+$","\\1",V9)) %>%
    dplyr::mutate(V9 = NULL)
  
  saveRDS(gencode.31, 'tmp/gencode.31.Rds')
} else {
  gencode.31 <- readRDS('tmp/gencode.31.Rds')
}




# full dataset ----

blacklist <- gsam.rna.metadata %>%
  dplyr::filter(blacklist.pca == T | pat.with.IDH == T) %>%
  dplyr::pull('sid') %>%
  c("CAO1-replicate", "GAS2-replicate","FAB2") %>% # replicates
  c("KAC2-new") %>% # dna contaminiation
  c("BAI2") %>% # Has BAI2-new, which likely wasn't in the PCA blacklist
  unique()


gsam.rnaseq.expression <- "data/gsam/output/tables/gsam_featureCounts_readcounts_new.txt" %>%
  read.delim(stringsAsFactors = F,comment="#") %>%
  `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1", colnames(.) ,fixed=F)) %>%
  `colnames<-`(gsub(".","-", colnames(.) ,fixed=T)) %>%
  dplyr::left_join(gencode.31, by=c('Geneid' = 'ENSG') ) %>%
  dplyr::filter(gene_type %in% c('protein_coding')) %>%  # ,'lncRNA'
  dplyr::filter(V1 != "chrM") %>%
  dplyr::mutate(gene_type = NULL) %>%
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
  dplyr::select(-all_of(paste0("V",1:8))) %>%
  dplyr::filter(grepl("_PAR_Y",rn) == F) %>% # 42 duplicate ENSEMBL IDs that are located on both X and Y
  tibble::column_to_rownames("rn") %>%
  dplyr::select( -all_of(blacklist))



#colnames(gsam.rnaseq.expression)[colnames(gsam.rnaseq.expression) %in% gsam.viii.rnaseq$sample == F]


stopifnot(colnames(gsam.rnaseq.expression) %in% gsam.viii.rnaseq$sample)
stopifnot(sum(duplicated(gsub("^([A-Z0-9]+).+$","\\1",rownames(gsam.rnaseq.expression)))) == 0) # no duplicate ens ids

rm(blacklist)


# Santoesha:
# - CAO1.replicate, GAS2.replicate, FAB2 want dat waren replicated
# - blacklistp.pca
# - KAC2-new << van nieuwe samples

# BAI2, KAC1-new en KAE1-new hebben geen subtypes?
# BAI2: should be replaced by BAI2-new
#
# in featureCounts;
# KAC1-new: ja - zat ook in Santoesha's classi
# "KAC1-new" %in% colnames(gsam.rnaseq.expression.vst)
# KAE1-new: ja - zat ook in Santoesha's classi
# "KAE1-new" %in% colnames(gsam.rnaseq.expression.vst)

# VST transform ----


## expression data ----


sel <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |>  # replicates
  dplyr::filter(tumour.percentage.dna >= 15) |> 
  dplyr::pull(sid)
#n.gsam.samples <- length(sel)


gsam.gene.expression.all <- gsam.rnaseq.expression |> 
  dplyr::select(sel) |> 
  dplyr::filter(rowSums(dplyr::across()) > ncol(dplyr::across()) * 3)
stopifnot(colnames(gsam.gene.expression.all) == sel)


gsam.gene.expression.all.vst <- gsam.gene.expression.all %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay()


#rm(gsam.gene.expression.all)


# if(!exists('gsam.rnaseq.expression.vst')) {
#   tmp.1 <- gsam.rnaseq.expression %>%
#       dplyr::filter(rowSums(.) >= ncol(.) * 3)
# 
#   tmp.2 <- gsam.viii.rnaseq %>%
#       dplyr::mutate(wt.reads = NULL) %>%
#       dplyr::mutate(n.reads = NULL) %>%
#       dplyr::mutate(egfrviii.pct = NULL) %>%
#       dplyr::mutate(sid = NULL) %>%
#       dplyr::rename(EGFRvIII = vIII.reads) %>%
#       tibble::column_to_rownames('sample') %>%
#       t() %>%
#       as.data.frame(stringsAsFactors = F) %>%
#       dplyr::select(colnames(tmp.1))
# 
#     tmp <- rbind(tmp.1, tmp.2)
#   
#   rm(tmp.1, tmp.2)
#   
#   cond <- as.factor(paste0('r',sample(c(1,2),ncol(tmp), T)))
#   tmp <- DESeq2::DESeqDataSetFromMatrix(tmp, S4Vectors::DataFrame(cond), ~cond)
#   gsam.rnaseq.expression.vst <- SummarizedExperiment::assay(DESeq2::vst(tmp,blind=T))
#   rm(cond, tmp)
# 
# }


# ## egfr high/low 
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
# ## MDM2 high/low 
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
# ## ---- CDK4 high/low 
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



