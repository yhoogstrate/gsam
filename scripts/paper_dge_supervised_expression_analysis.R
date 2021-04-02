#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(DESeq2)


library(ggplot2)
library(ggrepel)

#library(pheatmap)
#library(fgsea)
#library(limma)

library(EnhancedVolcano)
library(patchwork)



# load data ----

#source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
#ensembl_genes <- get_ensembl_hsapiens_gene_ids()

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


source("scripts/R/ligands.R")
source("scripts/R/subtype_genes.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

source('scripts/R/wang_glioma_intrinsic_genes.R')

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set


# prepare data ----

## GSAM ----


# gsam.metadata.all
# gsam.metadata.r1
# gsam.metadata.r2

gsam.metadata.all <- gsam.rna.metadata %>%
  dplyr::filter(blacklist.pca == F) %>%
  dplyr::filter(pat.with.IDH == F) %>%
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
  dplyr::filter(tumour.percentage.dna >= 15) %>%
  dplyr::mutate(tpc = 1 - (tumour.percentage.dna / 100))


gsam.metadata.all.paired <- gsam.metadata.all %>%
  dplyr::filter(pid %in% 
                  (gsam.metadata.all %>%
                     dplyr::group_by(pid) %>%
                     dplyr::tally() %>%
                     dplyr::filter(n == 2) %>% 
                     dplyr::ungroup() %>%
                     dplyr::filter(!duplicated(pid)) %>%
                     dplyr::pull(pid))
                  ) %>%
  dplyr::mutate(pid = as.factor(as.character(pid))) # re-factor?




gsam.metadata.r1 <- gsam.metadata.all %>%
  dplyr::filter(resection == "r1")

gsam.metadata.r2 <- gsam.metadata.all %>%
  dplyr::filter(resection == "r2")


# uncomment out the minimum TPC
wilcox.test(gsam.metadata.all %>%
              dplyr::filter(resection == "r1") %>%
              dplyr::pull(tumour.percentage.dna) ,
            gsam.metadata.all %>%
              dplyr::filter(resection == "r2") %>%
              dplyr::pull(tumour.percentage.dna)
            , alternative = "two.sided")



# gsam.gene.expression.all
# gsam.bfg.expression.all

gsam.gene.expression.all <- gsam.rnaseq.expression %>%
  dplyr::select(gsam.metadata.all$sid)
stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)


gsam.gene.expression.all.paired <- gsam.gene.expression.all %>%
  dplyr::select(gsam.metadata.all.paired$sid)
stopifnot(colnames(gsam.gene.expression.all.paired) == gsam.metadata.all.paired$sid)



# gsam.bfg.expression.all <- gsam.rnaseq.expression %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   dplyr::filter( # bona fide glioma genes (BFGs)
#     (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
#       (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
#   ) %>%
#   tibble::column_to_rownames('gid') %>%
#   dplyr::select(gsam.metadata.all$sid)
# 
# stopifnot(colnames(gsam.bfg.expression.all) == gsam.metadata.all$sid)


# gsam.gene.expression.all.vst
# TODO gsam.gene.expression.R1.vst
# TODO gsam.gene.expression.R2.vst




gsam.gene.expression.all.vst <- gsam.gene.expression.all %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay()




## GLASS ----



glass.metadata.all  <- glass.gbm.rnaseq.metadata %>%
  dplyr::mutate(condition = ifelse(resection == "TP","Primary","NotPrimary")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Primary","NotPrimary") )) %>%
  dplyr::filter(idh_status == "IDHwt") %>% # exclude the three IDH mutants according to Synapse WGS/WES VCF files
  dplyr::filter(pid %in%  (glass.gbm.rnaseq.metadata %>%
                             dplyr::filter(grade %in% c('II','III')) %>%
                             dplyr::pull(pid)) == F )



glass.gene.expression.all <- glass.gbm.rnaseq.expression %>%
  dplyr::select(glass.metadata.all$sid)


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)




## per-gene results table ----


results.out <- dplyr::full_join(
    gsam.gene.expression.all %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::select(gid) %>%
    dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
    dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
    dplyr::mutate(in.gsam = T) %>%
    dplyr::mutate(chr = as.factor(gsub("^.+(chr[^:]+):.+$","\\1",gid)))
    ,
    glass.gene.expression.all %>%
      tibble::rownames_to_column('ensembl_id') %>%
      dplyr::select(ensembl_id) %>% 
      dplyr::mutate(in.glass = T)
    , by = c('ensembl_id'='ensembl_id')) %>%
  dplyr::left_join(
    glass.gencode.v19 %>%
      dplyr::arrange( transcript_id) %>%
      dplyr::select(c('gene_symbol','gene_id')) %>%
      dplyr::filter(!duplicated(gene_id)) %>%
      dplyr::rename(ensembl_id = gene_id) %>%
      dplyr::rename(hugo_symbol.gencode = gene_symbol) ,
    by = c ('ensembl_id'='ensembl_id')) %>%
  dplyr::mutate(hugo_symbol = ifelse(is.na(hugo_symbol) , hugo_symbol.gencode , hugo_symbol )) %>%
  dplyr::mutate(hugo_symbol.gencode = NULL) %>%
  dplyr::mutate(is.bfg = (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) | (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol))
    


stopifnot(sum(duplicated(results.out$ensembl_id)) == 0)




# look for transporter genes / influx / efflux genes within a volcano like plot [hedgehog signalling?]
# PharmGKB: ABCB1, OPRM1, INHBA, ITGBL1, FCGR2B
r1 <- read.delim('data/pharmgkb/relationships.tsv',sep="\t",stringsAsFactors = F) %>%
  dplyr::filter((Entity1_type == "Gene" & Entity2_type == "Chemical")) %>%
  dplyr::filter(Association == "not associated") %>%
  dplyr::pull(Entity1_name) %>%
  unique()

r2<- read.delim('data/pharmgkb/relationships.tsv',sep="\t",stringsAsFactors = F) %>%
  dplyr::filter((Entity2_type == "Gene" & Entity1_type == "Chemical")) %>%
  dplyr::filter(Association == "not associated") %>%
  dplyr::pull(Entity2_name) %>%
  unique()


results.out <- results.out %>%
  dplyr::mutate(pharma.relation = hugo_symbol %in% c(r1, r2))

rm(r1, r2)





# corr TPC [GSAM] ----


stopifnot(colnames(gsam.gene.expression.all.vst) == gsam.metadata.all$sid)
gsam.gene.expression.all.cor.estimate <- data.frame(apply(gsam.gene.expression.all.vst, 1, function (x)  cor.test(x, gsam.metadata.all$tumour.percentage.dna) %>% purrr::pluck('estimate') ), stringsAsFactors = F) %>%
  `colnames<-`("estimate") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor.statistic <- data.frame(apply(gsam.gene.expression.all.vst, 1, function (x)  cor.test(x, gsam.metadata.all$tumour.percentage.dna) %>% purrr::pluck('statistic') ), stringsAsFactors = F) %>%
  `colnames<-`("statistic") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor.p.value <- data.frame(apply(gsam.gene.expression.all.vst, 1, function (x)  cor.test(x, gsam.metadata.all$tumour.percentage.dna) %>% purrr::pluck('p.value') ), stringsAsFactors = F) %>%
  `colnames<-`("p.value") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor <- gsam.gene.expression.all.cor.estimate %>%
  dplyr::left_join(gsam.gene.expression.all.cor.statistic , by=c('gid' = 'gid') ) %>%
  dplyr::left_join(gsam.gene.expression.all.cor.p.value , by=c('gid' = 'gid') ) %>%
  `colnames<-`(paste0(colnames(.), ".gsam.cor.tpc")) %>%
  dplyr::rename(gid = gid.gsam.cor.tpc)

rm(gsam.gene.expression.all.cor.estimate , gsam.gene.expression.all.cor.statistic , gsam.gene.expression.all.cor.p.value)



results.out <- results.out %>%
  dplyr::left_join(gsam.gene.expression.all.cor , by = c('gid' = 'gid')) 

stopifnot("statistic.gsam.cor.tpc" %in% colnames(gsam.gene.expression.all.cor))
stopifnot("statistic.gsam.cor.tpc" %in% colnames(results.out))




# 
# 
# plot(gsam.gene.expression.all.cor$statistic,abs(gsam.gene.expression.all.cor$estimate) + runif(nrow(gsam.gene.expression.all.cor), 0, 0.2) , pch=19,cex=0.05)
# 
# 
# 
# stopifnot(gsam.res.tpc.all$gid == gsam.gene.expression.all.cor$gid)
# plot(gsam.res.tpc.all$log2FoldChange , gsam.gene.expression.all.cor$statistic , pch=19,cex=0.05)
# 
# 
# 
# a = gsam.gene.expression.all.cor %>%
#   dplyr::arrange(-p.value) %>%
#   dplyr::top_n(500, -p.value) %>%
#   dplyr::pull(hugo_symbol)
# 
# b = gsam.gene.expression.all.cor %>%
#   dplyr::arrange(-p.value) %>%
#   dplyr::top_n(500, -p.value) %>%
#   dplyr::pull(ensembl_id)






# DE unpaired all [G-SAM] ----



gsam.gene.res.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                              dplyr::filter(rowSums(.) > ncol(.) * 3)
                                            , gsam.metadata.all, ~resection ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>%
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  `colnames<-`(paste0(colnames(.),".gsam.res")) %>%
  dplyr::rename(gid = gid.gsam.res)

results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res , by = c('gid' = 'gid'))



# TODO remove .paired twice:
gsam.gene.res.tpc.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                                  dplyr::filter(rowSums(.) > ncol(.) * 3),
                                                gsam.metadata.all, ~tpc + resection ) %>% # + resection; corrected for tpc
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  #dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
  #dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
  `colnames<-`(paste0(colnames(.),".gsam.tpc.res")) %>%
  dplyr::rename(gid = gid.gsam.tpc.res)

results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.tpc.res , by = c('gid' = 'gid'))



# DE unpaired all [GLASS] ----

# no per-sample tumor percentage available in Synapse portal


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)

glass.gene.res.res <- DESeqDataSetFromMatrix(glass.gene.expression.all %>%
                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
                                             , glass.metadata.all, ~condition ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('ensembl_id') %>% 
  `colnames<-`(paste0(colnames(.),".glass.res")) %>%
  dplyr::rename(ensembl_id = ensembl_id.glass.res)


results.out <- results.out %>%
  dplyr::left_join(glass.gene.res.res, by = c('ensembl_id' = 'ensembl_id'))





# < : : : deprecated code : : : > -----

# [x] DE unpaired ~ TPC as numeric [GSAM] ----

# quite sketchy test - sensitive to outliers and many NA corrected padj's
# 
# stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)
# 
# gsam.res.tpc.all <- DESeqDataSetFromMatrix(gsam.gene.expression.all, gsam.metadata.all, ~tpc ) %>% # + tpc
#   DESeq(parallel=T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   `colnames<-`(paste0(colnames(.),".gsam.tpc"))  %>%
#   dplyr::rename(gid = gid.gsam.tpc)
#   
# #sum(is.na ( gsam.res.tpc.all$padj.tpc ))
# #sum(is.na ( gsam.res.tpc.all$p.value.tpc ))
# 
# results.out <- results.out %>%
#   dplyr::left_join(gsam.res.tpc.all , by = c('gid' = 'gid')) 
# 


# 
# gsam.res.tpc.all %>%
#   dplyr::top_n(500, padj.tpc) %>% dplyr::pull(hugo_symbol.tpc)




# [x] zelfde maar dan voor alleen R1 ----

# 
# 
# expression <- gsam.rnaseq.expression %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   #  dplyr::filter( # bona fide glioma genes (BFGs)
#   #    (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
#   #      (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
#   #  ) %>%
#   tibble::column_to_rownames('gid') %>%
#   dplyr::select(metadata$sid)
# 
# stopifnot(colnames(expression) == metadata$sid)
# 
# 
# res.tpc.r1 <- DESeqDataSetFromMatrix(expression, metadata, ~tpc ) %>% # + tpc
#   DESeq(parallel=T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   #dplyr::arrange(padj) %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".r1"))
# 
# 
# 
# [x] zelfde maar dan voor alleen R2 ----
# 
# rm(metadata, expression)
# 
# metadata <- gsam.rna.metadata %>%
#   dplyr::filter(blacklist.pca == F) %>%
#   dplyr::filter(pat.with.IDH == F) %>%
#   dplyr::filter(resection == 'r2') %>%
#   dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
#   dplyr::filter(tumour.percentage.dna >= 25) %>%
#   dplyr::mutate(tpc = 1 - (tumour.percentage.dna / 100))
# 
# 
# expression <- gsam.rnaseq.expression %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   #  dplyr::filter( # bona fide glioma genes (BFGs)
#   #    (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
#   #      (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
#   #  ) %>%
#   tibble::column_to_rownames('gid') %>%
#   dplyr::select(metadata$sid)
# 
# stopifnot(colnames(expression) == metadata$sid)
# 
# 
# res.tpc.r2 <- DESeqDataSetFromMatrix(expression, metadata, ~tpc ) %>% # + tpc
#   DESeq(parallel=T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   #dplyr::arrange(padj) %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".r2"))
# 
# 
# 
# plt <- res.tpc.r1 %>%
#   dplyr::left_join( res.tpc.r2 , by = c('gid.r1'='gid.r2'))
# 
# ggplot(plt, aes(x =log2FoldChange.r1 , y = log2FoldChange.r2, label =  hugo_symbol.r2 )) +
#   geom_text_repel(data = subset(plt, abs(log2FoldChange.r1 - log2FoldChange.r2) > 8.5 ) , size=2) +
#   geom_point(cex=0.5) 
# 



# [x] DE unpaired BFGs ----
# 
# # TODO exclude IDH mutants [check]
# # TODO x-check replicates/duplicates? [check] << taken out using 'gsam.rnaseq.expression'
# # TODO exclude low tumour percentage ? << not if BFGs are explitily used & double test?
# 
# 
# gsam.bfg.res.tpc <- DESeqDataSetFromMatrix(gsam.bfg.expression.all, gsam.metadata.all, ~tpc ) %>% # + tpc
#   DESeq(parallel = T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::arrange(padj) %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".tpc"))
# 
# 
# gsam.bfg.res.res <- DESeqDataSetFromMatrix(gsam.bfg.expression.all, gsam.metadata.all, ~resection ) %>% # + resection
#   DESeq(parallel = T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   dplyr::arrange(padj) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".res"))
# 
# 
# gsam.bfg.res.tpc.res <- DESeqDataSetFromMatrix(gsam.bfg.expression.all, gsam.metadata.all, ~tpc + resection ) %>% # + resection; corrected for tpc
#   DESeq(parallel = T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   dplyr::arrange(padj) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   `colnames<-`(paste0(colnames(.),".tpc.res"))
# 
# 
# 
# gsam.bfg.res.combined <- gsam.bfg.res.tpc %>%
#   dplyr::full_join(gsam.bfg.res.res, by=c('gid.tpc' = 'gid.res')) %>%
#   dplyr::full_join(gsam.bfg.res.tpc.res, by=c('gid.tpc' = 'gid.tpc.res')) %>%
#   dplyr::full_join(gsam.gene.expression.all.cor, by=c('gid.tpc' = 'gid')) %>%
#   dplyr::filter(!is.na(padj.res))%>%
#   dplyr::filter(!is.na(padj.tpc))%>%
#   dplyr::filter(!is.na(padj.tpc.res))
# 
# 
# 
# 
# p1 <- ggplot(gsam.bfg.res.combined, aes(x=log2FoldChange.res ,
#                                         y=statistic,
#                                         col=significant.res,
#                                         label=hugo_symbol.tpc.res ) ) + 
#   geom_point(data=subset(gsam.bfg.res.combined, significant.res == F ),pch=19,cex=0.05) +
#   geom_point(data=subset(gsam.bfg.res.combined, significant.res == T ),pch=19,cex=0.5) +
#   scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
#   geom_text_repel(data = subset(gsam.bfg.res.combined, significant.res == T & abs(log2FoldChange.res) > 1), size = 2.4 )  +
#   labs(x = "log2FC R1 vs. R2 (unpaired)",
#        y="Correlation t-statistic with tumour percentage",
#        col="Difference significant (R1 ~ R2)"
#   ) +
#   youri_gg_theme
# 
# p2 <- ggplot(gsam.bfg.res.combined, aes(x=log2FoldChange.tpc.res ,
#                                         y= statistic,
#                                         col=significant.tpc.res,
#                                         label=hugo_symbol.tpc.res ) ) + 
#   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == F ),pch=19,cex=0.05) +
#   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == T ),pch=19,cex=0.5) +
#   scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
#   geom_text_repel(data = subset(gsam.bfg.res.combined, significant.tpc.res == T & abs(log2FoldChange.tpc.res) > 1), size = 2.4 )  +
#   labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage",
#        col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
#   youri_gg_theme
# 
# p1 + p2
# 
# # 
# # 
# # p3 <- ggplot(gsam.bfg.res.combined, aes(x=log2FoldChange.tpc.res, y= log2FoldChange.tpc, col=significant.tpc.res, label=hugo_symbol.tpc.res ) ) + 
# #   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == F ),pch=19,cex=0.05) +
# #   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == T ),pch=19,cex=0.3) +
# #   scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
# #   #geom_text_repel(data = subset(gsam.bfg.res.combined, significant.tpc.res == T & abs(log2FoldChange.tpc.res) > 1), size = 2.4 )  +
# #   youri_gg_theme
# # 


# [x] DE paired all [G-SAM] ----

# https://support.bioconductor.org/p/59481/
# fitType='local' or 'mean' 

if(!file.exists('tmp/gsam.gene.res.res.paired.Rds')) {
  
  stopifnot(F)
  
  
  # gsam.gene.res.res.paired <- DESeqDataSetFromMatrix(gsam.gene.expression.all.paired %>%
  #                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
  #                                             , gsam.metadata.all.paired , ~pid + resection ) %>% # + resection
  #   DESeq(parallel = T, fitType="local") %>%
  #   results() %>% 
  #   as.data.frame() %>%
  #   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 )  %>%
  #   dplyr::arrange(padj) %>% 
  #   tibble::rownames_to_column('gid') %>%
  #   `colnames<-`(paste0(colnames(.),".res.paired")) %>%
  #   dplyr::rename(gid = gid.res.paired)
  # 
  # saveRDS(gsam.gene.res.res.paired, file = "tmp/gsam.gene.res.res.paired.Rds")
  
  
  
  stopifnot(gsam.metadata.all.paired$sid == colnames(gsam.gene.expression.all.paired))
  # gsam.gene.res.tpc.res.paired <- DESeqDataSetFromMatrix(gsam.gene.expression.all.paired %>%
  #                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
  #                                             , gsam.metadata.all.paired , ~tpc + pid + resection ) %>% # + resection
  #   DESeq(parallel = T, fitType="local") %>%
  #   results() %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  #   dplyr::arrange(padj) %>%
  #   tibble::rownames_to_column('gid') %>%
  #   `colnames<-`(paste0(colnames(.),".tpc.res.paired")) %>%
  #   dplyr::rename(gid = gid.tpc.res.paired)
  # 
  # saveRDS(gsam.gene.res.tpc.res.paired, file = "tmp/gsam.gene.res.tpc.res.paired.Rds")
  
  stopifnot(gsam.metadata.all$sid == colnames(gsam.gene.expression.all))
  test.a <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , gsam.metadata.all , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.a")) %>%
    dplyr::rename(gid = gid.test.a)
  
  test.a.ntpc <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                          dplyr::filter(rowSums(.) > ncol(.) * 3)
                                        , gsam.metadata.all , ~ resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.a.ntpc")) %>%
    dplyr::rename(gid = gid.test.a.ntpc)
  
  
  test.b <- DESeqDataSetFromMatrix(gsam.gene.expression.all.paired %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , gsam.metadata.all.paired , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.b")) %>%
    dplyr::rename(gid = gid.test.b)
  
  
  m <- gsam.metadata.all %>%
    dplyr::filter() %>% 
    dplyr::filter(sid %in% c('FAF2', 'AIA2', 'HAF2-new') == F )  
  e <- gsam.gene.expression.all %>% dplyr::select(m$sid)
  stopifnot(m$sid == colnames(e))
  test.c <- DESeqDataSetFromMatrix(e %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , m , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.c")) %>%
    dplyr::rename(gid = gid.test.c)
  
  
  m <- gsam.metadata.all %>%
    dplyr::filter() %>% 
    dplyr::filter(sid %in% c('FAF2', 'AIA2', 'HAF2-new', "AOA2","BAE2","AAX1","EBB2") == F )  
  e <- gsam.gene.expression.all %>% dplyr::select(m$sid)
  stopifnot(m$sid == colnames(e))
  test.d <- DESeqDataSetFromMatrix(e %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , m , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.d")) %>%
    dplyr::rename(gid = gid.test.d)
  
  
  
  
  
  subset(gsam.gene.res.tpc.res.paired, gid %in% a)$log2FoldChange.tpc.res.paired %>% length()
  subset(test.a, gid %in% a)$log2FoldChange.test.a %>% length()
  subset(test.a.ntpc, gid %in% a)$log2FoldChange.test.a.ntpc %>% length()
  subset(test.b, gid %in% a)$log2FoldChange.test.b %>% length()
  subset(test.c, gid %in% a)$log2FoldChange.test.c %>% length()
  subset(test.d, gid %in% a)$log2FoldChange.test.d %>% length()
  #%>% length()
  
  plt <- gsam.gene.res.tpc.res.paired %>%
    dplyr::left_join(test.a , by  = c('gid'='gid')) %>%
    dplyr::mutate(col = ifelse(lfcSE.test.a  > 0.3 , "outlier", "regular") )
  ggplot(plt , aes(x=stat.tpc.res.paired,
                   y=stat.test.a , 
                   col=col )) +
    geom_point(data = subset(plt, col == "regular"), pch=19,cex=0.5) +
    geom_point(data = subset(plt, col == "outlier"), pch=19,cex=0.85)
  
  
  aa = gsam.gene.res.tpc.res%>% dplyr::filter( padj.tpc.res < 0.2 ) %>% pull(lfcSE.tpc.res)
  ab = gsam.gene.res.tpc.res %>% dplyr::filter( gid %in% a ) %>% pull(lfcSE.tpc.res)
  plot(density(aa))
  abline(v=ab)
  abline(v=0.3, col="red")
  
  
  # a = 
  # [1] "ENSG00000162571.13_5|TTLL10|chr1:1109264-1133315(+)"       "ENSG00000137975.8_3|CLCA2|chr1:86889854-86922236(+)"      
  # [3] "ENSG00000143196.5_4|DPT|chr1:168664706-168698444(-)"       "ENSG00000159173.19_4|TNNI1|chr1:201372896-201398994(-)"   
  # [5] "ENSG00000122180.5_3|MYOG|chr1:203052257-203055140(-)"      "ENSG00000168530.16_4|MYL1|chr2:211154874-211179898(-)"    
  # [7] "ENSG00000287059.1_1|AC090004.2|chr3:14080147-14088681(-)"  "ENSG00000185290.4_4|NUPR2|chr7:56182374-56184110(-)"      
  # [9] "ENSG00000285670.1_2|AC006970.3|chr7:56282672-56288440(+)"  "ENSG00000147573.17_5|TRIM55|chr8:67039131-67087720(+)"    
  # [11] "ENSG00000215182.6|MUC5AC|chr11:1151580-1222364(+)"         "ENSG00000129152.4_3|MYOD1|chr11:17741118-17743683(+)"     
  # [13] "ENSG00000230657.6_3|PRB4|chr12:11460017-11463369(-)"       "ENSG00000251655.6_3|PRB1|chr12:11504757-11548500(-)"      
  # [15] "ENSG00000121335.12_7|PRB2|chr12:11544474-11653975(-)"      "ENSG00000279134.1_5|AC090643.1|chr12:58488697-58491796(-)"
  # [17] "ENSG00000283361.2_7|CFAP97D2|chr13:114920166-114988559(+)" "ENSG00000226777.7_5|FAM30A|chr14:106383838-106398502(+)"  
  # [19] "ENSG00000260496.3_6|AC009041.1|chr16:1041151-1050926(-)"   "ENSG00000262152.7_5|LINC00514|chr16:3038257-3052017(+)"   
  # [21] "ENSG00000260034.1_6|LCMT1-AS2|chr16:25151898-25160353(-)"  "ENSG00000133020.4_3|MYH8|chr17:10293639-10325267(-)"      
  # [23] "ENSG00000128422.17_6|KRT17|chr17:39775694-39781094(-)"     "ENSG00000175894.18_7|TSPEAR|chr21:45917776-46131487(-)"   
  # [25] "ENSG00000187268.12_5|FAM9C|chrX:13053736-13062801(-)"     
  
  
} else {
  
  gsam.gene.res.res.paired <- readRDS("tmp/gsam.gene.res.res.paired.Rds")
  gsam.gene.res.tpc.res.paired <- readRDS("tmp/gsam.gene.res.tpc.res.paired.Rds")
  
}


results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res.paired, by = c('gid' = 'gid')) %>%
  dplyr::left_join(gsam.gene.res.tpc.res.paired, by = c('gid' = 'gid'))






plt <- results.out  %>%
  dplyr::mutate(show.label = !is.na(log2FoldChange.tpc.res.paired) & 
                  !is.na(log2FoldChange.tpc.res) & 
                  abs(log2FoldChange.tpc.res.paired - log2FoldChange.tpc.res) > 2
  )


#%>%
#  dplyr::filter(!is.na(log2FoldChange.res) & !is.na(statistic.cor.tpc)) %>%
#dplyr::mutate(is.limited.res = as.character(log2FoldChange.res > 3)) %>% # change pch to something that is limited
#dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3 , log2FoldChange.res)) %>%
#dplyr::mutate(is.limited.tpc.res = as.character(log2FoldChange.tpc.res > 3)) %>% # change pch to something that is limited
#dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3 , log2FoldChange.tpc.res))


# MYL1 & OPALIN & MYOD1 ??? 


p1 <- ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc , label=hugo_symbol
                      #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                      #shape = is.limited.res ,
                      #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  #geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  #scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  #scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  geom_text_repel(data =subset(plt, show.label == T),col="red") +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage") +
  xlim(-2.5,4)


p2 <- ggplot(plt, aes(x = log2FoldChange.tpc.res.paired , y = statistic.cor.tpc , label=hugo_symbol
                      #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                      #shape = is.limited.res ,
                      #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  #geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  #scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  #scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  geom_text_repel(data =subset(plt, show.label == T),col='red') +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage") +
  xlim(-2.5,4)


p1 + p2 




results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res.paired, by = c('gid' = 'gid')) %>%
  dplyr::left_join(gsam.gene.res.tpc.res.paired, by = c('gid' = 'gid'))
plt <- results.out  
p.corr <- ggplot(plt, aes(x = log2FoldChange.tpc.res.paired , y = log2FoldChange.tpc.res , label=hugo_symbol
                          #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                          #shape = is.limited.res ,
                          #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  geom_point(pch=19,cex=0.25, data =subset(plt, gid %in% a),col="red") +
  #geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  #scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  #scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  geom_text_repel(data =subset(plt, gid %in% a),col="red") +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage") + 
  xlim(-2,4) +
  ylim(-2,4) 

p.uncorr + p.corr
subset(plt, gid %in% a) $log2FoldChange.tpc.res
subset(plt, gid %in% a) $log2FoldChange.tpc.res.paired

# .paired
# > subset(plt, gid %in% a) $log2FoldChange.tpc.res
# [1]  2.616752  2.705314  2.839108        NA        NA        NA  4.413780 -2.573126  2.801220        NA  5.968508        NA  2.413724  8.257744  8.114305
# [16]  3.549533  3.056017  3.897042  4.550233 -2.336730  2.838127        NA  3.528260  2.586812        NA





lm(  log2FoldChange.tpc.res.paired ~ log2FoldChange.tpc.res  , data=plt)
lm(  log2FoldChange.tpc.res ~ log2FoldChange.tpc.res.paired   , data=plt)

c <- plt %>%
  dplyr::filter(!is.na(log2FoldChange.tpc.res) & !is.na(log2FoldChange.tpc.res.paired))
#& abs(log2FoldChange.tpc.res) < 2 & abs(log2FoldChange.tpc.res.paired) < 2 )
cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired  )
cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired, method = "spearman"  )

cor(  c$stat.tpc.res , c$stat.tpc.res.paired  )
cor(  c$stat.tpc.res , c$stat.tpc.res.paired, method = "spearman"  )



cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired  )
cor(  c$stat.tpc.res , c$stat.tpc.res.paired  )
c <- c %>% dplyr::filter(!is.na(padj.tpc.res) & !is.na (padj.tpc.res.paired  ) )
cor(  c$padj.tpc.res , c$padj.tpc.res.paired , method="" )
cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired, method = "spearman"  )
#cor(  c$stat.tpc.res , c$stat.tpc.res.paired, method = "spearman"  )





#ggplot(plt, aes(x = stat.tpc.res.paired , y = stat.tpc.res , label=hugo_symbol
ggplot(plt, aes(x = stat.tpc.res.paired , y = stat.tpc.res , label=hugo_symbol
                #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                #shape = is.limited.res ,
                #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(data=subset(plt, show.label == F), pch=19,cex=0.05) +
  geom_point(data=subset(plt, show.label == T), pch=19,cex=0.25,col="red") +
  geom_text_repel(data =subset(plt, show.label == T),col="red") +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage")
#+
#xlim(-3,10)+
#ylim(-3,10)



a = plt %>% dplyr::filter(show.label == T ) %>% pull(gid)


b = as.data.frame(gsam.gene.expression.all.vst) %>% 
  dplyr::filter(rownames(.) %in% a)

c = prcomp(t(b))
screeplot(c)


d = as.data.frame(c$x) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(in.paired.data = gid  %in% colnames(gsam.gene.expression.all.paired) ) %>%
  dplyr::mutate(res = gsub("^...(.).*$","R\\1",gid))


ggplot(d, aes(x=PC1, y=PC3, label=gid, col=in.paired.data) ) + 
  geom_point() + 
  geom_text_repel(data = subset(d, in.paired.data == F))



#AOA2
#BAE2
#AAX1
#EBB2

b = data.frame(
  gen = as.numeric(gsam.gene.expression.all.vst %>% as.data.frame() %>% dplyr::filter(rownames(.) == a[9]) ),
  sid=colnames(gsam.gene.expression.all.vst) ) %>%
  dplyr::mutate(in.paired.data = sid  %in% colnames(gsam.gene.expression.all.paired) ) %>%
  dplyr::mutate(res = gsub("^...(.).*$","R\\1",sid)) %>%
  dplyr::left_join(gsam.metadata.all %>% dplyr::select(c('sid', 'tpc')) , by = c('sid' = 'sid')) %>%
  dplyr::mutate(pid = as.factor(gsub("^(...).*$","\\1",sid)))


# AOA2, BAE2
ggplot(b, aes(x = tpc, y = gen, col=sid %in% c("AOA2","BAE2"), label=sid, shape = in.paired.data, group = pid)) +
  #geom_line(lwd=0.25,size=0.4,col="gray60") +
  geom_point(size=2)  +
  geom_text_repel(data = subset(b, gen > 5))



b1 <- b %>% dplyr::filter(in.paired.data) %>% dplyr::filter(res == 'R1')
b2 <- b %>% dplyr::filter(in.paired.data) %>% dplyr::filter(res == 'R2')
colnames(b1) <- paste0(colnames(b1), ".R1")
colnames(b2) <- paste0(colnames(b2), ".R2")
c = dplyr::left_join(b1, b2, by=c('pid.R1' = 'pid.R2')) %>%
  dplyr::mutate(gen.avg = (gen.R2 + gen.R1) / 2) %>%
  dplyr::mutate(dist = gen.R2 - gen.R1) %>%
  dplyr::mutate(dist.tpc = tpc.R2 - tpc.R1) %>%
  dplyr::arrange(dist) %>%
  dplyr::mutate(col = ifelse(sid.R1 %in% c('FAF1', 'AIA1','HAF1-new') , 1 , 2) )

## strong diff FAF2 / AIA2 / HAF2-new ?

plot(c$gen.R1 - c$gen.avg, y = 1:nrow(c), pch=19,col="red")
points(c$gen.R2 - c$gen.avg, y = 1:nrow(c), pch=19,col="blue")


plot(c$dist.tpc , c$dist, col=c$col, pch=19) # deze zou het meeste moeten zeggen
d <- b %>% dplyr::filter(in.paired.data == F)
points(d$tpc - median(c(c$tpc.R1 , c$tpc.R2)), d$gen - median(c(c$gen.R1 , c$gen.R2)) ,col=4, pch=19)
abline(h=0, col="gray30",lwd=0.3)
#plot((c$tpc.R1 + c$tpc.R2) / 2 , c$dist)
#plot(c$tpc.R1 , c$dist)
#plot(c$tpc.R2 , c$dist)


plot(b$tpc, b$gen, col=ifelse(b$res == "R1",2,4), pch=19)


# < : / : deprecated code : : : > -----

# plots ----

## 2.1 :: TPC violin ----


m1 <- gsam.metadata.all %>%
  dplyr::filter(!is.na(tumour.percentage.dna) & resection == "r1") %>%
  dplyr::pull(tumour.percentage.dna) %>%
  median()

m2 <- gsam.metadata.all %>%
  dplyr::filter(!is.na(tumour.percentage.dna) & resection == "r2") %>%
  dplyr::pull(tumour.percentage.dna) %>%
  median()


ggplot(gsam.metadata.all %>% dplyr::mutate(resection = ifelse(resection == 'r1', "Resection 1", "Resection 2"))
       , aes(x = resection, y = tumour.percentage.dna, col=resection)) +
  geom_violin() +
  geom_segment(aes(x=0.55, y=m1, xend=1.45, yend=m1),lty=1,lwd=0.15, col="gray40") +
  geom_segment(aes(x=1.55, y=m2, xend=2.45, yend=m2),lty=1,lwd=0.15, col="gray40") +
  geom_jitter( position=position_jitter(0.2), size=0.9) +
  ylim(0, 100) +
  labs(x = NULL, col = NULL, y = "WES estimate tumour cell percentage" ) +
  job_gg_theme
rm(m1,m2)


ggsave('output/figures/tpc_estimate_wes.png', width = 5.7, height = 4.3)



## 2.3 :: FC res x corr TPC G-SAM  [diag; chr7 ; chr10 ] ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(log2FoldChange.gsam.res > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2, 2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(col.chr7 = case_when(is.limited.gsam.res == T ~ "truncated" , chr == "chr7" ~ "at chr7", chr != "chr7" ~ "not at chr7")) %>%
  dplyr::mutate(col.chr10 = case_when(is.limited.gsam.res == T ~ "truncated" , chr == "chr10" ~ "at chr10", chr != "chr10" ~ "not at chr10")) 


p1 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = is.limited.gsam.res ,
                      size = is.limited.gsam.res  ,
                      col = is.limited.gsam.res   ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.gsam.res > 0.05),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.35))    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2,2) +
  ylim(-16.5,16.5)

p2 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = col.chr7,
                      size = col.chr7  ,
                      col = col.chr7  ,
                      fill =  col.chr7  
) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(data = subset(plt, chr != "chr7")) +
  geom_point(data = subset(plt, chr == "chr7", stroke=0.001)) +
  scale_shape_manual(values = c('truncated'= 4, 'not at chr7' = 19 , 'at chr7'= 21 )    ) +
  scale_size_manual(values = c('truncated'= 0.85, 'not at chr7' = 0.1 , 'at chr7'= 0.8 )    ) +
  scale_color_manual(values = c('truncated'= 'black', 'not at chr7' = rgb(0,0,0,0.35), 'at chr7'= rgb(0,0,0,0.5)  )    ) +
  scale_fill_manual(values = c('truncated'= 'black', 'not at chr7' = rgb(0,0,0,0.35), 'at chr7'= rgb(1,0,0,1)  )    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2,2) +
  ylim(-16.5,16.5)

p3 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = col.chr10,
                      size = col.chr10  ,
                      col = col.chr10  ,
                      fill =  col.chr10  
) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(data = subset(plt, chr != "chr10")) +
  geom_point(data = subset(plt, chr == "chr10", stroke=0.001)) +
  scale_shape_manual(values = c('truncated'= 4, 'not at chr10' = 19 , 'at chr10'= 21 )    ) +
  scale_size_manual(values = c('truncated'= 0.85, 'not at chr10' = 0.1 , 'at chr10'= 0.8 )    ) +
  scale_color_manual(values = c('truncated'= 'black', 'not at chr10' = rgb(0,0,0,0.35), 'at chr10'= rgb(0,0,0,0.5)  )    ) +
  scale_fill_manual(values = c('truncated'= 'black', 'not at chr10' = rgb(0,0,0,0.35), 'at chr10'= rgb(1,0,0,1)  )    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2,2) +
  ylim(-16.5,16.5)

p1 / p2 / p3


ggsave('output/figures/paper_dge_log2FoldChange.res_x_statistic.cor.tpc.png', width = 5.7 * 1.4 , height = 4 * 1.4 * 3)




## 2.5 :: G-SAM.res & GLASS.res ----


### a :: GLASS res x G-SAM TPC corr ----

plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res) & !is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) 

p1 <- ggplot(plt, aes(x = log2FoldChange.glass.res ,
                      y =  statistic.gsam.cor.tpc,
                    shape = is.limited.glass.res ,
                    size = is.limited.glass.res  ,
                    col = is.limited.glass.res  ) ) +
      geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
      geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
      geom_point() +
      geom_smooth(data = subset(plt, padj.glass.res > 0.05 &  is.limited.glass.res == "FALSE"),method="lm",
                  se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
      scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
      scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
      scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
      youri_gg_theme + 
      labs(x = "log2FC R1 vs. R2 (unpaired)",
           y="Correlation t-statistic with tumour percentage",
           shape = "Truncated at x-axis",
           size="Truncated at x-axis",
           col="Truncated at x-axis") +
      xlim(-2.5,2.5) +
      ylim(-16.5,16.5)




### b :: GLASS res x G-SAM res ----


plt <- results.out %>%
  dplyr::filter(!is.na(stat.glass.res) & !is.na(stat.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(stat.gsam.res) > 8)) %>% # change pch to something that is limited
  dplyr::mutate(stat.gsam.res = ifelse(stat.gsam.res > 8, 8 , stat.gsam.res)) %>%
  dplyr::mutate(stat.gsam.res = ifelse(stat.gsam.res < -8, -8 , stat.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(stat.glass.res) > 8)) %>% # change pch to something that is limited
  dplyr::mutate(stat.glass.res = ifelse(stat.glass.res > 8, 8 , stat.glass.res))  %>%
  dplyr::mutate(stat.glass.res = ifelse(stat.glass.res < -8, -8 , stat.glass.res))  %>%
  dplyr::mutate(is.limited = is.limited.gsam.res == "TRUE" | is.limited.glass.res == "TRUE")



p2 <- ggplot(plt, aes(x = stat.gsam.res ,
                y =  stat.glass.res ,
                shape = is.limited ,
                size = is.limited  ,
                col = is.limited  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, is.limited == F), method="lm", se = FALSE,  formula=y ~ x, col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
  youri_gg_theme + 
  labs(x = "DESeq2 stat G-SAM",
       y="DESeq2 stat GLASS",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") + 
  xlim(-8,8) +
  ylim(-8,8)




p1 / p2



ggsave('output/figures/paper_dge_glass_x_gsam_uncorercted.png', width = 5.7 * 1.4 , height = 4 * 1.4 * 2)



# ggplot(plt, aes(x = stat.gsam.res ,
#                 y =  stat.glass.res ,
#                 shape = is.limited ,
#                 size = is.limited  ,
#                 col = is.limited  ) ) +
#   geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
#   geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
#   geom_point() +
#   #geom_smooth(data = subset(plt, is.limited == F), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
#   scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
#   scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
#   scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
#   youri_gg_theme + 
#   labs(x = "log2FC R1 vs. R2 G-SAM",
#        y="log2FC R1 vs. R2 GLASS",
#        shape = "Truncated at x-axis",
#        size="Truncated at x-axis",
#        col="Truncated at x-axis") +
#   xlim(-2.5,2.5) +
#   ylim(-2.5, 2.5)



## 2.6 :: GLASS & ch19 ----

# wees zeker dat beide plots exact dezelfde genen bevatten


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res) & !is.na(log2FoldChange.gsam.res) & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res))  %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res))  %>%
  dplyr::mutate(is.chr19 = as.character(chr == "chr19")) %>%
  dplyr::mutate(size.type.glass = case_when(
    is.limited.glass.res == "TRUE" ~ "limited",
    is.chr19 == "TRUE" ~ "chr19",
    TRUE ~ "not chr19")) %>%
  dplyr::mutate(size.type.gsam = case_when(
    is.limited.gsam.res == "TRUE" ~ "limited",
    is.chr19 == "TRUE" ~ "chr19",
    TRUE ~ "not chr19"))



### a :: GLASS ----



plt$is.limited.gsam.res %>% as.factor %>% summary

plt$is.limited %>% as.factor %>% summary
plt$size.type %>% as.factor %>% summary


p1 <- ggplot(plt, aes(x = log2FoldChange.glass.res ,
                      y =  statistic.gsam.cor.tpc,
                      shape = is.limited.glass.res ,
                      size = size.type.glass  ,
                      col = is.chr19))   +
  geom_point(data=subset(plt, is.chr19 == 'FALSE')) +
  geom_point(data=subset(plt, is.chr19 == 'TRUE')) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19))  +
  scale_size_manual(values = c('chr19'= 0.65,
                               'not chr19' = 0.1,
                               'limited' = 0.86)) +
  scale_color_manual(values = c('FALSE'= 'black', 'TRUE' = rgb(1,0,0,0.85)))  +
  youri_gg_theme + 
  labs(x = "log2FC GLASS",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis")  +
  xlim(-2.5,2.5)



### b :: GLASS ----


p2 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                y =  statistic.gsam.cor.tpc,
                shape = is.limited.gsam.res ,
                size = size.type.gsam  ,
                col = is.chr19))   +
  geom_point(data=subset(plt, is.chr19 == 'FALSE')) +
  geom_point(data=subset(plt, is.chr19 == 'TRUE')) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19))  +
  scale_size_manual(values = c('chr19'= 0.65,
                               'not chr19' = 0.1,
                               'limited' = 0.86)) +
  scale_color_manual(values = c('FALSE'= 'black', 'TRUE' = rgb(1,0,0,0.85)))  +
  youri_gg_theme + 
  labs(x = "log2FC G-SAM",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis")  +
  xlim(-2.5,2.5)



p1 + p2

ggsave('output/figures/paper_dge_glass_x_gsam_chr19.png', width = 5.7 * 1.1 * 2 , height = 4 * 1.6 )




## 2.7 :: G-SAM FC corrected  ----



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(label.chr7 = case_when(chr == "chr7" ~ "chr7", TRUE ~ "remaining" )) %>%
  dplyr::mutate(label.chr10 = case_when(chr == "chr10" ~ "chr10", TRUE ~ "remaining" ))




### A :: Corrected ----


p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res ,
                y=statistic.gsam.cor.tpc ) ) + 
  geom_point(cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.05 &  is.limited.gsam.tpc.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  
  scale_color_manual(values = c('remaining'='black','chr7'='red') ) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       #,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-2.5,2.5)



### B :: Un-corrected ----


p2 <- ggplot(plt, aes(x=log2FoldChange.gsam.res ,
                y=statistic.gsam.cor.tpc ) ) + 
  geom_point(cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.res > 0.05 &  is.limited.gsam.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  
  scale_color_manual(values = c('remaining'='black','chr7'='red') ) +
  labs(x = "log2FC R1 vs. R2",
       y="Correlation t-statistic with tumour percentage"
       #,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-2.5,2.5)


p1 + p2


ggsave('output/figures/paper_dge_gsam_corrected_uncorrected.png', width = 5.7 * 1.1 * 2 , height = 4 * 1.6 )





## 2.10 Corrected LFC + labels ----

# signi in both corrected and uncorrected, but lfcSE outlier in corrected
# "TTLL10"     
# "MYO3A"     
# "AC090791.1" "AC090124.2" "AC090643.1" "AC009041.1" "LINC00514"  "KRT9"       "CKM"        "PI3"    


plt <- results.out %>%
  dplyr::filter(is.na(padj.glass.res) | (!is.na(padj.glass.res) & padj.glass.res  < 0.01) ) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c('CDK4','MDM2','GLI1','GLIS1','EGFR','MHMT','CD4', "CD55" ))
  #dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res > 0.3 & log2FoldChange.gsam.res < 0.01 )
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5)



ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res ,
                y=statistic.gsam.cor.tpc ,
                col=significant,
                label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.35) +
  geom_text_repel(data=subset(plt, show.label == T), col="blue", size=2.5 ,
                  nudge_x = 3.1, direction = "y", hjust = "left" #, lwd=0.5
                  ) +
  scale_color_manual(values = c('TRUE'='red','FALSE'=rgb(0,0,0,0.35))) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)



plot(density(results.out %>% filter(!is.na(lfcSE.gsam.tpc.res) ) %>% pull(lfcSE.gsam.tpc.res)))


## supervised clustering / reconstruction ----
# kijk of R1 & R2 scheiden; en de invloed van tpc








## g:profiler upregulated ----

# significantly DOWN regulated genes are related to blood vessels / angiogenesis
gsam.gene.res.combined %>%
  dplyr::filter(log2FoldChange.tpc.res < 0) %>%
  dplyr::filter(padj.tpc.res < 0.01) %>%
  #dplyr::top_n(500, -padj.tpc.res) %>%
  dplyr::pull(hugo_symbol.tpc.res)


## g:profiler downregulated top1000 ----

# significantly UP regulated genes are related to synaptic signalling and neurons / neuronal development / cell junctions
gsam.gene.res.combined %>%
  dplyr::filter(log2FoldChange.tpc.res > 0) %>%
  dplyr::filter(padj.tpc.res < 0.01) %>%
  dplyr::top_n(1000, -padj.tpc.res) %>%
  #dplyr::arrange(padj.tpc.res) %>%
  dplyr::pull(hugo_symbol.tpc.res)


## g:profiler all top1000 ----

# DE regulated genes are related to
gsam.gene.res.combined %>%
  dplyr::filter(padj.tpc.res < 0.01 & p.value >= 0.05) %>%
  dplyr::top_n(1000, -padj.tpc.res) %>%
  #dplyr::arrange(padj.tpc.res) %>%
  dplyr::pull(hugo_symbol.tpc.res)



## plot pharmaco genes ----


p1 <- ggplot(gsam.gene.res.combined, aes(x=log2FoldChange.res ,
                                         y=statistic,
                                         #shape=significant.tpc.res,
                                         col = pharma.relation,
                                         label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == F ), cex=0.05, pch=19) +
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #scale_shape_manual(values = c('TRUE'=1, 'FALSE' = 19) ) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.res == T & abs(log2FoldChange.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme

p2 <- ggplot(gsam.gene.res.combined, aes(x=log2FoldChange.tpc.res ,
                                         y= statistic,
                                         #shape=significant.tpc.res,
                                         col = pharma.relation,
                                         label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == F ), cex=0.05, pch=19) +
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #scale_shape_manual(values = c('TRUE'=1, 'FALSE' = 19) ) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.tpc.res == T & abs(log2FoldChange.tpc.res) > 1.5), size = 2.4 )  +
  geom_text_repel(data = subset(gsam.gene.res.combined, significant.tpc.res == T & pharma.relation), size = 5,
                  #box.padding = 1.5,
                  nudge_x = 3.1, direction = "y", hjust = "left" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


## plot synaptic signalling GO:0099537 ----

tss <- read.delim('data/gProfiler_trans-synaptic_signaling_GO.0099537.tsv',header=T, sep=",")
tss <- read.delim('data/gProfiler_Morphine_addiction_KEGG.05032.tsv',header=T, sep=",")


#summary(as.factor(gsub("^.+(chr[^:]+):.+$","\\1",plt %>% dplyr::filter(significant.tpc.res) %>% pull(gid.res) )))



plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = ensembl_id.res %in% tss$converted_alias ) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label == F & significant.res == T ), cex=0.15, pch=19, alpha=1) +
  geom_point(data=subset(plt, show.label == F & significant.res == F ), cex=0.15, pch=19, alpha=0.4) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label == F & significant.tpc.res == T ), cex=0.15, pch=19, alpha=1) +
  geom_point(data=subset(plt, show.label == F & significant.tpc.res == F ), cex=0.15, pch=19, alpha=0.4) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot glioma stem cells (GSC) markers ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("PROM1|CD44|FUT4|CD70|S100A4|ALDH1A3|NANOG$|^SOX2$|NES", hugo_symbol.tpc.res)) %>%
  #grepl("MGMT|(TP53)$|(EGFR)$|ABCB(1|11|4|5)$|CCN2|(HTRA3)$|GLI1|MTOR$", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot PI3K/Akt pathway  ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("AKT[0-9]$|TSC[1-9]$|^PTEN$|CDKN2[AB]$|^PIP[0-9]$|OAF|^MTOR|IRS1$|NFKB1", hugo_symbol.tpc.res)) %>%
  # grepl("MPG$|HMGA2|^POLB|PARP1$", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


#
gid <- 'FLT1|KDR'
rownames(gsam.gene.expression.all)[grepl(gid,rownames(gsam.gene.expression.all))]
plt$gid.res[grepl(gid,plt$gid.res )]



## plot BER mechanism  ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("MPG$|HMGA2|^POLB|PARP1$|HDAC", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


## plot MMR mechanism  ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("MLH1|MLH3|MSH2|MSH3|MSH6|PMS1|^PMS2$|^POLD[1-4]$|PCNA|^RPA$|HMGB1$|RFC|LIG1", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


## plot MGMT ----

# https://www.sciencedirect.com/science/article/pii/S2352304216300162?via%3Dihub


plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
  #                grepl("PROM1|CD44|FUT4|CD70|S100A4|ALDH1A3|NANOG$|^SOX2$|NES", hugo_symbol.tpc.res)) %>%
  grepl("MGMT|(TP53)$|(EGFR)$|ABCB(1|11|4|5)$|CCN2|(HTRA3)$|GLI1|MTOR$", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


gid <- 'SOX2'
rownames(gsam.gene.expression.all)[grepl(gid,rownames(gsam.gene.expression.all))]
plt$gid.res[grepl(gid,plt$gid.res )]



## plot transporter genes


## hh signalling ? -----

## subtype genes ----
### subtype genes [CL] ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Classical") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Classical") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2

### subtype genes [PN] ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Proneural") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Proneural") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2




### subtype genes [MES] ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot ABCB1 ----


plt <- data.frame(ABCB1 = as.numeric(gsam.gene.expression.all.vst %>%
                                       as.data.frame() %>%
                                       dplyr::filter(grepl("ABCB1\\|", rownames(.)))) ,
                  res = gsub("^...(.).*$","R\\1",colnames( gsam.gene.expression.all.vst)) )

ggplot(plt, aes(x = res, y = ABCB1)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y = "ABCB1 relative expression")

rm(plt)

## plot MGMT ----


plt <- data.frame(MGMT = as.numeric(gsam.gene.expression.all.vst %>%
                                      as.data.frame() %>%
                                      dplyr::filter(grepl("TP53\\|", rownames(.)))) ,
                  res = gsub("^...(.).*$","R\\1",colnames( gsam.gene.expression.all.vst)))

ggplot(plt, aes(x = res, y = MGMT)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y = "MGMT relative expression")

rm(plt)



## plot MGMT + TP53 + p21 + Chk1 + Chk2 + ATM ----
## "Increased levels of MGMT or loss of the mismatch repair capacity confer resistance to temozolomide (3)" - doi:10.1158/1535-7163.MCT-05-0428


## plot GLI1 ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot Akt pathway [TMZ resist?] ----


## plot [BMPR1B, CTGF, CYP4F2, EDNRB, ELL, EPHA3, FOS, GJA1, GPM6A, HTRA3, IGFBP2, IGFBP7, LEF1, RAI, RGL, SEC3L1, SRP72, SSB1, ZNF436] doi:10.1158/1535-7163.MCT-05-0428 ----






# - - - - -

# 
# EnhancedVolcano(res,
#                 lab = res$hugo_symbol,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.01)





#save.image(file = "dge_supervised_expression_analysis.RData")
#load(file = "dge_supervised_expression_analysis.RData")

# TP53 is interesting

ggid <- 'SAA3'
gid <- res %>% dplyr::filter(hugo_symbol == ggid) %>% dplyr::pull('gid')
plt <- data.frame(sid = colnames(vst),  expr = as.numeric(as.data.frame(vst) %>% dplyr::filter(rownames(.) == gid))) %>%
  dplyr::left_join(metadata, by = c ('sid' = 'sid') )
ggplot(plt, aes(x = tumour.percentage.dna, y= expr,col=resection)) + #)) + # , 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  labs(y = paste0("Expression: ", ggid))




EnhancedVolcano(glass.gene.res.res,
                lab = glass.gene.res.res$gene_symbol ,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01)


## dplyr inner join (!)

combined.gene.res.res <- dplyr::inner_join(
  gsam.gene.res.res , glass.gene.res.res , by=c('ensembl_id.res'='ensembl_id.glass.res') ) %>%
  dplyr::inner_join( gsam.gene.res.tpc.res, by=c('ensembl_id.res'='ensembl_id.tpc.res') ) %>%
  dplyr::left_join(gsam.gene.expression.all.cor, by=c('ensembl_id.res' = 'ensembl_id')) %>%
  dplyr::mutate(chr = gsub("^.+(chr[^:]+):.+$","\\1",gid))



head(gsam.gene.res.res)
head(glass.gene.res.res)


dim(gsam.gene.res.res)

dim(combined.gene.res.res)


## plot GSAM x GLASS uncorrected ----

plt <- combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) # plotting limit


p1 <- ggplot(plt, aes(x=log2FoldChange.res,
                                        y= log2FoldChange.glass.res
                                        #,col=significant.tpc.res
                      )) +
  #geom_smooth (method = "lm",lty=2,lwd=0.7,se = F) +
  geom_vline(xintercept = 0, col="red", lwd=0.2) +
  geom_hline(yintercept = 0, col="red", lwd=0.2) +
  geom_point(pch=19,cex=0.3) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "log2FC R1 vs. R2 G-SAM (unpaired)",
       y= "log2FC R1 vs. R2 GLSS (unpaired)",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-2,3) +
  ylim(-2,3)

p2 <- ggplot(plt,  aes(x=log2FoldChange.tpc.res,
                       y= log2FoldChange.glass.res,
                       col = significant.tpc.res,
                       label=hugo_symbol.tpc.res)) +
  geom_point(pch=19,cex=0.3) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(combined.gene.res.res, significant.tpc.res == T & (log2FoldChange.tpc.res > 1.5 | log2FoldChange.tpc.res < 0.85)), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme


mx <- plt %>% dplyr::filter(chr == "chr19") %>% dplyr::pull(log2FoldChange.tpc.res) %>% median()
my <- plt %>% dplyr::filter(chr == "chr19") %>% dplyr::pull(log2FoldChange.glass.res) %>% median()

p3 <- ggplot(plt,  aes(x= log2FoldChange.tpc.res,
                       y= log2FoldChange.glass.res,
                       col = chr == "chr19",
                       label=hugo_symbol.tpc.res)) +
  geom_point(data = dplyr::filter(plt , chr != "chr19") , pch=19,cex=0.3) +
  geom_point(data = dplyr::filter(plt , chr == "chr19") , pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "FC GLASS",
       y="FC G-SAM",
       col="Difference significant (R1 ~ R2)"
  ) +
  geom_hline(yintercept=my) +
  geom_vline(xintercept=mx) +
  youri_gg_theme


p1 + p2 + p3





## de COL genen van het MES type gaan aanzienlijk omlaag


## plot GLASS x GSAM signi ----

plt <-  a = combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) %>% # plotting limit
  dplyr::mutate(significant = significant.tpc.res & padj.glass.res < 0.05 ) %>%
  dplyr::mutate(show.label = significant ) %>%
  dplyr::filter(significant == T) %>%
  dplyr::pull(hugo_symbol.tpc.res) 
  
wang.glioma.intrinsic.genes%>% dplyr::filter (Gene_Symbol %in% a )


ggplot(plt,  aes(x= log2FoldChange.tpc.res,
                 y= log2FoldChange.glass.res,
                 col = show.label,
                 label=hugo_symbol.tpc.res)) +
  geom_point(data = dplyr::filter(plt , show.label == F) , pch=19,cex=0.3) +
  geom_point(data = dplyr::filter(plt , show.label == T) , pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "FC GLASS",
       y="FC G-SAM",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme



## plot chr [chr19]


combined.gene.res.res %>%
  dplyr::filter(significant.glass.res) %>%
  dplyr::pull(chr) %>%
  as.factor %>%
  summary()

plt <- combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) # plotting limit


# clearly chr19 is more to the right
chr.candidate = "chr19"
ggplot(combined.gene.res.res, aes(x= log2FoldChange.glass.res,
                                  y= statistic,
                                  label=hugo_symbol.tpc.res,
                                  #col=significant.glass.res,
                                  #col = ensembl_id.res %in% (wang.glioma.intrinsic.genes %>% dplyr::filter( Subtyping_Signature_Gene. == "Mesenchymal") %>% dplyr::pull(ENSG.short))
                                  col = chr == chr.candidate
)) +
  geom_point(data = subset(combined.gene.res.res, chr != chr.candidate), pch=19,cex=0.3) +
  geom_point(data = subset(combined.gene.res.res, chr == chr.candidate), pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(combined.gene.res.res, significant.glass.res == T & abs(log2FoldChange.glass.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme



## plot GLI1 ----



plt <- combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) %>% # plotting limit
  dplyr::mutate(show.label = grepl("GLI1", hugo_symbol.tpc.res))


ggplot(plt, aes(x= log2FoldChange.glass.res,
                                  y= statistic,
                                  label=hugo_symbol.tpc.res,
                                  col = show.label )) +
  geom_point(data = subset(plt, show.label == F), pch=19,cex=0.3) +
  geom_point(data = subset(plt, show.label == T), pch=19,cex=1.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label), size = 3.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme

