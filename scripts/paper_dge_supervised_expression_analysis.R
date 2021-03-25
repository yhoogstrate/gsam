# Set wd

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


gsam.bfg.expression.all <- gsam.rnaseq.expression %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
  dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
  dplyr::filter( # bona fide glioma genes (BFGs)
    (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
      (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
  ) %>%
  tibble::column_to_rownames('gid') %>%
  dplyr::select(gsam.metadata.all$sid)

stopifnot(colnames(gsam.bfg.expression.all) == gsam.metadata.all$sid)


# gsam.gene.expression.all.vst
# TODO gsam.gene.expression.R1.vst
# TODO gsam.gene.expression.R2.vst

# TODO gsam.bfg.expression.all.vst
# TODO gsam.bfg.expression.R1.vst
# TODO gsam.bfg.expression.R2.vst



gsam.gene.expression.all.vst <- gsam.gene.expression.all %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = paste0("r",round(runif( ncol(.) , 1, 2)))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay()



ggplot(gsam.metadata.all %>% dplyr::mutate(resection = ifelse(resection == 'r1', "Resection 1", "Resection 2"))
       , aes(x = resection, y = tumour.percentage.dna, col=resection)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylim(0, 100) +
  labs(x = "Resection", col = NULL, y = "WES estimate tumour cell percentage" ) +
  job_gg_theme


## GLASS ----


# GLSS-SM-R099

glass.metadata.all  <- glass.gbm.rnaseq.metadata %>%
  dplyr::mutate(condition = ifelse(resection == "TP","Primary","NotPrimary")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Primary","NotPrimary") ))
#dplyr::filter(pid != "GLSS-SM-R099")  %>%
#dplyr::pull(condition) %>% as.factor %>% summary()
#  dplyr::arrange(pid, condition, sid)

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
  `colnames<-`(paste0(colnames(.), ".cor.tpc")) %>%
  dplyr::rename(gid = gid.cor.tpc)

rm(gsam.gene.expression.all.cor.estimate , gsam.gene.expression.all.cor.statistic , gsam.gene.expression.all.cor.p.value)



results.out <- results.out %>%
  dplyr::left_join(gsam.gene.expression.all.cor , by = c('gid' = 'gid')) 



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




# DE unpaired ~ TPC [GSAM] ----

# quite sketchy test - sensitive to outliers and many NA corrected padj's

stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)

gsam.res.tpc.all <- DESeqDataSetFromMatrix(gsam.gene.expression.all, gsam.metadata.all, ~tpc ) %>% # + tpc
  DESeq(parallel=T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  tibble::rownames_to_column('gid') %>%
  `colnames<-`(paste0(colnames(.),".tpc"))  %>%
  dplyr::rename(gid = gid.tpc)
  
#sum(is.na ( gsam.res.tpc.all$padj.tpc ))
#sum(is.na ( gsam.res.tpc.all$p.value.tpc ))

results.out <- results.out %>%
  dplyr::left_join(gsam.res.tpc.all , by = c('gid' = 'gid')) 



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


# DE unpaired all [GSAM] ----



gsam.gene.res.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                              dplyr::filter(rowSums(.) > ncol(.) * 3)
                                            , gsam.metadata.all, ~resection ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  #dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
  #dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
  `colnames<-`(paste0(colnames(.),".res")) %>%
  dplyr::rename(gid = gid.res)


gsam.gene.res.tpc.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                                  dplyr::filter(rowSums(.) > ncol(.) * 3), gsam.metadata.all, ~tpc + resection ) %>% # + resection; corrected for tpc
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  #dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
  #dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
  `colnames<-`(paste0(colnames(.),".tpc.res")) %>%
  dplyr::rename(gid = gid.tpc.res)



results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res , by = c('gid' = 'gid')) %>%
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




# plots ----

## FC un-corrected x correlation TPC G-SAM ----



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.res) & !is.na(statistic.cor.tpc)) %>%
  dplyr::mutate(is.limited = as.character(log2FoldChange.res > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3 , log2FoldChange.res))


p1 <- ggplot(plt, aes(x = log2FoldChange.res , y = statistic.cor.tpc, pch=is.limited,
                      shape = is.limited,
                      size = is.limited   ) ) +
#p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",  se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=1) +
  scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  youri_gg_theme



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.res) & !is.na(statistic.cor.tpc)) %>%
  dplyr::mutate(is.limited = as.character(log2FoldChange.tpc.res > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3 , log2FoldChange.tpc.res))


p2 <- ggplot(plt, aes(x = log2FoldChange.tpc.res ,
                      y =  statistic.cor.tpc,
                      shape = is.limited,
                      size = is.limited ) ) +
#p2 <- ggplot(plt, aes(x = stat.tpc.res , y =  statistic.cor.tpc ) ) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.tpc.res > 0.01),method="lm",  se = FALSE, formula=y ~ x, orientation="y", col="red" , size=1) +
  scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  youri_gg_theme

p1 + p2





## FC corrected x correlation TPC ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3 ,log2FoldChange.tpc.res)) %>%
  dplyr::mutate(show.label = log2FoldChange.tpc.res > 0 & padj.tpc.res < 0.000383)

a = plt %>% dplyr::filter(show.label) %>% dplyr::pull(ensembl_id)


ggplot(plt, aes(x=log2FoldChange.tpc.res ,  y=statistic, 
                col=show.label
) ) + 
  geom_point(data=subset( plt, show.label == F ),pch=19,cex=0.05) +
  geom_point(data=subset( plt, show.label == T ),pch=19,cex=0.7) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Top 500 over-expressed"
  ) +
  youri_gg_theme



plt <- gsam.gene.res.combined %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3 ,log2FoldChange.tpc.res)) %>%
  dplyr::mutate(show.label = log2FoldChange.tpc.res < 0 & padj.tpc.res < 0.0089)

a = plt %>% dplyr::filter(show.label) %>% dplyr::pull(ensembl_id) 



ggplot(plt, aes(x=log2FoldChange.tpc.res ,  y=statistic, 
                col=show.label
) ) + 
  geom_point(data=subset( plt, show.label == F ),pch=19,cex=0.05) +
  geom_point(data=subset( plt, show.label == T ),pch=19,cex=0.7) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Top 500 under-expressed"
  ) +
  youri_gg_theme





p2 <- ggplot(gsam.gene.res.combined, aes(x=log2FoldChange.res ,
                                         y=statistic,
                                         #col=significant.res,
                                         col = is.bfg,
                                         label=hugo_symbol.tpc.res
) ) + 
  geom_point(data=subset(gsam.gene.res.combined, is.bfg == T ),pch=19,cex=0.05) +
  #geom_point(data=subset(gsam.gene.res.combined, significant.res == T ),pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.res == T & abs(log2FoldChange.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme

p1 + p2


## plot Fc corrected x correlation TPC +chr7 + chr10 ----


plt <- gsam.gene.res.combined %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 2.5, 2.5, log2FoldChange.res))%>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))

p1 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y=statistic, label=hugo_symbol.tpc.res ) ) + 
  geom_point( pch=19,cex=0.05) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #geom_point(data=subset(gsam.gene.res.combined, significant.res == T ),pch=19,cex=0.5) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.res == T & abs(log2FoldChange.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired; tpc corrected)",
       y="Correlation t-statistic with tumour percentage"
       #,col="Difference significant (R1 ~ R2)"
  ) +
  #geom_smooth(method="lm") + 
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col=chr == "chr7", label=hugo_symbol.tpc.res ) ) + 
  geom_point(data=subset(plt, chr != "chr7"), pch=19,cex=0.05) +
  geom_point(data=subset(plt, chr == "chr7"), pch=19,cex=0.65) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #geom_point(data=subset(gsam.gene.res.combined, significant.res == T ),pch=19,cex=0.5) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.res == T & abs(log2FoldChange.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired; uncorrected)",
       y="Correlation t-statistic with tumour percentage",
       col="Is at chr7"
  ) +
  #geom_smooth(method="lm") + 
  youri_gg_theme

p3 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col=chr == "chr10", label=hugo_symbol.tpc.res ) ) + 
  geom_point(data=subset(plt, chr != "chr10"), pch=19,cex=0.05) +
  geom_point(data=subset(plt, chr == "chr10"), pch=19,cex=0.65) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #geom_point(data=subset(gsam.gene.res.combined, significant.res == T ),pch=19,cex=0.5) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.res == T & abs(log2FoldChange.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired; uncorrected)",
       y="Correlation t-statistic with tumour percentage",
       col="Is at chr10"
  ) +
  #geom_smooth(method="lm") + 
  youri_gg_theme

p1 + p2 + p3



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


## plot chr [ch19] ----



plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = F ) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  %>%
  dplyr::filter(abs(padj.tpc.res) > 0.1)

m1 <- plt %>%
  dplyr::filter(chr == "chr19") %>%
  dplyr::pull(log2FoldChange.tpc.res) %>% median()

m2 <- plt %>%
  dplyr::filter(chr != "chr19") %>%
  dplyr::pull(log2FoldChange.tpc.res) %>% median()




p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = chr == "chr19" , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, chr != "chr19" ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, chr == "chr19" ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  geom_smooth(data = subset(plt, abs(log2FoldChange.res) <= 0.5), method='lm', formula= y~x) +
  youri_gg_theme

p2 <- ggplot(plt ,
             aes(x=log2FoldChange.tpc.res , y= statistic, col = chr == "chr19"  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, chr != "chr19" ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, chr == "chr19" ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  geom_vline(xintercept = m1) + 
  geom_vline(xintercept = m2) +
  youri_gg_theme

p1 + p2



## supervised clustering / reconstruction ----
# kijk of R1 & R2 scheiden; en de invloed van tpc






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



# PAIRED analysis ----


