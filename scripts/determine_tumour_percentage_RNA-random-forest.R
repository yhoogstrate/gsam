#!/usr/bin/env R


# load libs ----


library(tidyverse)
library(randomForest)
library(caret)
library(ranger)


# load data ----

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


#source("scripts/R/ligands.R")
#source("scripts/R/subtype_genes.R") # Verhaak/Wang/TCGA + GliTS Redux
#source("scripts/R/patel_scRNA_clusters.R")
#source("scripts/R/neftel_meta_modules.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

#source('scripts/R/wang_glioma_intrinsic_genes.R')

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set + metedata


# load data ----


## GSAM ----



gsam.metadata.all <- gsam.rna.metadata %>%
  dplyr::filter(blacklist.pca == F) %>%
  dplyr::filter(pat.with.IDH == F) %>%
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
  dplyr::filter(tumour.percentage.dna >= 15) %>%
  #dplyr::mutate(tpc = 1 - (tumour.percentage.dna / 100)) %>%
  dplyr::mutate(sid = paste0('GSAM-', sid)) %>%
  dplyr::select(c('sid', tumour.percentage.dna))

gsam.gene.expression.all <- gsam.rnaseq.expression %>%
  `colnames<-`(paste0('GSAM-', colnames(.))) %>%
  dplyr::select(gsam.metadata.all$sid)

stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)



## GLASS ----



glass.metadata.all  <- glass.gbm.rnaseq.metadata

glass.gene.expression.all <- glass.gbm.rnaseq.expression %>%
  dplyr::select(glass.metadata.all$sid)

stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)




## per-gene results table ----



gene.metadata <- dplyr::full_join(
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
  dplyr::filter(in.gsam == T & in.glass == T)

stopifnot(sum(duplicated(gene.metadata$ensembl_id)) == 0)



all.gene.expression <- gene.metadata %>%
  dplyr::select(c('gid','ensembl_id')) %>%
  dplyr::left_join(gsam.gene.expression.all %>% tibble::rownames_to_column('gid'), by=c('gid' = 'gid')) %>%
  dplyr::left_join(glass.gene.expression.all %>% tibble::rownames_to_column('ensembl_id'), by=c('ensembl_id' = 'ensembl_id')) %>%
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('ensembl_id')


all.gene.expression.vst <- all.gene.expression %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  limma::removeBatchEffect(as.factor(gsub("^(..).*$","\\1",colnames(.)) == "GS")) %>% # remove batch effect :)
  as.data.frame()


# features ----

# top100 cor w/ tpc in glass (50pos & 50neg)

features1 <- c("ENSG00000168137", "ENSG00000167785", "ENSG00000171649", "ENSG00000198521", "ENSG00000204514", "ENSG00000175414", "ENSG00000066117", "ENSG00000106344",
               "ENSG00000188321", "ENSG00000078177", "ENSG00000197299", "ENSG00000008311", "ENSG00000197128", "ENSG00000138443", "ENSG00000105708", "ENSG00000197928",
               "ENSG00000198799", "ENSG00000131115", "ENSG00000103037", "ENSG00000204519", "ENSG00000198551", "ENSG00000189079", "ENSG00000007392", "ENSG00000135164",
               "ENSG00000105866", "ENSG00000164828", "ENSG00000105486", "ENSG00000167637", "ENSG00000254004", "ENSG00000177839", "ENSG00000162086", "ENSG00000134744",
               "ENSG00000120784", "ENSG00000147274", "ENSG00000129351", "ENSG00000151612", "ENSG00000166704", "ENSG00000106443", "ENSG00000004139", "ENSG00000122779",
               "ENSG00000160961", "ENSG00000167380", "ENSG00000122970", "ENSG00000071575", "ENSG00000263001", "ENSG00000236609", "ENSG00000197647", "ENSG00000104885",
               "ENSG00000113387", "ENSG00000167384")

features2 <- c("ENSG00000111885", "ENSG00000141506", "ENSG00000166272", "ENSG00000198624", "ENSG00000197746", "ENSG00000151726", "ENSG00000107968", "ENSG00000266412",
               "ENSG00000131370", "ENSG00000100365", "ENSG00000082397", "ENSG00000204161", "ENSG00000108639", "ENSG00000248905", "ENSG00000084070", "ENSG00000110079",
               "ENSG00000155252", "ENSG00000134996", "ENSG00000130775", "ENSG00000122359", "ENSG00000110324", "ENSG00000165457", "ENSG00000155926", "ENSG00000134516",
               "ENSG00000173372", "ENSG00000163131", "ENSG00000153071", "ENSG00000101336", "ENSG00000175155", "ENSG00000128805", "ENSG00000148180", "ENSG00000142185",
               "ENSG00000111912", "ENSG00000155659", "ENSG00000137462", "ENSG00000180353", "ENSG00000147459", "ENSG00000136250", "ENSG00000235568", "ENSG00000198879",
               "ENSG00000197142", "ENSG00000143119", "ENSG00000138964", "ENSG00000167613", "ENSG00000183484", "ENSG00000130830", "ENSG00000101160", "ENSG00000205744",
               "ENSG00000100368", "ENSG00000012779")




# Re-predict TPC for G-SAM ----


metadata <- gsam.metadata.all %>%
  dplyr::select(c(sid, tumour.percentage.dna))

expression.data <- all.gene.expression.vst %>%
  dplyr::select(metadata$sid)


## randomForest [n.5000 s.1+3+3+7] ----



test.out.all <- data.frame()

for(i in 1:10) {
  print(i)
  
  train.metadata <- metadata[1:ncol(expression.data) %% 10 != (i-1),]
  test.metadata <- metadata[1:ncol(expression.data) %% 10 == (i-1),]
    
  train.expression.data <- expression.data %>% dplyr::select(train.metadata$sid)
  test.expression.data <- expression.data %>% dplyr::select(test.metadata$sid)
  
  train.data <- train.expression.data[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(tumour.percentage.dna = train.metadata$tumour.percentage.dna)
  
  test.data <- test.expression.data[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(tumour.percentage.dna = test.metadata$tumour.percentage.dna)
  

  set.seed(1+3+3+7)
  
  test.rf <- randomForest(tumour.percentage.dna ~ .,
                          data=train.data,
                          ntree = 5000,
                          #mtry=3,
                          importance=TRUE
                          #na.action=na.omit
                          )

  test.data$tumour.percentage.dna.predicted <- predict(test.rf, test.data)
  
  test.out.all <- rbind(test.data, test.out.all)
}

# voeg toe aan metadata
gsam.metadata.all <- gsam.metadata.all %>%
  dplyr::left_join(
    test.out.all %>%
      dplyr::select('tumour.percentage.dna.predicted') %>%
      tibble::rownames_to_column('sid') %>%
      dplyr::rename(tumour.percentage.dna.predicted.rf.n5000.s14.imp = tumour.percentage.dna.predicted),
    by = c('sid'='sid'))


plt <- gsam.metadata.all %>%
  dplyr::mutate(tumour.percentage.dna.pred = tumour.percentage.dna.predicted.rf.n5000.s14.imp) %>%
  dplyr::filter(!is.na(tumour.percentage.dna.pred) & !is.na(tumour.percentage.dna.pred))

ggplot(data=plt, aes(x=tumour.percentage.dna, y=tumour.percentage.dna.pred)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100)
cor(plt$tumour.percentage.dna , plt$tumour.percentage.dna.pred)



#plot(test.out.all$tumour.percentage.dna , test.out.all$tumour.percentage.dna.predicted)
#cor(test.out.all$tumour.percentage.dna , test.out.all$tumour.percentage.dna.predicted)
##n=5000, seed=1+3+3+7: 0.7617365


## caret ----


test.out.all <- data.frame()

for(i in 1:10) {
  print(i)
  
  train.metadata <- metadata[1:ncol(expression.data) %% 10 != (i-1),]
  test.metadata <- metadata[1:ncol(expression.data) %% 10 == (i-1),]
  
  train.expression.data <- expression.data %>% dplyr::select(train.metadata$sid)
  test.expression.data <- expression.data %>% dplyr::select(test.metadata$sid)
  
  train.data <- train.expression.data[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(tumour.percentage.dna = train.metadata$tumour.percentage.dna)
  
  test.data <- test.expression.data[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(tumour.percentage.dna = test.metadata$tumour.percentage.dna)
  
  
  set.seed(1+3+3+7)
  
  test.caret <- caret::train(
    x = train.data %>% dplyr::mutate(tumour.percentage.dna = NULL),
    y = train.data$tumour.percentage.dna ,
    trControl = caret::trainControl(
      returnData = TRUE,
      method = "repeatedcv",
      number = 5,
      repeats= 15
    ),
    method = "ranger",
    importance = "permutation"
  )
  
  
  test.data$tumour.percentage.dna.predicted <- predict(test.caret, test.data)
  
  test.out.all <- rbind(test.data, test.out.all)
  
}


# voeg toe aan metadata
gsam.metadata.all <- gsam.metadata.all %>%
  dplyr::left_join(
    test.out.all %>%
      dplyr::select('tumour.percentage.dna.predicted') %>%
      tibble::rownames_to_column('sid') %>%
      dplyr::rename(tumour.percentage.dna.predicted.caret.5.5 = tumour.percentage.dna.predicted),
    by = c('sid'='sid'))


plt <- gsam.metadata.all %>%
  dplyr::mutate(tumour.percentage.dna.pred = tumour.percentage.dna.predicted.caret.5.5) %>%
  dplyr::filter(!is.na(tumour.percentage.dna.pred) & !is.na(tumour.percentage.dna.pred))

ggplot(data=plt, aes(x=tumour.percentage.dna, y=tumour.percentage.dna.pred)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100)
cor(plt$tumour.percentage.dna , plt$tumour.percentage.dna.pred)
#0.7581211 [caret is echt beter dan RF]




# Predict GLASS ----

## randomForest ----


train.data <- all.gene.expression.vst %>%
  dplyr::select(gsam.metadata.all$sid) %>%
  dplyr::filter(rownames(.) %in% c(features1, features2)) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(metadata %>% dplyr::select(c('sid','tumour.percentage.dna')) ,
    by = c('sid'='sid')) %>%
  tibble::column_to_rownames('sid')


test.data <- all.gene.expression.vst %>%
  dplyr::select(glass.metadata.all$sid) %>%
  dplyr::filter(rownames(.) %in% c(features1, features2)) %>%
  t() %>%
  as.data.frame()


test.rf <- randomForest(tumour.percentage.dna ~ .,
                        data=train.data,
                        ntree = 5000 )

test.data$tumour.percentage.dna.predicted <- predict(test.rf, test.data)



## caret ----


train.data <- all.gene.expression.vst %>%
  dplyr::select(gsam.metadata.all$sid) %>%
  dplyr::filter(rownames(.) %in% c(features1, features2)) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(metadata %>% dplyr::select(c('sid','tumour.percentage.dna')) ,
                   by = c('sid'='sid')) %>%
  tibble::column_to_rownames('sid')


test.data <- all.gene.expression.vst %>%
  dplyr::select(glass.metadata.all$sid) %>%
  dplyr::filter(rownames(.) %in% c(features1, features2)) %>%
  t() %>%
  as.data.frame()


set.seed(1+3+3+7)

test.caret <- caret::train(
  x = train.data %>% dplyr::mutate(tumour.percentage.dna = NULL),
  y = train.data$tumour.percentage.dna ,
  trControl = caret::trainControl(
    returnData = TRUE,
    method = "repeatedcv",
    number = 5,
    repeats= 15
  ),
  method = "ranger",
  importance = "permutation"
)


p <- data.frame(sid = rownames(test.data),
                tumour.percentage.dna.imputed.caret = predict(test.caret, test.data) )

write.table(p, file="output/tables/GLASS.tumour.percentage.dna.imputed.caret.txt")


#  From 15 runs it's pretty consistent


# set.seed(1+3+3+8)
# 
# test.caret <- caret::train(
#   x = train.data %>% dplyr::mutate(tumour.percentage.dna = NULL),
#   y = train.data$tumour.percentage.dna ,
#   trControl = caret::trainControl(
#     returnData = TRUE,
#     method = "repeatedcv",
#     number = 5,
#     repeats= 15
#   ),
#   method = "ranger",
#   importance = "permutation"
# )
# 
# 
# q <- data.frame(sid = rownames(test.data),
#                 tumour.percentage.dna.imputed.caret = predict(test.caret, test.data) )
# 
# 
# plot(p$tumour.percentage.dna.imputed.caret , q$tumour.percentage.dna.imputed.caret)
# 

