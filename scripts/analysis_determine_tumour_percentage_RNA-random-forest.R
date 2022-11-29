#!/usr/bin/env R


# load libs ----


#library(tidyverse)
#library(randomForest)
#library(caret)
#library(ranger)


# load data ----

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set + metadata




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

# # top100 cor w/ tpc in glass (50pos & 50neg)
# 
# features1 <- c("ENSG00000168137", "ENSG00000167785", "ENSG00000171649", "ENSG00000198521", "ENSG00000204514", "ENSG00000175414", "ENSG00000066117", "ENSG00000106344",
#                "ENSG00000188321", "ENSG00000078177", "ENSG00000197299", "ENSG00000008311", "ENSG00000197128", "ENSG00000138443", "ENSG00000105708", "ENSG00000197928",
#                "ENSG00000198799", "ENSG00000131115", "ENSG00000103037", "ENSG00000204519", "ENSG00000198551", "ENSG00000189079", "ENSG00000007392", "ENSG00000135164",
#                "ENSG00000105866", "ENSG00000164828", "ENSG00000105486", "ENSG00000167637", "ENSG00000254004", "ENSG00000177839", "ENSG00000162086", "ENSG00000134744",
#                "ENSG00000120784", "ENSG00000147274", "ENSG00000129351", "ENSG00000151612", "ENSG00000166704", "ENSG00000106443", "ENSG00000004139", "ENSG00000122779",
#                "ENSG00000160961", "ENSG00000167380", "ENSG00000122970", "ENSG00000071575", "ENSG00000263001", "ENSG00000236609", "ENSG00000197647", "ENSG00000104885",
#                "ENSG00000113387", "ENSG00000167384")
# 
# features2 <- c("ENSG00000111885", "ENSG00000141506", "ENSG00000166272", "ENSG00000198624", "ENSG00000197746", "ENSG00000151726", "ENSG00000107968", "ENSG00000266412",
#                "ENSG00000131370", "ENSG00000100365", "ENSG00000082397", "ENSG00000204161", "ENSG00000108639", "ENSG00000248905", "ENSG00000084070", "ENSG00000110079",
#                "ENSG00000155252", "ENSG00000134996", "ENSG00000130775", "ENSG00000122359", "ENSG00000110324", "ENSG00000165457", "ENSG00000155926", "ENSG00000134516",
#                "ENSG00000173372", "ENSG00000163131", "ENSG00000153071", "ENSG00000101336", "ENSG00000175155", "ENSG00000128805", "ENSG00000148180", "ENSG00000142185",
#                "ENSG00000111912", "ENSG00000155659", "ENSG00000137462", "ENSG00000180353", "ENSG00000147459", "ENSG00000136250", "ENSG00000235568", "ENSG00000198879",
#                "ENSG00000197142", "ENSG00000143119", "ENSG00000138964", "ENSG00000167613", "ENSG00000183484", "ENSG00000130830", "ENSG00000101160", "ENSG00000205744",
#                "ENSG00000100368", "ENSG00000012779")
# 






# Re-predict TPC for G-SAM ----


metadata <- gsam.metadata.all %>%
  dplyr::select(c(sid, tumour.percentage.dna))

expression.data <- all.gene.expression.vst %>%
  dplyr::select(metadata$sid)


## randomForest [n.5000 s.1+3+3+7] ----


test.out.all <- data.frame()

for(i in 1:10) {
  print(i)
  #i = 3

  train.metadata <- metadata[1:ncol(expression.data) %% 10 != (i-1),]
  test.metadata <- metadata[1:ncol(expression.data) %% 10 == (i-1),]

  train.expression.data <- expression.data %>% dplyr::select(train.metadata$sid)
  test.expression.data <- expression.data %>% dplyr::select(test.metadata$sid)

  train.data <- train.expression.data %>% # [rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(tumour.percentage.dna = train.metadata$tumour.percentage.dna)

  test.data <- test.expression.data %>% #[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(tumour.percentage.dna = test.metadata$tumour.percentage.dna)


  set.seed(1+3+3+7)
  k = 150
  candidate.features <- readRDS(file = 'tmp/results.out.Rds') %>%
    dplyr::select(c(ensembl_id, statistic.gsam.cor.tpc)) %>%
    dplyr::filter(ensembl_id %in% rownames(expression.data)) %>%
    dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>%
    dplyr::arrange(statistic.gsam.cor.tpc) %>%
    dplyr::mutate(top.low = c(rep(T, k), rep(F, nrow(.) - k))) %>%
    dplyr::arrange(-statistic.gsam.cor.tpc) %>%
    dplyr::mutate(top.high = c(rep(T, k), rep(F, nrow(.) - k))) %>%
    dplyr::filter(top.high | top.low)

  features <- Boruta::Boruta(train.data$tumour.percentage.dna~. ,
                             data = train.data %>%
                               dplyr::select(
                                 candidate.features %>% dplyr::pull(ensembl_id)
                               )) %>%
    purrr::pluck('finalDecision') %>%
    as.data.frame() %>%
    dplyr::rename(boruta.status = '.') %>%
    tibble::rownames_to_column('ensembl_id') %>%
    dplyr::filter(boruta.status %in% c("Confirmed","Tentative")) %>%
    dplyr::left_join(candidate.features, by=c('ensembl_id' = 'ensembl_id'))

  dim(features)
  dim(train.data)
  dim(test.data)

  train.data <- train.data %>%
    dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna'))

  test.data <- test.data %>%
    dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna'))

  dim(train.data)
  dim(test.data)


  test.rf <- randomForest(tumour.percentage.dna ~ .,
                          data=train.data,
                          ntree = 5000,
                          #mtry=3,
                          
                          importance=TRUE
                          
                          #na.action=na.omit
                          )

  test.data$tumour.percentage.dna.predicted <- predict(test.rf, test.data)

  test.out.all <- rbind(test.data %>% dplyr::select('tumour.percentage.dna', 'tumour.percentage.dna.predicted'), test.out.all)

}




# voeg toe aan metadata
gsam.metadata.all <- gsam.metadata.all %>%
  dplyr::mutate(tumour.percentage.dna.predicted = NULL) %>% # force overwrite
  dplyr::left_join(
    test.out.all %>%
      dplyr::select('tumour.percentage.dna.predicted') %>%
      tibble::rownames_to_column('sid'),
    by = c('sid'='sid'))

#gsam.metadata.all$tumour.percentage.dna.predicted


plt <- gsam.metadata.all %>%
  dplyr::filter(!is.na(tumour.percentage.dna.predicted)) %>%
  dplyr::select(sid, tumour.percentage.dna, tumour.percentage.dna.predicted)
#write.table(plt, file="output/tables/GSAM.tumour.percentage.dna.imputed.rf.txt") 




plt <- read.table('output/tables/GSAM.tumour.percentage.dna.imputed.rf.txt')

err <- ModelMetrics::rmse(plt$tumour.percentage.dna , plt$tumour.percentage.dna.predicted)
ggplot(data=plt, aes(x=tumour.percentage.dna, y=tumour.percentage.dna.predicted)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100) +
  youri_gg_theme +
  geom_smooth(method="lm") + 
  labs(x = "Tumor cell percentage (WES)", y="Tumor cell percentage (Random forest re-fit on expression data)")
cor(plt$tumour.percentage.dna , plt$tumour.percentage.dna.predicted)

ggsave("output/figures/gsam_tpc_vs_imputed_tpc_rf5000.svg", width=7 , height=6)
ggsave("output/figures/gsam_tpc_vs_imputed_tpc_rf5000.pdf", width=7 , height=6)




# ## caret ----
# 
# 
# test.out.all <- data.frame()
# 
# for(i in 1:10) {
#   print(i)
#   
#   train.metadata <- metadata[1:ncol(expression.data) %% 10 != (i-1),]
#   test.metadata <- metadata[1:ncol(expression.data) %% 10 == (i-1),]
#   
#   train.expression.data <- expression.data %>% dplyr::select(train.metadata$sid)
#   test.expression.data <- expression.data %>% dplyr::select(test.metadata$sid)
#   
#   train.data <- train.expression.data %>% # [rownames(train.expression.data) %in% c(features1, features2),] %>%
#     t() %>%
#     as.data.frame() %>%
#     dplyr::mutate(tumour.percentage.dna = train.metadata$tumour.percentage.dna)
#   
#   test.data <- test.expression.data %>% # [rownames(train.expression.data) %in% c(features1, features2),] %>%
#     t() %>%
#     as.data.frame() %>%
#     dplyr::mutate(tumour.percentage.dna = test.metadata$tumour.percentage.dna)
#   
#   
#   set.seed(1+3+3+7)
#   k = 150
#   candidate.features <- readRDS(file = 'tmp/results.out.Rds') %>%
#     dplyr::select(c(ensembl_id, statistic.gsam.cor.tpc)) %>%
#     dplyr::filter(ensembl_id %in% rownames(expression.data)) %>%
#     dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>%
#     dplyr::arrange(statistic.gsam.cor.tpc) %>%
#     dplyr::mutate(top.low = c(rep(T, k), rep(F, nrow(.) - k))) %>%
#     dplyr::arrange(-statistic.gsam.cor.tpc) %>%
#     dplyr::mutate(top.high = c(rep(T, k), rep(F, nrow(.) - k))) %>%
#     dplyr::filter(top.low | top.high)
#     
#   
#   features <- Boruta::Boruta(train.data$tumour.percentage.dna~. ,
#                              data = train.data %>% dplyr::select( candidate.features %>% dplyr::pull(ensembl_id) )) %>%
#     purrr::pluck('finalDecision') %>%
#     as.data.frame() %>%
#     dplyr::rename(boruta.status = '.') %>%
#     tibble::rownames_to_column('ensembl_id') %>%
#     dplyr::filter(boruta.status %in% c("Confirmed","Tentative")) %>%
#     dplyr::left_join(candidate.features, by=c('ensembl_id' = 'ensembl_id'))
#   
#   dim(features)
#   dim(train.data)
#   dim(test.data)
#   
#   train.data <- train.data %>%
#     dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna'))
#   
#   test.data <- test.data %>%
#     dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna'))
#   
#   dim(train.data)
#   dim(test.data)
#   
#   
#   
#   test.caret <- caret::train(
#     x = train.data %>% dplyr::mutate(tumour.percentage.dna = NULL),
#     y = train.data$tumour.percentage.dna ,
#     trControl = caret::trainControl(
#       returnData = TRUE,
#       method = "repeatedcv",
#       number = 10,
#       repeats= 8
#     ),
#     method = "ranger",
#     importance = "permutation"
#   )
#   
#   
#   test.data$tumour.percentage.dna.predicted <- predict(test.caret, test.data)
#   
#   test.out.all <- rbind(test.data %>% dplyr::select('tumour.percentage.dna', 'tumour.percentage.dna.predicted'), test.out.all)
#   
# }
# 
# 
# 
# 
# 
# # voeg toe aan metadata
# gsam.metadata.all <- gsam.metadata.all %>%
#   dplyr::mutate(tumour.percentage.dna.predicted = NULL) %>% # force overwrite
#   dplyr::left_join(
#     test.out.all %>%
#       dplyr::select('tumour.percentage.dna.predicted') %>%
#       tibble::rownames_to_column('sid'),
#     by = c('sid'='sid'))
# 
# gsam.metadata.all$tumour.percentage.dna.predicted
# 
# plt <- gsam.metadata.all %>%
#   dplyr::filter(!is.na(tumour.percentage.dna.predicted) & !is.na(tumour.percentage.dna.predicted))
# 
# 
# 
# plt <- read.table('output/tables/GSAM.tumour.percentage.dna.imputed.caret.txt')
# sqrt(sum((plt$tumour.percentage.dna/100 - plt$tumour.percentage.dna.predicted/100)^2))
# ggplot(data=plt, aes(x=tumour.percentage.dna, y=tumour.percentage.dna.predicted)) +
#   geom_point() +
#   xlim(0,100) +
#   ylim(0,100) +
#   youri_gg_theme +
#   geom_smooth(method="lm") + 
#   labs(x = "Tumor cell percentage (WES)", y="Tumor cell percentage (Random forest re-fit on expression data)")
# cor(plt$tumour.percentage.dna , plt$tumour.percentage.dna.predicted)




# Predict GLASS ----

## randomForest ----


train.data <- all.gene.expression.vst %>%
  dplyr::filter(rownames(.) %in% rownames( all.gene.expression.vst)) %>%
  dplyr::select(gsam.metadata.all$sid) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(metadata %>% dplyr::select(c('sid','tumour.percentage.dna')) ,
    by = c('sid'='sid')) %>%
  tibble::column_to_rownames('sid')


test.data <- all.gene.expression.vst %>%
  dplyr::filter(rownames(.) %in% rownames( all.gene.expression.vst)) %>%
  dplyr::select(glass.metadata.all$sid) %>%
  t() %>%
  as.data.frame()


stopifnot(ncol(train.data) == (ncol(test.data) + 1))



set.seed(1+3+3+7)
k = 150
candidate.features <- readRDS(file = 'tmp/results.out.Rds') %>%
  dplyr::select(c(ensembl_id, statistic.gsam.cor.tpc)) %>%
  dplyr::filter(ensembl_id %in% rownames(all.gene.expression.vst)) %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::arrange(statistic.gsam.cor.tpc) %>%
  dplyr::mutate(top.low = c(rep(T, k), rep(F, nrow(.) - k))) %>%
  dplyr::arrange(-statistic.gsam.cor.tpc) %>%
  dplyr::mutate(top.high = c(rep(T, k), rep(F, nrow(.) - k))) %>%
  dplyr::filter(top.high | top.low)


dim(features)
dim(train.data)
dim(test.data)

train.data <- train.data %>%
  dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna'))

test.data <- test.data %>%
  dplyr::select(features$ensembl_id)

dim(train.data)
dim(test.data)


test.rf <- randomForest(tumour.percentage.dna ~ .,
                        data=train.data,
                        ntree = 5000 )

test.data$tumour.percentage.dna.predicted <- predict(test.rf, test.data)

write.table(test.data %>%
              dplyr::select(tumour.percentage.dna.predicted) %>%
              dplyr::rename(tumour.percentage.dna.imputed.rf = tumour.percentage.dna.predicted) %>%
              tibble::rownames_to_column('sid'), file="output/tables/GLASS.tumour.percentage.dna.imputed.rf.txt")




## caret - - - -
# 
# 
# train.data <- all.gene.expression.vst %>%
#   dplyr::select(gsam.metadata.all$sid) %>%
#   dplyr::filter(rownames(.) %in% c(features1, features2)) %>%
#   t() %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column('sid') %>%
#   dplyr::left_join(metadata %>% dplyr::select(c('sid','tumour.percentage.dna')) ,
#                    by = c('sid'='sid')) %>%
#   tibble::column_to_rownames('sid')
# 
# 
# test.data <- all.gene.expression.vst %>%
#   dplyr::select(glass.metadata.all$sid) %>%
#   dplyr::filter(rownames(.) %in% c(features1, features2)) %>%
#   t() %>%
#   as.data.frame()
# 
# 
# set.seed(1+3+3+7)
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
# p <- data.frame(sid = rownames(test.data),
#                 tumour.percentage.dna.imputed.caret = predict(test.caret, test.data) )
# 
# write.table(p, file="output/tables/GLASS.tumour.percentage.dna.imputed.caret.txt")
# 

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




plt <- read.table('output/tables/GLASS.tumour.percentage.dna.imputed.caret.txt') %>% 
  dplyr::left_join(read.table('output/tables/GLASS.tumour.percentage.dna.imputed.rf.txt') , by=c('sid'='sid'))


# ggplot(plt, aes(
#   x = tumour.percentage.dna.imputed.caret,
#   y = tumour.percentage.dna.predicted.rf
# )) + geom_point()


# 2022 RF imputed purities ----


## prep and normalise expression data ----

gene.metadata.all.patients <- dplyr::full_join(
  gsam.rnaseq.expression %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::select(gid) %>%
    dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
    dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
    dplyr::mutate(in.gsam = T) %>%
    dplyr::mutate(chr = as.factor(gsub("^.+(chr[^:]+):.+$","\\1",gid)))
  ,
  glass.gbm.rnaseq.expression.all.samples %>%
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



stopifnot(sum(duplicated(gene.metadata.all.patients$ensembl_id)) == 0)



all.gene.expression.all.patients <- gene.metadata.all.patients |> 
  dplyr::select(c('gid','ensembl_id')) |> 
  dplyr::left_join(gsam.rnaseq.expression |> tibble::rownames_to_column('gid'),by=c('gid' = 'gid')) |> 
  dplyr::left_join(glass.gbm.rnaseq.expression.all.samples %>% tibble::rownames_to_column('ensembl_id'), by=c('ensembl_id' = 'ensembl_id')) |> 
  dplyr::mutate(gid = NULL) |> 
  tibble::column_to_rownames('ensembl_id')


tmp.metadata <- data.frame(sid = colnames(all.gene.expression.all.patients)) |> 
  dplyr::mutate(cond = runif(n(),1,2))  |> 
  dplyr::mutate(cond = round(cond))  |> 
  dplyr::mutate(cond = paste0("r", cond)) |> 
  dplyr::mutate(cond = as.factor(cond)) |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |> 
      dplyr::select(aliquot_barcode, aliquot_batch_synapse, tumour.percentage.dna.cnv.2022),
    by=c('sid'='aliquot_barcode')
  ) |> 
  dplyr::mutate(aliquot_batch_synapse = ifelse(is.na(aliquot_batch_synapse),"G-SAM", aliquot_batch_synapse)) |> 
  dplyr::left_join(
    gsam.rna.metadata |>
      dplyr::select(sid, tumour.percentage.dna)
    ,by=('sid'='sid')
  ) |> 
  dplyr::mutate(tumour.percentage.dna.cnv.2022 = ifelse(is.na(tumour.percentage.dna.cnv.2022), tumour.percentage.dna, tumour.percentage.dna.cnv.2022)) |> 
  dplyr::mutate(tumour.percentage.dna = NULL) |> 
  tibble::column_to_rownames('sid')



stopifnot(rownames(tmp.metadata) == colnames(all.gene.expression.all.patients))
all.gene.expression.vst.all.patients <- all.gene.expression.all.patients |> 
  DESeq2::DESeqDataSetFromMatrix(tmp.metadata, ~cond) |> 
  DESeq2::vst(blind=T) |> 
  SummarizedExperiment::assay() |> 
  limma::removeBatchEffect(as.factor(tmp.metadata$aliquot_batch_synapse)) |> # remove batch effects: aliquot_batch_synapse
  as.data.frame(stringsAsFactors=F)


## A: ~ gsam purities ----



k <- 150
candidate.features <- readRDS(file = 'tmp/results.out.Rds') |> 
  dplyr::select(c(ensembl_id, statistic.gsam.cor.tpc)) |> 
  dplyr::filter(ensembl_id %in% rownames(all.gene.expression.vst.all.patients)) |> 
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) |> 
  dplyr::arrange(statistic.gsam.cor.tpc) |> 
  dplyr::mutate(top.low = c(rep(T, k), rep(F, n() - k))) |> 
  dplyr::arrange(-statistic.gsam.cor.tpc) |> 
  dplyr::mutate(top.high = c(rep(T, k), rep(F, n() - k))) |> 
  dplyr::filter(top.high | top.low)



train.data <- all.gene.expression.vst.all.patients |> 
  dplyr::select(
    tmp.metadata |>
      dplyr::filter(aliquot_batch_synapse == "G-SAM" & !is.na(tumour.percentage.dna.cnv.2022)) |>
      tibble::rownames_to_column('sid') |>
      dplyr::pull('sid')) |> 
  t() |> 
  as.data.frame(stringsAsFactors = F) |> 
  dplyr::select(candidate.features$ensembl_id) |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::left_join(tmp.metadata |> 
    tibble::rownames_to_column('sid') |> 
    dplyr::select(c('sid','tumour.percentage.dna.cnv.2022')) ,
    by = c('sid'='sid')) %>%
  tibble::column_to_rownames('sid') |> 
  dplyr::relocate(tumour.percentage.dna.cnv.2022, .before = tidyselect::everything()) # move to front, for convenience



test.data <- all.gene.expression.vst.all.patients %>%
  dplyr::select(
    tmp.metadata |>
      dplyr::filter(aliquot_batch_synapse != "G-SAM") |>
      tibble::rownames_to_column('sid') |>
      dplyr::pull('sid')) |> 
  t() |> 
  as.data.frame(stringsAsFactors = F) |> 
  dplyr::select(candidate.features$ensembl_id)


stopifnot(ncol(train.data) == (ncol(test.data) + 1))





set.seed(1+3+3+7)
model.rf.all.patients <- randomForest::randomForest(tumour.percentage.dna.cnv.2022 ~ .,
                        data=train.data,
                        ntree = 5000 )

test.data$tumour.percentage.dna.imputed.rf.2022.all.patients.A <- predict(model.rf.all.patients, test.data)

write.table(test.data |> 
    dplyr::select(tumour.percentage.dna.imputed.rf.2022.all.patients.A), file="output/tables/GLASS.tumour.percentage.dna.imputed.rf.A.2022.all.patients.A.txt")



rm(candidate.features, train.data,test.data, model.rf.all.patients)


## B: ~ glass purities only [10xCV*] ----

# use CNV calls from glass only

# 35 zonder CNV purity -> alsnog imputen
tmp.metadata |>
  dplyr::filter(aliquot_batch_synapse != "G-SAM" & is.na(tumour.percentage.dna.cnv.2022)) |> 
  dim()



test.out.all <- data.frame()
for(i in 1:10) {
  print(i)
  
  train.metadata <- tmp.metadata |> 
    dplyr::filter(aliquot_batch_synapse != "G-SAM") |> 
    dplyr::filter(1:n() %% 10 != (i-1) ) |> 
    dplyr::filter(!is.na(tumour.percentage.dna.cnv.2022)) |> 
    tibble::rownames_to_column('aliquot_barcode')
  
  test.metadata <- tmp.metadata |> 
    dplyr::filter( aliquot_batch_synapse != "G-SAM") |> 
    dplyr::filter(1:n() %% 10 == (i-1) )|> 
    tibble::rownames_to_column('aliquot_barcode')
  
  stopifnot(train.metadata$aliquot_barcode %in% test.metadata$aliquot_barcode == F)
  
  
  train.expression.data <- all.gene.expression.vst.all.patients |> 
    dplyr::select(train.metadata$aliquot_barcode)
  test.expression.data <- all.gene.expression.vst.all.patients |> 
    dplyr::select(test.metadata$aliquot_barcode)
  
  #'@ correlate train.expression.data genes with purity
  k <- 1500
  candidate.features <- train.expression.data |> 
    t() |> 
    as.data.frame(stringsAsFactors=F) |> 
    pbapply::pblapply(cor,y=train.metadata$tumour.percentage.dna.cnv.2022) |> 
    unlist() |> 
    data.frame() |> 
    `colnames<-`('cor') |> 
    dplyr::filter(!is.na(cor)) |> 
    dplyr::arrange(cor) |> 
    dplyr::mutate(top.low = c(rep(T, k), rep(F, n() - k))) |> 
    dplyr::arrange(-cor) |> 
    dplyr::mutate(top.high = c(rep(T, k), rep(F, n() - k))) |> 
    dplyr::filter(top.high | top.low) |> 
    tibble::rownames_to_column('ensembl_id')
  
  
  train.data <- train.expression.data |> # [rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() |> 
    as.data.frame(stringsAsFactors=F) |> 
    dplyr::select(candidate.features$ensembl_id) |> 
    tibble::rownames_to_column('aliquot_barcode') |> 
    dplyr::left_join(train.metadata |>
                       dplyr::select(aliquot_barcode, tumour.percentage.dna.cnv.2022)
                     ,by=c('aliquot_barcode'='aliquot_barcode')) |> 
    dplyr::relocate(tumour.percentage.dna.cnv.2022, .before = tidyselect::everything()) |>  # move to front, easier
    tibble::column_to_rownames('aliquot_barcode')
  
  
  test.data <- test.expression.data |>#[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() |> 
    as.data.frame(stringsAsFactors=F) |> 
    dplyr::select(candidate.features$ensembl_id)
  
  
  set.seed(1+3+3+7)
  features <- Boruta::Boruta(tumour.percentage.dna.cnv.2022 ~. ,
                             data = train.data ) %>%
    purrr::pluck('finalDecision') %>%
    as.data.frame(stringsAsFactor = F) %>%
    dplyr::rename(boruta.status = '.') %>%
    tibble::rownames_to_column('ensembl_id') %>%
    dplyr::filter(boruta.status %in% c("Confirmed","Tentative")) 
  
  print(dim(features))
  dim(train.data)
  dim(test.data)
  
  train.data <- train.data |> 
    dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna.cnv.2022'))
  
  test.data <- test.data |> 
    dplyr::select(c(features$ensembl_id))
  
  dim(train.data)
  dim(test.data)
  stopifnot((ncol(test.data) + 1) == ncol(train.data))
  
  model.rf <- randomForest::randomForest(tumour.percentage.dna.cnv.2022 ~ .,
                                         data=train.data,
                                         ntree = 5000,
                                         #mtry=3,
                                         
                                         importance=TRUE
                                         
                                         #na.action=na.omit
  )
  
  test.data$tumour.percentage.dna.imputed.rf.2022.all.patients.B <- predict(model.rf, test.data)
  
  
  test.out.all <- test.out.all |> 
    rbind(test.data |>  dplyr::select('tumour.percentage.dna.imputed.rf.2022.all.patients.B'))
}




#write.table(test.out.all |> 
#              dplyr::select(tumour.percentage.dna.imputed.rf.2022.all.patients.B), file="output/tables/GLASS.tumour.percentage.dna.imputed.rf.A.2022.all.patients.B.txt")


## C: ~ g-sam + glass purities only [10xCV*] ----

# 35 zonder CNV purity -> alsnog imputen
tmp.metadata |>
  dplyr::filter(aliquot_batch_synapse != "G-SAM" & is.na(tumour.percentage.dna.cnv.2022)) |> 
  dim()

train.gsam <- tmp.metadata |> 
  dplyr::filter(aliquot_batch_synapse == "G-SAM") |> 
  dplyr::filter(!is.na(tumour.percentage.dna.cnv.2022)) |> 
  tibble::rownames_to_column('aliquot_barcode')


test.out.all.C <- data.frame()
for(i in 1:10) {
  print(i)
  
  train.metadata <- tmp.metadata |> 
    dplyr::filter(aliquot_batch_synapse != "G-SAM") |> 
    dplyr::filter(1:n() %% 10 != (i-1) ) |> 
    dplyr::filter(!is.na(tumour.percentage.dna.cnv.2022)) |> 
    tibble::rownames_to_column('aliquot_barcode') |> 
    rbind(train.gsam)
  
  test.metadata <- tmp.metadata |> 
    dplyr::filter( aliquot_batch_synapse != "G-SAM") |> 
    dplyr::filter(1:n() %% 10 == (i-1) )|> 
    tibble::rownames_to_column('aliquot_barcode')
  
  stopifnot(train.metadata$aliquot_barcode %in% test.metadata$aliquot_barcode == F)
  
  
  train.expression.data <- all.gene.expression.vst.all.patients |> 
    dplyr::select(train.metadata$aliquot_barcode)
  test.expression.data <- all.gene.expression.vst.all.patients |> 
    dplyr::select(test.metadata$aliquot_barcode)
  
  #'@ correlate train.expression.data genes with purity
  k <- 250
  candidate.features <- train.expression.data |> 
    t() |> 
    as.data.frame(stringsAsFactors=F) |> 
    pbapply::pblapply(cor,y=train.metadata$tumour.percentage.dna.cnv.2022) |> 
    unlist() |> 
    data.frame() |> 
    `colnames<-`('cor') |> 
    dplyr::filter(!is.na(cor)) |> 
    dplyr::arrange(cor) |> 
    dplyr::mutate(top.low = c(rep(T, k), rep(F, n() - k))) |> 
    dplyr::arrange(-cor) |> 
    dplyr::mutate(top.high = c(rep(T, k), rep(F, n() - k))) |> 
    dplyr::filter(top.high | top.low) |> 
    tibble::rownames_to_column('ensembl_id')
  
  
  train.data <- train.expression.data |> # [rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() |> 
    as.data.frame(stringsAsFactors=F) |> 
    dplyr::select(candidate.features$ensembl_id) |> 
    tibble::rownames_to_column('aliquot_barcode') |> 
    dplyr::left_join(train.metadata |>
                       dplyr::select(aliquot_barcode, tumour.percentage.dna.cnv.2022)
                     ,by=c('aliquot_barcode'='aliquot_barcode')) |> 
    dplyr::relocate(tumour.percentage.dna.cnv.2022, .before = tidyselect::everything()) |>  # move to front, easier
    tibble::column_to_rownames('aliquot_barcode')
  
  
  test.data <- test.expression.data |>#[rownames(train.expression.data) %in% c(features1, features2),] %>%
    t() |> 
    as.data.frame(stringsAsFactors=F) |> 
    dplyr::select(candidate.features$ensembl_id)
  
  
  set.seed(1+3+3+7)
  features <- Boruta::Boruta(tumour.percentage.dna.cnv.2022 ~. ,
                             data = train.data ) %>%
    purrr::pluck('finalDecision') %>%
    as.data.frame(stringsAsFactor = F) %>%
    dplyr::rename(boruta.status = '.') %>%
    tibble::rownames_to_column('ensembl_id') %>%
    dplyr::filter(boruta.status %in% c("Confirmed","Tentative")) 
  
  print(dim(features))
  dim(train.data)
  dim(test.data)
  
  train.data <- train.data |> 
    dplyr::select(c(features$ensembl_id, 'tumour.percentage.dna.cnv.2022'))
  
  test.data <- test.data |> 
    dplyr::select(c(features$ensembl_id))
  
  dim(train.data)
  dim(test.data)
  stopifnot((ncol(test.data) + 1) == ncol(train.data))
  
  set.seed(1+3+3+7)
  model.rf <- randomForest::randomForest(tumour.percentage.dna.cnv.2022 ~ .,
                                         data=train.data,
                                         ntree = 5000,
                                         #mtry=3,
                                         
                                         importance=TRUE
                                         
                                         #na.action=na.omit
  )
  
  test.data$tumour.percentage.dna.imputed.rf.2022.all.patients.C <- predict(model.rf, test.data)
  
  
  test.out.all.C <- test.out.all.C |> 
    rbind(test.data |>  dplyr::select('tumour.percentage.dna.imputed.rf.2022.all.patients.C'))
}




stopifnot(tmp.metadata |> dplyr::filter(aliquot_batch_synapse != "G-SAM") |> rownames() %in% rownames(test.out.all.C)) # ensure all have been predicted


write.table(test.out.all.C |> 
           dplyr::select(tumour.percentage.dna.imputed.rf.2022.all.patients.C),
           file="output/tables/GLASS.tumour.percentage.dna.imputed.rf.A.2022.all.patients.C.txt")
 




## plots and stats ----



plt <- test.out.all |> 
  tibble::rownames_to_column('aliquot_barcode') |> 
  dplyr::left_join(
    tmp.metadata |> 
      dplyr::select(tumour.percentage.dna.cnv.2022) |> 
      tibble::rownames_to_column('aliquot_barcode')
    ,by=c('aliquot_barcode'='aliquot_barcode')
  ) |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |> 
      dplyr::select(aliquot_barcode,purity.synapse.rna,is.primary,aliquot_batch_synapse)
    ,by=c('aliquot_barcode'='aliquot_barcode')
  ) |> 
  dplyr::left_join(
    read.table("output/tables/GLASS.tumour.percentage.dna.imputed.rf.A.2022.all.patients.A.txt") |> 
      tibble::rownames_to_column('aliquot_barcode')
    ,by=c('aliquot_barcode'='aliquot_barcode')
  ) |> 
  dplyr::left_join(
    test.out.all.C |> 
      tibble::rownames_to_column('aliquot_barcode')
    ,by=c('aliquot_barcode'='aliquot_barcode')
  ) |> 
  dplyr::mutate(purity.synapse.rna = purity.synapse.rna * 100) |> 
  dplyr::mutate(dist.C = tumour.percentage.dna.imputed.rf.2022.all.patients.C - tumour.percentage.dna.cnv.2022) |> 
  dplyr::mutate(dist.B = tumour.percentage.dna.imputed.rf.2022.all.patients.B - tumour.percentage.dna.cnv.2022) |> 
  dplyr::mutate(dist.A = tumour.percentage.dna.imputed.rf.2022.all.patients.A - tumour.percentage.dna.cnv.2022) |> 
  dplyr::mutate(dist.EST = purity.synapse.rna - tumour.percentage.dna.cnv.2022) 


stopifnot(sum(is.na(plt$tumour.percentage.dna.imputed.rf.2022.all.patients.A)) == 0)
stopifnot(sum(is.na(plt$tumour.percentage.dna.imputed.rf.2022.all.patients.B)) == 0)
stopifnot(sum(is.na(plt$tumour.percentage.dna.imputed.rf.2022.all.patients.C)) == 0)



### Find most appropriate fit ----



lm.A = lm(tumour.percentage.dna.imputed.rf.2022.all.patients.A ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))
lm.B = lm(tumour.percentage.dna.imputed.rf.2022.all.patients.B ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))
lm.C = lm(tumour.percentage.dna.imputed.rf.2022.all.patients.C ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))

test <- plt |> dplyr::filter(!is.na(tumour.percentage.dna.cnv.2022))
sqrt(sum((test$tumour.percentage.dna.cnv.2022 - test$tumour.percentage.dna.imputed.rf.2022.all.patients.A)^2))
sqrt(sum((test$tumour.percentage.dna.cnv.2022 - test$tumour.percentage.dna.imputed.rf.2022.all.patients.B)^2)) # best fit, lowest squared error
sqrt(sum((test$tumour.percentage.dna.cnv.2022 - test$tumour.percentage.dna.imputed.rf.2022.all.patients.C)^2))



### plots ----


plt <- plt |> 
  dplyr::left_join(
    data.frame(residuals.EST = lm(purity.synapse.rna ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))$residuals) |> 
      tibble::rownames_to_column('aliquot_barcode'),
    by=c('aliquot_barcode'='aliquot_barcode'),suffix=c('','')
  )

plt <- plt |> 
  dplyr::left_join(
    data.frame(residuals.A = lm(tumour.percentage.dna.imputed.rf.2022.all.patients.A ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))$residuals) |> 
      tibble::rownames_to_column('aliquot_barcode'),
    by=c('aliquot_barcode'='aliquot_barcode'),suffix=c('','')
  )

plt <- plt |> 
  dplyr::left_join(
    data.frame(residuals.B = lm(tumour.percentage.dna.imputed.rf.2022.all.patients.B ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))$residuals) |> 
      tibble::rownames_to_column('aliquot_barcode'),
    by=c('aliquot_barcode'='aliquot_barcode'),suffix=c('','')
  )

plt <- plt |> 
  dplyr::left_join(
    data.frame(residuals.C = lm(tumour.percentage.dna.imputed.rf.2022.all.patients.C ~ tumour.percentage.dna.cnv.2022, plt |>  tibble::column_to_rownames('aliquot_barcode'))$residuals) |> 
      tibble::rownames_to_column('aliquot_barcode'),
    by=c('aliquot_barcode'='aliquot_barcode'),suffix=c('','')
  )



# mean(lm.A$residuals^2) fewer data points
mean(lm.B$residuals^2)
mean(lm.C$residuals^2)




plt <- plt |> 
  dplyr::mutate(outlier.A = abs(residuals.A) > 20) |> 
  dplyr::mutate(outlier.B = abs(residuals.B) > 20) |> 
  dplyr::mutate(outlier.C = abs(residuals.C) > 20) |> 
  dplyr::mutate(outlier.EST = abs(residuals.EST) > 20)


plot(plt$residuals.B,plt$residuals.A)
plot(plt$residuals.B,plt$residuals.EST)





ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022, purity.synapse.rna, col=outlier.B)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100) + 
  geom_smooth(aes(col = NULL),method="lm",se=F,col="black")


# ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022, y=tumour.percentage.dna.imputed.rf.2022.all.patients.A,col=outlier.EST)) +
#   geom_point() +
#   xlim(0,100) +
#   ylim(0,100)


ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022, 
               y=tumour.percentage.dna.imputed.rf.2022.all.patients.B,
               label=aliquot_barcode,
               col=aliquot_batch_synapse)) +
  geom_point() +
  ggrepel::geom_text_repel(data = plt |> dplyr::filter(outlier.B == T),cex=3) +
  xlim(0,100) +
  ylim(0,100)


# ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022, y=tumour.percentage.dna.imputed.rf.2022.all.patients.C,col=outlier.EST)) +
#   geom_point() +
#   xlim(0,100) +
#   ylim(0,100)


plot(plt$dist.A, plt$dist.B)
plot(plt$dist.A, plt$dist.EST)



ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022, y=purity.synapse.rna)) +
  geom_point()


ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022, y=tumour.percentage.dna.imputed.rf.2022.all.patients)) +
  geom_point()


ggplot(plt,aes(x=tumour.percentage.dna.cnv.2022,
               y=tumour.percentage.dna.imputed.rf.2022.all.patients.B,
               label=aliquot_barcode,
               col=aliquot_batch_synapse)) + # predicted.GLASS.batch,aliquot_batch_synapse
  geom_point() +
  #facet_grid(cols = vars(aliquot_batch_synapse)) +
  ggrepel::geom_text_repel(data = plt |> dplyr::filter(dist.B > 30 ),cex=3)



ggplot(plt, aes(x=aliquot_batch_synapse, y=dist, col=predicted.GLASS.batch)) +
  geom_violin(alpha=0.1) + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  youri_gg_theme



corrplot(cor(as.matrix(plt |>  dplyr::filter(!is.na(tumour.percentage.dna.cnv.2022)) |> tibble::column_to_rownames('aliquot_barcode'))),tl.cex=0.4)




ggplot(plt, aes(x=is.primary, y=tumour.percentage.dna.imputed.rf.2022.all.patients)) +
  geom_violin(alpha=0.1) + 
  ggbeeswarm::geom_beeswarm(pch=21,cex=1.2) +
  youri_gg_theme


wilcox.test(
  plt |> dplyr::filter(is.primary == T) |> dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.A),
  plt |> dplyr::filter(is.primary != T) |> dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.A)
) |> purrr::pluck('p.value')


wilcox.test(
  plt |> dplyr::filter(is.primary == T) |> dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.B),
  plt |> dplyr::filter(is.primary != T) |> dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.B)
) |> purrr::pluck('p.value')


wilcox.test(
  plt |> dplyr::filter(is.primary == T) |> dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.C),
  plt |> dplyr::filter(is.primary != T) |> dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.C)
) |> purrr::pluck('p.value')


wilcox.test(
  plt |> dplyr::filter(is.primary == T) |> dplyr::pull(tumour.percentage.dna.cnv.2022),
  plt |> dplyr::filter(is.primary != T) |> dplyr::pull(tumour.percentage.dna.cnv.2022)
) |> purrr::pluck('p.value')



tmp <- plt |>
  dplyr::select(is.primary, tumour.percentage.dna.cnv.2022,
                             tumour.percentage.dna.imputed.rf.2022.all.patients.A,
                             tumour.percentage.dna.imputed.rf.2022.all.patients.B,
                             tumour.percentage.dna.imputed.rf.2022.all.patients.C) |> 
  dplyr::mutate(purity = ifelse(is.na(tumour.percentage.dna.cnv.2022),tumour.percentage.dna.imputed.rf.2022.all.patients.B ,tumour.percentage.dna.cnv.2022))



wilcox.test(
  tmp |> dplyr::filter(is.primary == T) |> dplyr::pull(purity),
  tmp |> dplyr::filter(is.primary != T) |> dplyr::pull(purity)
) |> purrr::pluck('p.value')


## plot again, from metadata ----


#glass.gbm.rnaseq.metadata.all.samples$tumour.percentage.2022.source

ggplot(glass.gbm.rnaseq.metadata.all.samples, 
       aes(x=tumor.purity.cnv.pct.2022,
           y=tumour.percentage.dna.imputed.rf.2022.all.patients.B,
           label=sid.label)) +
  geom_point() +
  ggrepel::geom_text_repel(data = glass.gbm.rnaseq.metadata.all.samples |>  dplyr::filter(
    abs(
      tumor.purity.cnv.pct.2022 - tumour.percentage.dna.imputed.rf.2022.all.patients.B
    ) > 25
    
  ))


#glass.gbm.rnaseq.metadata.all.samples$tumor.purity.cnv.pct.2022
#glass.gbm.rnaseq.metadata.all.samples$tumour.percentage.dna.imputed.rf.2022.all.patients.B


# 〰 © Dr. Youri Hoogstrate 〰 ----


