#!/usr/bin/env R 

# load libs ----


library(EPIC)



# load data ----


source("scripts/R/gsam_rna-seq_expression.R")
source("scripts/R/gsam_metadata.R")


# counts to tpm ----
# use stricly protein coding genes to avoid biasses by rRNA etc.


# [x] read counts
gsam.rnaseq.expression

# [x] feature / gene lengths
feature.lengths <- gsam.rnaseq.expression %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::select(c('gid')) %>%
  dplyr::mutate(ENSG = gsub('\\|.+$','',gid) ) %>%
  dplyr::left_join(
    read.delim("data/gsam/output/tables/gsam_featureCounts_readcounts_new.txt", stringsAsFactors = F,comment="#") %>%
      dplyr::select(c("Geneid", "Length"))
    , by=c('ENSG'='Geneid'))


stopifnot(rownames(gsam.rnaseq.expression) == feature.lengths$gid)


# [x] fragment lengths
fragment.lengths <- data.frame(sid = colnames(gsam.rnaseq.expression)) %>%
  dplyr::left_join(read.delim("output/tables/picard_insertsizemetrics_combined.txt", header=T, stringsAsFactors = F), by = c('sid'='sample'))

stopifnot(rownames(gsam.rnaseq.expression) == feature.lengths$sid)


# slowkow/counts_to_tpm.R ----
#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}


counts.tpm <- counts_to_tpm(gsam.rnaseq.expression , feature.lengths$Length , fragment.lengths$MEDIAN_INSERT_SIZE) %>%
  as.data.frame %>%
  tibble::rownames_to_column('gid') %>%
  #dplyr::mutate(gid.full = gid) %>%
  dplyr::filter(grepl('ENSG00000284917.1_2', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000248835.2', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000203812.2', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000269620.1', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000285053.1_3', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000284741.1_3', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000283706.2_3', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000285258.1_2', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000268791.1', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000063438.17_5', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000280987.4_6', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000285441.1_2', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000168255.20_5', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000285292.1_3', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000258724.1_6', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000187522.16_6', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000284934.1_2', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000205147.3', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000267110.1_7', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000267706.3_8', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000123584.7', gid) == F ) %>%
  dplyr::filter(grepl('ENSG00000286070.1_1', gid) == F ) %>%
  dplyr::mutate(gid =  gsub("^.+\\|([^\\]+)\\|.+$","\\1", gid ) ) %>%
  tibble::column_to_rownames('gid')
  

# make out ----

out <- gsam.patient.metadata %>%
  dplyr::select(c('studyID')) %>%
  dplyr::rename(pid = studyID) %>%
  dplyr::left_join(
    data.frame(sid.R1 = colnames(counts.tpm)) %>%
      dplyr::mutate(resection = gsub('^...(.).*$','r\\1',sid.R1) ) %>%
      dplyr::filter(resection == 'r1') %>%
      dplyr::mutate(pid = gsub('^(...).+$','\\1',sid.R1), resection = NULL)
    ,  by = c('pid'='pid') ) %>%
  dplyr::left_join(
    data.frame(sid.R2 = colnames(counts.tpm)) %>%
      dplyr::mutate(resection = gsub('^...(.).*$','r\\1',sid.R2) ) %>%
      dplyr::filter(resection == 'r2') %>%
      dplyr::mutate(pid = gsub('^(...).+$','\\1',sid.R2), resection = NULL)
    ,  by = c('pid'='pid') )





rbfox3 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("RBFOX3", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(RBFOX3 = `ENSG00000167281.19_4|RBFOX3|chr17:77085427-77512230(-)`)


cd163 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000177575", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(CD163 = `ENSG00000177575.12_3|CD163|chr12:7623409-7656489(-)`)



out <- out %>%
  dplyr::left_join(rbfox3 %>% dplyr::rename(RBFOX3.R1 = RBFOX3), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(rbfox3 %>% dplyr::rename(RBFOX3.R2 = RBFOX3), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(cd163 %>% dplyr::rename(CD163.R1 = CD163), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(cd163 %>% dplyr::rename(CD163.R2 = CD163), by = c('sid.R2'='sid'))




# epic ----




out.epic <- EPIC(bulk = counts.tpm)
out.epic$mRNAProportions <- as.data.frame(out.epic$mRNAProportions) %>%
  `colnames<-`(paste0('EPIC.',colnames(.))) %>%
  tibble::rownames_to_column('sid')


out <- out %>%
  dplyr::left_join(out.epic$mRNAProportions %>% `colnames<-`(paste0(colnames(.),'.R1')) , by = c('sid.R1'='sid.R1')) %>%
  dplyr::left_join(out.epic$mRNAProportions %>% `colnames<-`(paste0(colnames(.),'.R2')) , by = c('sid.R2'='sid.R2'))


# @todo volgende:
out.epic.paired <- dplyr::inner_join(out.epic.r1, out.epic.r2, by=c('pid'='pid'))  %>%
  #dplyr::mutate(Macrophages.log.odds = log(Macrophages.r2) - log(Macrophages.r1) ) %>%
  #dplyr::mutate(Macrophages.log.odds = -log( log(Macrophages.r2) / log(Macrophages.r1) )) %>% 
  #dplyr::mutate(tpc.log.odds =         log( tumour.percentage.dna.r2 / tumour.percentage.dna.r1)) %>%
  dplyr::mutate(Bcells.log.odds =      -log( log(Bcells.r2) / log(Bcells.r1) )) %>% 
  dplyr::mutate(CAFs.log.odds =        -log( log(CAFs.r2) / log(CAFs.r1) )) %>% 
  dplyr::mutate(CD4_Tcells.log.odds =  -log( log(CD4_Tcells.r2) / log(CD4_Tcells.r1) )) %>% 
  dplyr::mutate(CD8_Tcells.log.odds =  -log( log(CD8_Tcells.r2) / log(CD8_Tcells.r1) )) %>% 
  dplyr::mutate(Endothelial.log.odds = -log( log(Endothelial.r2) / log(Endothelial.r1) )) %>% 
  dplyr::mutate(NKcells.log.odds =     -log( log(NKcells.r2) / log(NKcells.r1) )) %>% 
  dplyr::mutate(otherCells.log.odds =  -log( log(otherCells.r2) / log(otherCells.r1) )) %>%
  dplyr::mutate(Macrophages.log.odds =  log(Macrophages.r2 / Macrophages.r1)) %>% 
  dplyr::mutate(tpc.log.odds =          tumour.percentage.dna.r2 - tumour.percentage.dna.r1) %>%
  dplyr::filter(!is.na(tpc.log.odds))



# quantiseq ----

# https://icbi.i-med.ac.at/software/quantiseq/doc/


data.quantiseq <- immunedeconv::deconvolute(counts.tpm ,
                                            method = 'quantiseq',
                                            tumor = T,
                                            scale_mrna = T)
out.quantiseq <- data.quantiseq %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('cell_type') %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(paste0('QS.',colnames(.))) %>%
  tibble::rownames_to_column('sid')


out <- out %>%
  dplyr::left_join(out.quantiseq %>% `colnames<-`(paste0(colnames(.),'.R1')) , by = c('sid.R1'='sid.R1')) %>%
  dplyr::left_join(out.quantiseq %>% `colnames<-`(paste0(colnames(.),'.R2')) , by = c('sid.R2'='sid.R2'))


# todo

out.paired <- dplyr::inner_join(out.r1, out.r2, by=c('pid'='pid'))  %>%
  #dplyr::mutate(Macrophages.log.odds.2 =  Macrophages.r2 - Macrophages.r1 ) %>% 
  
  dplyr::mutate(Bcells.log.odds =      log(Bcells.r2 / Bcells.r1)) %>% 
  dplyr::mutate(CAFs.log.odds =        log(CAFs.r2 / CAFs.r1)) %>% 
  dplyr::mutate(CD4_Tcells.log.odds =  log(CD4_Tcells.r2 / CD4_Tcells.r1)) %>% 
  dplyr::mutate(CD8_Tcells.log.odds =  log(CD8_Tcells.r2 / CD8_Tcells.r1)) %>% 
  dplyr::mutate(Endothelial.log.odds = log(Endothelial.r2 / Endothelial.r1)) %>% 
  
  dplyr::mutate(NKcells.log.odds =     log(NKcells.r2 / NKcells.r1)) %>% 
  dplyr::mutate(otherCells.log.odds =  log(otherCells.r2 / otherCells.r1)) %>%
  
  dplyr::mutate(RBFOX3.log.odds =      log( RBFOX3.r2 / RBFOX3.r1 )) %>%
  dplyr::mutate(CD163.log.odds =       log( CD163.r2 / CD163.r1 )) %>%
  
  dplyr::mutate(Macrophages.log.odds.1 =  log(Macrophages.r2 / Macrophages.r1 )) %>% 
  #dplyr::mutate(tpc.log.odds.1 =         log( tumour.percentage.dna.r2 / tumour.percentage.dna.r1)) %>%
  dplyr::mutate(tpc.log.odds.2 =          tumour.percentage.dna.r2 - tumour.percentage.dna.r1) %>%
  dplyr::filter(!is.na(tpc.log.odds.2) & pid != 'ECD')


# todo make checks on NA's and duplicate id's


# @ done @ make plots

# & pid != 'ECD'

rm(out.epic, out.epic.r1, out.epic.r2)  


out %>% dplyr::left_join(out.epic.paired , by = c('pid'='pid') )



corrplot::corrplot(
  cor(
    out.epic.paired %>%
      dplyr::filter(!is.na(tpc.log.odds)) %>%
      dplyr::select(contains("log.odds")) )
)


plot(out.paired$tpc.log.odds,  out.paired$Macrophages.log.odds)#, ylim = c(2.8,-2)) 
cor(out.paired$tpc.log.odds,  out.paired$Macrophages.log.odds)
abline(h=0)
abline(v=0)


plot(out.paired$tpc.log.odds,  out.paired$CD163.log.odds)
cor(out.paired$tpc.log.odds,  out.paired$CD163.log.odds)



plot(out.paired$CD163.log.odds,  out.paired$Macrophages.log.odds)





plot(sort(out.paired$Bcells.log.odds))
abline(h=0)

plot(sort(out.paired$CAFs.log.odds))
abline(h=0)

plot(sort(out.paired$CD4_Tcells.log.odds))
abline(h=0)

plot(sort(out.paired$CD8_Tcells.log.odds))
abline(h=0)

plot(sort(out.paired$Endothelial.log.odds))
abline(h=0)

plot(sort(out.paired$Macrophages.log.odds))
abline(h=0)

plot(sort(out.paired$NKcells.log.odds))
abline(h=0)

plot(sort(out.paired$otherCells.log.odds))
abline(h=0)

plot(sort(out.paired$tpc.log.odds))
abline(h=0)




# check correlation any of cell types w/ CD163 vst?



plt <- out$mRNAProportions %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(rbfox3, by=c('sid'='sid') ) %>%
  dplyr::left_join(cd163, by=c('sid'='sid') ) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::select(c('sid','tumour.percentage.dna')) , by= c('sid'='sid'))  %>%
  dplyr::filter(!is.na(tumour.percentage.dna)) %>%
  tibble::column_to_rownames('sid')


corrplot::corrplot(cor(plt))


plot(plt$CD163,  plt$Macrophages)
plot(plt$RBFOX3,  plt$CD8_Tcells)






plot(log(plt$Macrophages) , plt$CD163)





plot(out$cellFractions[,1] * 100,  out.2$`B cell` * 100, xlim=c(0,8) , ylim=c(0,8))


plot(out$cellFractions[,6] * 100,  out.2$`Macrophage M1` * 100, xlim=c(0,18) , ylim=c(0,18))
plot(out$cellFractions[,6] * 100,  out.2$`Macrophage M2` * 100, xlim=c(0,28) , ylim=c(0,28))








plot(out.paired$Macrophages.log.odds.1 , out.paired$Macrophages.log.odds.2)
plot(out.paired$tpc.log.odds.1 , out.paired$tpc.log.odds.2)

plot(out.paired$tpc.log.odds , out.paired$Macrophages.log.odds)
#text(out.paired$tpc.log.odds , out.paired$Macrophages.log.odds , out.paired$pid)
cor(out.paired$tpc.log.odds , out.paired$Macrophages.log.odds)


corrplot::corrplot(cor(out.paired %>% dplyr::select(contains('odds'))))


plot(out.paired$tpc.log.odds.2 ,  out.paired$Macrophages.log.odds.1)#, ylim = c(2.8,-2)) 
abline(h=0)
abline(v=0)



