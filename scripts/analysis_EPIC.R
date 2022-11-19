#!/usr/bin/env R 

# load libs ----


#library(EPIC)


# load data ----


source("scripts/load_G-SAM_metadata.R")
source("scripts/load_G-SAM_expression_data.R")

# counts to tpm ----


# [x] read counts
# gsam.rnaseq.expression

# [x] feature / gene lengths
tmp.feature.lengths <- gsam.rnaseq.expression %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::select(c('gid')) %>%
  dplyr::mutate(ENSG = gsub('\\|.+$','',gid) ) %>%
  dplyr::left_join(
    read.delim("data/gsam/output/tables/gsam_featureCounts_readcounts_new.txt", stringsAsFactors = F,comment="#") %>%
      dplyr::select(c("Geneid", "Length"))
    , by=c('ENSG'='Geneid'))


stopifnot(rownames(gsam.rnaseq.expression) == tmp.feature.lengths$gid)


# [x] fragment lengths
tmp.fragment.lengths <- data.frame(sid = colnames(gsam.rnaseq.expression)) %>%
  dplyr::left_join(
    read.delim("data/gsam/output/tables/picard_insertsizemetrics_combined.txt", header=T, stringsAsFactors = F),
    by = c('sid'='sample'),suffix=c('',''))

stopifnot(rownames(gsam.rnaseq.expression) == tmp.feature.lengths$sid)


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


counts.tpm <- counts_to_tpm(gsam.rnaseq.expression , tmp.feature.lengths$Length , tmp.fragment.lengths$MEDIAN_INSERT_SIZE) %>%
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

rm(tmp.feature.lengths, tmp.fragment.lengths)  



# add epic from file to out(s) ----



if(file.exists('tmp/out.epic.Rds')) {
  
  print("EPIC output already present")
  
} else {
  
  out.epic <- EPIC(bulk = counts.tpm)
  out.epic$mRNAProportions <- as.data.frame(out.epic$mRNAProportions) %>%
    `colnames<-`(paste0('EPIC.',colnames(.))) %>%
    tibble::rownames_to_column('sid')
  
  saveRDS(out.epic, file = "tmp/out.epic.Rds")
  
}


# 
# out.patient <- out.patient %>%
#   dplyr::left_join(out.epic$mRNAProportions %>% `colnames<-`(paste0(colnames(.),'.R1')) , by = c('sid.R1'='sid.R1')) %>%
#   dplyr::left_join(out.epic$mRNAProportions %>% `colnames<-`(paste0(colnames(.),'.R2')) , by = c('sid.R2'='sid.R2'))
# 
# out.patient <- out.patient %>%
#   dplyr::mutate( `EPIC.Bcells.log.odds` = log( ( `EPIC.Bcells.R2` + 1 ) /  ( `EPIC.Bcells.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.CAFs.log.odds` = log( ( `EPIC.CAFs.R2` + 1 ) /  ( `EPIC.CAFs.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.CD4_Tcells.log.odds` = log( ( `EPIC.CD4_Tcells.R2` + 1 ) /  ( `EPIC.CD4_Tcells.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.CD8_Tcells.log.odds` = log( ( `EPIC.CD8_Tcells.R2` + 1 ) /  ( `EPIC.CD8_Tcells.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.Endothelial.log.odds` = log( ( `EPIC.Endothelial.R2` + 1 ) /  ( `EPIC.Endothelial.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.Macrophages.log.odds` = log( ( `EPIC.Macrophages.R2` + 1 ) /  ( `EPIC.Macrophages.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.NKcells.log.odds` = log( ( `EPIC.NKcells.R2` + 1 ) /  ( `EPIC.NKcells.R1` + 1 ))) %>%
#   dplyr::mutate( `EPIC.otherCells.log.odds` = log( ( `EPIC.otherCells.R2` + 1 ) /  ( `EPIC.otherCells.R1` + 1 )))
# 
#   
# 
# 
# out.sample <- out.sample %>%
#   dplyr::left_join( out.epic$mRNAProportions , by=c('sid' = 'sid') )
# 
# 
# 
# 
# rm(out.epic)



# quantiseq - - - -

## results are not as good as EPICs [https://doi.org/10.1093/bioinformatics/btz363]
## discrepancies w/ EPIC, requiring to chooose one
# 
# 
# 
# if (file.exists('tmp/out.quantiseq.Rds') &
#     file.exists('tmp/data.quantiseq.Rds')) {
#   data.quantiseq <- readRDS(file = 'tmp/data.quantiseq.Rds')
#   out.quantiseq <- readRDS(file = 'tmp/out.quantiseq.Rds')
#   print("Reading QUANTISEQ data from cache")
#   
# } else {
#   # https://icbi.i-med.ac.at/software/quantiseq/doc/
#   data.quantiseq <- immunedeconv::deconvolute(
#     counts.tpm ,
#     method = 'quantiseq',
#     tumor = T,
#     scale_mrna = T
#   )
#   out.quantiseq <- data.quantiseq %>%
#     as.data.frame() %>%
#     tibble::column_to_rownames('cell_type') %>%
#     t() %>%
#     as.data.frame() %>%
#     `colnames<-`(paste0('QS.', colnames(.))) %>%
#     #`colnames<-`(make.names(colnames(.))) %>%
#     `colnames<-`(gsub(' ', '.', colnames(.))) %>%
#     tibble::rownames_to_column('sid')
#   
#   saveRDS(data.quantiseq, file = "tmp/data.quantiseq.Rds")
#   saveRDS(out.quantiseq, file = "tmp/out.quantiseq.Rds")
# }
# 
# 
# 
# 
# out.patient <- out.patient %>%
#   dplyr::left_join(out.quantiseq %>% `colnames<-`(paste0(colnames(.), '.R1')) ,
#                    by = c('sid.R1' = 'sid.R1')) %>%
#   dplyr::left_join(out.quantiseq %>% `colnames<-`(paste0(colnames(.), '.R2')) ,
#                    by = c('sid.R2' = 'sid.R2')) %>%
#   `colnames<-`(gsub(' ', '.', colnames(.)))
# 
# out.patient <- out.patient %>%
#   dplyr::mutate(`QS.B.cell.log.odds` = log((`QS.B.cell.R2` + 1) /  (`QS.B.cell.R1` + 1))) %>%
#   
#   #dplyr::mutate( `QS.Macrophage.M1.log.odds` = log( ( `QS.Macrophage.M1.R2` + 1 ) /  ( `QS.Macrophage.M1.R1` + 1 ))) %>%
#   #dplyr::mutate( `QS.Macrophage.M2.log.odds` = log( ( `QS.Macrophage.M2.R2` + 1 ) /  ( `QS.Macrophage.M2.R1` + 1 ))) %>%
#   
#   dplyr::mutate(`QS.Macrophage.M1.log.odds` =  (`QS.Macrophage.M1.R2` - `QS.Macrophage.M1.R1`)) %>%
#   dplyr::mutate(`QS.Macrophage.M2.log.odds` =  (`QS.Macrophage.M2.R2` - `QS.Macrophage.M2.R1`)) %>%
#   
#   dplyr::mutate(`QS.Monocyte.log.odds` = log((`QS.Monocyte.R2` + 1) /  (`QS.Monocyte.R1` + 1))) %>%
#   dplyr::mutate(`QS.Neutrophil.log.odds` = log((`QS.Neutrophil.R2` + 1) /  (`QS.Neutrophil.R1` + 1))) %>%
#   dplyr::mutate(`QS.NK.cell.log.odds` = log((`QS.NK.cell.R2` + 1) /  (`QS.NK.cell.R1` + 1))) %>%
#   dplyr::mutate(`QS.T.cell.CD4+.(non-regulatory).log.odds` = log(
#     (`QS.T.cell.CD4+.(non-regulatory).R2` + 1) /  (`QS.T.cell.CD4+.(non-regulatory).R1` + 1)
#   )) %>%
#   dplyr::mutate(`QS.T.cell.CD8+.log.odds` = log((`QS.T.cell.CD8+.R2` + 1) /  (`QS.T.cell.CD8+.R1` + 1))) %>%
#   dplyr::mutate(`QS.T.cell.regulatory.(Tregs).log.odds` = log(
#     (`QS.T.cell.regulatory.(Tregs).R2` + 1) /  (`QS.T.cell.regulatory.(Tregs).R1` + 1)
#   )) %>%
#   dplyr::mutate(`QS.Myeloid.dendritic.cell.log.odds` = log(
#     (`QS.Myeloid.dendritic.cell.R2` + 1) /  (`QS.Myeloid.dendritic.cell.R1` + 1)
#   )) %>%
#   dplyr::mutate(`QS.uncharacterized.cell.log.odds` = log(
#     (`QS.uncharacterized.cell.R2` + 1) /  (`QS.uncharacterized.cell.R1` + 1)
#   ))
# 
# 
# 
# 
# 
# out.sample <- out.sample %>%
#   dplyr::left_join(out.quantiseq , by = c('sid' = 'sid'))
# 
# 
# 
# 
# rm(out.quantiseq)



rm(counts_to_tpm, counts.tpm)



