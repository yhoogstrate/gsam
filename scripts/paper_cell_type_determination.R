#!/usr/bin/env R 

# load libs ----


library(EPIC)
library(patchwork)



# load data ----


source("scripts/R/gsam_rna-seq_expression.R")
source("scripts/R/gsam_metadata.R")

source("scripts/R/youri_gg_theme.R")
source("scripts/R/job_gg_theme.R")

source("scripts/R/palette.R")



# counts to tpm ----
# use stricly protein coding genes to avoid biasses by rRNA etc.


# [x] read counts
# gsam.rnaseq.expression

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
  dplyr::left_join(read.delim("data/gsam/output/tables/picard_insertsizemetrics_combined.txt", header=T, stringsAsFactors = F), by = c('sid'='sample'))

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
  


# make out.sample and out.patients ----


plt.ids <- gsam.rna.metadata %>%
  dplyr::filter(blacklist.pca == F) %>%
  dplyr::filter(pat.with.IDH == F) %>%
  #dplyr::filter(is.na(tumour.percentage.dna) | tumour.percentage.dna >= 15) %>%  # refactor
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
  dplyr::select(c('pid', 'sid')) %>%
  dplyr::mutate(resection = gsub("^...(.).*$","R\\1",sid)) %>%
  dplyr::mutate(pid = as.factor(as.character(pid))) 




out.sample <- gsam.rna.metadata %>%
  dplyr::select(
    c(
      'sid',
      'pid',
      'resection',
      'tumour.percentage.dna' ,
      'NMF.123456.PCA.SVM.class',
      
      'NMF.123456.PCA.SVM.Classical.p',
      'NMF.123456.PCA.SVM.Proneural.p',
      'NMF.123456.PCA.SVM.Mesenchymal.p',
      
      'NMF:123456.1',
      'NMF:123456.2',
      'NMF:123456.3',
      
      'pat.with.IDH',
      'old.batch',
      'blacklist.pca'
    )
  ) %>%
  dplyr::filter(sid %in% plt.ids$sid)



dim(out.sample)

# only full pairs ~ 125
tmp <- out.sample %>%
  dplyr::select(resection, pid) %>% 
  dplyr::mutate(pid = as.factor(as.character(pid))) %>%
  dplyr::group_by(pid) %>%
  tally() %>%
  dplyr::filter(n >= 2) %>%
  dplyr::pull(pid)

# zelfde 125 als hoofd tabel
out.sample <- out.sample %>%
  dplyr::filter(pid %in% tmp)



out.patient <- gsam.patient.metadata %>%
  dplyr::select(c('studyID')) %>%
  dplyr::rename(pid = studyID) %>%
  dplyr::filter(pid %in% tmp) %>% 
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

length(unique(as.character(out.patient$pid)))




stopifnot( gsam.rnaseq.expression.vst > 1) # log fold changes need log( >=1 / >= 1) 

rbfox3 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("RBFOX3", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`RBFOX3 [neur]` = `ENSG00000167281.19_4|RBFOX3|chr17:77085427-77512230(-)`)


cd163 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000177575", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`CD163 [TAM]` = `ENSG00000177575.12_3|CD163|chr12:7623409-7656489(-)`)


nos2 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000007171", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(NOS2 = `ENSG00000007171.17_5|NOS2|chr17:26083792-26127525(-)`)

vav3 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000134215", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(VAV3 = `ENSG00000134215.16_4|VAV3|chr1:108113783-108507802(-)`)

egfr <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000146648", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`EGFR [CL]` = `ENSG00000146648.18_5|EGFR|chr7:55086710-55279321(+)`)

egfrviii <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("EGFRvIII", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid')

lox <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000113083", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`LOX [MES]` = `ENSG00000113083.14_5|LOX|chr5:121398890-121414055(-)`)

erbb3 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000065361", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`ERBB3 [PN]` = `ENSG00000065361.15_4|ERBB3|chr12:56470583-56497289(+)`)

tmem144<- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000164124", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`TMEM144 [olig]` = `ENSG00000164124.10_3|TMEM144|chr4:159122756-159176563(+)`)

nostrin<- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000163072", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`NOSTRIN [endo]` = `ENSG00000163072.15_5|NOSTRIN|chr2:169643049-169722024(+)`)

cachd1<- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(grepl("ENSG00000158966", gid)) %>% 
  column_to_rownames('gid') %>%
  t() %>%
  as.data.frame %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(`CACHD1 [astr]` = `ENSG00000158966.15_4|CACHD1|chr1:64935812-65158741(+)`)




out.patient <- out.patient %>%
  dplyr::left_join(rbfox3 %>% dplyr::rename(RBFOX3.R1 = `RBFOX3 [neur]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(rbfox3 %>% dplyr::rename(RBFOX3.R2 = `RBFOX3 [neur]`), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(cd163 %>% dplyr::rename(CD163.R1 = `CD163 [TAM]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(cd163 %>% dplyr::rename(CD163.R2 = `CD163 [TAM]`), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(nos2 %>% dplyr::rename(NOS2.R1 = NOS2), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(nos2 %>% dplyr::rename(NOS2.R2 = NOS2), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(vav3 %>% dplyr::rename(VAV3.R1 = VAV3), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(vav3 %>% dplyr::rename(VAV3.R2 = VAV3), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(egfr %>% dplyr::rename(EGFR.R1 = `EGFR [CL]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(egfr %>% dplyr::rename(EGFR.R2 = `EGFR [CL]`), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(egfrviii %>% dplyr::rename(EGFRvIII.R1 = EGFRvIII), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(egfrviii %>% dplyr::rename(EGFRvIII.R2 = EGFRvIII), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(lox %>% dplyr::rename(LOX.R1 = `LOX [MES]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(lox %>% dplyr::rename(LOX.R2 = `LOX [MES]`), by = c('sid.R2'='sid')) %>%
  dplyr::left_join(erbb3 %>% dplyr::rename(ERBB3.R1 = `ERBB3 [PN]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(erbb3 %>% dplyr::rename(ERBB3.R2 = `ERBB3 [PN]`), by = c('sid.R2'='sid')) %>%
  
  dplyr::left_join(tmem144 %>% dplyr::rename(TMEM144.R1 = `TMEM144 [olig]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(tmem144 %>% dplyr::rename(TMEM144.R2 = `TMEM144 [olig]`), by = c('sid.R2'='sid')) %>%
  
  dplyr::left_join(nostrin %>% dplyr::rename(NOSTRIN.R1 = `NOSTRIN [endo]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(nostrin %>% dplyr::rename(NOSTRIN.R2 = `NOSTRIN [endo]`), by = c('sid.R2'='sid')) %>%
  
  dplyr::left_join(cachd1 %>% dplyr::rename(CACHD1.R1 = `CACHD1 [astr]`), by = c('sid.R1'='sid')) %>%
  dplyr::left_join(cachd1 %>% dplyr::rename(CACHD1.R2 = `CACHD1 [astr]`), by = c('sid.R2'='sid')) %>%
  
  
  dplyr::left_join(
    gsam.rna.metadata %>%
      dplyr::select(c('sid','tumour.percentage.dna')) %>%
      dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna)
    ,by = c('sid.R1'='sid')) %>% 
  dplyr::left_join(
    gsam.rna.metadata %>%
      dplyr::select(c('sid','tumour.percentage.dna')) %>%
      dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna)
    ,by = c('sid.R2'='sid')) %>%
  dplyr::filter(!is.na(tumour.percentage.dna.R1) & !is.na(tumour.percentage.dna.R2) ) %>%
  dplyr::mutate(tumour.percentage.dna.log.odds = log(tumour.percentage.dna.R2 / tumour.percentage.dna.R1 )) %>%
  
  dplyr::mutate(`RBFOX3.lfc [neur]` = log(RBFOX3.R2 / RBFOX3.R1 )) %>%
  dplyr::mutate(`CD163.lfc [TAM]` =  log(CD163.R2 / CD163.R1 )) %>%
  dplyr::mutate(`NOS2.lfc` =  log(NOS2.R2 / NOS2.R1 )) %>%
  dplyr::mutate(`VAV3.lfc` =  log(VAV3.R2 / VAV3.R1 )) %>%
  dplyr::mutate(`EGFR.lfc [CL]` =  log(EGFR.R2 / EGFR.R1 )) %>%
  dplyr::mutate(`EGFRvIII.lfc` =  log(EGFRvIII.R2 / EGFRvIII.R1 )) %>%
  dplyr::mutate(`LOX.lfc [MES]` =  log(LOX.R2 / LOX.R1 )) %>%
  dplyr::mutate(`ERBB3.lfc [PN]` =  log(ERBB3.R2 / ERBB3.R1 ))  %>%
  dplyr::mutate(`TMEM144.lfc [olig]` =  log(TMEM144.R2 / TMEM144.R1 ))  %>%
  dplyr::mutate(`NOSTRIN.lfc [endo]` =  log(NOSTRIN.R2 / NOSTRIN.R1 ))  %>%
  dplyr::mutate(`CACHD1.lfc [astr]` =  log(CACHD1.R2 / CACHD1.R1 )) 

length(unique(as.character(out.patient$pid)))





stopifnot(is.na(out.patient$tumour.percentage.dna.log.odds) == F)


out.sample <- out.sample %>%
  dplyr::left_join(rbfox3 , by = c('sid'='sid')) %>%
  dplyr::left_join(cd163 , by = c('sid'='sid')) %>%
  dplyr::left_join(nos2, by = c('sid'='sid')) %>%
  dplyr::left_join(vav3 , by = c('sid'='sid')) %>% 
  dplyr::left_join(egfr, by = c('sid'='sid')) %>%
  dplyr::left_join(egfrviii , by = c('sid'='sid')) %>%
  dplyr::left_join(lox , by = c('sid'='sid')) %>%
  dplyr::left_join(erbb3, by = c('sid'='sid')) %>%
  dplyr::left_join(tmem144, by = c('sid'='sid')) %>%
  dplyr::left_join(nostrin, by = c('sid'='sid')) %>%
  dplyr::left_join(cachd1, by = c('sid'='sid')) %>%
  dplyr::filter(!is.na(NOS2)) 



  
# add epic from file to out(s) ----



if(file.exists('tmp/out.epic.Rds')) {
  
  out.epic <- readRDS(file ='tmp/out.epic.Rds')
  print("Reading EPIC data from cache")
  
} else {
  
  out.epic <- EPIC(bulk = counts.tpm)
  out.epic$mRNAProportions <- as.data.frame(out.epic$mRNAProportions) %>%
    `colnames<-`(paste0('EPIC.',colnames(.))) %>%
    tibble::rownames_to_column('sid')
  
  #saveRDS(out.epic, file = "tmp/out.epic.Rds")
  
}



out.patient <- out.patient %>%
  dplyr::left_join(out.epic$mRNAProportions %>% `colnames<-`(paste0(colnames(.),'.R1')) , by = c('sid.R1'='sid.R1')) %>%
  dplyr::left_join(out.epic$mRNAProportions %>% `colnames<-`(paste0(colnames(.),'.R2')) , by = c('sid.R2'='sid.R2'))

out.patient <- out.patient %>%
  dplyr::mutate( `EPIC.Bcells.log.odds` = log( ( `EPIC.Bcells.R2` + 1 ) /  ( `EPIC.Bcells.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.CAFs.log.odds` = log( ( `EPIC.CAFs.R2` + 1 ) /  ( `EPIC.CAFs.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.CD4_Tcells.log.odds` = log( ( `EPIC.CD4_Tcells.R2` + 1 ) /  ( `EPIC.CD4_Tcells.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.CD8_Tcells.log.odds` = log( ( `EPIC.CD8_Tcells.R2` + 1 ) /  ( `EPIC.CD8_Tcells.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.Endothelial.log.odds` = log( ( `EPIC.Endothelial.R2` + 1 ) /  ( `EPIC.Endothelial.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.Macrophages.log.odds` = log( ( `EPIC.Macrophages.R2` + 1 ) /  ( `EPIC.Macrophages.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.NKcells.log.odds` = log( ( `EPIC.NKcells.R2` + 1 ) /  ( `EPIC.NKcells.R1` + 1 ))) %>%
  dplyr::mutate( `EPIC.otherCells.log.odds` = log( ( `EPIC.otherCells.R2` + 1 ) /  ( `EPIC.otherCells.R1` + 1 )))

  


out.sample <- out.sample %>%
  dplyr::left_join( out.epic$mRNAProportions , by=c('sid' = 'sid') )




rm(out.epic)



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




# Additional changes  ----






out.patient <- out.patient %>%  
  dplyr::filter(pid != 'ECD') %>% # excessive outlier, R1 or R2 RNA odd?
  dplyr::left_join(
    gsam.rna.metadata %>%
      dplyr::filter(resection == 'r1') %>%
      dplyr::select(c('pid','NMF.123456.PCA.SVM.class', 'NMF.123456.PCA.SVM.Classical.p','NMF.123456.PCA.SVM.Proneural.p','NMF.123456.PCA.SVM.Mesenchymal.p','NMF:123456.1')) %>%
      dplyr::filter(!is.na(NMF.123456.PCA.SVM.class) ) %>%
      dplyr::rename(NMF.123456.PCA.SVM.class.R1 = NMF.123456.PCA.SVM.class) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Classical.p.R1 = NMF.123456.PCA.SVM.Classical.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Proneural.p.R1 = NMF.123456.PCA.SVM.Proneural.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Mesenchymal.p.R1 = NMF.123456.PCA.SVM.Mesenchymal.p) %>% 
      dplyr::rename('NMF:123456.1.R1' = 'NMF:123456.1')
      
    , by=c('pid'='pid')) %>%
  dplyr::left_join(
    gsam.rna.metadata %>%
      dplyr::filter(resection == 'r2') %>%
      dplyr::select(c('pid','NMF.123456.PCA.SVM.class', 'NMF.123456.PCA.SVM.Classical.p','NMF.123456.PCA.SVM.Proneural.p','NMF.123456.PCA.SVM.Mesenchymal.p','NMF:123456.1')) %>%
      dplyr::filter(!is.na(NMF.123456.PCA.SVM.class) ) %>%
      dplyr::rename(NMF.123456.PCA.SVM.class.R2 = NMF.123456.PCA.SVM.class) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Classical.p.R2 = NMF.123456.PCA.SVM.Classical.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Proneural.p.R2 = NMF.123456.PCA.SVM.Proneural.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Mesenchymal.p.R2 = NMF.123456.PCA.SVM.Mesenchymal.p) %>% 
      dplyr::rename('NMF:123456.1.R2' = 'NMF:123456.1')
    , by=c('pid'='pid')) %>%
  dplyr::mutate(NMF.123456.PCA.SVM.Classical.p.log.odds = log( (NMF.123456.PCA.SVM.Classical.p.R2+1) / (NMF.123456.PCA.SVM.Classical.p.R1+1) )) %>%
  dplyr::mutate(NMF.123456.PCA.SVM.Mesenchymal.p.log.odds = log( (NMF.123456.PCA.SVM.Mesenchymal.p.R2+1) / (NMF.123456.PCA.SVM.Mesenchymal.p.R1+1) )) %>%
  dplyr::mutate(NMF.123456.PCA.SVM.Proneural.p.log.odds = log( (NMF.123456.PCA.SVM.Proneural.p.R2+1) / (NMF.123456.PCA.SVM.Proneural.p.R1+1 ))) %>% 
  dplyr::mutate(`NMF:123456.1.log.odds` = log( (`NMF:123456.1.R1` + 1) / (`NMF:123456.1.R2` +1)))
  
  
  #dplyr::mutate(NMF.123456.PCA.SVM.Classical.p.log.odds2 = NMF.123456.PCA.SVM.Classical.p.R2 - NMF.123456.PCA.SVM.Classical.p.R1 ) %>%
  #dplyr::mutate(NMF.123456.PCA.SVM.Mesenchymal.p.log.odds2 = NMF.123456.PCA.SVM.Mesenchymal.p.R2 - NMF.123456.PCA.SVM.Mesenchymal.p.R1 ) %>%
  #dplyr::mutate(NMF.123456.PCA.SVM.Proneural.p.log.odds2 = NMF.123456.PCA.SVM.Proneural.p.R2 - NMF.123456.PCA.SVM.Proneural.p.R1 )







# plots ----

## 1: corplot per sample ----



plt <- out.sample %>%
  dplyr::rename(`--- [   NMF.SVM.Classical.p   ] ---` = NMF.123456.PCA.SVM.Classical.p) %>%
  dplyr::rename(`--- [   NMF.SVM.Mesenchymal.p   ] ---` = NMF.123456.PCA.SVM.Mesenchymal.p) %>%
  dplyr::rename(`--- [   NMF.SVM.Proneural.p   ] ---` = NMF.123456.PCA.SVM.Proneural.p) %>%
  #dplyr::filter(resection == 'r1') %>%
  dplyr::select(-c('pid','resection','NMF.123456.PCA.SVM.class')) %>%
  tibble::column_to_rownames('sid') %>%
  cor(   )
  #as.data.frame %>%
  #dplyr::select(contains("EPIC")) %>%
  #tibble::rownames_to_column('type') %>%
  #dplyr::filter(grepl("^QS", type)) %>% 
  #tibble::column_to_rownames('type') %>%
  #as.matrix


corrplot::corrplot(plt, order='hclust' )






## 2: corrplot per patient (change/time) ----



plt <- out.patient %>%
  dplyr::rename(SVM.CL.p.log.odds = NMF.123456.PCA.SVM.Classical.p.log.odds) %>%
  dplyr::rename(SVM.MES.p.log.odds = NMF.123456.PCA.SVM.Mesenchymal.p.log.odds) %>%
  dplyr::rename(SVM.PN.p.log.odds = NMF.123456.PCA.SVM.Proneural.p.log.odds) %>%
  #dplyr::filter(NMF.123456.PCA.SVM.class.R1 == "Mesenchymal") %>%
  tibble::column_to_rownames('pid') %>%
  dplyr::select(contains('log.odds')) %>%
  `colnames<-`(gsub('.log.odds','.lfc',colnames(.))) 


corrplot::corrplot(  cor( plt , method="pearson" ) , order="hclust")
#corrplot::corrplot(  cor( plt ) %>% `colnames<-`(NULL) )







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





## 3: per patient stats ----




out <- out.patient

plt <- out %>%
  dplyr::mutate(order = match(1:nrow(.), order(`tumour.percentage.dna.log.odds`))) 



plt <- out %>%
  dplyr::mutate(order = match(1:nrow(.), order(`QS.Macrophage.M2.log.odds`)))



plt <- out %>%
  dplyr::mutate(order = match(1:nrow(.), order(`RBFOX3.log.odds`)))




plt <- out %>%
  dplyr::mutate(order = match(1:nrow(.), order(2 * `VAV3.log.odds` + `NOS2.log.odds`)))





p1 <- ggplot(plt, aes(x =  reorder(pid, order), y = `tumour.percentage.dna.log.odds` ))  +
  geom_bar(stat="identity") + labs(x=NULL)

p2 <- ggplot(plt, aes(x =  reorder(pid, order), y = `CD163.log.odds` ))  +
  geom_bar(stat="identity") + labs(x=NULL)

p3 <- ggplot(plt, aes(x =  reorder(pid, order), y = `EPIC.Macrophages.log.odds` ))  +
  geom_bar(stat="identity") + labs(x=NULL)

#p4 <- ggplot(plt, aes(x =  reorder(pid, order), y = `QS.Macrophage.M1.log.odds` ))  +
#  geom_bar(stat="identity")

p5 <- ggplot(plt, aes(x =  reorder(pid, order), y = `QS.Macrophage.M2.log.odds` ))  +
  geom_bar(stat="identity") + labs(x=NULL)

p6a <- ggplot(plt, aes(x =  reorder(pid, order), y = `NOS2.log.odds` ))  +
  geom_bar(stat="identity") + labs(x=NULL)
p6b <- ggplot(plt, aes(x =  reorder(pid, order), y = `VAV3.log.odds` ))  +
  geom_bar(stat="identity") + labs(x=NULL)


plt.class.R1 <- 
  rbind(
    plt %>%
      dplyr::select(c(pid, order, p.SVM.Classical.R1)) %>%
      dplyr::mutate(class.R1 = "Classical") %>%
      dplyr::rename(probability.R1 = p.SVM.Classical.R1),
    plt %>%
      dplyr::select(c(pid, order, p.SVM.Proneural.R1)) %>%
      dplyr::mutate(class.R1 = "Proneural") %>%
      dplyr::rename(probability.R1 = p.SVM.Proneural.R1) ,
    plt %>%
      dplyr::select(c(pid, order, p.SVM.Mesenchymal.R1)) %>%
      dplyr::mutate(class.R1 = "Mesenchymal") %>%
      dplyr::rename(probability.R1 = p.SVM.Mesenchymal.R1)
  )

p7 <- ggplot(plt.class.R1, aes(x =  reorder(pid, order), y = `probability.R1` , fill=class.R1 ))  +
  geom_bar(stat="identity") + labs(x=NULL)
p7 <- ggplot(plt, aes(x =  reorder(pid, order), y = 1 , fill=NMF.123456.PCA.SVM.class.R1 ))  +
  geom_bar(stat="identity") + labs(x=NULL)




plt.class.R2 <- 
  rbind(
    plt %>%
      dplyr::select(c(pid, order, p.SVM.Classical.R2)) %>%
      dplyr::mutate(class.R2 = "Classical") %>%
      dplyr::rename(probability.R2 = p.SVM.Classical.R2),
    plt %>%
      dplyr::select(c(pid, order, p.SVM.Proneural.R2)) %>%
      dplyr::mutate(class.R2 = "Proneural") %>%
      dplyr::rename(probability.R2 = p.SVM.Proneural.R2) ,
    plt %>%
      dplyr::select(c(pid, order, p.SVM.Mesenchymal.R2)) %>%
      dplyr::mutate(class.R2 = "Mesenchymal") %>%
      dplyr::rename(probability.R2 = p.SVM.Mesenchymal.R2)
  )

p8 <- ggplot(plt.class.R2, aes(x =  reorder(pid, order), y = `probability.R2` , fill=class.R2 ))  +
  geom_bar(stat="identity") + labs(x=NULL)

p8 <- ggplot(plt, aes(x =  reorder(pid, order), y = 1 , fill=NMF.123456.PCA.SVM.class.R2 ))  +
  geom_bar(stat="identity") + labs(x=NULL)




p1 / p2 / p3 / p5 / p7 / p8 / p6a / p6b





plot( out$p.SVM.Classical.R1 , out$tumour.percentage.dna.R1 )
plot( out$p.SVM.Classical.R2 , out$tumour.percentage.dna.R2 )

plot( out$p.SVM.Mesenchymal.R1 , out$tumour.percentage.dna.R1 )
plot( out$p.SVM.Mesenchymal.R2 , out$tumour.percentage.dna.R2 )

plot( out$p.SVM.Proneural.R1 , out$tumour.percentage.dna.R1 )
plot( out$p.SVM.Proneural.R2 , out$tumour.percentage.dna.R2 )




plot( out$tumour.percentage.dna.R1  , out$p.SVM.Classical.R1 )


corrplot::corrplot(
  out %>%
    #tibble::rownames_to_column('type') %>%
    dplyr::mutate(`QS.T.cell.CD8+.R1` = NULL) %>%
    #dplyr::filter(NMF.123456.PCA.SVM.class.R1 == "Classical") %>%
    dplyr::select(- c( NMF.123456.PCA.SVM.class.R1, NMF.123456.PCA.SVM.class.R2 ,sid.R1  , sid.R2)) %>%
    dplyr::select(contains("R1")) %>%
    cor ,
  order='hclust'
)


corrplot::corrplot(
  out %>%
    dplyr::select(- c( NMF.123456.PCA.SVM.class.R1, NMF.123456.PCA.SVM.class.R2 ,sid.R1  , sid.R2)) %>%
    dplyr::select(contains("R2")) %>%
    cor ,
  order="hclust"
)



EPIC.Bcells.log.odds # poor corr
EPIC.CAFs.log.odds # poor corr
EPIC.CD4_Tcells.log.odds # nice corr !!
EPIC.CD8_Tcells.log.odds # poor corr
EPIC.Endothelial.log.odds # some corr !
EPIC.Macrophages.log.odds # nice corr !!
EPIC.NKcells.log.odds # poor corr
EPIC.otherCells.log.odds # some corr ! [tumor cells itself ofc.]

ggplot(out.patient, aes(x = tumour.percentage.dna.log.odds, 
                        y = EPIC.Macrophages.log.odds
                        #, col=NMF.123456.PCA.SVM.class.R1
                        )) +
  geom_point() +
  geom_smooth(method="lm")





plt <- out.sample %>%
  dplyr::filter(pid %in%  out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds,
                                                   EPIC.Macrophages.log.odds,
                                                   EPIC.CD4_Tcells.log.odds
                                                   )), by=c('pid'='pid')) %>%
  dplyr::mutate(tumour.percentage.status = ifelse(tumour.percentage.dna.log.odds > 0, "increase", "decrease"))


plt <- plt %>%
  dplyr::select(pid, resection, tumour.percentage.dna, tumour.percentage.status, EPIC.Macrophages, EPIC.CD4_Tcells) %>% 
  reshape2::melt(id = c('pid','resection', 'tumour.percentage.dna', 'tumour.percentage.status'))


ggplot(plt, aes(x=tumour.percentage.dna, y=value, group=pid, col= tumour.percentage.status)) + 
  facet_grid(cols = vars(variable), scales = "free",space="free") +
  geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
  geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
  labs(x = "Tumor cell percentage (DNA-seq estimate)",y='EPIC deconvolution score') +
  scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1" )) +
  xlim(0, 100) +
  ylim(0, 0.5) +
  youri_gg_theme


# p1 <- ggplot(plt, aes(x=tumour.percentage.dna, y=EPIC.Macrophages, group=pid, col= tumour.percentage.status)) + 
#   geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
#   geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
#   labs(x = "Tumor cell percentage (DNA-seq estimate)",y="EPIC Macrophages score") +
#   scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1" )) +
#   xlim(0, 100) +
#   ylim(0, 0.5) +
#   youri_gg_theme
# 
# p2 <- ggplot(plt, aes(x=tumour.percentage.dna, y=EPIC.CD4_Tcells, group=pid, col= tumour.percentage.status)) + 
#   geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
#   geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
#   labs(x = "Tumor cell percentage (DNA-seq estimate)",y="EPIC CD4 T-cells score") +
#   scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1" )) +
#   xlim(0, 100) +
#   ylim(0, 0.5) +
#   youri_gg_theme
# 
# 
# p1 + p2
# 


ggsave("output/figures/epic_tumor_percentage_traversal.pdf", width = 12, height=4.8)
ggsave("output/figures/epic_tumor_percentage_traversal_4.2.pdf", width = 12, height=4.2)





plt <- out.sample %>%
  dplyr::filter(pid %in% out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds,
                                                   tumour.percentage.dna.R1,
                                                   tumour.percentage.dna.R2,
                                                   EPIC.Macrophages.log.odds,EPIC.CD4_Tcells.log.odds
                                                   ,`NMF:123456.1.log.odds`
                                                   #,`NMF:123456.2.log.odds`
                                                   #,`NMF:123456.3.log.odds`
                                                   )), by=c('pid'='pid')) %>%
  dplyr::mutate(tumour.percentage.dna.delta = tumour.percentage.dna.R2 - tumour.percentage.dna.R1 ) %>%
  dplyr::mutate(panel = ifelse(
    pid %in% (out.sample %>% dplyr::filter(!is.na(tumour.percentage.dna) & tumour.percentage.dna < 15) %>% dplyr::pull(pid)),
    "< 15%", ">= 15%"
  )) %>% 
  dplyr::mutate(panel = factor(panel, levels=c("< 15%", ">= 15%"))) # fix order



plt <- rbind(plt %>% dplyr::mutate(y = tumour.percentage.dna, type="Tumor purity",
                                   sign = ifelse(tumour.percentage.dna.log.odds > 0, "increase", "decrease")),
             plt %>% dplyr::mutate(y = EPIC.Macrophages, type="EPIC Macrophages score",
                                   sign = ifelse(EPIC.Macrophages.log.odds > 0, "increase", "decrease")),
             plt %>% dplyr::mutate(y = EPIC.CD4_Tcells, type="EPIC CD4 T-cells score",
                                   sign = ifelse(EPIC.CD4_Tcells.log.odds > 0, "increase", "decrease"))
             ,plt %>% dplyr::mutate(y = `NMF:123456.1`, type="NMF W1 (MES)",
                                    sign = ifelse(`NMF:123456.1.log.odds` > 0, "increase", "decrease")) 
             # ,plt %>% dplyr::mutate(y = `NMF:123456.2`, type="NMF W2 (MES)",
             #                        sign = ifelse(`NMF:123456.2.log.odds` > 0, "increase", "decrease"))
             # ,plt %>% dplyr::mutate(y = `NMF:123456.3`, type="NMF W3 (MES)",
             #                        sign = ifelse(`NMF:123456.3.log.odds` > 0, "increase", "decrease"))
             )


p1 <- ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=y, col=sign)) +
#p1 <- ggplot(plt, aes(x = reorder(pid, `NMF:123456.1.log.odds`), y=y, col=sign)) +
  geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2, alpha=0.8) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8, lwd=1.05 )  +
  scale_color_manual(values = c('increase'='#bb5f6c', 'decrease'='#79b1b1')) +
  facet_grid(rows=vars(type),cols=vars(panel), scales="free", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  labs(y=NULL,x="G-SAM patient",col="Longitudinal direction")

p1




plt <- out.sample %>%
  dplyr::filter(pid %in% out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds)), by=c('pid'='pid')) %>%
  dplyr::mutate(panel = ifelse(
    pid %in% (out.sample %>% dplyr::filter(!is.na(tumour.percentage.dna) & tumour.percentage.dna < 15) %>% dplyr::pull(pid)),
    "< 15%", ">= 15%"
  )) %>% 
  dplyr::select(pid, panel, tumour.percentage.dna.log.odds) %>%
  dplyr::distinct() %>%
  dplyr::left_join(out.patient %>% dplyr::select(pid, NMF.123456.PCA.SVM.class.R1, NMF.123456.PCA.SVM.class.R2), by = c('pid' = 'pid')) %>%
  dplyr::rename(`Sub-type R1` = NMF.123456.PCA.SVM.class.R1) %>%
  dplyr::rename(`Sub-type R2` = NMF.123456.PCA.SVM.class.R2) %>%
  dplyr::mutate(panel = factor(panel, levels=c("< 15%", ">= 15%"))) %>%
  dplyr::left_join(gsam.patient.metadata %>%
      dplyr::select(studyID, 
                    initialMGMT,HM,
                    AXIN2,APC,JAK2,
                    RB1,MSH2,BRCA1,
                    BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,
                    NF1,ERBB3,
                    EGFR,TP53BP1,TP53,PIK3R1,PIK3CA,TSC2,
                    # PI3K, TP53Signaling, Wnt, Telomere, RTK, RAS, DNADamage, primary7, primary10, recurrent7, recurrent10, cnStatusEGFRs,
                    
                    cnStatusRB1s,
                    cnStatusNF1s,
                    cnStatusCDK4s,
                    cnStatusMDM2s,
                    
                    SETD2,PDGFRA
                    ),by = c('pid' = 'studyID')
  ) %>%
  reshape2::melt(id = c('pid','panel', 'tumour.percentage.dna.log.odds')) %>% 
  dplyr::mutate(value = gsub('Methylated','Yes',value)) %>% 
  dplyr::mutate(value = gsub('Unmethylated','No',value)) %>%
  dplyr::mutate(value = gsub('Normal','No',value)) %>%
  dplyr::mutate(value = gsub('Loss','Yes',value)) %>%
  dplyr::mutate(value = gsub('^Gain$','Yes',value)) %>%
  dplyr::mutate(value = gsub('^Gained$','Gained/increased',value)) %>%
  dplyr::mutate(value = gsub('^Lost$','Lost/decreased',value)) %>%
  dplyr::mutate(ypanel = case_when(
    grepl("NMF|Sub-type",variable) ~ "A",
    variable %in% c("EGFR",  "cnStatusEGFRs", "NF1") ~ "B",
    T ~ "C"
    )) %>% 
  dplyr::mutate(ypanel = factor(ypanel, levels=c('A','B','C')))



p2 <- ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=variable, fill = value)) +
  facet_grid(cols=vars(panel), rows=vars(ypanel), scales="free", space="free") + 
  geom_tile(colour = "black", size = 0.3) +
  theme_bw() +
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  scale_fill_manual(values = c( subtype_colors ,
"Wildtype"="white",
"No"="white",

"Gained/increased"="#bb5f6c", # rood #bb5f6c
"Lost/decreased"="#79b1b1", # lichtlauw #79b1b1
"Stable"="#2e415e", # donker blauw #2e415e
"Yes" = "#2e415e", # zelfde als stable #2e415e
  
"NA"="grey"  )) + 
  labs(y=NULL,x=NULL)




p1 / p2 +  plot_layout(heights = c(2, 1))





ggsave("output/figures/epic_tumor_percentage_traversal_vertical.pdf", width = 16 * 1.05, height=11.5 * 1.05)




## VAV3 en NOS2 staan in de GLITS/REDUX set


