#!/usr/bin/env R


#'@description This scripts classifies Bulk RNA samples according to the
#' Wang/Verhaak ssgsea.GBM.classification. It is a rather old script making
#' use of plyr, which should detach dplyr for safety


# load expression data ----


source('scripts/load_G-SAM_expression_data.R')
source('scripts/load_GLASS_data.R')
source('scripts/load_results.out.R')

detach("dplyr")
detach("tidyverse")


# Select and merge samples ----


# make data.frame with selected patient samples and batches
tmp.combined.metadata <- rbind(
  gsam.rna.metadata %>%
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::select(sid) |>
    dplyr::mutate(dataset = "G-SAM"),
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::select(sid) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |>
      dplyr::select(aliquot_barcode, aliquot_batch_synapse),
    by=c('sid'='aliquot_barcode')
  ) |> 
  dplyr::mutate(batch = ifelse(dataset == "G-SAM","G-SAM",aliquot_batch_synapse)) |> 
  dplyr::mutate(aliquot_batch_synapse = NULL ) |> 
  dplyr::mutate(batch = gsub("[",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = gsub("]",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = gsub(" ",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = gsub("-",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = as.factor(batch))



# 150 run ----

## Select raw expression data ----


# get read count tables
tmp.gsam.gene.expression.all <- gsam.rnaseq.expression %>%
  dplyr::select(
    tmp.combined.metadata %>%
      dplyr::filter(dataset == "G-SAM") |>
      dplyr::pull(sid)
  ) %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3)
stopifnot(colnames(tmp.gsam.gene.expression.all) == 
            tmp.combined.metadata %>%
            dplyr::filter(dataset == "G-SAM") |>
            dplyr::pull(sid)
)


tmp.glass.gene.expression.all <- glass.gbm.rnaseq.expression.all.samples %>%
  dplyr::select(
    tmp.combined.metadata %>%
      dplyr::filter(dataset == "GLASS") |>
      dplyr::pull(sid)
  ) %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3)
stopifnot(colnames(tmp.glass.gene.expression.all) ==     tmp.combined.metadata %>%
            dplyr::filter(dataset == "GLASS") |>
            dplyr::pull(sid)
)



### combine ----

# DACH1 has distinct ENSEMBL ID's across the datasets
#"ENSG00000165659" %in% (results.out |> dplyr::filter(!is.na(TCGA.subtype.marker ))  |> dplyr::pull(ensembl_id))
#"ENSG00000276644" %in% (results.out |> dplyr::filter(!is.na(TCGA.subtype.marker ))  |> dplyr::pull(ensembl_id))

# Merge count tables, perform VST transformation and perform batch effect correction
tmp.combined.gene.expression <- dplyr::inner_join(
  tmp.gsam.gene.expression.all %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::mutate(gid = gsub('^(ENSG[0-9]+).+$','\\1', gid)) # change to ENS id's
  ,
  tmp.glass.gene.expression.all %>%
    tibble::rownames_to_column('gid') |> 
    dplyr::mutate(gid = ifelse(gid == "ENSG00000165659", "ENSG00000276644", gid)) # DACH1 equivalent
  ,
  by=c('gid'='gid') ) %>%
  tibble::column_to_rownames('gid') %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% (results.out |> dplyr::filter(!is.na(TCGA.subtype.marker)) |> dplyr::pull(ensembl_id))) |> 
  tibble::column_to_rownames('gid') %>%
  limma::removeBatchEffect(tmp.combined.metadata$batch) # remove batch effects
stopifnot(tmp.combined.metadata$sid == colnames(tmp.combined.gene.expression))
stopifnot(nrow(tmp.combined.gene.expression) == 150)

rm(tmp.gsam.gene.expression.all, tmp.glass.gene.expression.all)



# export ----

# Export to 'GCT' format
tmp = tmp.combined.gene.expression |> 
  as.data.frame() |> 
  tibble::rownames_to_column('ensembl_id') |> 
  dplyr::left_join(results.out |> dplyr::select(ensembl_id, hugo_symbol), by=c('ensembl_id'='ensembl_id'),suffix=c('','')) |> 
  dplyr::rename(Name = hugo_symbol) |> 
  dplyr::mutate(Description = "description") |> 
  dplyr::mutate(ensembl_id = NULL) |> 
  dplyr::relocate(Name, Description, .before = tidyselect::everything()) # move columns to front


# there seems no header line argument in write.table so do it manually
write.table("#comment-line1", file="/tmp/tmp.txt",quote=F,col.name=F,row.name=F) 
write.table("#comment-line2", file="/tmp/tmp.txt",quote=F,col.name=F,row.name=F,append=T) 
write.table(stringi::stri_paste(colnames(tmp), collapse='\t'), file="/tmp/tmp.txt",quote=F,col.name=F,row.name=F,append=T)
write.table(tmp, file="/tmp/tmp.gct",append=T, quote=F, row.name=F, col.name=F,sep="\t")


# x-check with these lines:
#read.table("/tmp/tmp.gct",header=T,row.names=1,skip=2)[-1][1:10,1:10]



# proceed ----


# get and unzip; ssgsea.GBM.classification_1.0.tar.gz - attached to wang supplement
# R CMD INSTALL data/wang/ssgsea.GBM.classification

detach("dplyr")
library(ssgsea.GBM.classification)


# gives strange rm() warnings of variables that don't exist, but still completes
runSsGSEAwithPermutation("/tmp/tmp.gct",100000)


file.rename(from='/tmp/p_result_tmp.gct.txt',
           to='output/tables/ssgsea.GBM.classification_p_result_tmp.gct.txt')


