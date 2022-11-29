#!/usr/bin/env R


# generates input table for GlioVis/SubtypeME


# load expression data ----


source('scripts/load_G-SAM_expression_data.R')
source('scripts/load_GLASS_data.R')
source('scripts/load_results.out.R')


# select and merge samples ----


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
    dplyr::filter(tumour.percentage.2022 >= 15) |>
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



stopifnot(nrow(tmp.combined.metadata) == (287+216))


# 150 run ----

## export raw expression data ----


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
    dplyr::pull(sid))


tmp.glass.gene.expression.all <- glass.gbm.rnaseq.expression.all.samples %>%
  dplyr::select(
    tmp.combined.metadata %>%
      dplyr::filter(dataset == "GLASS") |>
      dplyr::pull(sid)
  ) %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3)
stopifnot(colnames(tmp.glass.gene.expression.all) == tmp.combined.metadata %>%
  dplyr::filter(dataset == "GLASS") |>
  dplyr::pull(sid))



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
  
  dplyr::left_join(results.out |> dplyr::select(hugo_symbol, TCGA.subtype.marker, ensembl_id), 
                   by=c('gid'='ensembl_id'),
                   suffix=c('','')) |> 
  dplyr::filter(!is.na(TCGA.subtype.marker)) %>%
  dplyr::mutate(TCGA.subtype.marker = NULL) %>%
  dplyr::mutate(ensembl_id = NULL) %>%
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  limma::removeBatchEffect(tmp.combined.metadata$batch) # remove batch effects
stopifnot(tmp.combined.metadata$sid == colnames(tmp.combined.gene.expression))
stopifnot(nrow(tmp.combined.gene.expression) == 150)

rm(tmp.gsam.gene.expression.all, tmp.glass.gene.expression.all)


### fix names ----


export <- tmp.combined.gene.expression %>%
  as.data.frame() %>%
  t() %>%
  `colnames<-`(gsub("^.+\\|([^\\|]+)\\|.+$","\\1",colnames(.))) %>%
  `colnames<-`(gsub('^LHFPL6$', 'LHFP',colnames(.))) %>%
  `colnames<-`(gsub('^JPT1$', 'HN1',colnames(.))) %>%
  `colnames<-`(gsub('^PLAAT1$', 'HRASLS',colnames(.))) %>%
  `colnames<-`(gsub('^PAK5$', 'PAK7',colnames(.))) %>%
  `colnames<-`(gsub('^ZFP96B$', 'ZNF643',colnames(.))) %>%
  `colnames<-`(gsub('^DGLUCY$', 'C14orf159',colnames(.))) %>%
  `colnames<-`(gsub('^EFCAB14$', 'KIAA0494',colnames(.))) %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('Sample')

export[] <- lapply(export, as.character) # quotes..

write.csv(export,"output/tables/analysis_SubtypeME_export_2022.csv", row.names=F, quote=T)

# web results saved to: output/tables/analysis_SubtypeME_GlioVis_output_2022.csv


