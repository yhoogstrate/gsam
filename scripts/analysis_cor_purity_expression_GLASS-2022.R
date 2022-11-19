#!/usr/bin/env R 


# load data and metadata ----


if(!exists('glass.gbm.rnaseq.metadata.all.samples') | !exists('glass.gbm.rnaseq.expression.all.samples')) {
  source('scripts/load_glass_expression_data.R')
}


## shape data ----


stopifnot(glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode == colnames(glass.gbm.rnaseq.expression.all.samples))

tmp.metadata <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(cond = runif(n(),1,2))  |> 
  dplyr::mutate(cond = round(cond))  |> 
  dplyr::mutate(cond = paste0("r", cond)) |> 
  dplyr::mutate(cond = as.factor(cond))



expression.vst <- glass.gbm.rnaseq.expression.all.samples |> 
  dplyr::select(tmp.metadata$aliquot_barcode) |> 
  DESeq2::DESeqDataSetFromMatrix(tmp.metadata, ~cond) |> 
  DESeq2::vst(blind=T) |> 
  SummarizedExperiment::assay() |> 
  limma::removeBatchEffect(as.factor(tmp.metadata$aliquot_batch_synapse))




## perform cor ----


stopifnot(colnames(expression.vst) == tmp.metadata$aliquot_barcode)

cor.stats <- expression.vst |> # remove batch effects: aliquot_batch_synapse
  t() |> 
  as.data.frame(stringsAsFactors=F) |> 
  pbapply::pblapply(cor.test , y=tmp.metadata$tumour.percentage.2022, method = "pearson") |> 
  #pbapply::pblapply(cor.test , y=tmp.metadata$tumour.percentage.dna.imputed.rf.2022.all.patients.B, method = "pearson") |>
  unlist() |> 
  data.frame() |> 
  `colnames<-`('value') |> 
  tibble::rownames_to_column('key') |> 
  dplyr::mutate(ensembl_id = gsub("^([^\\.]+)\\..+$","\\1",key))|> 
  dplyr::mutate(key = gsub("^[^\\.]+\\.(.+)$","\\1",key)) |> 
  tidyr::pivot_wider(names_from = key, values_from = value) |> 
  dplyr::select(ensembl_id, estimate.cor, statistic.t, p.value) |> 
  dplyr::rename_with( ~ paste0(.x, ".glass-2022.cor.tpc")) |> 
  dplyr::rename(`ensembl_id` = `ensembl_id.glass-2022.cor.tpc`) |> 
  dplyr::mutate(`estimate.cor.glass-2022.cor.tpc` = as.numeric( `estimate.cor.glass-2022.cor.tpc` )) |> 
  dplyr::mutate(`statistic.t.glass-2022.cor.tpc` = as.numeric( `statistic.t.glass-2022.cor.tpc` )) |> 
  dplyr::mutate(`p.value.glass-2022.cor.tpc` = as.numeric( `p.value.glass-2022.cor.tpc` )) |> 
  as.data.frame()


## export ----


saveRDS(cor.stats, file='tmp/analysis_cor_purity_expression_GLASS-2022.Rds')


## cleanup ----


rm(expression.vst, tmp.metadata, cor.stats)




