#!/usr/bin/env R 


# load data and metadata ----


if(!exists('glass.gbm.rnaseq.metadata.all.samples') | !exists('glass.gbm.rnaseq.expression.all.samples')) {
  source('scripts/load_glass_expression_data.R')
}



metadata <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(condition = factor(ifelse(is.primary,"primary","recurrence"), levels=c('primary','recurrence')) )


expression.data <- glass.gbm.rnaseq.expression.all.samples |> 
  dplyr::select(metadata$aliquot_barcode)



stopifnot(colnames(expression.data) == metadata$aliquot_barcode)


# DE unpaired all [GLASS 2022] ----

#'@todo x-check batches w/ bam's

analysis.dge.glass.2022 <- expression.data %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3) %>%
  DESeq2::DESeqDataSetFromMatrix( metadata, ~condition ) %>% # + resection
  DESeq2::DESeq(parallel = T) %>%
  DESeq2::results() %>% 
  as.data.frame(stringsAsFactors=F) %>%
  dplyr::arrange(padj, log2FoldChange) %>%
  tibble::rownames_to_column('ensembl_id') %>% 
  dplyr::rename_with( ~ paste0(.x, ".glass-2022.res")) %>%
  dplyr::rename(ensembl_id = `ensembl_id.glass-2022.res`)


# export ----


saveRDS(analysis.dge.glass.2022,"cache/analysis_DGE_GLASS-2020.Rds")




# cleanup ----


rm(expression.data, metadata)





