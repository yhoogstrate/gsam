#!/usr/bin/env R 


# load data and metadata ----


if(!exists('glass.gbm.rnaseq.metadata.all.samples') | !exists('glass.gbm.rnaseq.expression.all.samples')) {
  source('scripts/load_GLASS_data.R')
}



metadata <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(condition = factor(ifelse(is.primary,"primary","recurrence"), levels=c('primary','recurrence'))) |> 
  dplyr::mutate(aliquot_batch_synapse = gsub("[",".",aliquot_batch_synapse,fixed=T)) |> 
  dplyr::mutate(aliquot_batch_synapse = gsub("]",".",aliquot_batch_synapse,fixed=T)) |> 
  dplyr::mutate(aliquot_batch_synapse = gsub(" ",".",aliquot_batch_synapse,fixed=T)) |> 
  dplyr::mutate(aliquot_batch_synapse = gsub("-",".",aliquot_batch_synapse,fixed=T)) |>  
  dplyr::mutate(aliquot_batch_synapse = as.factor(aliquot_batch_synapse)) |> 
  dplyr::mutate(tpc = 1 - (tumour.percentage.2022 / 100))
  




expression.data <- glass.gbm.rnaseq.expression.all.samples |> 
  dplyr::select(metadata$aliquot_barcode)



stopifnot(colnames(expression.data) == metadata$aliquot_barcode)

dim(expression.data)


# DE unpaired all [GLASS 2022] ----

# - have no elegant substitution for rowsums and ncol filter in 1 cmd
#   #dplyr::rowwise() |> 
#   #dplyr::mutate(rowsum = sum(dplyr::across())) |> 
#   #dplyr::ungroup() |> 
#   #dplyr::mutate(ncol = length(dplyr::across()) - 1) |> 
#   #dplyr::filter(rowsum > ncol * 3) |> 
#   #dplyr::mutate(ncol = NULL, rowsum = NULL) |> 



stopifnot(sum(is.na(expression.data)) == 0)
stopifnot(sum(is.na(metadata$tumour.percentage.2022)) == 0)
stopifnot(metadata$aliquot_barcode == colnames(expression.data))

analysis.dge.glass.2022 <- expression.data %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3) %>%
  DESeq2::DESeqDataSetFromMatrix(metadata, ~aliquot_batch_synapse + condition ) %>%  # + resection
  DESeq2::DESeq(parallel = T) %>%
  DESeq2::results() %>%
  as.data.frame(stringsAsFactors=F) %>%
  dplyr::arrange(padj, log2FoldChange) %>%
  tibble::rownames_to_column('ensembl_id') %>%
  dplyr::rename_with( ~ paste0(.x, ".glass-2022.res")) %>%
  dplyr::rename(ensembl_id = `ensembl_id.glass-2022.res`)



saveRDS(analysis.dge.glass.2022, "tmp/analysis_DGE_GLASS-2022.Rds")



# DE unpaired all ~ tpc [GLASS 2022] ----


stopifnot(sum(is.na(expression.data)) == 0)
stopifnot(sum(is.na(metadata$tumour.percentage.2022)) == 0)
stopifnot(metadata$aliquot_barcode == colnames(expression.data))

analysis.dge.glass.tpc.2022 <- expression.data %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3) %>%
  DESeq2::DESeqDataSetFromMatrix(metadata, ~ aliquot_batch_synapse + tpc + condition) %>% # + resection
  DESeq2::DESeq(parallel = F) %>%
  DESeq2::results() %>%
  as.data.frame(stringsAsFactors=F) %>%
  dplyr::arrange(padj, log2FoldChange) %>%
  tibble::rownames_to_column('ensembl_id') %>%
  dplyr::rename_with( ~ paste0(.x, ".glass-2022.tpc.res")) %>%
  dplyr::rename(ensembl_id = `ensembl_id.glass-2022.tpc.res`)


analysis.dge.glass.tpc.2022 |>
  dplyr::filter(ensembl_id == "ENSG00000164692") |> 
  dplyr::mutate(`baseMean.glass-2022.tpc.res` = NULL)




saveRDS(analysis.dge.glass.tpc.2022, "tmp/analysis_DGE_GLASS-2022.tpc.Rds")







# sandbox ----



## view pwr ----

readRDS("tmp/DEpower.Rds") |> 
  dplyr::arrange(ENSG00000164692) |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |> 
      dplyr::select(aliquot_barcode, 
                    #tumour.percentage.dna.cnv.2022,
                    tumour.percentage.2022),
    by=c('excl'='aliquot_barcode')
  ) |> 
  head(n=10)


## more ----

analysis.dge.glass.tpc.2022 |> 
  dplyr::filter(`padj.glass-2022.tpc.res` < 0.01 & abs(`log2FoldChange.glass-2022.tpc.res`) > 0.5) |> 
  dim()


results.out <- readRDS(file = 'tmp/results.out.Rds')


sig.2021 <- results.out %>%
  dplyr::filter( !is.na(log2FoldChange.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(log2FoldChange.glass.tpc.res)  ) %>%
  dplyr::filter( !is.na(padj.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(padj.glass.tpc.res)  ) %>%
  
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") ) %>%
  
  dplyr::mutate(significant.2021 = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 &
                  abs(log2FoldChange.glass.tpc.res) > 0.5 & 
                  direction.gsam.tpc.res == direction.glass.tpc.res) |> 
  
  dplyr::filter(significant.2021) |> 
  dplyr::pull(ensembl_id)





sig.2022 <- results.out |>
  dplyr::select(-contains("glass.tpc")) |> 
  
  dplyr::left_join(analysis.dge.glass.tpc.2022, by=c('ensembl_id'='ensembl_id'),suffix=c('','') ) %>%
  
  dplyr::filter( !is.na(log2FoldChange.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(`log2FoldChange.glass-2022.tpc.res`)  ) %>%
  dplyr::filter( !is.na(padj.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(`padj.glass-2022.tpc.res`)  ) %>%
  
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(`log2FoldChange.glass-2022.tpc.res` > 0 , "up", "down") ) %>%
  
  dplyr::mutate(significant.2022 = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 &
                  abs(`log2FoldChange.glass-2022.tpc.res`) > 0.5 & 
                  direction.gsam.tpc.res == direction.glass.tpc.res) |> 
  
  dplyr::filter(significant.2022 == T) |> 
  dplyr::pull(ensembl_id)


dim(plt)


length(intersect(sig.2021,sig.2022))



results.out |> 
  dplyr::filter(ensembl_id %in% sig.2021) |> 
  dplyr::filter(ensembl_id %in% sig.2022 == F) |> 
  View()


results.out |> 
  dplyr::filter(ensembl_id %in% sig.2021) |> 
  dplyr::filter(ensembl_id %in% sig.2022) |> 
  View()



# cleanup ----


rm(expression.data, metadata)

rm(analysis.dge.glass.2022)
rm(analysis.dge.glass.tpc.2022)



