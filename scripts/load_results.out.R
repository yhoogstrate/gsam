#!/usr/bin/env R

# load generic file ----

results.out <- readRDS(file = 'tmp/results.out.Rds') 



# DGE GLASS 2022 ----


tmp <- readRDS("tmp/analysis_DGE_GLASS-2022.tpc.Rds")


results.out <- results.out |> 
  dplyr::left_join(
    tmp,
    by=c('ensembl_id'='ensembl_id'), suffix=c('','')
  )

rm(tmp)


# cor GLASS 2022 ----

tmp <- readRDS('tmp/analysis_cor_purity_expression_GLASS-2022.Rds')

results.out <- results.out |> 
  dplyr::left_join(
    tmp,
    by=c('ensembl_id'='ensembl_id'), suffix=c('','')
  )

rm(tmp)


