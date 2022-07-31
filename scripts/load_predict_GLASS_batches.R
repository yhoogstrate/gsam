#!/usr/bin/env R


stopifnot(file.exists('cache/analysis_predict_GLASS_batches.R'))

if(!exists('glass.gbm.rnaseq.metadata.all.samples')) {
  source('scripts/load_glass_expression_data.R') # w/ metadata
}


# load ----



tmp <- readRDS("cache/analysis_predict_GLASS_batches.R")


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::left_join(tmp, by=c('aliquot_barcode'='aliquot_barcode'), suffix=c('',''))


rm(tmp)



plot(table(glass.gbm.rnaseq.metadata.all.samples$predicted.GLASS.batch,glass.gbm.rnaseq.metadata.all.samples$aliquot_batch_synapse))


plt <- table(glass.gbm.rnaseq.metadata.all.samples$predicted.GLASS.batch,glass.gbm.rnaseq.metadata.all.samples$aliquot_batch_synapse)
plt <- as.data.frame(table(plt))



ggplot(plt, aes(x=Var1))+
  geom_bar(y = Freq)
