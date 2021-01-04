#!/usr/bin/env R

# ---- load libs ----

library(tidyverse)

# ---- load cfg ----

source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')

# ---- Load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/expression_matrix.R')
source('scripts/R/expression_matrix_VST.R')

# ---- PCA only the subtype genes ----

plt <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
  t()

plt <- prcomp(plt)$x %>%
  data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
  , by = c('sid' = 'sid')) %>%
  dplyr::mutate(tumour.percentage.dna = ifelse(is.na(tumour.percentage.dna), -1 , tumour.percentage.dna) ) %>%
  dplyr::mutate(tumour.percentage.dna.bin = cut(tumour.percentage.dna, breaks=3) )





ggplot2::ggplot(data=plt, mapping=aes(x=PC1, y=PC2,
                                      #col=gliovis.majority_call,
                                      #shape =  tumour.percentage.dna.bin
                                      col =  tumour.percentage.dna.bin
                                      )) +
  geom_point() +
  youri_gg_theme +
  labs(colour="GBM subtype")





ggplot2::ggplot(data=plt, mapping=aes(x=gliovis.majority_call ,
                                      y  = tumour.percentage.dna
)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  #youri_gg_theme +
  labs(colour="GBM subtype")




