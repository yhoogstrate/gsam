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

# ---- boxplot + jitter ----


ggplot2::ggplot(data=plt, mapping=aes(x=gliovis.majority_call , y  = tumour.percentage.dna )) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  youri_gg_theme +
  labs(colour="GBM subtype")


# ---- PCA on subtype genes ----

plt.data <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
  t()

plt <- prcomp(plt.data)$x %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
  , by = c('sid' = 'sid')) %>%
  dplyr::mutate(tumour.percentage.dna = ifelse(is.na(tumour.percentage.dna), -1 , tumour.percentage.dna) ) %>%
  dplyr::mutate(tumour.percentage.dna.bin = cut(tumour.percentage.dna, breaks=3) ) %>%
  dplyr::mutate(tumour.percentage.dna = ifelse(tumour.percentage.dna == -1, NA , tumour.percentage.dna) )



# ---- t-SNE on subtype genes ----


plt <- Rtsne::Rtsne(plt.data)$Y %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(c('tsne1', 'tsne2')) %>%
  dplyr::mutate(sid = rownames(plt.data)) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
    , by = c('sid' = 'sid')) %>%
  dplyr::mutate(pid = gsub("^(...).+$","\\1", sid))


ggplot2::ggplot(plt, aes(x = tsne1, y=tsne2, col = gliovis.majority_call, shape =  tumour.percentage.dna >= 10, group=pid)) + 
  geom_line(alpha=0.2) +
  geom_point() +
  youri_gg_theme +
  labs(colour="GBM subtype") +
  scale_shape_manual(values=c(19, 1))



# ---- UMAP on subtype genes ----


plt <- uwot::umap(plt.data) %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(c('umap1', 'umap2')) %>%
  dplyr::mutate(sid = rownames(plt.data)) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
    , by = c('sid' = 'sid')) %>%
  dplyr::mutate(pid = gsub("^(...).+$","\\1", sid))



ggplot2::ggplot(plt, aes(x = umap1, y=umap2, col = gliovis.majority_call, shape =  tumour.percentage.dna >= 10, group=pid)) + 
  geom_line(alpha=0.2) +
  geom_point() +
  youri_gg_theme +
  labs(colour="GBM subtype") +
  scale_shape_manual(values=c(19, 1))



