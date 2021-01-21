#!/usr/bin/env R

# ---- load libs ----

library(tidyverse)

# ---- load cfg ----

source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')

# ---- Load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/gsam_expression_matrix.R')
source('scripts/R/gsam_expression_matrix_VST.R')
source('scripts/R/subtype_genes.R')


# ---- export subtype data ----

plt.data <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
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

plt.data[] <- lapply(plt.data, as.character) # quotes..

#write.csv(plt.data,"output/tables/gliovis-subtype-input-table.csv", row.names=F, quote=T)



#rm(plt.data)

  
# ---- boxplot + jitter ----


plt <- gsam.rna.metadata %>%
  dplyr::filter(!is.na(gliovis.majority_call))

ggplot2::ggplot(data=plt, mapping=aes(x=gliovis.majority_call , y  = tumour.percentage.dna )) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  youri_gg_theme +
  labs(colour="GBM subtype")


# ---- PCA on subtype genes ----

plt.data <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
  t()

plt <- prcomp(plt.data, scale=F)$x %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
  , by = c('sid' = 'sid')) %>%
  dplyr::mutate(tumour.percentage.dna = ifelse(is.na(tumour.percentage.dna), -1 , tumour.percentage.dna) ) %>%
  dplyr::mutate(tumour.percentage.dna.bin = cut(tumour.percentage.dna, breaks=3) ) %>%
  dplyr::mutate(tumour.percentage.dna = ifelse(tumour.percentage.dna == -1, NA , tumour.percentage.dna) ) %>%
  dplyr::mutate(pid = gsub("^(...).+$","\\1", sid)) %>%
  dplyr::mutate(resection = as.factor(gsub("^...(.).*$","R\\1", sid)) )


# PC2 correleert met pecentage
# maar dit kan komen omdat die percentages anders zijn per subtype
# binnen een subtype moet deze correlatie eigenlijk best laag zijn... checken dus


# PC3 = tumour purity...
ggplot2::ggplot(plt, aes(x = PC1, y=PC2, col = gliovis.majority_call, shape =  tumour.percentage.dna >= 20, group=pid)) + 
  geom_line(alpha=0.2, col="gray60") +
  geom_point(size = 3) +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype") +
  scale_shape_manual(values=c(19, 1))
ggsave("output/figures/paper_subtype_pca.png",width=10 * 1.3, height=6 * 1.3)

# PC3 = tumour purity...
ggplot2::ggplot(plt, aes(x = PC1, y=PC2, col = gliovis.majority_call, shape = resection, group=pid)) + 
  geom_line(alpha=0.2, col="gray60") +
  geom_point(size = 3) +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype") +
  scale_shape_manual(values=c(19, 1))
ggsave("output/figures/paper_subtype_pca_res.png",width=10 * 1.3, height=6 * 1.3)


#plot(plt$tumour.percentage.dna , plt$PC2)
#plot(plt$tumour.percentage.dna , plt$PC3)

ggplot(plt , aes(x = PC2 , y = PC3, col = tumour.percentage.dna.bin) ) +
  geom_point()

# ---- t-SNE on subtype genes ----


plt.data <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
  t()

plt <- Rtsne::Rtsne(plt.data)$Y %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(c('tsne1', 'tsne2')) %>%
  dplyr::mutate(sid = rownames(plt.data)) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
    , by = c('sid' = 'sid')) %>%
  dplyr::mutate(pid = gsub("^(...).+$","\\1", sid)) %>%
  dplyr::mutate(resection = as.factor(gsub("^...(.).*$","R\\1", sid)) )


ggplot2::ggplot(plt, aes(x = tsne1, y=tsne2, col = gliovis.majority_call, shape =  tumour.percentage.dna >= 10, group=pid)) + 
  geom_line(alpha=0.2, col="gray60") +
  geom_point(size = 3) +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype",x="t-SNE 1", y="t-SNE 2") +
  scale_shape_manual(values=c(19, 1))

ggsave("output/figures/paper_subtype_tSNE.png",width=10 * 1.3, height=6 * 1.3)


ggplot2::ggplot(plt, aes(x = tsne1, y=tsne2, col = gliovis.majority_call, shape =  resection , group=pid)) + 
  geom_line(alpha=0.2, col="gray60") +
  geom_point(size = 3) +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype",x="t-SNE 1", y="t-SNE 2") +
  scale_shape_manual(values=c(19, 1))

ggsave("output/figures/paper_subtype_tSNE_res.png",width=10 * 1.3, height=6 * 1.3)



# ---- UMAP on subtype genes ----

plt.data <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
  t()

plt <- uwot::umap(plt.data) %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(c('umap1', 'umap2')) %>%
  dplyr::mutate(sid = rownames(plt.data)) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::select(c('sid', 'gliovis.majority_call', 'tumour.percentage.dna'))
    , by = c('sid' = 'sid')) %>%
  dplyr::mutate(pid = gsub("^(...).+$","\\1", sid)) %>%
  dplyr::mutate(resection = as.factor(gsub("^...(.).*$","R\\1", sid)) )



ggplot2::ggplot(plt, aes(x = umap1, y=umap2, col = gliovis.majority_call, shape =  tumour.percentage.dna >= 10, group=pid)) + 
  geom_line(alpha=0.2, col="gray60") +
  geom_point(size = 3) +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype", x="UMAP 1" , y = "UMAP 2") +
  scale_shape_manual(values=c(19, 1))

ggsave("output/figures/paper_subtype_UMAP.png",width=10 * 1.3, height=6 * 1.3)




ggplot2::ggplot(plt, aes(x = umap1, y=umap2, col = gliovis.majority_call, shape =  resection , group=pid)) + 
  geom_line(alpha=0.2, col="gray60") +
  geom_point(size = 3) +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype", x="UMAP 1" , y = "UMAP 2") +
  scale_shape_manual(values=c(19, 1))

ggsave("output/figures/paper_subtype_UMAP_res.png",width=10 * 1.3, height=6 * 1.3)




