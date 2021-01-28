#!/usr/bin/env R

# load libs ----

library(tidyverse)
library(NMF)
library(scales)


# load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/gsam_rna-seq_expression.R')

source('data/wang/msig.library.12.R') # license is incomplete regarding sharing, can't include it in source tree

source('scripts/R/subtype_genes.R')
source('scripts/R/wang_glioma_intrinsic_genes.R')

source('scripts/R/youri_gg_theme.R')


source('scripts/R/palette.R')


# NMF ----

seeds <- c(123456, 12345, 1337, 909)

gsam_nmf_150 <- {{}}
gsam_nmf_15k <- {{}}
gsam_nmf_7k <- {{}}

## Normalize for NMF ----

### normalize the VST values ----
tmp <- gsam.rnaseq.expression.vst

if(min(tmp) < 0) { # not the case
  tmp <- tmp - min(tmp) + .Machine$double.eps
}


### make 15k set ----
### make 7k set ----
### make 150 set ----


metadata <- gsam.rna.metadata %>%
  dplyr::filter(sid %in% colnames(tmp))


gsam.rnaseq.expression.vst.150 <- tmp %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(ensg = gsub("\\..+$","",gid)) %>%
  dplyr::filter(ensg %in%
                  (wang.glioma.intrinsic.genes %>%
                    dplyr::filter(Subtyping_Signature_Gene. != "") %>%
                    dplyr::pull(ENSG.short))
                  ) %>%
  dplyr::mutate(ensg = NULL) %>%
  tibble::column_to_rownames('gid')


stopifnot(ncol(gsam.rnaseq.expression.vst.150) == nrow(metadata))
stopifnot(colnames(gsam.rnaseq.expression.vst.150) %in% metadata$sid)
stopifnot( metadata$sid %in% colnames(gsam.rnaseq.expression.vst.150))


gsam.rnaseq.expression.vst.150 <- gsam.rnaseq.expression.vst.150 %>%
  dplyr::select(metadata$sid)


stopifnot(colnames(gsam.rnaseq.expression.vst.150) == metadata$sid)



### exclude IDH muts [after normalisation, so they can be mapped back into the coordinates later] ----

## Perform NMF by Wang-code [seeds: ] ----

# TEST REGULAR PCA - SHOULD WORK(!)
p <- prcomp(t(gsam.rnaseq.expression.vst.150))#screeplot(p)


plt.pca.150   <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], gliovis = as.factor(metadata$gliovis.majority_call))
plt.pca.150.p <- ggplot(plt.pca.150, aes(x = PC1, y = PC2, col = gliovis) ) + 
  geom_point(size=1.5) +
  youri_gg_theme + 
  scale_color_manual(name = NULL, values =  c('Classical'='#6ba6e5',# blue
                                              'Mesenchymal'='#eab509',#mustard
                                              'Proneural'='#ff5f68'),#red/pink
                     guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))




for(seed in seeds) {
  gsam_nmf_150[[as.character(seed)]] <- NMF(as.matrix(gsam.rnaseq.expression.vst.150), 3, seed = seed)
}

gsam_nmf_150[["nmf"]] <- nmf(as.data.frame(gsam.rnaseq.expression.vst.150), 3, seed=12345)


plt.nmf.150   <- as.data.frame(t(gsam_nmf_150$`123456`$H)) %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::rename(NMF1 = V1, NMF2 = V2, NMF3 = V3) %>%
  dplyr::left_join(metadata %>% dplyr::select(c('sid','gliovis.majority_call','pat.with.IDH')), by=c('sid' = 'sid')) %>%
  dplyr::rename(gliovis = gliovis.majority_call)

plt.nmf.150.p <- ggplot(plt.nmf.150, aes(x = NMF1, y = NMF2, col = gliovis) ) + 
  geom_point(size=1.5) +
  youri_gg_theme + 
  scale_color_manual(name = NULL, values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_x_continuous(limits= c(min(plt.nmf.150$NMF1), max(plt.nmf.150$NMF1)) * 1.1  ) +
  scale_y_continuous(limits= c(min(plt.nmf.150$NMF2), max(plt.nmf.150$NMF2)) * 1.1  )



#gsam_nmf_150$`909`$membership
#gsam_nmf_150$`12345`$membership
#gsam_nmf_150$`123456`$membership
#gsam_nmf_150$`1337`$membership


# comparison w/ gliovis
#plt <- data.frame(nmf = as.factor(gsam_nmf_150$`1337`$membership) ,
#                  gliovis = as.character(metadata$gliovis.majority_call ) ) %>%
#  dplyr::mutate(gliovis = ifelse(is.na(gliovis), 'NA', as.character(gliovis)  )) %>%
#  table() %>%
#  as.data.frame()

#ggplot(plt, aes(x = nmf, y = Freq, col=gliovis)  ) +
#  geom_bar(stat="identity")




# clearly the nicest fit (!!!!!)
# PCA to also fit NMF 3 into two dimensions (is a better fit)
p <- prcomp(t(gsam_nmf_150$`123456`$H))
plt.nmf.pca.150   <- data.frame(NMF.PC1 = p$x[,1], NMF.PC2 = p$x[,2], gliovis = as.factor(metadata$gliovis.majority_call))
plt.nmf.pca.150.p <- ggplot(plt.nmf.pca.150 , aes(x = NMF.PC1, y = NMF.PC2, col = gliovis) ) + 
  geom_point(size=1.5) +
  youri_gg_theme + 
  scale_color_manual(name = NULL, values =  subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_x_continuous(limits= c(min(plt.nmf.pca.150$NMF.PC1), max(plt.nmf.pca.150$NMF.PC1)) * 1.1  ) +
  scale_y_continuous(limits= c(min(plt.nmf.pca.150$NMF.PC2), max(plt.nmf.pca.150$NMF.PC2)) * 1.1  )
  



# perform linear classifier on the PCA
train.in <- p$x %>%
  as.data.frame() %>%
  dplyr::mutate(gliovis = metadata$gliovis.majority_call)

test.lda <- MASS::lda(gliovis ~ PC1 + PC2, data = train.in)


#plot(c(-6,6), c(-6,6), type="n")
#points(
#  p$x[,1],
#  p$x[,2],
#  col = as.factor(metadata$gliovis.majority_call)
#  #col = as.factor(gsam_nmf_150$`123456`$membership)
#  ,pch=19, cex=0.6
#)# screeplot(p)

# ---- brute force class labels ----

resolution <- 250 # 1000 x 1000 data points

range_pc1 = seq(from = min(p$x[,1]), to = max(p$x[,1]), length.out = resolution) * 1.1
range_pc2 = seq(from = min(p$x[,2]), to = max(p$x[,2]), length.out = resolution) * 1.1

range_df = expand.grid(PC1 = range_pc1, PC2 = range_pc2)
range_predictions = predict(test.lda, newdata = range_df)



plt.nmf.pca.lda.150 <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], class = metadata$gliovis.majority_call, type = "GlioVis") %>%
  rbind( data.frame(PC1 = range_df$PC1 ,PC2 = range_df$PC2 ,class = range_predictions$class,type = "LDA on NVM+PCA"))

plt.nmf.pca.lda.150.p <- ggplot(plt.nmf.pca.lda.150, aes(x = PC1, y = PC2, col = class)) + 
  geom_raster(data = subset(plt.nmf.pca.lda.150, type != "GlioVis"), aes(fill = factor(class)),alpha=0.1) +
  geom_contour(data= subset(plt.nmf.pca.lda.150, type != "GlioVis"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2
               ,breaks=c(1.5,2.5)
               ) +
  geom_point(data = subset(plt.nmf.pca.lda.150, type == "GlioVis"), size=1.5) +
  youri_gg_theme +
  labs(x = "PC1 on NMF[k=3] (150 S-Genes)",
       y = "PC2 on NMF[k=3] (150 S-Genes)",
       col = "GlioVis Majority subtype",
       fill = "NMF/PCA/LDA subtype") + 
  scale_color_manual(name = NULL, values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)







plt.nmf.pca.lda.reclass.150 <- plt.nmf.pca.lda.150 %>%
  dplyr::filter(type == 'GlioVis') %>%
  dplyr::select(c('PC1','PC2'))

range_predictions = predict(test.lda, newdata = plt.nmf.pca.lda.reclass.150)

plt.nmf.pca.lda.reclass.150$class <- range_predictions$class
plt.nmf.pca.lda.reclass.150$type <- "LDA on NVM+PCA data"
#plt.nmf.pca.lda.reclass.150 <- rbind(plt.nmf.pca.lda.reclass.150, plt.nmf.pca.lda.150 %>%
#                                       dplyr::filter(type != 'GlioVis'))

rm(range_predictions)




plt.nmf.pca.lda.reclass.150.p <- ggplot(plt.nmf.pca.lda.reclass.150, aes(x = PC1, y = PC2, col = class)) + 
  geom_contour(data= subset(plt.nmf.pca.lda.150, type != "GlioVis"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2
               ,breaks=c(1.5,2.5)
  ) +
  geom_point(data = subset(plt.nmf.pca.lda.reclass.150, type == "LDA on NVM+PCA data"), size=1.5) +
  youri_gg_theme +
  labs(x = "PC1 on NMF[k=3] (150 S-Genes)",
       y = "PC2 on NMF[k=3] (150 S-Genes)",
       col = "GlioVis Majority subtype",
       fill = "NMF/PCA/LDA subtype") + 
  scale_color_manual(name = NULL, values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



library(patchwork)

#plt.pca.150.p + 
  
plt.nmf.150.p + plt.nmf.pca.150.p + plt.nmf.pca.lda.150.p + plt.nmf.pca.lda.reclass.150.p

## Perform PCA on NMF coordinates ? ----

## Export ----



