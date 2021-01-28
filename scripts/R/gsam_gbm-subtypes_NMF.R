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
p <- prcomp(t(gsam.rnaseq.expression.vst.150))
plot(
  p$x[,1],
  p$x[,2],
  
  col = as.factor(metadata$gliovis.majority_call)
)

#screeplot(p)




for(seed in seeds) {
  gsam_nmf_150[[as.character(seed)]] <- NMF(as.matrix(gsam.rnaseq.expression.vst.150), 3, seed = seed)
}

gsam_nmf_150[["nmf"]] <- nmf(as.data.frame(gsam.rnaseq.expression.vst.150), 3, seed=12345)


gsam_nmf_150$`909`$membership
gsam_nmf_150$`12345`$membership
gsam_nmf_150$`123456`$membership
gsam_nmf_150$`1337`$membership


# comparison w/ gliovis
plt <- data.frame(nmf = as.factor(gsam_nmf_150$`1337`$membership) ,
                  gliovis = as.character(metadata$gliovis.majority_call ) ) %>%
  dplyr::mutate(gliovis = ifelse(is.na(gliovis), 'NA', as.character(gliovis)  )) %>%
  table() %>%
  as.data.frame()

ggplot(plt, aes(x = nmf, y = Freq, col=gliovis)  ) +
  geom_bar(stat="identity")


# # PCA to also fit NMF 3 into two dimensions (is a better fit)
# p <- prcomp(t(gsam_nmf_150$`909`$H))
# plot(
#   p$x[,1],
#   p$x[,2],
#   
#   #col = as.factor(metadata$gliovis.majority_call)
#   col = as.factor(gsam_nmf_150$`909`$membership)
# )
# 
# # PCA to also fit NMF 3 into two dimensions (is a better fit)
# p <- prcomp(t(gsam_nmf_150$`1337`$H))
# plot(
#   p$x[,1],
#   p$x[,2],
#   
#   #col = as.factor(metadata$gliovis.majority_call)
#   col = as.factor(gsam_nmf_150$`1337`$membership)
# )
# 
# 
# # PCA to also fit NMF 3 into two dimensions (is a better fit)
# p <- prcomp(t(gsam_nmf_150$`1337`$H))
# plot(
#   p$x[,1],
#   p$x[,2],
#   
#   #col = as.factor(metadata$gliovis.majority_call)
#   col = as.factor(gsam_nmf_150$`1337`$membership)
# )
# 
# 
# # PCA to also fit NMF 3 into two dimensions (is a better fit)
# p <- prcomp(t(gsam_nmf_150$`12345`$H))
# plot(
#   p$x[,1],
#   p$x[,2],
#   
#   #col = as.factor(metadata$gliovis.majority_call)
#   col = as.factor(gsam_nmf_150$`12345`$membership)
# )



# clearly the nicest fit (!!!!!)
# PCA to also fit NMF 3 into two dimensions (is a better fit)
p <- prcomp(t(gsam_nmf_150$`123456`$H))
# perform linear classifier on the PCA
train.in <- p$x %>%
  as.data.frame() %>%
  dplyr::mutate(gliovis = metadata$gliovis.majority_call)

test.lda <- MASS::lda(gliovis ~ PC1 + PC2, data = train.in)


plot(c(-6,6), c(-6,6), type="n")
points(
  p$x[,1],
  p$x[,2],
  
  col = as.factor(metadata$gliovis.majority_call)
  #col = as.factor(gsam_nmf_150$`123456`$membership)
  
  ,pch=19, cex=0.6
)
# screeplot(p)


points( test.lda$means[,1] ,  test.lda$means[,2] , col = rgb(0.5,0.5,0.5,0.05) , pch=19, cex = 25)

# ---- brute force class labels ----

resolution <- 250 # 1000 x 1000 data points

range_pc1 = seq(from = min(p$x[,1]), to = max(p$x[,1]), length.out = resolution) * 1.2
range_pc2 = seq(from = min(p$x[,2]), to = max(p$x[,2]), length.out = resolution) * 1.2

range_df = expand.grid(PC1 = range_pc1, PC2 = range_pc2)
range_predictions = predict(test.lda, newdata = range_df)

# range_df[1:10,1:2]
points(
  range_df$PC1,
  range_df$PC2,
  
  col=as.factor(range_predictions$class), pch=19, cex=0.01 , alpha=0.1 )




plt <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], class = metadata$gliovis.majority_call, type = "GlioVis") %>%
  rbind( 
    data.frame(
      PC1 = range_df$PC1 ,
      PC2 = range_df$PC2 ,
      class = range_predictions$class,
      type = "LDA on NVM+PCA"
    )
    )


ggplot(plt, aes(x = PC1, y = PC2, col = class)) + 
  #geom_point(data = subset(plt, type != "GlioVis"), alpha = 0.15, size = 0.4) +
  geom_raster(data = subset(plt, type != "GlioVis"), aes(fill = factor(class)),alpha=0.2) +
  geom_contour(data= subset(plt, type != "GlioVis"), aes(z=as.numeric(class)), colour="gray", alpha=0.85,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "GlioVis")) +
  youri_gg_theme +
  labs(x = "PC1 on NMF[k=3] (150 C-genes)",
       y = "PC2 on NMF[k=3] (150 C-genes)",
       col = "GlioVis Majority subtype",
       fill = "NMF/PCA/LDA subtype")
  



## Perform PCA on NMF coordinates ? ----

## Export ----



