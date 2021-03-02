#!/usr/bin/env R

# ---- load libs ----

library(tidyverse)
library(pheatmap)
library(NMF)
library(corrplot)

# ---- load cfg ----

#source('scripts/R/job_gg_theme.R')
#source('scripts/R/youri_gg_theme.R')

# ---- Load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/expression_matrix_VST.R')

# ---- load NMF ----


expression.data <- expression_matrix_full_new.vst %>%
  dplyr::filter(rownames(.) %in% c(subtype.mesenchymal$id, subtype.proneural$id, subtype.classical$id)) %>%
  t() %>%
  scale() %>%
  - min(df) %>%
  `colnames<-`(  gsub("^.+\\|([^\\]+)\\|.+$","\\1",colnames(.) )   )  %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(  gsub(".","-",colnames(.) , fixed=T ) )


metadata <- data.frame(sid = rownames(df)) %>%
  dplyr::left_join(gsam.rna.metadata , by = c('sid' = 'sid' ) )



# NMF moet goed gescaled zijn, net als pheatmap i./c.m. dist() ...
tmp <- expression.data %>%
  t() %>%
  scale() %>%
   - min(.)


set.seed(12345) # same seed as Verhaak / Wang
#nmf.fit.2 <- nmf(df, 2 , seed = 12345)
nmf.fit.3 <- nmf(df, 3 , seed = 12345)
nmf.fit.4 <- nmf(df, 4 , seed = 12345)
nmf.fit.5 <- nmf(df, 5 , seed = 12345)
nmf.fit.6 <- nmf(df, 6 , seed = 12345)




# ---- pheatmap of 150 vst genes ----

pheatmap(t(expression.data), scale="column", 
         ar = metadata %>%
                dplyr::select('sid', 'gliovis.majority_call','gliovis.knn_call','gliovis.svm_call') %>%
                tibble::column_to_rownames('sid'))



#, annotation_row = metadata %>% dplyr::select('sid', 'gliovis.majority_call','gliovis.knn_call','gliovis.svm_call') %>%
#           tibble::column_to_rownames('sid') )




coefmap(minfit(fit), info = TRUE, annCol = ar, hclustfun='complete', distfun='correlation')


basismap(fit, info = TRUE)
consensusmap(fit)

#basiscor(fit)
#profcor(fit)


fit <- fit.3
p <- data.frame(sid = colnames(coef(fit)),
                nmf.1 = coef(fit)[1,],
                nmf.2 = coef(fit)[2,]
                #,nmf.3 = fit.3@fit@H[3,]
                ) %>%
  dplyr::left_join(metadata , by=c('sid' = 'sid'))

ggplot(p, aes(x = nmf.1 , y = nmf.2 , col= gliovis.majority_call) ) + 
  geom_point()




# basismap(res, main="Metagenes", annRow=list(d, e), tracks=c(Metagene=':basis'))
