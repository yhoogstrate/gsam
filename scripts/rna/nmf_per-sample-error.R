#!/usr/bin/env R

# load libs ----

library(tidyverse)
library(NMF)

# load cfg ----

#source('scripts/R/job_gg_theme.R')
#source('scripts/R/youri_gg_theme.R')

# Load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/gsam_expression_matrix_VST.R')

source('scripts/R/wang_glioma_intrinsic_genes.R')


# NMF ----
## prepare data ----

# NMF moet goed gescaled zijn, net als pheatmap i./c.m. dist() ...
tmp <- expression_matrix_full_new.vst %>%
  dplyr::filter(gsub("^([^\\.]+).+$","\\1",rownames(.))  %in% wang.glioma.intrinsic.genes$ENSG.short) %>%
  `rownames<-`( gsub("^([^\\.]+).+$","\\1",rownames(.))  ) %>%
  t() %>%
  scale() %>%
   - min(.)


## run NMF ----


set.seed(12345) # same seed as Verhaak / Wang

###nmf.fit.2 <- nmf(df, 2 , seed = 12345)

#nmf.fit.3 <- nmf(tmp, 3 , seed = 12345)
#saveRDS(nmf.fit.3 , file="tmp/nmf_per-sample-error_nmf.fit.3.Rds")

#nmf.fit.3 <- readRDS("tmp/nmf_per-sample-error_nmf.fit.3.Rds")

###nmf.fit.4 <- nmf(df, 4 , seed = 12345)
###nmf.fit.5 <- nmf(df, 5 , seed = 12345)
###nmf.fit.6 <- nmf(df, 6 , seed = 12345)


## reconstruct and calculate error ----

# https://stats.stackexchange.com/questions/354611/pattern-of-out-of-sample-reconstruction-error-in-nmf-cross-validation-why-is-it

#reconstructed <- nmf.fit.3@fit@W %*% nmf.fit.3@fit@H
#err <- reconstructed - tmp
#err.sq <- err * err
#frobenius.norm <- sum(err.sq) # I guess this should equal the residuals value from NMF?
#nmf_per.sample.error <- sqrt( rowSums(err.sq) / nrow(err.sq))
#barplot(sort(sample.rmse))
#saveRDS(nmf_per.sample.error , file="tmp/nmf_per-sample-error.Rds")


## cache ----

nmf_per.sample.error <- readRDS("tmp/nmf_per-sample-error.Rds")






