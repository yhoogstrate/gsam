#!/usr/bin/env R

# -------------------------------------------------------------- #
# THIS FILE IS OLD AND OBSOLETE
# 
# Please use: 'scripts/combined_gbm-subtypes_NMF.R' instead
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# -------------------------------------------------------------- #


# load libs ----


library(tidyverse)
library(NMF)
library(scales)

#library(MASS)
library(e1071)
library(fitdistrplus)
library(patchwork)
library(ggplot2)
library(circlize)
library(ggrepel)


# load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/gsam_rna-seq_expression.R')

source('data/wang/msig.library.12.R') # license is incomplete regarding sharing, can't include it in source tree

source('scripts/R/subtype_genes.R')
source('scripts/R/wang_glioma_intrinsic_genes.R')

source('scripts/R/youri_gg_theme.R')


source('scripts/R/palette.R')


# Per sample stats ----

## Normalize for NMF ----

tmp <- gsam.rnaseq.expression.vst

if(min(tmp) < 0) { # not the case
  tmp <- tmp - min(tmp) + .Machine$double.eps
}


## S-150 genes subset ----


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




plt.single <- metadata %>%
  dplyr::select(c('sid', 'gliovis.majority_call', 'pat.with.IDH', 'resection', 'gliovis.knn_call', 'gliovis.svm_call', 'gliovis.gsea_call','Gravendeel.Centroid.Class.Full','Gravendeel.Centroid.Class.Subset')) %>%
  dplyr::mutate(resection = gsub('r','R',resection)) %>%
  dplyr::mutate(pid = gsub('^(...).*$','\\1',sid))



## Exclude IDH muts [after normalisation, so they can be mapped back into the coordinates later] ----
## not really possible with NMF

## Add regular PCA ---- 


p <- prcomp(t(gsam.rnaseq.expression.vst.150))#screeplot(p)

plt.single <- plt.single %>%
  dplyr::left_join(
    p %>%
      purrr::pluck('x') %>%
      as.data.frame() %>%
      tibble::rownames_to_column('sid') %>%
      dplyr::select(c('sid','PC1','PC2','PC3','PC4'))
    , by=c('sid' = 'sid') )

rm(p)



## NMF:1-3 using Wang-code ----


if(!file.exists("tmp/gsam_nmf_150.Rds")) {
  gsam_nmf_150 <- {{}}

  seeds <- c(123456 , 12345, 1337, 909) # 123456 is the one used by their paper
  
  for(seed in seeds) {
    gsam_nmf_150[[as.character(seed)]] <- NMF(as.matrix(gsam.rnaseq.expression.vst.150), 3, seed = seed)
  }
  
  gsam_nmf_150[["nmf"]] <- nmf(as.data.frame(gsam.rnaseq.expression.vst.150), 3, seed=12345)
  
  #saveRDS(gsam_nmf_150, "tmp/gsam_nmf_150.Rds")
} else {
  gsam_nmf_150 <- readRDS("tmp/gsam_nmf_150.Rds")
}




plt.single <- plt.single %>%
  dplyr::left_join(
    gsam_nmf_150 %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:123456.1','NMF:123456.2','NMF:123456.3')) %>%
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(`NMF:123456.membership` = as.factor (gsam_nmf_150 %>% purrr::pluck('123456') %>% purrr::pluck('membership') %>% gsub('^(.+)$','NMF-cluster:\\1',.) ) )
    , by=c('sid' = 'sid') )




## PCA:1-2 over NMF:1-3 ----
# clearly the nicest fit (!!!!!) compared w/ primary pca & NMF:1+2
# PCA to also fit NMF 3 into two dimensions (is a better fit)



p <- plt.single %>%
  dplyr::select(c('sid',"NMF:123456.1","NMF:123456.2","NMF:123456.3")) %>%
  tibble::column_to_rownames('sid') %>%
  #t() %>%
  prcomp()


plt.single <- plt.single %>%
  dplyr::left_join(
    p %>%
      purrr::pluck('x') %>%
      as.data.frame() %>%
      `colnames<-`(gsub('^PC','NMF:123456.PC',colnames(.))) %>%
      tibble::rownames_to_column('sid') %>%
      dplyr::select(c('sid','NMF:123456.PC1','NMF:123456.PC2','NMF:123456.PC3'))
    , by=c('sid' = 'sid') ) %>%
  dplyr::mutate(`NMF:123456.PC1.n` = as.numeric(scale(`NMF:123456.PC1`))) %>% # scaled PC's, to be used for Eucledian distances
  dplyr::mutate(`NMF:123456.PC2.n` = as.numeric(scale(`NMF:123456.PC2`))) %>%
  dplyr::mutate(`NMF:123456.PC3.n` = as.numeric(scale(`NMF:123456.PC3`)))

#png("output/figures/tmp-screeplot.png")
#screeplot(p)
#dev.off()


rm(p)



## LDA (re-)classification ----

# # SVM > LDA
# build classifier

s150.pca.nmf.subtype.classifier.lda <-
  MASS::lda(
    `gliovis.majority_call` ~ `NMF:123456.PC1` + `NMF:123456.PC2` ,
    data = train.in <- plt.single %>%
      dplyr::select(
        c(
          'sid',
          'NMF:123456.PC1',
          'NMF:123456.PC2',
          'gliovis.majority_call'
        )
      ) %>%
      tibble::column_to_rownames('sid') %>%
      as.data.frame()
  )
# re-fit against new classifier


plt.single <- plt.single %>% dplyr::left_join(
  predict(s150.pca.nmf.subtype.classifier.lda, newdata = plt.single %>%
            dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2', 'gliovis.majority_call')) %>%
            tibble::column_to_rownames('sid') %>%
            as.data.frame() ) %>%
    data.frame() %>%
    dplyr::select(c('class', 'posterior.Classical', 'posterior.Mesenchymal', 'posterior.Proneural')) %>%
    `colnames<-`( gsub('^(.+)$','NMF:123456.PCA.LDA.\\1', colnames(.))  ) %>%
    tibble::rownames_to_column('sid')
  , by=c('sid'='sid') )




### calc LDA Contour ----

resolution <- 250 # 1000 x 1000 data points

off_x <- (max(plt.single$`NMF:123456.PC1`) - min(plt.single$`NMF:123456.PC1`)) * 0.05
off_y <- (max(plt.single$`NMF:123456.PC2`) - min(plt.single$`NMF:123456.PC2`)) * 0.05


range_pc1 = seq(from = min(plt.single$`NMF:123456.PC1`) - off_x, to = max(plt.single$`NMF:123456.PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(plt.single$`NMF:123456.PC2`) - off_y, to = max(plt.single$`NMF:123456.PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:123456.PC1' = range_pc1, 'NMF:123456.PC2' = range_pc2)
nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.lda, newdata = range_df) %>%
  as.data.frame() %>%
  cbind(range_df) %>%
  dplyr::select(c('class', 'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
  dplyr::mutate(type="Contour")

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)



## SVM (re-classification) ----

s150.pca.nmf.subtype.classifier.svm <- svm(
  x =  plt.single %>%
    dplyr::select(c(
      'sid',
      'NMF:123456.PC1',
      'NMF:123456.PC2'
    )) %>%
    tibble::column_to_rownames('sid') %>%
    as.matrix() ,
  y = as.factor(plt.single$gliovis.majority_call),
  kernel = 'linear',
  tolerance = 0.0001,
  cost = 3
  )


plt.single <- plt.single %>% dplyr::left_join(
    predict(s150.pca.nmf.subtype.classifier.svm, newdata = plt.single %>% dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2')) %>% tibble::column_to_rownames('sid') %>% as.data.frame() ) %>%
    as.data.frame() %>%
    `colnames<-`( c('NMF:123456.PCA.SVM.class' )) %>%
    tibble::rownames_to_column('sid')
  , by=c('sid'='sid'))

sum(plt.single$gliovis.majority_call != plt.single$`NMF:123456.PCA.SVM.class`)


### calc SVM contour ----


resolution <- 250 # 1000 x 1000 data points

off_x <- (max(plt.single$`NMF:123456.PC1`) - min(plt.single$`NMF:123456.PC1`)) * 0.05
off_y <- (max(plt.single$`NMF:123456.PC2`) - min(plt.single$`NMF:123456.PC2`)) * 0.05


range_pc1 = seq(from = min(plt.single$`NMF:123456.PC1`) - off_x, to = max(plt.single$`NMF:123456.PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(plt.single$`NMF:123456.PC2`) - off_y, to = max(plt.single$`NMF:123456.PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:123456.PC1' = range_pc1, 'NMF:123456.PC2' = range_pc2)
nmf.pca.svm.countours <- predict(s150.pca.nmf.subtype.classifier.svm, newdata = range_df) %>%
  as.data.frame() %>%
  `colnames<-`(c('class')) %>%
  cbind(range_df) %>%
  dplyr::select(c('class', 'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
  dplyr::mutate(type="Contour")

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)





# Per pair stats ----
## Split table into pairs + eucl dist ----

plt.paired <- plt.single %>%
  dplyr::group_by(pid) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::top_n(1, sid) %>%
  dplyr::select(all_of('pid')) %>%
  as.data.frame() %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R1') %>% `colnames<-`( paste0(colnames(.) , ".R1") ) , by=c('pid'='pid.R1') ) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R2') %>% `colnames<-`( paste0(colnames(.) , ".R2") ) , by=c('pid'='pid.R2') ) %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                      (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2)) %>%
  dplyr::mutate(gliovis.status = as.factor(ifelse(gliovis.majority_call.R1 == gliovis.majority_call.R2, "Stable", "Transition"))) %>%
  dplyr::mutate(`NMF:123456.membership.status` = as.factor(ifelse(`NMF:123456.membership.R1`== `NMF:123456.membership.R2`, "Stable", "Transition"))) %>%
  dplyr::mutate(`NMF:123456.PCA.LDA.status` = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition")))
  

  


plot(density(plt.paired$eucledian.dist))
hist(plt.paired$eucledian.dist,breaks=15)


## table w/ non-matching tumour pairs (control) ----


plt.anti.paired <- expand.grid(
  plt.single %>%
    dplyr::filter(resection == 'R1') %>%
    dplyr::pull(pid)
  ,
  plt.single %>%
    dplyr::filter(resection == 'R2') %>%
    dplyr::pull(pid)
  ) %>%
  dplyr::rename(pid.R1 = Var1) %>%
  dplyr::rename(pid.R2 = Var2) %>%
  dplyr::mutate(pid.R1 = as.character(pid.R1)) %>%
  dplyr::mutate(pid.R2 = as.character(pid.R2)) %>%
  dplyr::filter(`pid.R1` != `pid.R2`) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R1') %>% `colnames<-`( paste0(colnames(.) , ".R1") ) , by=c('pid.R1'='pid.R1') ) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R2') %>% `colnames<-`( paste0(colnames(.) , ".R2") ) , by=c('pid.R2'='pid.R2') ) %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                      (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2)) %>%
  #dplyr::mutate(gliovis.status = as.factor(ifelse(gliovis.majority_call.R1 == gliovis.majority_call.R2, "Stable", "Transition"))) %>%
  #dplyr::mutate(`NMF:123456.membership.status` = as.factor(ifelse(`NMF:123456.membership.R1`== `NMF:123456.membership.R2`, "Stable", "Transition"))) %>%
  dplyr::mutate(`NMF:123456.PCA.LDA.status` = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition")))




plot(density(plt.anti.paired %>% dplyr::pull(eucledian.dist)))
# classical heeft kleinere distances, by chance (!)
plot(density(plt.anti.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Classical') %>% dplyr::pull(eucledian.dist)))
plot(density(plt.anti.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal') %>% dplyr::pull(eucledian.dist)))
plot(density(plt.anti.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Proneural') %>% dplyr::pull(eucledian.dist)))



wilcox.test(plt.anti.paired$eucledian.dist, plt.paired$eucledian.dist)
median(plt.paired$eucledian.dist)
median(plt.anti.paired$eucledian.dist)
median(plt.anti.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Classical') %>% dplyr::pull(eucledian.dist))
median(plt.anti.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal') %>% dplyr::pull(eucledian.dist))
median(plt.anti.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Proneural') %>% dplyr::pull(eucledian.dist))


## Fit gamma distr on eucl distances ----

library(fitdistrplus)
gamma.paired.fit <- fitdist(plt.paired$eucledian.dist, "gamma")
gamma.anti.paired.fit <- fitdist(plt.anti.paired$eucledian.dist, "gamma")


plot(gamma.paired.fit)
median(plt.paired$eucledian.dist)

pgamma(c( median(plt.paired$eucledian.dist) ), shape = gamma.paired.fit$estimate[1], rate = gamma.paired.fit$estimate[2])
qgamma(c( 0.5098686 ), shape = gamma.paired.fit$estimate[1], rate = gamma.paired.fit$estimate[2])




plot(gamma.paired.fit)
plot(gamma.anti.paired.fit)




## Determine transitions by chance ----
### + 1 gamma CDS on eucledian ----


#plt.paired <- plt.paired %>%
#dplyr::mutate(p.gamma.paired = 1 - pgamma(c( eucledian.dist ), shape = gamma.paired.fit$estimate[1], rate = gamma.paired.fit$estimate[2]) )%>%
#dplyr::mutate(p.gamma.anti.paired = 1 - pgamma(c( eucledian.dist ), shape = gamma.anti.paired.fit$estimate[1], rate = gamma.anti.paired.fit$estimate[2]) )
# 
# 
# plt.anti.paired <- plt.anti.paired %>% 
#   dplyr::mutate(p.gamma.paired = 1 - pgamma(c( eucledian.dist ), shape = gamma.paired.fit$estimate[1], rate = gamma.paired.fit$estimate[2]) )
# 
# 
# 
# 
# 
# plt.paired <- plt.paired %>%
#   dplyr::left_join(
#     
#     plt.anti.paired %>%
#       dplyr::filter(sid.R1 %in% plt.paired$sid.R1) %>%
#       dplyr::select(c('sid.R1', 'sid.R2',"eucledian.dist", "p.gamma.paired",
#                       'NMF:123456.PCA.LDA.class.R1','NMF:123456.PCA.LDA.class.R2')) %>%
#       dplyr::group_by(sid.R1 , `NMF:123456.PCA.LDA.class.R2`) %>%
#       dplyr::summarise(sum.p.gamma.paired = sum(p.gamma.paired), .groups='drop' ) %>%
#       dplyr::ungroup() %>% 
#       tidyr::pivot_wider(names_from = `NMF:123456.PCA.LDA.class.R2`, values_from = `sum.p.gamma.paired`) %>%
#       dplyr::mutate(tmp.sum = Classical + Mesenchymal + Proneural,
#                     p.gamma.to.cl = Classical / tmp.sum,
#                     p.gamma.to.mes = Mesenchymal / tmp.sum,
#                     p.gamma.to.pr = Proneural / tmp.sum,
#                     tmp.sum = NULL, Classical = NULL, Mesenchymal = NULL,  Proneural = NULL)
#     
#     ,by=c('sid.R1'='sid.R1')) %>%
#   dplyr::mutate(
#     p.gamma.stable = case_when(
#       `NMF:123456.PCA.LDA.class.R1` == 'Classical' ~ p.gamma.to.cl,
#       `NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal' ~ p.gamma.to.mes,
#       `NMF:123456.PCA.LDA.class.R1` == 'Proneural' ~ p.gamma.to.pr)
#   ) %>%
#   dplyr::mutate(p.gamma.transition = 1 - p.gamma.stable)
# 
# 
# 
# plt.paired %>% dplyr::pull(p.gamma.transition) %>% median()
# plt.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == "Classical")  %>% dplyr::pull(p.gamma.transition) %>% median()
# plt.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == "Proneural")  %>% dplyr::pull(p.gamma.transition) %>% median()
# plt.paired %>% dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == "Mesenchymal")  %>% dplyr::pull(p.gamma.transition) %>% median()
# 
# 
# # check in plt.anti.paired:
# # mean + median eucledian distance per patient in mesenchymal & see if this is close to median in the fit
# 
# plt.anti.paired %>%
#   dplyr::group_by(sid.R1) %>%
#   dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>%
#   dplyr::pull(eucledian.dist) %>% median()
# 
# 
# pgamma(c( 1.8 ), shape = gamma.paired.fit$estimate[1], rate = gamma.paired.fit$estimate[2])

# 
# plt.paired <- plt.paired %>% dplyr::left_join(
#              rbind(
#               plt.paired %>%
#                 dplyr::select(-c('pid') ),
#               plt.anti.paired %>%
#                 dplyr::select(-c("pid.R1", "pid.R2", "p.gamma.paired"))) %>%
#               dplyr::left_join(
#                 plt.paired %>%
#                   dplyr::select(c('sid.R1', 'eucledian.dist')) %>%
#                   dplyr::rename(eucledian.dist.max = eucledian.dist)
#                 , by=c('sid.R1'='sid.R1')
#               ) %>%
#               dplyr::arrange('sid.R1') %>%
#               dplyr::select(c('sid.R1', 'sid.R2', 'NMF:123456.PCA.LDA.class.R1', 'NMF:123456.PCA.LDA.class.R2', 'eucledian.dist', 'eucledian.dist.max')) %>%
#               dplyr::filter(eucledian.dist <= eucledian.dist.max) %>%
#               dplyr::group_by(sid.R1) %>%
#               dplyr::summarise(n.eucledian.stable = sum(`NMF:123456.PCA.LDA.class.R1`==`NMF:123456.PCA.LDA.class.R2`),
#                                n.eucledian.transition = sum(`NMF:123456.PCA.LDA.class.R1`!=`NMF:123456.PCA.LDA.class.R2`), .groups='drop'
#                                ) %>%
#               dplyr::ungroup() %>%
#               as.data.frame() %>%
#               dplyr::mutate(n.eucledian.sum = n.eucledian.stable + n.eucledian.transition) %>%
#               dplyr::mutate(p.eucledian.stable = n.eucledian.stable / n.eucledian.sum) %>%
#               dplyr::mutate(p.eucledian.transition = 1 - p.eucledian.stable)
#              ,by=c('sid.R1'='sid.R1')
#             )
# 
# 
# 
# p.stable <- plt.paired$p.eucledian.stable
# events <- plt.paired$`NMF:123456.PCA.LDA.class.R1` == plt.paired$`NMF:123456.PCA.LDA.class.R2`
# 
# observed = c(stable=sum(events), trans=sum(events == F))
# expected = c(stable=sum(p.stable) , (sum(1 - p.stable) ))
# 
# chisq = sum((expected - observed)^2 / expected)
# #chi = sqrt(chisq)
# 
# # kans dat een transitie plaats vind
# expected / sum(expected) * 100
# 
# # kans dat het aantal gevonden transities plaatsvind op basis van toeval
# pchisq(chisq, df=1, lower.tail=FALSE)



# ---- + 2 non-parameteric & perimeter ----


plt.paired <- plt.paired %>% 
  dplyr::left_join(
    rbind(
      plt.paired %>%
        dplyr::select(-c('pid','gliovis.status','NMF:123456.membership.status','NMF:123456.PCA.LDA.status') ),
      plt.anti.paired %>%
        dplyr::select(-c("pid.R1", "pid.R2",'NMF:123456.PCA.LDA.status'))) %>%
      dplyr::left_join(
        plt.paired %>%
          dplyr::select(c('sid.R1', 'eucledian.dist')) %>%
          dplyr::rename(eucledian.dist.max = eucledian.dist)
        , by=c('sid.R1'='sid.R1')
      ) %>%
      dplyr::arrange('sid.R1') %>%
      dplyr::select(c('sid.R1', 'sid.R2', 'NMF:123456.PCA.LDA.class.R1', 'NMF:123456.PCA.LDA.class.R2', 'eucledian.dist', 'eucledian.dist.max')) %>%
      dplyr::filter(eucledian.dist <= eucledian.dist.max) %>%
      dplyr::group_by(sid.R1) %>%
      dplyr::summarise(n.eucledian.stable = sum(`NMF:123456.PCA.LDA.class.R1`==`NMF:123456.PCA.LDA.class.R2`),
                       n.eucledian.transition = sum(`NMF:123456.PCA.LDA.class.R1`!=`NMF:123456.PCA.LDA.class.R2`),
                       n.eucledian.cl = sum(`NMF:123456.PCA.LDA.class.R2` == 'Classical'),
                       n.eucledian.mes = sum(`NMF:123456.PCA.LDA.class.R2` == 'Mesenchymal'),
                       n.eucledian.pn = sum(`NMF:123456.PCA.LDA.class.R2` == 'Proneural'),
                       .groups='drop'
      ) %>%
      dplyr::ungroup() %>%
      as.data.frame() %>%
      dplyr::mutate(n.eucledian.sum = n.eucledian.stable + n.eucledian.transition) %>%
      dplyr::mutate(p.eucledian.stable = n.eucledian.stable / n.eucledian.sum) %>%
      dplyr::mutate(p.eucledian.transition = 1 - p.eucledian.stable) %>%
      dplyr::mutate(p.eucledian.cl = n.eucledian.cl / n.eucledian.sum) %>%
      dplyr::mutate(p.eucledian.mes = n.eucledian.mes / n.eucledian.sum) %>%
      dplyr::mutate(p.eucledian.pn = n.eucledian.pn / n.eucledian.sum)
    ,by=c('sid.R1'='sid.R1'))




# a = cc[1:7]
# b = cs[1:7]
# probs <- c()
# for(k in 0:3) {
#   
#   com = t(combn(length(a), k))
#   for(i in 1:nrow(com)) {
#     combi <- com[i,]
#     
#     p <- a / 100
#     p[combi] <- 1 - p[combi]
#     
#     probs <- c(probs, prod(p))
#   }
# 
# }
# print(sum(probs))



# Perform tests ----

## find overall transition rate(s)) ----


plt.paired %>%
  dplyr::pull(gliovis.status) %>%
  summary
# Stable Transition 
# 99         61


plt.paired %>%
  dplyr::pull('NMF:123456.PCA.LDA.status') %>%
  summary()
# Stable Transition 
# 103         57


plt.paired %>%
  dplyr::pull('NMF:123456.membership.status') %>%
  summary()
# Stable Transition 
# 92         68 




events <- plt.paired %>% dplyr::pull('NMF:123456.PCA.LDA.status')
observed = c(stable=sum(events == "Stable"), trans=sum(events == "Transition"))
expected = c(stable=plt.paired$p.eucledian.stable %>% sum() , trans=(1 - plt.paired$p.eucledian.stable) %>% sum()  )

chisq = sum((expected - observed)^2 / expected) #chi = sqrt(chisq)
# kans dat een transitie plaats vind #expected / sum(expected) * 100
pchisq(chisq, df=1, lower.tail=FALSE)# kans dat het aantal gevonden transities plaatsvind op basis van toeval


observed
expected
(observed[2] / (observed[1] + observed[2])) / 
  (expected[2] / (expected[1] + expected[2]))




## find classical transition rate(s)) ----

plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Classical') %>%
  dplyr::pull('NMF:123456.PCA.LDA.status') %>%
  summary()

events <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Classical') %>%
  dplyr::pull('NMF:123456.PCA.LDA.status')


observed = c(stable=sum(events == "Stable"), trans=sum(events == "Transition"))
expected = c(stable=sum(  (plt.paired %>%
                             dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Classical') %>%
                             dplyr::pull(p.eucledian.stable)) ) ,
             trans=sum( 1 - (plt.paired %>%
                               dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Classical') %>%
                               dplyr::pull(p.eucledian.stable)) )   )



chisq = sum((expected - observed)^2 / expected)
#expected / sum(expected) * 100 # kans dat een transitie plaats vind
pchisq(chisq, df=1, lower.tail=FALSE)# kans dat het aantal gevonden transities plaatsvind op basis van toeval

observed
expected
(observed[2] / (observed[1] + observed[2])) / 
  (expected[2] / (expected[1] + expected[2]))



## find Mesenchymal transition rate(s)) ----

plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal') %>%
  dplyr::mutate(status.gliovis = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition"))) %>%
  dplyr::pull(status.gliovis) %>%
  summary()

events <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal') %>%
  dplyr::mutate(status.gliovis = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition"))) %>%
  dplyr::pull(status.gliovis)

observed = c(stable=sum(events == "Stable"), trans=sum(events == "Transition"))
expected = c(stable=sum(  (plt.paired %>%
                             dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal') %>%
                             dplyr::pull(p.eucledian.stable)) ) ,
             trans=sum( 1 - (plt.paired %>%
                               dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Mesenchymal') %>%
                               dplyr::pull(p.eucledian.stable)) )   )

chisq = sum((expected - observed)^2 / expected)


#expected / sum(expected) * 100 # kans dat een transitie plaats vind
pchisq(chisq, df=1, lower.tail=FALSE)# kans dat het aantal gevonden transities plaatsvind op basis van toeval


observed
expected
(observed[2] / (observed[1] + observed[2])) / 
  (expected[2] / (expected[1] + expected[2]))




## find Proneural transition rate(s)) ----


plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Proneural') %>%
  dplyr::mutate(status.gliovis = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition"))) %>%
  dplyr::pull(status.gliovis) %>%
  summary()

events <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Proneural') %>%
  dplyr::mutate(status.gliovis = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition"))) %>%
  dplyr::pull(status.gliovis)

observed = c(stable=sum(events == "Stable"), trans=sum(events == "Transition"))
expected = c(stable=sum(  (plt.paired %>%
                             dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Proneural') %>%
                             dplyr::pull(p.eucledian.stable)) ) ,
             trans=sum( 1 - (plt.paired %>%
                               dplyr::filter(`NMF:123456.PCA.LDA.class.R1` == 'Proneural') %>%
                               dplyr::pull(p.eucledian.stable)) )   )

chisq = sum((expected - observed)^2 / expected)


#expected / sum(expected) * 100 # kans dat een transitie plaats vind
pchisq(chisq, df=1, lower.tail=FALSE)# kans dat het aantal gevonden transities plaatsvind op basis van toeval


observed
expected
(observed[2] / (observed[1] + observed[2])) / 
  (expected[2] / (expected[1] + expected[2]))


## eucledian distances diff ----


plt <- plt.paired %>%
  dplyr::select(c('pid','eucledian.dist','NMF:123456.PCA.LDA.status','NMF:123456.PCA.LDA.class.R1'))


ggplot(plt, aes(x = `NMF:123456.PCA.LDA.status` , y = eucledian.dist , col =  `NMF:123456.PCA.LDA.class.R1`  )) +
  geom_point() +
  facet_grid(. ~ `NMF:123456.PCA.LDA.class.R1` ) +
  youri_gg_theme



# predict random transition and see if this is indeed closer ----


rnd <- plt.anti.paired[sample(1:nrow(plt.anti.paired) ),] %>%
  dplyr::left_join(
    plt.paired %>%
      dplyr::select(c('sid.R1', 'eucledian.dist')) %>%
      dplyr::rename(eucledian.dist.max = eucledian.dist)
    , by=c('sid.R1'='sid.R1')
  ) %>%
  dplyr::filter(eucledian.dist <= eucledian.dist.max) %>%
  dplyr::group_by(sid.R1) %>%
  filter(row_number()==1) %>%
  as.data.frame() %>%
  dplyr::left_join(
    rbind(
      plt.paired %>%
        dplyr::select(-c('pid') ),
      plt.anti.paired %>%
        dplyr::select(-c("pid.R1", "pid.R2", "p.gamma.paired"))) %>%
      dplyr::left_join(
        plt.paired %>%
          dplyr::select(c('sid.R1', 'eucledian.dist')) %>%
          dplyr::rename(eucledian.dist.max = eucledian.dist)
        , by=c('sid.R1'='sid.R1')
      ) %>%
      dplyr::arrange('sid.R1') %>%
      dplyr::select(c('sid.R1', 'sid.R2', 'NMF:123456.PCA.LDA.class.R1', 'NMF:123456.PCA.LDA.class.R2', 'eucledian.dist', 'eucledian.dist.max')) %>%
      dplyr::filter(eucledian.dist <= eucledian.dist.max) %>%
      dplyr::group_by(sid.R1) %>%
      dplyr::summarise(n.eucledian.stable = sum(`NMF:123456.PCA.LDA.class.R1`==`NMF:123456.PCA.LDA.class.R2`),
                       n.eucledian.transition = sum(`NMF:123456.PCA.LDA.class.R1`!=`NMF:123456.PCA.LDA.class.R2`), .groups='drop'
      ) %>%
      dplyr::ungroup() %>%
      as.data.frame() %>%
      dplyr::mutate(n.eucledian.sum = n.eucledian.stable + n.eucledian.transition) %>%
      dplyr::mutate(p.eucledian.stable = n.eucledian.stable / n.eucledian.sum) %>%
      dplyr::mutate(p.eucledian.transition = 1 - p.eucledian.stable)
    ,by=c('sid.R1'='sid.R1'))


rnd %>%
  dplyr::mutate(status.gliovis = as.factor(ifelse(`NMF:123456.PCA.LDA.class.R1` == `NMF:123456.PCA.LDA.class.R2`, "Stable", "Transition"))) %>%
  dplyr::pull(status.gliovis) %>%
  summary()

sum(rnd$p.eucledian.stable)
sum(1 - rnd$p.eucledian.stable)




# example test chi-sq

# p.stable <- c(0.9,0.5,0.7,0.8,0.8,0.2,0.9)
# events <- c(1,1,0,1,1,0,1)
# 
# observed = c(stable=5, trans=2)
# expected = c(stable=sum(p.stable) , (sum(1 - p.stable) ))
# 
# chisq = sum((expected - observed)^2 / expected)
# #chi = sqrt(chisq)
# 
# # kans dat een transitie plaats vind
# expected / sum(expected) * 100
# 
# # kans dat het aantal gevonden transities plaatsvind op basis van toeval
# pchisq(chisq, df=1, lower.tail=FALSE)
# 
# 
# # test combinations
# 
# kk = 0:2 # first all 
# nn = 10
# 
# for(k in kk) {
#   t(combn(nn, k))
# }




# Plots ----

## Original S-150 PCA space ----


ggplot(plt.single, aes(x = PC1, y = PC2, col = gliovis.majority_call) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  scale_color_manual(name = NULL, values =  c('Classical'='#6ba6e5',# blue
                                              'Mesenchymal'='#eab509',#mustard
                                              'Proneural'='#ff5f68'),#red/pink
                     guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


ggsave('output/figures/paper_subtypes_pca_S150G_PC1_PC2_GlioVis.png',width=7,height=7.2)



## NMF:1 & NMF:2 + NMF membership ----


ggplot(plt.single, aes(x = `NMF:123456.1`, y = `NMF:123456.2`, col = `NMF:123456.membership`) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  scale_color_manual(
    values = c('NMF-cluster:1'="#6ba6e5", # Classical
               'NMF-cluster:3'= "#eab509", # Mesenchymal
               'NMF-cluster:2'= "#ff5f68" # Proneural
    ),
    guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
  ) +
  scale_x_continuous(limits= c(min( plt.single$`NMF:123456.1`), max( plt.single$`NMF:123456.1`)) ,
                     expand = expansion(mult = c(.1, .1)) ) +
  scale_y_continuous(limits= c(min( plt.single$`NMF:123456.2`), max( plt.single$`NMF:123456.2`)) ,
                     expand = expansion(mult = c(.1, .1)) ) +
  labs(x = "NMF meta-feature 1", y="NMF meta-feature 2", col='Subtype by NMF meta-feature')


ggsave('output/figures/paper_subtypes_nmf_S150G_nmf1_nmf2_nmf-membership.png',width=7,height=7.2)




## NMF:1 & NMF:2 + GlioVis labels ----



ggplot(plt.single, aes(x = `NMF:123456.1`, y = `NMF:123456.2`, col = gliovis.majority_call) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  scale_color_manual(
    values = subtype_colors,
    guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
  ) +
  scale_x_continuous(expand = expansion(mult = c(.05, .05)) ) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05)) ) +
  labs(x = "NMF meta-feature 1", y="NMF meta-feature 2", col='Subtype by GlioVis [majority call]')


ggsave('output/figures/paper_subtypes_nmf_S150G_nmf1_nmf2_GlioVis.png',width=7,height=7.2)




## PCA:1+2(NMF:1+2+3) + GlioVis labels ----


ggplot(plt.single , aes(x = `NMF:123456.PC1` , y = `NMF:123456.PC2`, col = gliovis.majority_call) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  scale_color_manual(
    values = subtype_colors,
    guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
  ) +
  scale_x_continuous(expand = expansion(mult = c(.05, .05)) ) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05)) ) +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]')

ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_GlioVis.png',width=7,height=7.2)



## PCA:1+2(NMF:1+2+3) + GravenDeel labels ----



p1 <- ggplot(plt.single , aes(x = `NMF:123456.PC1` , y = `NMF:123456.PC2`, col = gliovis.majority_call) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  scale_x_continuous(expand = expansion(mult = c(.05, .05)) ) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05)) ) +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='GlioVis')

p2 <- ggplot(plt.single , aes(x = `NMF:123456.PC1` , y = `NMF:123456.PC2`, col = `NMF:123456.PCA.LDA.class`) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  scale_x_continuous(expand = expansion(mult = c(.05, .05)) ) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05)) ) +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='LDA-re-class')


p3 <- ggplot(plt.single , aes(x = `NMF:123456.PC1` , y = `NMF:123456.PC2`, col = Gravendeel.Centroid.Class.Full) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  #scale_color_manual(
  #values = subtype_colors,
  #  guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
  #) +
  scale_x_continuous(expand = expansion(mult = c(.05, .05)) ) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05)) ) +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Gravendeel centr. full')


p4 <- ggplot(plt.single , aes(x = `NMF:123456.PC1` , y = `NMF:123456.PC2`, col = Gravendeel.Centroid.Class.Subset) ) +
  geom_point(size=1.5) +
  youri_gg_theme +
  #scale_color_manual(
  #values = subtype_colors,
  #  guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
  #) +
  scale_x_continuous(expand = expansion(mult = c(.05, .05)) ) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05)) ) +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Gravendeel centr. sub')


(p1 + p2) / (p3 + p4)



## LDA status circlize [observed] ----


grid.col <- c(subtype_colors, subtype_colors, subtype_colors)
names(grid.col ) <- c(names(subtype_colors) , paste0('R1: ',names(subtype_colors)), paste0('R2: ',names(subtype_colors)))


circos.clear()
par(mar = rep(0, 4), cex=1)
circos.par(start.degree = 180)


pdf("output/figures/paper_subtypes_nmf_S150G_observed_transitions_circos.pdf")
#circos.initialize(sectors=sectors)
chordDiagram(x=table(paste0("R1: ",plt.paired$`NMF:123456.PCA.LDA.class.R1`),
                     paste0("R2: ",plt.paired$`NMF:123456.PCA.LDA.class.R2`)),
             grid.col = grid.col, 
             transparency = 0.3,
             big.gap = 30,
             small.gap = 3, 
             annotationTrackHeight = c(0.01, 0.08),
             annotationTrack = c("axis","grid")#"name",
             #annotationTrackHeight = convert_height(c(20, 60), "mm")
             #order = order,
)

#circos.track( track.index = 1, panel.fun = function(x, y) {
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  print(xlim)
  print(ylim)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim),si , sector.index = si, track.index =1 , 
              facing = "bending.inside", niceFacing = TRUE, col = "white",cex=0.8)
}

circos.clear()
dev.off()


## LDA status circlize [expected] ----



#observed <- table(paste0("R1: ",plt.paired$`NMF:123456.PCA.LDA.class.R1`),
#                  paste0("R2: ",plt.paired$`NMF:123456.PCA.LDA.class.R2`))

expected["R1: Classical","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.cl) %>% sum()
expected["R1: Classical","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.mes) %>% sum()
expected["R1: Classical","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.pn) %>% sum()

expected["R1: Mesenchymal","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.cl) %>% sum()
expected["R1: Mesenchymal","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.mes) %>% sum()
expected["R1: Mesenchymal","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.pn) %>% sum()

expected["R1: Proneural","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.cl) %>% sum() 
expected["R1: Proneural","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.mes) %>% sum() 
expected["R1: Proneural","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.pn) %>% sum() 



grid.col <- c(subtype_colors, subtype_colors, subtype_colors)
names(grid.col ) <- c(names(subtype_colors) , paste0('R1: ',names(subtype_colors)), paste0('R2: ',names(subtype_colors)))




circos.clear()
par(mar = rep(0, 4), cex=1)
circos.par(start.degree = 180)






pdf("output/figures/paper_subtypes_nmf_S150G_expected_transitions_circos.pdf")
#circos.initialize(sectors=sectors)
chordDiagram(x=expected,
             grid.col = grid.col, 
             transparency = 0.3,
             big.gap = 30,
             small.gap = 3, 
             annotationTrackHeight = c(0.01, 0.08),
             annotationTrack = c("axis","grid")#"name",
             #annotationTrackHeight = convert_height(c(20, 60), "mm")
             #order = order,
)

#circos.track( track.index = 1, panel.fun = function(x, y) {
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  print(xlim)
  print(ylim)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim),si , sector.index = si, track.index =1 , 
              facing = "bending.inside", niceFacing = TRUE, col = "white",cex=0.8)
}

circos.clear()
dev.off()




## LDA status circlize [expected] ----



#observed <- table(paste0("R1: ",plt.paired$`NMF:123456.PCA.LDA.class.R1`),
#                  paste0("R2: ",plt.paired$`NMF:123456.PCA.LDA.class.R2`))

expected["R1: Classical","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.cl) %>% sum()
expected["R1: Classical","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.mes) %>% sum()
expected["R1: Classical","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.pn) %>% sum()

expected["R1: Mesenchymal","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.cl) %>% sum()
expected["R1: Mesenchymal","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.mes) %>% sum()
expected["R1: Mesenchymal","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.pn) %>% sum()

expected["R1: Proneural","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.cl) %>% sum() 
expected["R1: Proneural","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.mes) %>% sum() 
expected["R1: Proneural","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.pn) %>% sum() 



grid.col <- c(subtype_colors, subtype_colors, subtype_colors)
names(grid.col ) <- c(names(subtype_colors) , paste0('R1: ',names(subtype_colors)), paste0('R2: ',names(subtype_colors)))




circos.clear()
par(mar = rep(0, 4), cex=1)
circos.par(start.degree = 180)






pdf("output/figures/paper_subtypes_nmf_S150G_expected_transitions_circos.pdf")
#circos.initialize(sectors=sectors)
chordDiagram(x=expected,
             grid.col = grid.col, 
             transparency = 0.3,
             big.gap = 30,
             small.gap = 3, 
             annotationTrackHeight = c(0.01, 0.08),
             annotationTrack = c("axis","grid")#"name",
             #annotationTrackHeight = convert_height(c(20, 60), "mm")
             #order = order,
)

#circos.track( track.index = 1, panel.fun = function(x, y) {
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  print(xlim)
  print(ylim)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim),si , sector.index = si, track.index =1 , 
              facing = "bending.inside", niceFacing = TRUE, col = "white",cex=0.8)
}

circos.clear()
dev.off()




## LDA status circlize [lfc observed( / expected)] ----



#observed <- table(paste0("R1: ",plt.paired$`NMF:123456.PCA.LDA.class.R1`),
#                  paste0("R2: ",plt.paired$`NMF:123456.PCA.LDA.class.R2`))

expected["R1: Classical","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.cl) %>% sum()
expected["R1: Classical","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.mes) %>% sum()
expected["R1: Classical","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Classical") %>% dplyr::pull(p.eucledian.pn) %>% sum()

expected["R1: Mesenchymal","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.cl) %>% sum()
expected["R1: Mesenchymal","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.mes) %>% sum()
expected["R1: Mesenchymal","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Mesenchymal") %>% dplyr::pull(p.eucledian.pn) %>% sum()

expected["R1: Proneural","R2: Classical"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.cl) %>% sum() 
expected["R1: Proneural","R2: Mesenchymal"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.mes) %>% sum() 
expected["R1: Proneural","R2: Proneural"] <- plt.paired %>% dplyr::filter( `NMF:123456.PCA.LDA.class.R1` == "Proneural") %>% dplyr::pull(p.eucledian.pn) %>% sum() 



grid.col <- c(subtype_colors, subtype_colors, subtype_colors)
names(grid.col ) <- c(names(subtype_colors) , paste0('R1: ',names(subtype_colors)), paste0('R2: ',names(subtype_colors)))




circos.clear()
par(mar = rep(0, 4), cex=1)
circos.par(start.degree = 180)

pdf("output/figures/paper_subtypes_nmf_S150G_observed_div_expected_transitions_circos.pdf")
#circos.initialize(sectors=sectors)
chordDiagram(x=observed / expected,
             grid.col = grid.col, 
             transparency = 0.3,
             big.gap = 30,
             small.gap = 3, 
             annotationTrackHeight = c(0.01, 0.08),
             annotationTrack = c("axis","grid")#"name",
             #annotationTrackHeight = convert_height(c(20, 60), "mm")
             #order = order,
)

#circos.track( track.index = 1, panel.fun = function(x, y) {
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  print(xlim)
  print(ylim)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim),si , sector.index = si, track.index =1 , 
              facing = "bending.inside", niceFacing = TRUE, col = "white",cex=0.8)
}

circos.clear()
dev.off()








## PCA:1+2(NMF:1+2+3) + LDA contours + GlioVis labels ----


#sum(plt.single$gliovis.majority_call != plt.single$`NMF:123456.PCA.SVM.class`)
#sum(plt.single$gliovis.majority_call != plt.single$`NMF:123456.PCA.LDA.class`)


plt <- rbind(plt.single %>%
    dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'gliovis.majority_call')) %>%
    dplyr::mutate(type = "Patient Sample") %>%
    dplyr::rename(class = gliovis.majority_call),
  nmf.pca.svm.countours)


ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)




ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_GlioVis_LDA-countours.png',width=7,height=7.2)




## PCA:1+2(NMF:1+2+3) + LDA contours + LDA labels ----

plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.LDA.class')) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.LDA.class`),
             nmf.pca.lda.countours)

ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)

ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours.png',width=7,height=7.2)


## PCA:1+2(NMF:1+2+3) + LDA contours + LDA labels + EGFR ampli ----


plt <- rbind(plt.single %>%
               dplyr::select(c('sid','pid','resection',`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.LDA.class')) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.LDA.class`),
             nmf.pca.lda.countours %>% 
               dplyr::mutate(pid=NA, sid=NA,resection=NA)
) %>%
  dplyr::left_join(
    gsam.patient.metadata %>%
      dplyr::select(c('studyID','EGFR')) %>%
      dplyr::filter(!is.na(EGFR)) ,
    by = c('pid'='studyID')
  )



#as.factor(gsam.patient.metadata$cnStatusEGFRs)


plt.egfr <- gsam.patient.metadata %>%
  dplyr::select(c('studyID','EGFR')) %>%
  dplyr::filter(!is.na(EGFR)) %>%
  dplyr::mutate(EGFR.R1 = case_when(
    EGFR == "Gained" ~ "Wildtype",
    EGFR == "Lost"  ~ "Mutated",
    EGFR == "Stable" ~ "Mutated",
    EGFR == "Wildtype" ~ "Wildtype" )) %>%
  dplyr::mutate(EGFR.R2 = case_when(
    EGFR == "Gained" ~ "Mutated",
    EGFR == "Lost"  ~ "Wildtype",
    EGFR == "Stable" ~ "Mutated",
    EGFR == "Wildtype" ~ "Wildtype" )) %>%
  dplyr::mutate(EGFR = NULL)



plt <- plt %>%
  dplyr::mutate(tmp = ifelse(resection == "R1", pid , NA) ) %>%
  dplyr::left_join( plt.egfr %>%
                      dplyr::select(c('studyID', 'EGFR.R1')) ,
                    by=c('tmp'='studyID')) %>%
  dplyr::mutate(tmp = ifelse(resection == "R2", pid , NA) ) %>%
  dplyr::left_join(  plt.egfr %>%
                       dplyr::select(c('studyID', 'EGFR.R2')) ,
                     by=c('tmp'='studyID')) %>%
  dplyr::mutate(tmp = NULL) %>%
  dplyr::mutate(EGFR = ifelse(is.na(EGFR.R1), EGFR.R2, EGFR.R1) ) %>%
  dplyr::mutate(EGFR.R1 = NULL, EGFR.R2 = NULL)





ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = EGFR)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclass_EGFRampl.png',width=7,height=7.2)



## PCA:1+2(NMF:1+2+3) + LDA contours + LDA labels + TP53 ----


plt <- rbind(plt.single %>%
               dplyr::select(c('sid','pid','resection',`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.LDA.class')) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.LDA.class`),
             nmf.pca.lda.countours %>% 
               dplyr::mutate(pid=NA, sid=NA,resection=NA)
) %>%
  dplyr::left_join(
    gsam.patient.metadata %>%
      dplyr::select(c('studyID','TP53')) %>%
      dplyr::filter(!is.na(TP53)) ,
    by = c('pid'='studyID')
  )


plt.TP53 <- gsam.patient.metadata %>%
  dplyr::select(c('studyID','TP53')) %>%
  dplyr::filter(!is.na(TP53)) %>%
  dplyr::mutate(TP53.R1 = case_when(
    TP53 == "Gained" ~ "Wildtype",
    TP53 == "Lost"  ~ "Mutated",
    TP53 == "Stable" ~ "Mutated",
    TP53 == "Wildtype" ~ "Wildtype" )) %>%
  dplyr::mutate(TP53.R2 = case_when(
    TP53 == "Gained" ~ "Mutated",
    TP53 == "Lost"  ~ "Wildtype",
    TP53 == "Stable" ~ "Mutated",
    TP53 == "Wildtype" ~ "Wildtype" )) %>%
  dplyr::mutate(TP53 = NULL)



plt <- plt %>%
  dplyr::mutate(tmp = ifelse(resection == "R1", pid , NA) ) %>%
  dplyr::left_join( plt.TP53 %>%
                      dplyr::select(c('studyID', 'TP53.R1')) ,
                    by=c('tmp'='studyID')) %>%
  dplyr::mutate(tmp = ifelse(resection == "R2", pid , NA) ) %>%
  dplyr::left_join(  plt.TP53 %>%
                       dplyr::select(c('studyID', 'TP53.R2')) ,
                     by=c('tmp'='studyID')) %>%
  dplyr::mutate(tmp = NULL) %>%
  dplyr::mutate(TP53 = ifelse(is.na(TP53.R1), TP53.R2, TP53.R1) ) %>%
  dplyr::mutate(TP53.R1 = NULL, TP53.R2 = NULL)

rm(plt.TP)




ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = TP53)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  #scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



## PCA:1+2(NMF:1+2+3) + LDA contours + LDA labels + NF1 ----


plt <- rbind(plt.single %>%
               dplyr::select(c('sid','pid','resection',`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.LDA.class')) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.LDA.class`),
             nmf.pca.lda.countours %>% 
               dplyr::mutate(pid=NA, sid=NA,resection=NA)
) %>%
  dplyr::left_join(
    gsam.patient.metadata %>%
      dplyr::select(c('studyID','NF1')) %>%
      dplyr::filter(!is.na(NF1)) ,
    by = c('pid'='studyID')
  )


plt.NF1 <- gsam.patient.metadata %>%
  dplyr::select(c('studyID','NF1','cnStatusNF1s')) %>%
  
  dplyr::mutate(NF1 = NULL ) %>%
  dplyr::rename(NF1 = cnStatusNF1s) %>%
  
  dplyr::filter(!is.na(NF1)) %>%
  dplyr::mutate(NF1.R1 = case_when(
    NF1 == "Gained" ~ "Wildtype",
    NF1 == "Lost"  ~ "Mutated",
    NF1 == "Stable" ~ "Mutated",
    NF1 == "Wildtype" ~ "Wildtype" )) %>%
  dplyr::mutate(NF1.R2 = case_when(
    NF1 == "Gained" ~ "Mutated",
    NF1 == "Lost"  ~ "Wildtype",
    NF1 == "Stable" ~ "Mutated",
    NF1 == "Wildtype" ~ "Wildtype" )) %>%
  dplyr::mutate(NF1 = NULL)



plt <- plt %>%
  dplyr::mutate(tmp = ifelse(resection == "R1", pid , NA) ) %>%
  dplyr::left_join( plt.NF1 %>%
                      dplyr::select(c('studyID', 'NF1.R1')) ,
                    by=c('tmp'='studyID')) %>%
  dplyr::mutate(tmp = ifelse(resection == "R2", pid , NA) ) %>%
  dplyr::left_join(  plt.NF1 %>%
                       dplyr::select(c('studyID', 'NF1.R2')) ,
                     by=c('tmp'='studyID')) %>%
  dplyr::mutate(tmp = NULL) %>%
  dplyr::mutate(NF1 = ifelse(is.na(NF1.R1), NF1.R2, NF1.R1) ) %>%
  dplyr::mutate(NF1.R1 = NULL, NF1.R2 = NULL)

rm(plt.NF1)




ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = NF1)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_fill_manual(values = subtype_colors)





## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [classical] ----


plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.LDA.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(primary.classical = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.LDA.class` == "Classical" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::arrange(resection, sid) %>%
  dplyr::mutate(pid = as.factor(pid))


#plt <- plt %>% dplyr::filter(pid == "AAT")


ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.LDA.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA, pid=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA, pid=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(`NMF:123456.PCA.LDA.class` != "Classical" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(`NMF:123456.PCA.LDA.class` == "Classical" & resection == "R1"), size=1.5) +
  geom_path(
    data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.classical == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))
    ,alpha = 0.1
  ) + 
  geom_text_repel(data = plt %>% dplyr::filter(resection == "R2" & `NMF:123456.PCA.LDA.status` != "Stable" & primary.classical == T &
                           `NMF:123456.PCA.LDA.posterior.Classical` < 0.25), size = 3)  +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="LDA class") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours_classical.png',width=7,height=5.5)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [proneural] ----

plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.LDA.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(pid = as.factor(pid)) %>%
  dplyr::mutate(primary.proneural = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.LDA.class` == "Proneural" & resection == "R1") %>% dplyr::pull(pid)))


ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.LDA.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  
  geom_point(data = plt %>% dplyr::filter(`NMF:123456.PCA.LDA.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(`NMF:123456.PCA.LDA.class` == "Proneural" & resection == "R1"), size=1.5) +
  geom_path(
    data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))
  ) + 
  geom_text_repel(data = plt %>% dplyr::filter(resection == "R2" & `NMF:123456.PCA.LDA.status` != "Stable" & primary.proneural == T), size = 3)  +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="LDA class") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours_proneural.png',width=7,height=5.5)


## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [mesenchymal] ----

plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.LDA.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(pid = as.factor(pid)) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.LDA.class` == "Mesenchymal" & resection == "R1") %>% dplyr::pull(pid)))



ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.LDA.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(`NMF:123456.PCA.LDA.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(`NMF:123456.PCA.LDA.class` == "Mesenchymal" & resection == "R1"), size=1.5) +
  geom_path(
    data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))
  ) + 
  geom_text_repel(data = plt %>% dplyr::filter(resection == "R2" & `NMF:123456.PCA.LDA.status` != "Stable" & primary.mesenchymal == T), size = 3)  +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="LDA class") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours_mesenchymal.png',width=7,height=5.5)


## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [hypermutated] ----

plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.LDA.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(pid = as.factor(pid)) %>%
  dplyr::mutate(hypermutated = pid %in% c('FAN', 'HAC', 'AQA', 'BAB', 'JAG', 'GAR', 'AAT', 'BAK', 'JAC', 'FAG', 'JAF', 'BAD', 'EAZ', 'CAO', 'FAL') )



ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.LDA.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA, pid=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.LDA.status' = NA, pid=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(hypermutated == F), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(hypermutated == T & resection == "R1"), size=1.5) +
  geom_path(
    data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(hypermutated == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))
  ) + 
  geom_text_repel(data = plt %>% dplyr::filter(resection == "R2" & hypermutated == T), size = 3)  +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="LDA class") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours_hypermutated.png',width=7,height=5.5)


## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [tumour-%] ----

plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.LDA.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(pid = as.factor(pid)) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::select(c('sid', 'tumour.percentage.dna')) ,
                   by=c('sid'='sid')) %>%
  dplyr::filter(!is.na(tumour.percentage.dna)) %>%
  dplyr::mutate(tpc = ifelse(tumour.percentage.dna <= 15 , "TPC <= 15", "TPC > 15") )

tmp <- plt %>%
  dplyr::filter(grepl("<", tpc)) %>% 
  dplyr::pull(pid) %>%
  as.character() %>%
  unique()

plt <- plt %>%
  dplyr::mutate(low.tpc.pair = pid %in% tmp)



ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.LDA.status` , label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate(tpc = NA, pid=NA, `NMF:123456.PCA.LDA.status`=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate(tpc = NA, pid=NA , `NMF:123456.PCA.LDA.status`=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(grepl(">",tpc)), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(grepl("<",tpc ) & resection == "R1" ), size=1.5) +
  geom_path(
    data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter( low.tpc.pair  ),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches"))
  ) + 
  geom_text_repel(data = plt %>% dplyr::filter(resection == "R2" & hypermutated == T), size = 3)  +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="LDA class") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours_hypermutated.png',width=7,height=5.5)




# Per sample predicted class bar / grid plot (geom_tile) ----

# 1 LDA reclassification
# 2 GlioVis majority
#  2.a KNN 
#  2.b SVM
#  2.c GSEA 
# 3 NMF clustering

#df <- expand.grid(x=x, y=x)


plt <- plt.single %>%
  dplyr::select(c('sid',
                  `NMF:123456.PCA.LDA.class`,
                  gliovis.majority_call,
                  gliovis.knn_call,
                  gliovis.svm_call,
                  gliovis.gsea_call,
                  `NMF:123456.membership` )) %>%
  dplyr::mutate(`NMF:123456.membership` = as.character(`NMF:123456.membership`) ) %>%
  dplyr::mutate(`NMF:123456.membership` = as.factor(case_when(
    `NMF:123456.membership` == "NMF-cluster:1" ~ "Classical",
    `NMF:123456.membership` == "NMF-cluster:2" ~ "Proneural",
    `NMF:123456.membership` == "NMF-cluster:3" ~ "Mesenchymal")))  %>%
  dplyr::mutate(order = match(1:nrow(.), order(`NMF:123456.PCA.LDA.class`,
                                               gliovis.majority_call,
                                               `NMF:123456.membership`,
                                               gliovis.svm_call,
                                               gliovis.gsea_call,
                                               gliovis.knn_call,
                                               sid))) %>%
  reshape2::melt(id=c("sid","order")) %>%
  dplyr::rename(classifier = variable) %>%
  dplyr::rename(predicted.class = value) %>%
  dplyr::mutate(classifier = 
                  factor(classifier, levels = c(
                    "NMF:123456.membership", "gliovis.gsea_call", "gliovis.svm_call", "gliovis.knn_call", "gliovis.majority_call", "NMF:123456.PCA.LDA.class"
                  ))) %>%
  dplyr::mutate(subclassifier = case_when(
    classifier == "NMF:123456.membership" ~ "no",
    classifier == "gliovis.majority_call" ~ "no",
    classifier == "NMF:123456.PCA.LDA.class" ~ "no",
    classifier == "gliovis.gsea_call" ~ "yes",
    classifier == "gliovis.svm_call" ~ "yes",
    classifier == "gliovis.knn_call" ~ "yes")) %>%
  dplyr::mutate(xmin=order,xmax=order+1) %>%
  dplyr::mutate(ymin = case_when(
      classifier == "NMF:123456.membership" ~ 1,
      classifier == "gliovis.gsea_call" ~ 1 + 3,
      classifier == "gliovis.svm_call" ~ 1 + 3 + 1,
      classifier == "gliovis.knn_call" ~  1 + 3 + 1 + 1,
      classifier == "gliovis.majority_call" ~ 1 + 3 + 1 + 1 + 1,
      classifier == "NMF:123456.PCA.LDA.class" ~ 1 + 3 + 1 + 1 + 1 + 3
    )) %>%
  dplyr::mutate(ymax = case_when(
      classifier == "NMF:123456.membership" ~ 1 + 3,
      classifier == "gliovis.gsea_call" ~ 1 + 3 + 1,
      classifier == "gliovis.svm_call" ~ 1 + 3 + 1 +1 ,
      classifier == "gliovis.knn_call" ~  1 + 3 + 1 + 1 + 1,
      classifier == "gliovis.majority_call" ~ 1 + 3 + 1 + 1 + 1 + 3,
      classifier == "NMF:123456.PCA.LDA.class" ~ 1 + 3 + 1 + 1 + 1 + 3 + 3
    ))
  

  


ggplot(plt , aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill = predicted.class, alpha=subclassifier) )  +
  geom_rect(colour="black",size=0.15) +
  scale_fill_manual(values = subtype_colors) + 
  scale_alpha_manual(values = c('no' = 0.8, 'yes' = 0.4)) +
  scale_y_continuous(breaks = c((1+4) / 2,(4 + 5) /2, (5 + 6) / 2, (6 + 7) / 2,(7 + 10) / 2,(10 + 13) / 2),
                     labels=c("NMF clustering","GlioVis: ssGSEA","GlioVis: SVM","GlioVis: KNN",
                              "GlioVis: marjoity call","NMF -> PCA -> LDA reclass")) + 
  scale_x_continuous(breaks =NULL ,expand=expansion(add=c(2,5)) ) + 
  theme(
    legend.position = 'bottom',
    panel.background = element_blank(),
    #axis.title.y = element_text(angle=90,vjust =2),
    #axis.title.x = element_text(vjust = -0.2),
    axis.text.x = element_blank(), #element_text(),
    #axis.line = element_line(colour="black"),
    panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank() )


ggsave('output/figures/paper_subtypes_classifier_tileplot.png',width=14 * 0.8,height=7.2 * 0.4)



# CNV diff CL->Mes ----


source('scripts/R/cnv_matrix.R')
source('scripts/R/chrom_sizes.R')

tmp <- c('FAN', 'HAC', 'AQA', 'BAB', 'JAG', 'GAR', 'AAT', 'BAK', 'JAC', 'FAG', 'JAF', 'BAD', 'EAZ', 'CAO', 'FAL')
plt <- data.frame(pid = tmp) %>%
  dplyr::left_join(gsam.cnv.metadata %>%
                     dplyr::filter(resection == "R1") %>%
                     dplyr::select(c('PD_ID', 'pid')) %>%
                     dplyr::rename(PD_ID.R1 = PD_ID)
                   , by=c('pid' = 'pid')) %>%
  dplyr::left_join(gsam.cnv.metadata %>%
                     dplyr::filter(resection == "R2") %>%
                     dplyr::select(c('PD_ID', 'pid')) %>%
                     dplyr::rename(PD_ID.R2 = PD_ID)
                   , by=c('pid' = 'pid')) %>%
  dplyr::mutate(R1.in.CNV = PD_ID.R1 %in% colnames(cnv_matrix)) %>%
  dplyr::mutate(R2.in.CNV = PD_ID.R2 %in% colnames(cnv_matrix)) %>%
  dplyr::filter(R1.in.CNV == T & R2.in.CNV == T) %>%
  dplyr::right_join(data.frame(pid = tmp) , by = c('pid'='pid') )


#a = cnv_matrix %>% dplyr::select(plt[1,]$PD_ID.R1)

patient = 1
plt.2 <- data.frame(r1 = (cnv_matrix %>% dplyr::pull(plt[patient,]$PD_ID.R1)) ,
                    r2 = (cnv_matrix %>% dplyr::pull(plt[patient,]$PD_ID.R2)),
                    region = rownames(cnv_matrix)) %>%
  dplyr::mutate(chr = as.factor(gsub(":.+$","",rownames(cnv_matrix)))) %>%
  dplyr::mutate(diff = r1 - r2)  %>%
  dplyr::filter(chr %in% c("chrX", "chrY") == F)



ggplot(plt.2 , aes(x = r1, y = r2 , col = chr , label = region) ) +
  geom_point(size=0.1) + 
  geom_text(data = subset(plt.2, abs(diff) > 3))


# TODO: put in separate script
data <- read.delim('data/gsam/output/tables/cnv_copynumber-ratio.cns_all.txt', stringsAsFactors = F) %>%
  dplyr::mutate(gene = NULL)


overlap_segments <- function(d1, d2) {
  df <- data.frame()
  chrs <- unique(c(d1$chromosome , d2$chromosome))
  
  for(chr in chrs) {
    starts <- sort(unique(c(
      d1 %>% dplyr::filter(chromosome == chr) %>% dplyr::pull(start) ,
      d2 %>% dplyr::filter(chromosome == chr) %>% dplyr::pull(start)  
    )))
    
    ends <- sort(unique(c(
      d1 %>% dplyr::filter(chromosome == chr) %>% dplyr::pull(end) ,
      d2 %>% dplyr::filter(chromosome == chr) %>% dplyr::pull(end)  
    ))) 

    print(chr)
    stopifnot(length(starts) == length(ends))
    
    v1 <- c()
    v2 <- c()
    
    for ( i in 1:length(starts)) {
      v1 <- c(v1, d1 %>%
                dplyr::filter( chromosome == chr ) %>%
                dplyr::filter( starts[i] >=  start &  ends[i] <= end) %>%
                dplyr::pull(log2))
      
      v2 <- c(v2, d2 %>%
                dplyr::filter( chromosome == chr ) %>%
                dplyr::filter( starts[i] >=  start &  ends[i] <= end) %>%
                dplyr::pull(log2))
    }
    
    df <- rbind(df, 
                data.frame(
                  start = starts,
                  end = ends, 
                  v1 = v1,
                  v2 = v2,
                  chr = paste0("chr",chr)) %>%
                  dplyr::mutate(dv = v2 - v1)
              )
  }
  return(df)
}


ddf <- data.frame()
for(i in 1:nrow(plt )) {
  
  #i = 6
  pat <- plt[i,]
  if(pat$pid %in% c("GAR") == F) {
    print(pat)
    
    d1 <- data %>% 
      dplyr::filter(sample.id == pat$PD_ID.R1 )
    
    d2 <- data %>%
      dplyr::filter(sample.id == pat$PD_ID.R2)
    
    
    tmp = overlap_segments(d1, d2) %>%
      dplyr::mutate(chr_offset = chrs_hg19_s[.$chr]) %>%
      dplyr::mutate(pid = pat$pid)
    
    ddf <- rbind(ddf, tmp)
  }
}


ddf <- ddf %>%
  dplyr::filter(chr %in% c("chrX" , "chrY") == F) %>%
  dplyr::mutate(x    = start + chr_offset,
                xend = end + chr_offset,
                y = dv,
                yend = dv)


#ggplot(subset(ddf, chr=="chr9"), aes(x = x , xend = xend, y=y, yend =yend, group = pid , col=pid ) ) + 
ggplot(subset(ddf, T), aes(x = x , xend = xend, y=y, yend =yend, group = pid , col=chr == "chr7" | chr == "chr9" ) ) + 
  geom_segment(lwd=1.3) +
  ylim(-2.5, 2.5) + 
  labs(y = "Copynumber DIFFERENCE R1 ~ R2")


plot(c(0 , max (ddf$end) ), c (min(c(ddf$v1, df$v2)) , max(  c(ddf$v1, ddf$v2) ) ) , type="n")



abline(h=0, lwd=0.5)
#for(i in 1:nrow(df)) {
#  lines(c(df$start[i] , df$end[i] ) , c(df$v1[i],df$v1[i]) , col="blue") 
#}
#for(i in 1:nrow(df)) {
#  lines(c(df$start[i] , df$end[i] ) , c(df$v2[i],df$v2[i]) , col="red") 
#}
for(i in 1:nrow(df)) {
  lines(c(df$start[i] , df$end[i] ) , c(df$v2[i] - df$v1[i],df$v2[i] - df$v1[i]) , col="darkgray", lwd=3) 
}



d1g <- GenomicRanges::makeGRangesFromDataFrame(
  d1 %>%
    dplyr::mutate(chromosome = paste0( "chr", chromosome)) %>%
    dplyr::mutate(end = end - 1) 
  , keep.extra.columns = T)

d2g <- GenomicRanges::makeGRangesFromDataFrame(
  d2 %>%
    dplyr::mutate(chromosome = paste0( "chr", chromosome)) %>%
    dplyr::mutate(end = end - 1) 
  , keep.extra.columns = T)






asd1 = subsetByOverlaps(d1g, d2g)
asd2 = setdiff(d1g, d2g)

# setdiff

# co-expr alg ----


df<-data.frame(site.x=c("A","A","A","B","B","C"),   
               site.y=c("B","C","D","C","D","D"),
               Distance=c(67,57,64,60,67,60))
n <- max(table(df$site.x)) + 1  # +1,  so we have diagonal of 
res <- lapply(with(df, split(Distance, df$site.x)), function(x) c(rep(NA, n - length(x)), x))
res <- do.call("rbind", res)
res <- rbind(res, rep(NA, n))
res <- as.dist(t(res))


df <- 
  data.frame(
  site.x = c(rep('a',7),
             rep('b',6),
             rep('c',5),
             rep('d',4),
             rep('EGFR',3),
             rep('f',2),
             rep('SOCS2',1),
             rep('h',0))
,
  site.y = c(
    c('b','c','d','EGFR','f','SOCS2','h'),
    c('c','d','EGFR','f','SOCS2','h'),
    c('d','EGFR','f','SOCS2','h'),
    c('EGFR','f','SOCS2','h'),
    c('f','SOCS2','h'),
    c('SOCS2','h'),
    c('h')
  ),
Distance = c(
  rep(2.566018937, 14),
  2.198769158,
  2.566018937,
  2.377492545,
  rep(2.566018937, 6),
  1,
  rep(2.566018937, 4) )
)


n <- max(table(df$site.x)) + 1  # +1,  so we have diagonal of 
res <- lapply(with(df, split(Distance, df$site.x)), function(x) c(rep(NA, n - length(x)), x))
res <- do.call("rbind", res)
res <- rbind(res, rep(NA, n))
res <- as.dist(t(res))





fit <- cmdscale(res, eig=TRUE, k=2)
plt <- fit$points %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(gid = ifelse(gid == "", "h", gid))

ggplot(plt, aes(x=V1, y=V2, label=gid)) +
  geom_point() +
  ggrepel::geom_text_repel()


# 〰 © Dr. Youri Hoogstrate 〰 ----

