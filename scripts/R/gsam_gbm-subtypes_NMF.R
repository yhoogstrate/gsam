#!/usr/bin/env R

# load libs ----

library(tidyverse)
library(NMF)
library(scales)

#library(MASS)
library(fitdistrplus)
library(patchwork)
library(ggplot2)


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
  dplyr::select(c('sid', 'gliovis.majority_call', 'pat.with.IDH', 'resection')) %>%
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



## + NMF:1-3 using Wang-code ----


if(!file.exists("tmp/gsam_nmf_150.Rds")) {
  gsam_nmf_150 <- {{}}
  #gsam_nmf_15k <- {{}}
  #gsam_nmf_7k <- {{}}
  
  seeds <- c(123456 , 12345, 1337, 909) # 123456 is the one used by their paper
  
  for(seed in seeds) {
    gsam_nmf_150[[as.character(seed)]] <- NMF(as.matrix(gsam.rnaseq.expression.vst.150), 3, seed = seed)
  }
  
  gsam_nmf_150[["nmf"]] <- nmf(as.data.frame(gsam.rnaseq.expression.vst.150), 3, seed=12345)
  
  saveRDS(gsam_nmf_150, "tmp/gsam_nmf_150.Rds")
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




## + PCA:1-2 over NMF:1-3 ----
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



## + LDA classification ----


# build classifier

s150.pca.nmf.subtype.classifier.lda <- MASS::lda(`gliovis.majority_call` ~ `NMF:123456.PC1` + `NMF:123456.PC2` , data = train.in <- plt.single %>%
                              dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2', 'gliovis.majority_call')) %>%
                              tibble::column_to_rownames('sid') %>%
                              as.data.frame())


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
### 1 gamma CDS on eucledian ----


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



# ---- 2 non-parameteric & perimeter ----


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
                       n.eucledian.transition = sum(`NMF:123456.PCA.LDA.class.R1`!=`NMF:123456.PCA.LDA.class.R2`), .groups='drop'
      ) %>%
      dplyr::ungroup() %>%
      as.data.frame() %>%
      dplyr::mutate(n.eucledian.sum = n.eucledian.stable + n.eucledian.transition) %>%
      dplyr::mutate(p.eucledian.stable = n.eucledian.stable / n.eucledian.sum) %>%
      dplyr::mutate(p.eucledian.transition = 1 - p.eucledian.stable)
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


#plt.nmf.150.p <-
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
ggsave("output/figures/paper_subtypes_nmf_nmf_1_2_gliovis.png",width=7,height=7.2)




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

#ggsave("output/figures/paper_subtypes_nmf_pca_1_2_gliovis.png",width=7,height=7.2)



## contour ----

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


# nmf.pca.lda.countours <- data.frame(PC1 = plt.single$`NMF:123456.PC1`,
#                                   PC2 = plt.single$`NMF:123456.PC2`,
#                                   class = metadata$gliovis.majority_call,
#                                   type = "GlioVis") %>%
#   rbind( data.frame(PC1 = range_df$`NMF:123456.PC1`,
#                     PC2 = range_df$`NMF:123456.PC2`,
#                     class = range_predictions$class,
#                     type = "LDA on NVM + PCA"))



## PCA:1+2(NMF:1+2+3) + LDA countours + GlioVis labels ----

plt <- rbind(plt.single %>%
    dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'gliovis.majority_call')) %>%
    dplyr::mutate(type = "Patient Sample") %>%
    dplyr::rename(class = gliovis.majority_call),
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
ggsave("output/figures/paper_subtypes_nmf_pca_1_2_lda-contour.png",width=7,height=7.2)



## PCA:1+2(NMF:1+2+3) + LDA countours + LDA labels ----

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




