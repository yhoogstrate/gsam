#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(NMF)
library(scales)

#library(MASS)
library(fitdistrplus)
library(patchwork)
library(ggplot2)
library(circlize)
library(ggrepel)
library(e1071)

library(splancs) # https://stackoverflow.com/questions/17571602/r-filter-coordinates
library(combinat)


# load data ----


source('scripts/R/gsam_rna-seq_expression.R') # recursively calls metadata

source('scripts/R/glass_expression_matrix.R')


source('data/wang/msig.library.12.R')  # no license w/ code provided, can't include it in source tree

source('scripts/R/subtype_genes.R')
source('scripts/R/wang_glioma_intrinsic_genes.R')

source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')


source('scripts/R/palette.R')


# Per sample stats ----

## GLASS: vst.150 ----


glass.gbm.rnaseq.expression.vst.150 <- glass.gbm.rnaseq.expression.vst %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in%  (wang.glioma.intrinsic.genes %>%
                             dplyr::filter(Subtyping_Signature_Gene. != "") %>%
                             dplyr::pull(ENSG.short) %>%
                             gsub("ENSG00000276644","ENSG00000165659",.,fixed=T)) # DACH1
  ) %>%
  dplyr::mutate(ensg = NULL) %>%
  dplyr::mutate(gid = gsub("ENSG00000276644","DACH1",gid) ) %>% 
  dplyr::mutate(gid = gsub("ENSG00000165659","DACH1",gid) ) %>%
  tibble::column_to_rownames('gid')


stopifnot(min(glass.gbm.rnaseq.expression.vst.150) > 0)
stopifnot(nrow(glass.gbm.rnaseq.expression.vst.150) == 150)


## GLASS: metadata ----


glass.metadata <- glass.gbm.rnaseq.metadata %>%
  dplyr::slice(match(colnames(glass.gbm.rnaseq.expression.vst.150), sid))


stopifnot(colnames(glass.gbm.rnaseq.expression.vst.150) == glass.metadata$sid)
stopifnot(ncol(glass.gbm.rnaseq.expression.vst.150) == nrow(glass.metadata))
stopifnot(colnames(glass.gbm.rnaseq.expression.vst.150) %in% glass.metadata$sid)
stopifnot(glass.metadata$sid %in% colnames(glass.gbm.rnaseq.expression.vst.150))
stopifnot(colnames(glass.gbm.rnaseq.expression.vst.150) == glass.metadata$sid)



## G-SAM: vst.150 ----


gsam.rnaseq.expression.vst.150 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  `colnames<-`(paste0("GSAM-",colnames(.) )   ) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(ensg = gsub("\\..+$","",gid)) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(ensg %in%
                  (wang.glioma.intrinsic.genes %>%
                     dplyr::filter(Subtyping_Signature_Gene. != "") %>%
                     dplyr::pull(ENSG.short))
  ) %>%
  dplyr::mutate(ensg = gsub("ENSG00000276644","DACH1",ensg) ) %>% 
  dplyr::mutate(ensg = gsub("ENSG00000165659","DACH1",ensg) ) %>%
  tibble::column_to_rownames('ensg') 


stopifnot(min(gsam.rnaseq.expression.vst.150) > 0)
stopifnot(nrow(gsam.rnaseq.expression.vst.150) == 150)




## G-SAM: metadata ----


gsam.metadata <- gsam.rna.metadata %>%
  dplyr::mutate(sid = paste0("GSAM-",sid)) %>%
  dplyr::slice(match(colnames(gsam.rnaseq.expression.vst.150), sid))


stopifnot(ncol(gsam.rnaseq.expression.vst.150) > 0)
stopifnot(ncol(gsam.rnaseq.expression.vst.150) == nrow(gsam.metadata))
stopifnot(colnames(gsam.rnaseq.expression.vst.150) %in% gsam.metadata$sid)
stopifnot(gsam.metadata$sid %in% colnames(gsam.rnaseq.expression.vst.150))


gsam.rnaseq.expression.vst.150 <- gsam.rnaseq.expression.vst.150 %>%
  dplyr::select(gsam.metadata$sid)


stopifnot(colnames(gsam.rnaseq.expression.vst.150) == gsam.metadata$sid)


## Combine GLASS & G-SAM ----


#sva::ComBat(as.factor(gsub("^(..).*$","\\1",colnames(.))))


combi.rnaseq.expression.vst.150 <- dplyr::full_join(
  glass.gbm.rnaseq.expression.vst.150 %>%
    tibble::rownames_to_column('gid')
  ,
  gsam.rnaseq.expression.vst.150 %>%
    tibble::rownames_to_column('gid') ,
  by = c('gid'='gid') ) %>%
  dplyr::mutate(gid = NULL) %>%
  limma::removeBatchEffect(as.factor(gsub("^(..).*$","\\1",colnames(.))))



combi.metadata <- rbind(
  glass.metadata %>%
    dplyr::mutate(
      resection = case_when( resection == "TP" ~ "R1",
                             resection == "R1" ~ "R2",
                             resection == "R2" ~ "R3",
                             resection == "R3" ~ "R4" )) %>%
    dplyr::rename(subtype.public = GBM.transcriptional.subtype.Synapse) ,
  gsam.metadata %>%
    dplyr::select(c("sid","pid","resection","gliovis.majority_call")) %>%
    dplyr::rename(subtype.public = gliovis.majority_call ) %>%
    dplyr::mutate(sid.label = gsub("GSAM-","",sid) ) %>%
    dplyr::mutate(sid) %>%
    dplyr::mutate(resection = gsub("r","R",resection)) %>%
    dplyr::mutate(dataset = "GSAM") ) %>%
  dplyr::mutate(resection = as.factor(resection)) %>%
  dplyr::mutate(dataset = as.factor(dataset))



colnames(combi.rnaseq.expression.vst.150)
combi.metadata$sid

stopifnot(colnames(combi.rnaseq.expression.vst.150) == combi.metadata$sid)
stopifnot(sum(is.na(combi.rnaseq.expression.vst.150)) == 0)
stopifnot(nrow(combi.rnaseq.expression.vst.150) == 150)
stopifnot(ncol(combi.rnaseq.expression.vst.150) > 350)



plt.single <- combi.metadata


## Add regular PCA ---- 


p <- prcomp(t(combi.rnaseq.expression.vst.150 ))#screeplot(p)

p$x <- p$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join( combi.metadata , by = c('sid' = 'sid') )




ggplot(p$x, aes(x=PC1, y=PC2, col=dataset , label=sid.label ) ) +
  geom_point() +
  youri_gg_theme

ggplot(p$x, aes(x=PC1, y=PC2, col=subtype.public , label=sid.label ) ) +
  geom_point() +
  youri_gg_theme




## + NMF:1-3 using Wang-code ----


#if(!file.exists("tmp/combi.gbm_nmf_150.Rds") & F) {
  combi.gbm_nmf_150 <- {{}}
  seeds <- c(123456 ) #, 12345, 1337, 909) # 123456 is the one used by their paper
  
  for(seed in seeds) {
    combi.gbm_nmf_150[[as.character(seed)]] <- NMF(as.matrix(combi.rnaseq.expression.vst.150), 3, seed = seed)
  }
  
  #glass.gbm_nmf_150[["nmf"]] <- nmf(as.data.frame(glass.gbm.rnaseq.expression.vst.150), 3, seed=12345)
  
  # huidige mds is baas
  #saveRDS(glass.gbm_nmf_150, "tmp/combi.gbm_nmf_150.Rds")
#} else {
  #glass.gbm_nmf_150 <- readRDS("tmp/glass.gbm_nmf_150.Rds")
#}



plt.single <- plt.single %>%
  dplyr::left_join(
    combi.gbm_nmf_150 %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:123456.1','NMF:123456.2','NMF:123456.3')) %>%
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(`NMF:123456.membership` = as.factor (combi.gbm_nmf_150 %>% purrr::pluck('123456') %>% purrr::pluck('membership') %>% gsub('^(.+)$','NMF-cluster:\\1',.) ) )
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


#screeplot(p)


rm(p)



ggplot(plt.single, aes( x=`NMF:123456.PC1` , y = `NMF:123456.PC2` , col=`subtype.public` )) + 
  geom_point(size=2) +
  youri_gg_theme


# ## + LDA classification 
# # build classifier
# 
# s150.pca.nmf.subtype.classifier.lda <- MASS::lda(`subtype.public` ~ `NMF:123456.PC1` + `NMF:123456.PC2` , data = plt.single %>%
#                                                    dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2', 'subtype.public')) %>%
#                                                    tibble::column_to_rownames('sid') %>%
#                                                    as.data.frame())
# 
# 
# # re-fit against new classifier
# 
# 
# plt.single <- plt.single %>% dplyr::left_join(
#   predict(s150.pca.nmf.subtype.classifier.lda, newdata = plt.single %>%
#             dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2', 'subtype.public')) %>%
#             tibble::column_to_rownames('sid') %>%
#             as.data.frame() ) %>%
#     data.frame() %>%
#     dplyr::select(c('class', 'posterior.Classical', 'posterior.Mesenchymal', 'posterior.Proneural')) %>%
#     `colnames<-`( gsub('^(.+)$','NMF:123456.PCA.LDA.\\1', colnames(.))  ) %>%
#     tibble::rownames_to_column('sid')
#   , by=c('sid'='sid') )
# 


## + SVM classification ----

# Aanzienlijk lekkerdere fit dan LDA

s150.pca.nmf.subtype.classifier.svm <- svm(x = plt.single  %>%
                                             dplyr::select(c('sid' ,'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
                                             tibble::column_to_rownames('sid') ,
                                           y = plt.single$subtype.public, scale = F, tolerance = 0.0001, type = "C-classification", kernel = "linear", cost = 3, probability = T)



# re-fit samples

#attr(train.pred, 'decision.values')
plt.single <- plt.single %>% dplyr::left_join(
  data.frame(
    sid = plt.single$sid,
    
    'NMF:123456.PCA.SVM.class' = as.character(
      predict(object = s150.pca.nmf.subtype.classifier.svm , newdata =  (plt.single  %>%
                                                dplyr::select(c('sid' ,'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
                                                tibble::column_to_rownames('sid')) , decision.values = T, probability = T)),
    check.names = F
    ) , by = c('sid'='sid') )



#plt.single$`NMF:123456.PCA.SVM.class`



# ~~ Export stats ~~ ----

write.table(
  plt.single %>%
    dplyr::select(c("sid",
                    "NMF:123456.membership",
                    "NMF:123456.PC1", "NMF:123456.PC2", "NMF:123456.PC3",
                    "NMF:123456.PCA.SVM.class"
                    #"NMF:123456.PCA.LDA.class", "NMF:123456.PCA.LDA.posterior.Classical", "NMF:123456.PCA.LDA.posterior.Mesenchymal", "NMF:123456.PCA.LDA.posterior.Proneural"
    )),
  "output/tables/gsam_nmf_lda_data.txt")




# Per pair stats ----
## Split table into pairs + eucl dist ----

plt.paired <- plt.single %>%
  #dplyr::filter(resection %in% c("R1", "R2") ) %>%
  dplyr::group_by(pid) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::top_n(1, sid) %>%
  dplyr::select(all_of('pid')) %>%
  as.data.frame() %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R1') %>% `colnames<-`( paste0(colnames(.) , ".R1") ) , by=c('pid'='pid.R1') ) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection %in% c('R2', 'R3', 'R4') ) %>% `colnames<-`( paste0(colnames(.) , ".R2") ) , by=c('pid'='pid.R2') ) %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                        (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2)) %>%
  dplyr::mutate(subtype.public.status = as.factor(ifelse(subtype.public.R1 == subtype.public.R2, "Stable", "Transition"))) %>%
  #dplyr::mutate(`NMF:123456.membership.status` = as.factor(ifelse(`NMF:123456.membership.R1`== `NMF:123456.membership.R2`, "Stable", "Transition"))) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.factor(ifelse(`NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`, "Stable", "Transition")))



plt.paired %>%   dplyr::filter(pid %in% c('G-SM-R056-2', 'GLSS-HF-3081', 'GLSS-HF-2869')  )


plot(density(plt.paired$eucledian.dist))
hist(plt.paired$eucledian.dist,breaks=15)




# Plots ----


## Determine Contour ----


resolution <- 250 # 1000 x 1000 data points

off_x <- (max(plt.single$`NMF:123456.PC1`) - min(plt.single$`NMF:123456.PC1`)) * 0.05
off_y <- (max(plt.single$`NMF:123456.PC2`) - min(plt.single$`NMF:123456.PC2`)) * 0.05


range_pc1 = seq(from = min(plt.single$`NMF:123456.PC1`) - off_x, to = max(plt.single$`NMF:123456.PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(plt.single$`NMF:123456.PC2`) - off_y, to = max(plt.single$`NMF:123456.PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:123456.PC1' = range_pc1, 'NMF:123456.PC2' = range_pc2)
#nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.lda, newdata = range_df) %>%
nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = range_df) %>% data.frame(class = .) %>%
  #as.data.frame() %>%
  cbind(range_df) %>%
  dplyr::select(c('class', 'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
  dplyr::mutate(type="Contour")

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)



## PCA:1+2(NMF:1+2+3) + LDA contours + GlioVis labels ----


plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'subtype.public', sid.label, dataset)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = subtype.public),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA))



plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM")
       , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample" & dataset == "GSAM"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



plt.glass <- ggplot(plt %>% dplyr::filter(dataset != "GSAM")
       , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample" & dataset != "GSAM"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_GlioVis_LDA-countours.png',width=7,height=7.2)




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
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call] / GLASS consortium',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



#ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours.png',width=7,height=7.2)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [classical] ----


plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.classical = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::arrange(pid, resection)
  # %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Classical" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.classical == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Classical" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.classical == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass






ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_classical.png',width=15,height=6.5)





## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [proneural] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.proneural = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::arrange(pid, resection)
  # %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_proneural.png',width=15,height=6.5)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [mesenchymal] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::arrange(pid, resection)
  # %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_mesenchymal.png',width=15,height=6.5)



# KNN Bootstrapping ~ estimate by chance ----

k = 20
df = data.frame()
for(p in 1:nrow(plt.paired) ) {
  #p = 21

  # for pair
  target = plt.paired[p,]
  
  
  # find k closest R1's
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                               (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(k, -ed)
  
  
  #plot(plt.single$`NMF:123456.PC1` , plt.single$`NMF:123456.PC2`, pch=19 , col = as.numeric(plt.single $`NMF:123456.PCA.SVM.class`) + 1 , cex=0.5)
  
  
  #nodes <- data.frame(`NMF:123456.PC1` = target$`NMF:123456.PC1.R1`,
  #                    `NMF:123456.PC2` = target$`NMF:123456.PC2.R1`,
  #                    type="start", check.names=F)
  nodes <- data.frame()
  for(i in 1:k) {
    neighbour <- knn[i,] 
    
    #lines(c(neighbour$`NMF:123456.PC1.R1`,neighbour$`NMF:123456.PC1.R2`), c(neighbour$`NMF:123456.PC2.R1` , neighbour$`NMF:123456.PC2.R2`)  )
    #points(neighbour$`NMF:123456.PC1.R2` , neighbour$`NMF:123456.PC2.R2` , pch=8,cex=0.6  )
    
    delta_PC1 = target$`NMF:123456.PC1.R1` - neighbour$`NMF:123456.PC1.R1`
    delta_PC2 = target$`NMF:123456.PC2.R1` - neighbour$`NMF:123456.PC2.R1`
    
    
    nodes <- rbind(nodes,
              data.frame(`NMF:123456.PC1` = neighbour$`NMF:123456.PC1.R2` + delta_PC1 ,
                         `NMF:123456.PC2` = neighbour$`NMF:123456.PC2.R2` + delta_PC2 ,
                         type = "node", check.names=F))


  }

  nodes <- nodes %>%
    dplyr::mutate(class.svm = predict(s150.pca.nmf.subtype.classifier.svm , newdata = nodes %>% dplyr::mutate(type=NULL)) )
  
  df <- rbind(df, 
              data.frame(n.cl  = nodes %>% dplyr::filter(class.svm == "Classical") %>% nrow(),
                   n.mes = nodes %>% dplyr::filter(class.svm == "Mesenchymal") %>% nrow(),
                   n.pn  = nodes %>% dplyr::filter(class.svm == "Proneural") %>% nrow() ) %>%
    dplyr::mutate(p.cl = n.cl /  rowSums(.),
                  p.mes = n.mes /  rowSums(.),
                  p.pn = n.pn /  rowSums(.) ) %>%
    dplyr::mutate(pid = target$pid)
  )
  
}


df <- df %>%
  #group_by(pid) %>% 
  #dplyr::summarise( p.cl = mean(n.cl),     p.mes = mean(n.mes), p.pn = mean(n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))



df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()



#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()


df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()



#df.2 %>% filter( dataset.R1 == 'GSAM' & `NMF:123456.PCA.SVM.status`  == "Stable") %>% nrow()



# plot per-group eucledian distance between pairs?


# KNN Bootstrapping 2 ~ estimate by chance ----

# Bij deze bootstrapping zoeken we van de KNN de richtingen/angles
# En gebruiken bootstrappen we de eucledian distances van alle samples


# hoe berekenen we van 2 nodes de angle en verlengen/verkorten we deze
# R1 = ( 1,  1)
# R2 = (-2, -3)

change_length_line <- function(x1, y1, x2, y2, length_new) {
  # From the  line between points (x1, y1) , (x2 ,y2),
  # we want to create a new line with an identical angle
  # but of length `length_new`.
  
  dy <- y2 - y1
  dx <- x2 - x1
  
  #slope <- dy / dx
  angle <- atan2(dy , dx) # in rads
  
  length_x_new <- cos(angle) * length_new
  length_y_new <- sin(angle) * length_new
  
  x2_new <- x1 + length_x_new
  y2_new <- y1 + length_y_new
 
  return (data.frame(x = c(x1, x2_new) ,
             y = c(y1, y2_new) ,
             point = as.factor(c("initial start", "new end"))))
}


# plot(c(-4,4), c(-4,4), type="n")
# 
# points(1,1, pch=19)
# points(-2,-3, pch=19)
# lines(c(-2,1), c(-3,1))
# 
# df <- change_length_line(1,1,-2,-3, 1)
# points(df[2,1] , df[2,2], pch=19,col="red")




nn = 15
bootstrap_n = 25 # bootstrap iterations per sample-pair
df = data.frame()
for(p in 1:nrow(plt.paired) ) { # for pair
  #p = 21 
  target = plt.paired[p,]
  
  
  # find k closest R1's for angles in close proximity
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                                 (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(nn, -ed)
  

  #plot(plt.single$`NMF:123456.PC1` , plt.single$`NMF:123456.PC2`, pch=19 , col = as.numeric(plt.single $`NMF:123456.PCA.SVM.class`) + 1 , cex=0.5)
  #lines(c(target$`NMF:123456.PC1.R1`,target$`NMF:123456.PC1.R2`), c(target$`NMF:123456.PC2.R1` , target$`NMF:123456.PC2.R2`) , lwd=1.5 )
  #points(target$`NMF:123456.PC1.R1` , target$`NMF:123456.PC2.R1` , pch=4,cex=0.6   )
  #points(target$`NMF:123456.PC1.R2` , target$`NMF:123456.PC2.R2` , pch=4,cex=0.6   )
  
  
  for(i in 1:bootstrap_n) {
    
    # Take the angle from local subsampling
    neighbour <- knn %>% # RANDOM SHUFFLE `MET TERUGLEGGEN` want bootstrappen!!!
      dplyr::slice(sample(1 : dplyr::n() ) [1] )

    # take length from overall subsamples
    random_length <- plt.paired %>%
      dplyr::slice(sample(1 : dplyr::n() )) %>%
      dplyr::pull(eucledian.dist) %>%
      purrr::pluck(1)

    # calc new line length
    bootstrapped_line <- change_length_line(neighbour$`NMF:123456.PC1.n.R1`,
                                            neighbour$`NMF:123456.PC2.n.R1`,
                                            
                                            neighbour$`NMF:123456.PC1.n.R2`,
                                            neighbour$`NMF:123456.PC2.n.R2`,
                                            
                                            random_length)

    
    # fit to target's R1
    bootstrapped_line$x = bootstrapped_line$x + (target$`NMF:123456.PC1.n.R1` - neighbour$`NMF:123456.PC1.n.R1`)
    bootstrapped_line$y = bootstrapped_line$y + (target$`NMF:123456.PC2.n.R1` - neighbour$`NMF:123456.PC2.n.R1`)
    

    class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = bootstrapped_line %>%
              dplyr::filter(point == "new end") %>%
              dplyr::mutate(point = NULL) %>%
              dplyr::rename(`NMF:123456.PC1` = x) %>%
              dplyr::rename(`NMF:123456.PC2` = y) )
    
    df <- rbind(df, data.frame(x = bootstrapped_line$x[2],
               y = bootstrapped_line$y[2],
               pid = target$pid,
               neighbour = neighbour$pid,
               class = class ))
    
    
    # plot(plt.single$`NMF:123456.PC1.n` , plt.single$`NMF:123456.PC2.n`, pch=19 , col = as.numeric(plt.single $`NMF:123456.PCA.SVM.class`) + 1 , cex=0.5)
    # 
    # lines(
    #   c(target$`NMF:123456.PC1.n.R1` , target$`NMF:123456.PC1.n.R2`),
    #   c(target$`NMF:123456.PC2.n.R1` , target$`NMF:123456.PC2.n.R2`) , lwd=1.5 )
    # 
    # lines(c(neighbour$`NMF:123456.PC1.n.R1` , neighbour$`NMF:123456.PC1.n.R2`),
    #       c(neighbour$`NMF:123456.PC2.n.R1` , neighbour$`NMF:123456.PC2.n.R2`) , lwd=2  )
    # 
    # lines(c(bootstrapped_line[1,1], bootstrapped_line[2,1]),
    #       c(bootstrapped_line[1,2], bootstrapped_line[2,2]) , lwd=2  )
  }
}




df <- df %>%
  dplyr::group_by(pid) %>% 
  dplyr::summarise(n.cl  = sum(class == "Classical"),
                   n.mes = sum(class == "Mesenchymal"),
                   n.pn  = sum(class == "Proneural")) %>%
  dplyr::mutate(p.cl  = n.cl  / (n.cl + n.mes + n.pn) ,
                p.mes = n.mes / (n.cl + n.mes + n.pn) ,
                p.pn  = n.pn  / (n.cl + n.mes + n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))



df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()





#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()


df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()








# 〰 © Dr. Youri Hoogstrate 〰 ----













