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

## Normalize for NMF ----

tmp <- gsam.rnaseq.expression.vst

if(min(tmp) < 0) { # not the case
  tmp <- tmp - min(tmp) + .Machine$double.eps
}

## make 150 set ----


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


# # 
# 
# 
# plt.pca.150.p <- ggplot(plt.single, aes(x = PC1, y = PC2, col = gliovis.majority_call) ) + 
#   geom_point(size=1.5) +
#   youri_gg_theme + 
#   scale_color_manual(name = NULL, values =  c('Classical'='#6ba6e5',# blue
#                                               'Mesenchymal'='#eab509',#mustard
#                                               'Proneural'='#ff5f68'),#red/pink
#                      guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))
# 

## NMF using Wang-code ----

# actual nmf

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
      dplyr::mutate(`NMF:123456.membership` = as.factor (gsam_nmf_150 %>% purrr::pluck('123456') %>% purrr::pluck('membership') %>% gsub('^(.+)$','NMF:\\1',.) ) )
    , by=c('sid' = 'sid') )




# 



# ggplot(plt.nmf.150, aes(x = NMF1, y = NMF2, col = NMF.mem) ) + 
#   geom_point(size=1.5) +
#   youri_gg_theme + 
#   scale_color_manual(
#     name = "Subtype by NMF Meta-feature:",
#     values = c('1'="#6ba6e5", # Classical
#                '3'= "#eab509", # Mesenchymal 
#                '2'= "#ff5f68" # Proneural
#     ),
#     guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
#   ) +
#   scale_x_continuous(limits= c(min(plt.nmf.150$NMF1), max(plt.nmf.150$NMF1)) , expand = expansion(mult = c(.1, .1)) ) +
#   scale_y_continuous(limits= c(min(plt.nmf.150$NMF2), max(plt.nmf.150$NMF2)) , expand = expansion(mult = c(.1, .1)) ) 
# ggsave("output/figures/paper_subtypes_nmf_nmf_1_2_nmf.png",width=7,height=7.2)
# 
# 
# #plt.nmf.150.p <-
# ggplot(plt.nmf.150, aes(x = NMF1, y = NMF2, col = gliovis) ) + 
#   geom_point(size=1.5) +
#   youri_gg_theme + 
#   scale_color_manual(
#     name = "Subtype by GlioVis:",
#     values = subtype_colors,
#     guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
#   ) +
#   scale_x_continuous(limits= c(min(plt.nmf.150$NMF1), max(plt.nmf.150$NMF1)) , expand = expansion(mult = c(.1, .1)) ) +
#   scale_y_continuous(limits= c(min(plt.nmf.150$NMF2), max(plt.nmf.150$NMF2)) , expand = expansion(mult = c(.1, .1)) ) 
# ggsave("output/figures/paper_subtypes_nmf_nmf_1_2_gliovis.png",width=7,height=7.2)
# 





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


p <- plt.single %>%
  dplyr::select(c('sid',"NMF:123456.1","NMF:123456.2","NMF:123456.3")) %>%
  tibble::column_to_rownames('sid') %>%
  #t() %>%
  prcomp()


#png("output/figures/tmp-screeplot.png")
#screeplot(p)
#dev.off()


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

rm(p)



# 
# 
# p <- prcomp(t(gsam_nmf_150$`123456`$H))
# plt.nmf.pca.150   <- data.frame(NMF.PC1 = p$x[,1], NMF.PC2 = p$x[,2], gliovis = as.factor(metadata$gliovis.majority_call))
# plt.nmf.pca.150.p <- 
# 
#   ggplot(plt.nmf.pca.150 , aes(x = NMF.PC1, y = NMF.PC2, col = gliovis) ) + 
#   geom_point(size=1.5) +
#   youri_gg_theme + 
#   scale_color_manual(
#     name = "Subtype by GlioVis:",
#     values = subtype_colors,
#     guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
#   ) +
#   scale_x_continuous(limits= c(min(plt.nmf.pca.150$NMF.PC1), max(plt.nmf.pca.150$NMF.PC1)) * 1.1  ) +
#   scale_y_continuous(limits= c(min(plt.nmf.pca.150$NMF.PC2), max(plt.nmf.pca.150$NMF.PC2)) * 1.1  )
# ggsave("output/figures/paper_subtypes_nmf_pca_1_2_gliovis.png",width=7,height=7.2)
#   



classifier.lda <- MASS::lda(`gliovis.majority_call` ~ `NMF:123456.PC1` + `NMF:123456.PC2` , data = train.in <- plt.single %>%
                              dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2', 'gliovis.majority_call')) %>%
                              tibble::column_to_rownames('sid') %>%
                              as.data.frame())



plt.single <- plt.single %>% dplyr::left_join(
  predict(classifier.lda, newdata = plt.single %>%
            dplyr::select(c('sid','NMF:123456.PC1', 'NMF:123456.PC2', 'gliovis.majority_call')) %>%
            tibble::column_to_rownames('sid') %>%
            as.data.frame() ) %>%
    data.frame() %>%
    dplyr::select(c('class', 'posterior.Classical', 'posterior.Mesenchymal', 'posterior.Proneural')) %>%
    `colnames<-`( gsub('^(.+)$','NMF:123456.PCA.LDA.\\1', colnames(.))  ) %>%
    tibble::rownames_to_column('sid')
  ,  by=c('sid'='sid') )






#plot(c(-6,6), c(-6,6), type="n")
#points(
#  p$x[,1],
#  p$x[,2],
#  col = as.factor(metadata$gliovis.majority_call)
#  #col = as.factor(gsam_nmf_150$`123456`$membership)
#  ,pch=19, cex=0.6
#)# screeplot(p)

### brute force class labels for plotting ----

resolution <- 250 # 1000 x 1000 data points

range_pc1 = seq(from = min(p$x[,1]), to = max(p$x[,1]), length.out = resolution) * 1.1
range_pc2 = seq(from = min(p$x[,2]), to = max(p$x[,2]), length.out = resolution) * 1.1

range_df = expand.grid(PC1 = range_pc1, PC2 = range_pc2)
range_predictions = predict(test.lda, newdata = range_df)



plt.nmf.pca.lda.150 <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], class = metadata$gliovis.majority_call, type = "GlioVis") %>%
  rbind( data.frame(PC1 = range_df$PC1 ,PC2 = range_df$PC2 ,class = range_predictions$class,type = "LDA on NVM+PCA")) %>%
  dplyr::mutate(PC1.n = as.numeric(scale(PC1)),
                PC2.n = as.numeric(scale(PC2)))

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
ggsave("output/figures/paper_subtypes_nmf_pca_1_2_lda-contour.png",width=7,height=7.2)






plt.nmf.pca.lda.reclass.150 <- plt.nmf.pca.lda.150 %>%
  dplyr::filter(type == 'GlioVis') %>%
  dplyr::select(c('PC1','PC2', 'PC1.n', 'PC2.n'))

range_predictions = predict(test.lda, newdata = plt.nmf.pca.lda.reclass.150)

plt.nmf.pca.lda.reclass.150$class <- range_predictions$class
plt.nmf.pca.lda.reclass.150$type <- "LDA on NVM+PCA data"
#plt.nmf.pca.lda.reclass.150 <- rbind(plt.nmf.pca.lda.reclass.150, plt.nmf.pca.lda.150 %>%
#                                       dplyr::filter(type != 'GlioVis'))

rm(range_predictions)




plt.nmf.pca.lda.reclass.150.p <-

  ggplot(plt.nmf.pca.lda.reclass.150, aes(x = PC1, y = PC2, col = class)) + 
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
  scale_color_manual(
    name = "Subtype by LDA labels:",
    values = subtype_colors,
    guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)
  ) +
  scale_fill_manual(values = subtype_colors)
ggsave("output/figures/paper_subtypes_nmf_pca_1_2_lda-labels.png",width=7,height=7.2)



library(patchwork)

#plt.pca.150.p + 
  
plt.nmf.150.p + plt.nmf.pca.150.p + plt.nmf.pca.lda.150.p + plt.nmf.pca.lda.reclass.150.p

## Perform PCA on NMF coordinates ? ----

# Per pair stats ----


plt.paired <- plt.single %>%
  dplyr::group_by(pid) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::top_n(1, sid) %>%
  dplyr::select(all_of('pid')) %>%
  as.data.frame() %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R1') %>% `colnames<-`( paste0(colnames(.) , ".R1") ) , by=c('pid'='pid.R1') ) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R2') %>% `colnames<-`( paste0(colnames(.) , ".R2") ) , by=c('pid'='pid.R2') ) %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                      (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2))



plot(density(distances.paired$eucledian.dist))
hist(distances.paired$eucledian.dist,breaks=15)



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
                                      (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2))


wilcox.test(plt.anti.paired$eucledian.dist, plt.paired$eucledian.dist)
median(plt.anti.paired$eucledian.dist)
median(plt.paired$eucledian.dist)



## asd ----

# # determine the eucledian distance between any two pairs and make p-dist
# distances.paired <- plt.nmf.pca.lda.150 %>% 
#   dplyr::filter(type == 'GlioVis') %>%
#   tibble::rownames_to_column('sid') %>%
#   dplyr::mutate(pid = gsub('^(...).*','\\1',sid) ) %>%
#   dplyr::mutate(res = gsub('^...(.).*','R\\1',sid) ) %>%
#   dplyr::group_by(pid) %>%
#   dplyr::filter(n() == 2) %>%
#   dplyr::mutate(PC1 = NULL, PC2=NULL) # use normalised ones
# 
# distances.paired.R1 <- distances.paired %>%
#   dplyr::filter(res == "R1")
# 
# distances.paired.R2 <- distances.paired %>%
#   dplyr::filter(res == "R2")
# 
# stopifnot(distances.paired.R1$pid == distances.paired.R2$pid)
# 
# distances.paired <- data.frame(pid = distances.paired.R1$pid,
#                                eucledian.dist = sqrt((distances.paired.R1$PC1.n - distances.paired.R2$PC1.n)^2 +
#                                                      (distances.paired.R1$PC2.n - distances.paired.R2$PC2.n)^2),
#                                stringsAsFactors = F)
# 
# plot(density(distances.paired$eucledian.dist))
# hist(distances.paired$eucledian.dist,breaks=15)
# 

# calc distances of all R1 to any R2 from another tumour

t <- plt.nmf.pca.lda.150 %>% 
  dplyr::filter(type == 'GlioVis') %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(pid = gsub('^(...).*','\\1',sid) ) %>%
  dplyr::mutate(res = gsub('^...(.).*','R\\1',sid) ) %>%
  dplyr::mutate(PC1 = NULL, PC2=NULL) # use normaalised ones


z <- 0
distances.unpaired <- c()
for(i in 1:nrow(t)) {
  if( t[i,]$res == "R1" ) {
    for(j in 1:nrow(t)) {
      if(t[j,]$res == "R2" & t[i,]$pid != t[j,]$pid) {
        z <- z + 1
        
        distances.unpaired <- c(distances.unpaired, sqrt((t[i,]$PC1.n - t[j,]$PC1.n) ^ 2 + (t[i,]$PC2.n - t[j,]$PC2.n) ^ 2))
      }
    }
  }
}

wilcox.test(distances.unpaired, distances.paired$eucledian.dist)
median(distances.paired$eucledian.dist)
median(distances.unpaired)



## Voorspellen hoeveel swaps door toeval  ----
### 1 gamma CDS on eucledian ----


plt <- plt.nmf.pca.lda.reclass.150 %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(pid = gsub('^(...).*','\\1',sid) ) %>%
  dplyr::mutate(res = as.factor(gsub('^...(.).*','R\\1',sid)))%>%
  dplyr::mutate(i = 1:nrow(.) )



ggplot(subset(plt, i == 8 | res == "R2"), aes(x = PC1, y = PC2, col = class, shape=res, label=sid, group=pid)) + 
ggplot(plt, aes(x = PC1, y = PC2, col = class, shape=res, label=sid, group=pid)) + 
  geom_contour(data= subset(plt.nmf.pca.lda.150, type != "GlioVis") %>% dplyr::mutate(res='R1', sid=''), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,
               group=""
               ,breaks=c(1.5,2.5)
  ) +
  geom_point(size=1.5) +
  #geom_text() +
  youri_gg_theme +
  labs(x = "PC1 on NMF[k=3] (150 S-Genes)",
       y = "PC2 on NMF[k=3] (150 S-Genes)",
       col = "GlioVis Majority subtype",
       fill = "NMF/PCA/LDA subtype") + 
  scale_color_manual(name = NULL, values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors) + 
  scale_shape_manual(values = c('R1' = 19, 'R2' = 1)) +
  geom_line(alpha=0.2, col="gray60")


#fit
library(MASS)
library(fitdistrplus)
fit1 <- fitdist(distances.paired$eucledian.dist, "gamma")
plot(fit1)



u <- data.frame()
for(i in plt %>% dplyr::filter(res == "R1") %>% dplyr::pull(i)) {# pak de eerste R1:

  a = plt[i,]
  df <- data.frame(r1 = c(), r2 = c(), dist = c(), r2.subtype = c(), p = c())
  for(j in 1:nrow(plt) ) {
    if(i != j & plt[j,]$res == "R2" & plt[i,]$pid != plt[j,]$pid) {
      #print(plt[j,])
      
      d <- data.frame(r1 = plt[i,]$sid,
                      r2 = plt[j,]$sid,
                      dist = sqrt ( (plt[i,]$PC1 - plt[j,]$PC1)^2 + (plt[i,]$PC2 - plt[j,]$PC2)^2  ),
                      r2.subtype = plt[j,]$class ,
                      p = NA) %>%
        dplyr::mutate(p = 1 - pgamma(c( .$dist ), shape = fit1$estimate[1], rate = fit1$estimate[2]) )
      
      df <- rbind(df, d)
      rm(d)
    }
  }
  
  
  cl = df %>%
    dplyr::filter(r2.subtype == "Classical") %>%
    dplyr::pull('p') %>%
    sum()
  
  mes = df %>%
    dplyr::filter(r2.subtype == "Mesenchymal") %>%
    dplyr::pull('p') %>%
    sum()
  
  pr = df %>%
    dplyr::filter(r2.subtype == "Proneural") %>%
    dplyr::pull('p') %>%
    sum()
  
  
  trans = df %>%
    dplyr::filter(r2.subtype != a$class) %>%
    dplyr::pull('p') %>%
    sum() / sum(df$p) * 100
  
  ret = df %>%
    dplyr::filter(r2.subtype == a$class) %>%
    dplyr::pull('p') %>%
    sum() / sum(df$p) * 100
  
  
  u <- rbind(u, data.frame(
    sid = a$sid,
    class.r1 = a$class,
    i = a$i,
    
    p.cl = cl / (cl + mes + pr) * 100,
    p.mes = mes / (cl + mes + pr) * 100,
    p.pr = pr / (cl + mes + pr) * 100,
    
    p.retains = ret ,
    p.transion = trans
  ))

}


u %>% dplyr::pull(p.transion) %>% mean()
# [1] 12.50859

u %>% dplyr::filter(class.r1 == "Classical")  %>% dplyr::pull(p.transion) %>% mean()
# [1] 8.691213

u %>% dplyr::filter(class.r1 == "Proneural")  %>% dplyr::pull(p.transion) %>% mean()
# [1] 10.4027

u %>% dplyr::filter(class.r1 == "Mesenchymal")  %>% dplyr::pull(p.transion) %>% mean()
# [1] 25.66799


### 2 per pair, find swaps within same euc dist ----
### the regulat eucledian way does not incorporate that larger distances typicall happen in those furthest away..


# 1. find all r1's with matching r2


r1.paired <- plt.nmf.pca.lda.150 %>% 
  dplyr::filter(type == 'GlioVis') %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(pid = gsub('^(...).*','\\1',sid) ) %>%
  dplyr::mutate(res = gsub('^...(.).*','R\\1',sid) ) %>%
  dplyr::group_by(pid) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::ungroup() %>%
  dplyr::filter(res == "R1") %>%
  dplyr::pull(sid)



probable_transitions <- data.frame()
for(sid.r1 in r1.paired) {
  #sid.r1 <- 'CAD1'
  
  a <- plt %>% dplyr::filter(sid == sid.r1)
  a.r2 <- plt %>% dplyr::filter(pid == a$pid & res == "R2")
  
  d <- sqrt((a$PC1.n - a.r2$PC1.n)^2 + (a$PC2.n - a.r2$PC2.n)^2)
  
  to_cl <- 0
  to_pr <- 0
  to_mes <- 0
  
  for(sid.r2 in plt$sid) {
    b <- plt %>% dplyr::filter(sid == sid.r2)
    
    if(a$pid != b$pid & b$res == 'R2') {
    #if(b$res == 'R2') {
      d.b <- sqrt((a$PC1.n - b$PC1.n)^2 + (a$PC2.n - b$PC2.n)^2)
      
      if(d.b <= d) { # fits within 'circle'
        #print(b)
        
        probable_transitions <- rbind(
          probable_transitions ,
          data.frame(from = a$sid,
                     from.class = a$class,
                     to =  b$pid,
                     to.class = b$class,
                     dist = d.b )
        )
        
      }
    }
  }
}



probable_transitions.s <- probable_transitions %>%
  dplyr::group_by(from) %>%
  dplyr::summarize(to.cl = sum(to.class == 'Classical'),
                   to.mes = sum(to.class == 'Mesenchymal'),
                   to.pn = sum(to.class == 'Proneural')
                   ) %>%
  dplyr::mutate(to.n = to.cl + to.mes + to.pn) %>%
  dplyr::left_join(plt %>% dplyr::select(c('sid','class')) %>% dplyr::rename(class.r1 = class) , by = c('from' = 'sid')) %>%
  dplyr::mutate(p.retained = case_when(
    class.r1 == "Classical" ~ to.cl / to.n,
    class.r1 == "Mesenchymal" ~ to.mes / to.n,
    class.r1 == "Proneural" ~ to.pn / to.n
  )) %>%
  dplyr::mutate(p.transition = 1 - p.retained) %>%
  dplyr::full_join(
    data.frame(sid=r1.paired) %>%
      dplyr::left_join(
        plt %>%
          dplyr::select(c('sid','class')) %>%
          dplyr::rename(class.r1.tmp = class)
        , by=c('sid'='sid'))
    , by = c('from' = 'sid')) %>%
  dplyr::mutate(pid = gsub('^(...).*$','\\1', from) ) %>%
  dplyr::left_join(
    plt %>%
      dplyr::filter(res == 'R2') %>%
      dplyr::select(c('sid','class','pid') ) %>%
      dplyr::rename(class.r2 = class) %>%
      dplyr::rename(to = sid)
    ,
    by = c('pid' = 'pid')
  ) %>%
  dplyr::mutate(p.retained = ifelse(is.na(p.retained), 1, p.retained) ) %>%
  dplyr::mutate(p.transition = ifelse(is.na(p.transition), 1, p.transition) ) %>%
  dplyr::mutate(class.r1 = as.factor( ifelse( is.na(class.r1),as.character(class.r1.tmp), as.character(class.r1) ) ) ) %>%
  dplyr::mutate(class.r1.tmp = NULL) %>%
  dplyr::mutate(status = as.factor(ifelse(class.r1 == class.r2,"retained","transition")))


sum(probable_transitions.s$p.transition)
sum(probable_transitions.s$p.retained)
probable_transitions.s$status




# er zijn er 2 die zo dicht bij hun partner liggen dat er geen ttansities worden gevonden
probable_transitions.s %>%
  dplyr::pull(transition) %>%
  mean()
# [1] 0.1887275

probable_transitions.s %>%
  dplyr::filter(class.r1 == "Classical") %>%
  dplyr::pull(transition) %>% mean()
# [1] 0.133278

probable_transitions.s %>%
  dplyr::filter(class.r1 == "Proneural") %>%
  dplyr::pull(transition) %>% mean()
# [1] 0.2110109

probable_transitions.s %>%
  dplyr::filter(class.r1 == "Mesenchymal") %>%
  dplyr::pull(transition) %>% mean()
# [1] 0.1887275


## find actual transitions (distribution) ----


# TODO check aantal transities met gliovis gewone klassificate
# TODO check aantal transities met gliovis deze klassificate
# TODO check aantal transities met gliovis NMF klassificate


