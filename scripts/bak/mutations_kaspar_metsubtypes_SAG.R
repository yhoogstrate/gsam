#!/usr/bin/env R 

# load libs ----


library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(survival)
library(survminer)


# small check

# df.new <- read.table("output/tables/gsam_nmf_lda_data.new.txt") %>%
#   `colnames<-`(paste0(colnames(.), ".new"))
# df.old <- read.table("output/tables/gsam_nmf_lda_data.txt") %>%
#   `colnames<-`(paste0(colnames(.), ".old"))
# 
# df <- dplyr::inner_join(df.old, df.new, by=c('sid.old'='sid.new'))



#  load data ----


source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')
source('scripts/R/palette.R')

source('scripts/R/gsam_rna-seq_expression.R') # recursively calls metadata
source('scripts/R/gsam_metadata.R') # recursively calls metadata

#source('scripts/R/subtype_genes.R')
#source('scripts/R/glass_expression_matrix.R')
#source('data/wang/msig.library.12.R')  # no license w/ code provided, can't include it in source tree
#source('scripts/R/wang_glioma_intrinsic_genes.R')





#### ---- surv plot ---- ####

dat_surv <- dat %>%
            dplyr::mutate(survivalFromSecondSurgeryMonths = survivalFromSecondSurgeryDays * 0.032855) %>%
            dplyr::mutate(monthsToSecondSurgery = daysToSecondSurgery * 0.032855) %>%
            dplyr::select(pid,survivalFromSecondSurgeryMonths,monthsToSecondSurgery) %>%
            melt(., id.vars = "pid") 

surv_plot <- ggplot(dat_surv, aes(x=pid, y=value, fill=variable)) +
             geom_bar(position = "stack", stat='identity') +
             job_gg_theme +
             labs(y = "Survival (months)") +
             theme(axis.title.x=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             legend.title = element_blank()) +
             scale_fill_brewer(palette = "Reds") +
             coord_equal()

#### ---- treatment plot ---- ####
dat_patient <- dat %>%
               dplyr::select(pid,gender,HM,treatedWithTMZ,treatedWithRT) %>%
               melt(., "pid")

treatment_plot <- ggplot(dat_patient, aes(x = pid, y = variable, fill = value)) +
                  geom_tile(colour = "black", size = 0.3) +
                  theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.ticks.y =element_blank(),
                  legend.title = element_blank(),
                  legend.key.size = unit(0.3, 'cm')) +
                  scale_fill_brewer(palette = "PuOr") +
                  coord_equal() +
                  scale_y_discrete(labels= c("Sex","HM", "TMZ", "RT"))

#### ---- plot Kaspar ----
dat_kasp <- dat %>%
            dplyr::filter(IDH1 == "Wildtype") %>%
            dplyr::select(cnStatusEGFRs,AXIN2,APC,IDH2,IDH1,JAK2,BUB1B,RB1,MSH2,BRCA1,
                          MSH6,BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,KRAS,
                          MAP3K1,NF1,KIT,EPHA3,IGF1R,FGFR3,ERBB3,
                          EGFR,TP53BP1,TP53,PIK3R1,MTOR,TSC2,PIK3CA,
                          PIK3C2G,PTEN,pid,nMutationsRecurrent) 

ordering_kasp <- dat_kasp %>%
                 dplyr::select(pid,nMutationsRecurrent) 

dat_kasp <- dat_kasp %>%
            dplyr::select(-nMutationsRecurrent) %>%
            melt(., "pid") %>%
            left_join(.,ordering_kasp,by = "pid")

dat_kasp_plot <- ggplot(dat_kasp, aes(x = reorder(pid,-nMutationsRecurrent), y = variable, fill = value)) +
                 geom_tile(colour = "black", size = 0.3) +
                 theme(axis.title.x=element_blank(),
                 axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
                 axis.title.y=element_blank(),
                 axis.ticks.y =element_blank(),
                 legend.title = element_blank(),
                 axis.text=element_text(size=5),
                 legend.key.size = unit(0.3, 'cm')) +
                 scale_fill_manual(values = c("#bb5f6c", "#79b1b1", "#2e415e", "white","grey")) +
                 coord_equal() 

#ggsave("output/figures/figures_SAG/remake_kaspar_plot.pdf")

#### ---- cnv plot  ---- ####
pid_no_pair <- dat %>% 
               dplyr::filter(! is.na(NMF.123456.PCA.SVM.class)) %>%
               dplyr::group_by(pid) %>%
               dplyr::summarise(n = n()) %>%
               dplyr::filter(n == 1) %>%
               dplyr::pull(pid)

dat_cnv <- dat %>%
           dplyr::filter(! pid %in% pid_no_pair) %>%
           dplyr::filter(! is.na(NMF.123456.PCA.SVM.class)) %>%
           dplyr::group_by(pid) %>% 
           dplyr::mutate(transition=length(unique(NMF.123456.PCA.SVM.class))) %>%
           dplyr::mutate(transition = ifelse(transition == 2, T, F)) %>% 
           dplyr::mutate(CL_stable=ifelse(transition==F & NMF.123456.PCA.SVM.class == "Classical" & resection == "r2",T,F)) %>%
           dplyr::mutate(MES_stable=ifelse(transition==F & NMF.123456.PCA.SVM.class == "Mesenchymal" & resection == "r2",T,F)) %>%
           dplyr::mutate(PN_stable=ifelse(transition==F & NMF.123456.PCA.SVM.class == "Proneural" & resection == "r2",T,F)) %>%
           dplyr::mutate(to_CL=ifelse(transition==T & resection == "r2" & NMF.123456.PCA.SVM.class == "Classical",T,F)) %>%
           dplyr::mutate(to_MES=ifelse(transition==T & resection == "r2" & NMF.123456.PCA.SVM.class == "Mesenchymal",T,F)) %>%
           dplyr::mutate(to_PN=ifelse(transition==T & resection == "r2" & NMF.123456.PCA.SVM.class == "Proneural",T,F)) %>%
           dplyr::filter(resection == "r2") %>%
           dplyr::ungroup() %>%
           dplyr::select(pid,PTEN,PIK3C2G,PIK3CA,TSC2,MTOR,PIK3R1,TP53, TP53BP1, EGFR, ERBB3, FGFR3, IGF1R,
                         EPHA3, KIT, NF1, MAP3K1, KRAS, KMT2D, KMT2C, ARID2, SETD2, ATM, BRCA2, MSH6, BRCA1,
                         MSH2, RB1, BUB1B, JAK2, IDH1, IDH2, APC, AXIN2, cnStatusCDKN2ABs, transition,CL_stable,MES_stable,PN_stable,to_CL,to_MES,to_PN) %>%
           dplyr::mutate(order.x =  match(1:nrow(.), order(MES_stable, to_MES, EGFR)))

dat_cnv_plot <- melt(dat_cnv, id.vars=c("pid", "order.x"))

cnv_plot <- ggplot(dat_cnv_plot, aes(x = reorder(pid, order.x), y = variable, fill = value)) +
            geom_tile(colour = "black", size = 0.3) +
            theme(axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
                  axis.title.y=element_blank(),
                  axis.ticks.y =element_blank(),
                  legend.title = element_blank(),
                  axis.text=element_text(size=5),
                  legend.key.size = unit(0.3, 'cm')) +
                  scale_fill_manual(values = c("#f59d87", "grey", "orange", "#555b8a", "#eb4634", "white", "green")) +
                  coord_equal() 

cnv_plot

#ggsave("output/figures/figures_SAG/subtype_CNV_ordsubtype.pdf")

#### ---- cnv plot order on different vars ---- ####
ordering <- dat_cnv %>%
            dplyr::select(pid,TP53,transition)
ordering$TP53 <- as.factor(ordering$TP53)

dat_cnv_plot <- melt(dat_cnv, "pid") %>%
                left_join(.,ordering,by = "pid")

cnv_plot <- ggplot(dat_cnv_plot, aes(x = reorder(pid,TP53), y = variable, fill = value)) +
            geom_tile(colour = "black", size = 0.3) +
            theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y =element_blank(),
            legend.title = element_blank(),
            axis.text=element_text(size=4),
            legend.key.size = unit(0.3, 'cm')) +
            scale_fill_manual(values = c("#f59d87", "grey", "orange", "#555b8a", "#eb4634", "white", "green")) +
            coord_equal() 

#### ---- statistical test ----
tmp <- dat_cnv %>%
       dplyr::mutate(MES = ifelse(MES_stable == T | to_MES == T, TRUE, FALSE)) %>%
       dplyr::select(-CL_stable,-MES_stable,-PN_stable,-to_CL,-transition, -to_MES, -to_PN,-pid, -IDH1, -order.x) %>%
       type.convert()

chi_square <- lapply(tmp %>% dplyr::select(-MES), 
              function(x) chisq.test(tmp$MES, x))
chi_square <- as.data.frame(do.call(rbind, chi_square)[,c(1,3)])$p.value
chi_square <- t(as.data.frame(chi_square))

model <- glm(transition ~.,family=binomial(link='logit'),data=tmp)

tab_tp53 <- table(tmp$SETD2,tmp$transition) %>%
            sweep(., 2, colSums(.), FUN = "/")
tab_tp53

anova(model, test="Chisq")

#### ---- combined plot ---- ####
ordering <- dat_cnv %>%
  dplyr::select(pid,transition,CL_stable,PN_stable,MES_stable)

dat_cnv <- melt(dat_cnv, "pid") %>%
  left_join(.,ordering,by = "pid")

cnv_plot <- ggplot(dat_cnv, aes(x = reorder(reorder(reorder(pid,transition),to_MES),to_PN), y = variable, fill = value)) +
  geom_tile(colour = "black", size = 0.3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y =element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size=4),
        legend.key.size = unit(0.3, 'cm')) +
  scale_fill_manual(values = c("#f59d87", "grey", "orange", "#555b8a", "#eb4634", "white", "green")) +
  coord_equal() 


# :::::: [ new stuff ] :::::: ----
## load data ----


components <- gsam.rna.metadata %>%
  dplyr::select(contains(".component") | contains(".component") | 'sid' | 'pid') %>%
  dplyr::filter(!is.na(neuron.component))


components.paired <- dplyr::inner_join(
  components %>% 
    dplyr::filter(grepl("1", sid)) %>%
    `colnames<-`(paste0(colnames(.), ".R1")) %>%
    dplyr::rename(pid = pid.R1)
  ,
  components %>% 
    dplyr::filter(grepl("2", sid)) %>%
    `colnames<-`(paste0(colnames(.), ".R2")) %>%
    dplyr::rename(pid = pid.R2)
  , by=c('pid'='pid')) %>%
  dplyr::mutate(delta.neuron.component = neuron.component.R2 - neuron.component.R1) %>% 
  dplyr::mutate(delta.oligodendrocyte.component = oligodendrocyte.component.R2 - oligodendrocyte.component.R1) %>% 
  dplyr::mutate(delta.endothelial.component = endothelial.component.R2 - endothelial.component.R1) %>% 
  dplyr::mutate(delta.extracellular.matrix.component = extracellular.matrix.component.R2 - extracellular.matrix.component.R1)



##  mes  ----


##  neuronal  ----


rbfox3 <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::filter(grepl('RBFOX3', gid)) %>% 
  dplyr::mutate(gid=NULL) %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::rename(RBFOX3 = V1) %>% 
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(resection = paste0('R',gsub('^...(.).*$','\\1',sid))) %>%
  dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))

egfr <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::filter(grepl('EGFR', gid)) %>% 
  dplyr::mutate(gid=NULL) %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::rename(EGFR = V1) %>% 
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(resection = paste0('R',gsub('^...(.).*$','\\1',sid))) %>%
  dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
  dplyr::mutate(col = pid %in% c('CDD', 'CDA', 'JAK', 'HAE', 'AZG', 'AZF'))

ggplot(egfr, aes(x = reorder(sid, EGFR) , y=EGFR, col=col)) +
  geom_point()






plt.tight <- gsam.patient.metadata %>% 
  dplyr::filter(IDH1 == "Wildtype" & IDH2 == "Wildtype") %>%
  dplyr::rename(pid = studyID ) %>%
  dplyr::inner_join(components.paired , by=c('pid'='pid') ) %>%
  #dplyr::mutate(order.1 =  match(1:nrow(.), order(- (nMutationsRecurrent - nMutationsPrimary)))) %>%
  dplyr::mutate(order.1 =  match(1:nrow(.), order(delta.neuron.component))) %>%
  dplyr::left_join(
    rbfox3 %>% dplyr::filter(resection == "R1") %>%  dplyr::select(pid, RBFOX3) %>% dplyr::rename(RBFOX3.R1 = RBFOX3) , by = c('pid'='pid')
  ) %>%
  dplyr::left_join(
    rbfox3 %>% dplyr::filter(resection == "R2") %>%  dplyr::select(pid, RBFOX3) %>% dplyr::rename(RBFOX3.R2 = RBFOX3) , by = c('pid'='pid')
  ) %>% 
  dplyr::mutate(RBFOX3.delta = RBFOX3.R2 - RBFOX3.R1)

## is inderdaad negatieve component

# ggplot(rbfox3, aes(x = resection, y = V1, fill=resection)) +
#   geom_violin(draw_quantiles = c(0.6), col="black", alpha=0.2) +
#   geom_jitter( position=position_jitter(0.2), size=2.5, pch=21, col="black") +
#   labs(x = NULL, col = NA, y = "Tumor cell percentage" ) +
#   job_gg_theme


# plt.tight <- plt.tight %>%
#   dplyr::filter(grepl("^C", pid))



plt <- plt.tight %>%
  dplyr::select(c(pid, survivalDays , order.1 ))

p1 <- ggplot(plt, aes(x = reorder(pid, order.1), y=survivalDays )) +
  geom_col() + 
  labs(y="survival (Days)") +
  theme(axis.text.x = element_text(angle = 90))




plt <- plt.tight %>%
  dplyr::select(c(pid, neuron.component.R2 , neuron.component.R1 , order.1 ))  %>%
  melt(., id.vars=c("pid", "order.1")) %>%
  dplyr::mutate(variable = gsub('neuron.component.','',variable)) 

p2 <- ggplot(plt, aes(x = reorder(pid, order.1), y=value, col=variable )) +
  geom_point() + 
  labs(y="Neuronal component score") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))


plt <- plt.tight %>%
  dplyr::select(pid,PTEN,PIK3C2G,PIK3CA,TSC2,MTOR,PIK3R1,TP53, TP53BP1, EGFR, ERBB3, FGFR3, IGF1R,
                EPHA3, KIT, NF1, MAP3K1, KRAS, KMT2D, KMT2C, ARID2, SETD2, ATM, BRCA2, MSH6, BRCA1,
                MSH2, RB1, BUB1B, JAK2, IDH1, IDH2, APC, AXIN2, cnStatusCDKN2ABs, order.1 , mgmtStability,tumorLocation) %>%
  dplyr::mutate(mgmtStability = ifelse(mgmtStability == "Stable methylated","Stable",mgmtStability)) %>%
  dplyr::mutate(mgmtStability = ifelse(mgmtStability == "Stable unmethylated","Wildtype",mgmtStability)) %>%
  #dplyr::select(c(pid, NF1, EGFR, order.1)) %>%
  melt(., id.vars=c("pid", "order.1"))


p3 <- ggplot(plt, aes(x = reorder(pid, order.1), y = variable, fill = value)) +
  geom_tile(colour = "black", size = 0.3) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
        axis.title.y=element_blank(),
        axis.ticks.y =element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size=5),
        legend.key.size = unit(0.3, 'cm')) +
  #scale_fill_manual(values = c("#bb5f6c", "#79b1b1", "#2e415e", "white","grey")) +
  coord_equal() 




(p1 / p2 / p3)  + plot_layout( heights=c(1.5, 1.5, 4) )





# ggplot(plt.tight, aes(x = CDKN2AB, y=neuron.component.R2  )) + 
#   geom_violin(draw_quantiles = c(0.5), col="black") +
#   geom_jitter( position=position_jitter(0.2), size=0.9)
# 




##  oligodendrocyte scores etc. ----


plt.tight <- gsam.patient.metadata %>% 
  dplyr::filter(IDH1 == "Wildtype" & IDH2 == "Wildtype") %>%
  dplyr::rename(pid = studyID ) %>%
  dplyr::inner_join(components.paired , by=c('pid'='pid') ) %>%
  #dplyr::mutate(order.1 =  match(1:nrow(.), order(- (nMutationsRecurrent - nMutationsPrimary)))) %>%
  dplyr::mutate(order.1 =  match(1:nrow(.), order(delta.oligodendrocyte.component)))
  # dplyr::mutate(order.1 =  match(1:nrow(.), order(oligodendrocyte.component.R2))) 


plt <- plt.tight %>%
  dplyr::select(c(pid, survivalDays , order.1 ))

p1 <- ggplot(plt, aes(x = reorder(pid, order.1), y=survivalDays )) +
  geom_col() + 
  labs(y="survival (Days)")



plt <- plt.tight %>%
  dplyr::select(c(pid, oligodendrocyte.component.R2 , oligodendrocyte.component.R1 , order.1 ))  %>%
  melt(., id.vars=c("pid", "order.1")) %>%
  dplyr::mutate(variable = gsub('oligodendrocyte.component.','',variable)) 

p2 <- ggplot(plt, aes(x = reorder(pid, order.1), y=value, col=variable )) +
  geom_point() + 
  labs(y="Neuronal component score")



plt <- plt.tight %>%
  dplyr::select(pid,PTEN,PIK3C2G,PIK3CA,TSC2,MTOR,PIK3R1,TP53, TP53BP1, EGFR, ERBB3, FGFR3, IGF1R,
                 EPHA3, KIT, NF1, MAP3K1, KRAS, KMT2D, KMT2C, ARID2, SETD2, ATM, BRCA2, MSH6, BRCA1,
                 MSH2, RB1, BUB1B, JAK2, IDH1, IDH2, APC, AXIN2, cnStatusCDKN2ABs, order.1 , mgmtStability) %>%
  dplyr::mutate(mgmtStability = ifelse(mgmtStability == "Stable methylated","Stable",mgmtStability)) %>%
  dplyr::mutate(mgmtStability = ifelse(mgmtStability == "Stable unmethylated","Wildtype",mgmtStability)) %>%
  ###dplyr::select(c(pid, NF1, EGFR, order.1)) %>%
  melt(., id.vars=c("pid", "order.1"))


# plt <- plt.tight %>%
#   dplyr::select(pid, order.1 ,  tumorLocation) %>%
#   dplyr::mutate(tumorLocation = as.factor(tumorLocation)) %>% 
#   dplyr::mutate(tumorLocation.NA = is.na(tumorLocation)) %>% 
#   dplyr::mutate(tumorLocation.Fossa.posterior = tumorLocation == "Fossa posterior") %>% 
#   dplyr::mutate(tumorLocation.Frontal = tumorLocation == "Frontal") %>% 
#   dplyr::mutate(tumorLocation.Occipital = tumorLocation == "Occipital") %>% 
#   dplyr::mutate(tumorLocation.Perietal = tumorLocation == "Parietal") %>% 
#   dplyr::mutate(tumorLocation.Temporal = tumorLocation == "Temporal") %>%
#   dplyr::mutate(tumorLocation = NULL) %>% 
#   melt(., id.vars=c("pid", "order.1"))


p3 <- ggplot(plt, aes(x = reorder(pid, order.1), y = variable, fill = value)) +
  geom_tile(colour = "black", size = 0.3) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
        axis.title.y=element_blank(),
        axis.ticks.y =element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size=5),
        legend.key.size = unit(0.3, 'cm')) +
  ####scale_fill_manual(values = c("white", "black", "grey", "white","grey")) +
  scale_fill_manual(values = c("#bb5f6c", "#79b1b1", "#2e415e", "white","grey")) +
  coord_equal() 



(p1 / p2 / p3) + plot_layout( heights=c(1, 1.5, 4) )




## EM-genes & signature ----

# gsam.rna.metadata %>%
#   dplyr::filter(pid %in% c("CAV", "CBG", "CDF", "DAB", "EBD", "ECB", "FAL", "JAB", "JAD", "KAC", "BAW")) %>%
#   dplyr::select(pid, pat.with.IDH)



plt.tight <- gsam.patient.metadata %>% 
  dplyr::rename(pid = studyID ) %>%
  dplyr::filter(pid %in% (gsam.rna.metadata %>% dplyr::filter(pat.with.IDH) %>% dplyr::pull(pid) %>% unique()) == F) %>% 
  dplyr::inner_join(components.paired , by=c('pid'='pid') ) %>%
  #dplyr::mutate(order.1 =  match(1:nrow(.), order(- (nMutationsRecurrent - nMutationsPrimary)))) %>%
  #dplyr::mutate(order.1 =  match(1:nrow(.), order(extracellular.matrix.component.R2)))
  dplyr::mutate(order.1 =  match(1:nrow(.), order(delta.extracellular.matrix.component))) %>% 
  dplyr::left_join(
      gsam.rna.metadata %>%
        dplyr::filter(resection =='r1') %>% 
        dplyr::select(sid, NMF.123456.PCA.SVM.class) %>% 
        dplyr::rename(NMF.123456.PCA.SVM.class.R1 = NMF.123456.PCA.SVM.class)
      ,by = c('sid.R1' = 'sid')) %>%
  dplyr::left_join(
      gsam.rna.metadata %>%
        dplyr::filter(resection =='r2') %>% 
        dplyr::select(sid, NMF.123456.PCA.SVM.class) %>% 
        dplyr::rename(NMF.123456.PCA.SVM.class.R2 = NMF.123456.PCA.SVM.class),
      by = c('sid.R2' = 'sid'))


plt <- plt.tight %>%
  dplyr::select(c(pid, survivalDays, survivalFromSecondSurgeryDays , order.1 )) %>%
  dplyr::mutate(survivalBeforeSecondSurgery = survivalDays - survivalFromSecondSurgeryDays) %>%
  dplyr::mutate(survivalDays = NULL)


plt <- rbind(plt %>%
        dplyr::mutate(survivalFromSecondSurgeryDays = NULL) %>% 
        dplyr::rename(survivalDays = survivalBeforeSecondSurgery) %>%
        dplyr::mutate(survivalDays = -1 * survivalDays) %>% 
        dplyr::mutate(type = "Time to second surgery")
      ,
      plt %>%
        dplyr::mutate(survivalBeforeSecondSurgery = NULL) %>% 
        dplyr::rename(survivalDays = survivalFromSecondSurgeryDays) %>%
        dplyr::mutate(type = "Time from second surgery")
      )

p1 <- ggplot(plt, aes(x = reorder(pid, order.1), y=survivalDays, col=type )) +
  geom_col(fill="gray60", col="black", lwd=0.4) +
  labs(y="survival (Days)", x="G-SAM patient") + 
  job_gg_theme +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))



plt <- plt.tight %>%
  dplyr::select(c(pid, extracellular.matrix.component.R2 , extracellular.matrix.component.R1 , order.1 ))  %>%
  melt(., id.vars=c("pid", "order.1")) %>%
  dplyr::mutate(variable = as.factor(as.character(variable))) %>% # melt orders on R2, R1 ? weird, probably inital column order?
  dplyr::arrange(pid, variable) %>%
  dplyr::mutate(variable = gsub('extracellular.matrix.component.','',variable)) %>%
  dplyr::left_join(
    plt.tight %>%
      dplyr::mutate(em.pc.status = ifelse(extracellular.matrix.component.R2  > extracellular.matrix.component.R1 , "increase", "decrease")) %>%
      dplyr::select(c(pid, em.pc.status)),
    by=c('pid'='pid')
  )



p2 <- ggplot(plt, aes(x = reorder(pid, order.1), y=value, col=em.pc.status, group=pid))  +
  geom_hline(yintercept=2.5, lty=1, color = "#FFFFFF44",lwd=3) +
  geom_point(data = subset(plt, variable == "R1"), pch=19, cex=1.2, alpha=0.25) +
  geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 )  + 
  #labs(y="Extra-cellular Matrix principal component", x="G-SAM patient", col="Extra-cellular principal component status") +
  labs(x=NULL, col=NULL, y="Extra-cellular Matrix principal component") +
  youri_gg_theme +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c(
    "increase"="#bb5f6c", # rood #bb5f6c
    "decrease"="#79b1b1" # lichtlauw #79b1b
  )) +
  geom_hline(yintercept=2.5, lty=3, color = "red",lwd=0.5)



plt <- plt.tight %>%
  dplyr::select(pid, order.1,
                #PTEN,PIK3C2G,PIK3CA,TSC2,MTOR,PIK3R1,TP53, TP53BP1, EGFR, ERBB3, FGFR3, IGF1R,
                #EPHA3, KIT, NF1, MAP3K1, KRAS, KMT2D, KMT2C, ARID2, SETD2, ATM, BRCA2, MSH6, BRCA1,
                #MSH2, RB1, BUB1B, JAK2, IDH1, IDH2, APC, AXIN2, cnStatusCDKN2ABs, order.1 , mgmtStability,
                #HM, NMF.123456.PCA.SVM.class.R1 , NMF.123456.PCA.SVM.class.R2
                
                mgmtStability,HM,
                AXIN2,APC,JAK2,
                RB1,MSH2,BRCA1,
                BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,
                NF1,ERBB3,
                EGFR,TP53BP1,TP53,PIK3R1,PIK3CA,TSC2,
                # PI3K, TP53Signaling, Wnt, Telomere, RTK, RAS, DNADamage, primary7, primary10, recurrent7, recurrent10, cnStatusEGFRs,
                
                cnStatusCDKN2ABs, 
                cnStatusRB1s,
                cnStatusNF1s,
                cnStatusCDK4s,
                cnStatusMDM2s,
                
                SETD2,PDGFRA,
                
                NMF.123456.PCA.SVM.class.R1 , NMF.123456.PCA.SVM.class.R2
                ) %>%
    dplyr::mutate(mgmtStability = ifelse(mgmtStability == "Stable methylated","Stable",mgmtStability)) %>%
    dplyr::mutate(mgmtStability = ifelse(mgmtStability == "Stable unmethylated","Wildtype",mgmtStability)) %>%
  dplyr::rename(`Sub-type R1` = NMF.123456.PCA.SVM.class.R1) %>%
  dplyr::rename(`Sub-type R2` = NMF.123456.PCA.SVM.class.R2) %>%
  melt(., id.vars=c("pid", "order.1")) %>%
  dplyr::mutate(value = gsub('^Gained$','Gained/increased',value)) %>% 
  dplyr::mutate(value = gsub('^Lost$','Lost/decreased',value)) %>% 
  dplyr::mutate(panel = case_when(
    grepl("Sub-type", variable) ~ "A",
    variable %in% c("EGFR",  "cnStatusCDKN2ABs") ~ "B",
    T ~ "C"
  ))


p3 <- ggplot(plt, aes(x = reorder(pid, order.1), y = variable, fill = value)) +
  facet_grid(rows=vars(panel), scales="free", space="free") + 
  geom_tile(colour = "black", size = 0.3) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
        axis.title.y=element_blank(),
        axis.ticks.y =element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size=5),
        legend.key.size = unit(0.3, 'cm')) +
  scale_fill_manual(values = c( subtype_colors ,
                                "Wildtype"="white",
                                "No"="white",
                                
                                "Gained/increased"="#bb5f6c", # rood #bb5f6c
                                "Lost/decreased"="#79b1b1", # lichtlauw #79b1b1
                                "Stable"="#2e415e", # donker blauw #2e415e
                                "Yes" = "#2e415e", # zelfde als stable #2e415e
                                
                                "NA"="grey")) + 
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom')


(p1 / p2 / p3) + plot_layout( heights=c(2, 1.5, 2.5) )



ggsave("output/figures/em_pat_scores.pdf",heigh=15/1.5,width=15)





test <- plt.tight %>%
  dplyr::mutate(to.mes = NMF.123456.PCA.SVM.class.R1 != "Mesenchymal" & NMF.123456.PCA.SVM.class.R2 == "Mesenchymal") %>%
  dplyr::select(pid, to.mes, delta.extracellular.matrix.component)

wilcox.test(
  test %>%dplyr::filter(to.mes == F) %>% dplyr::pull(delta.extracellular.matrix.component) ,
  test %>%dplyr::filter(to.mes != F) %>% dplyr::pull(delta.extracellular.matrix.component) ,
)


plot(gsam.rna.metadata$extracellular.matrix.component,
     gsam.rna.metadata$`NMF:123456.1` )

plot(gsam.rna.metadata$extracellular.matrix.component,
     gsam.rna.metadata$`NMF:123456.2` )

plot(gsam.rna.metadata$extracellular.matrix.component,
     gsam.rna.metadata$`NMF:123456.3` )


tmp <- gsam.rna.metadata %>% 
  #dplyr::filter(resection == 'r1') %>% 
  dplyr::select(sid, pid, extracellular.matrix.component , `NMF:123456.1`,`NMF:123456.2`,`NMF:123456.3`, tumour.percentage.dna ) %>%
  dplyr::left_join(
    gsam.patient.metadata %>%dplyr::select(studyID, survivalDays, survivalFromSecondSurgeryDays),
    by=c('pid'='studyID')
  )

corrplot::corrplot(
cor(tmp %>%
    dplyr::mutate(pid = NULL, sid=NULL) %>%
    dplyr::filter(!is.na( `NMF:123456.1` )) %>%
    dplyr::filter(!is.na( tumour.percentage.dna )))
)


corrplot::corrplot(
  cor(tmp %>%
        dplyr::mutate(pid = NULL, sid=NULL) %>%
        dplyr::filter(!is.na( `NMF:123456.2` )) %>%
        dplyr::filter(!is.na( tumour.percentage.dna )))
)
corrplot::corrplot(
  cor(tmp %>%
        dplyr::mutate(pid = NULL, sid=NULL) %>%
        dplyr::filter(!is.na( `NMF:123456.3` )) %>%
        dplyr::filter(!is.na( tumour.percentage.dna )))
)


  

ggplot(plt.tight, aes(x=survivalDays , y = extracellular.matrix.component.R1)) + 
  geom_point()

ggplot(plt.tight, aes(x=survivalDays , y = extracellular.matrix.component.R2)) + 
  geom_point()

ggplot(plt.tight, aes(x=survivalFromSecondSurgeryDays , y = extracellular.matrix.component.R2)) + 
  geom_point()








plt <- rbind(plt.tight %>%
    dplyr::mutate(resection = "R1") %>%
    dplyr::mutate(extracellular.matrix.component = extracellular.matrix.component.R1) %>%
    dplyr::mutate(tts = survivalDays ) %>%
    dplyr::select(resection, extracellular.matrix.component, tts, pid,survivalFromSecondSurgeryDays,
                  survivalDays,extracellular.matrix.component.R1,extracellular.matrix.component.R2,
                  treatedWithTMZ,
                  tumorLocation
                  )
,
plt.tight %>%
  dplyr::mutate(resection = "R2") %>%
  dplyr::mutate(extracellular.matrix.component = extracellular.matrix.component.R2) %>%
  dplyr::mutate(tts = survivalFromSecondSurgeryDays ) %>%
  dplyr::select(resection, extracellular.matrix.component, tts, pid,survivalFromSecondSurgeryDays, survivalDays,extracellular.matrix.component.R1,extracellular.matrix.component.R2,
                treatedWithTMZ,
                tumorLocation
                )) %>%
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(survival.from.second.surgery = ifelse(survivalFromSecondSurgeryDays > 365, "> 1yr", "<= 1yr") ) %>%
  dplyr::mutate(survival.from.first.surgery = ifelse(survivalDays > 365* 2.5, "> 2.5yr", "<= 2.5yr") ) %>%
  dplyr::mutate(em.status.res2 = ifelse(extracellular.matrix.component.R2 > 2.5, "high", "low")) %>%
  dplyr::mutate(survival.from.second.surgery = factor(survival.from.second.surgery, levels=c("> 1yr", "<= 1yr")))



ggplot(plt, aes(y = extracellular.matrix.component , x=tts, group=pid, col=survival.from.second.surgery ))  +
  #ggplot(plt, aes(y = extracellular.matrix.component , x=tts, group=pid, col=treatedWithTMZ ))  +
  #ggplot(plt, aes(y = extracellular.matrix.component , x=tts, group=pid, col=extracellular.matrix.component.R2 ))  +
  geom_point(data = subset(plt, resection == "R1"), pch=19, cex=0.8, alpha=0.5) +
  geom_path(data = subset(plt,survival.from.second.surgery == "<= 1yr"),arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.6 )  + 
  geom_path(data = subset(plt,survival.from.second.surgery != "<= 1yr"),arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.25, lwd=3, col="white" )  + 
  geom_path(data = subset(plt,survival.from.second.surgery != "<= 1yr"),arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 1 )  + 
  labs(y = "Extracellular matrix signature", x = "Time until death") +
  geom_hline(yintercept = 1.5, lty=1, lwd=3.5, col="#FFFFFFBB") +
  geom_hline(yintercept = 1.5, lty=2, col="red", lwd=0.25) +
  geom_vline(xintercept = 0, lty=1, col="black", lwd=0.35) +
  geom_vline(xintercept = -6, lty=1, col="black", lwd=0.35) +
  scale_x_reverse() +
  scale_color_manual(values = c('#009E74', '#CB75A4') ) +
  youri_gg_theme +
  labs(col = 'Survival from 2nd surgery')


ggsave("output/figures/survival_x_em_res.png",width=8*1.5,height=5*1.1)



ggplot(plt, aes(y = extracellular.matrix.component , x=tts, group=pid, fill=em.status.res2 ))  +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.3, col="gray50" )  + 
  geom_point(data = subset(plt, resection == "R1"), pch=21, cex=1.2, alpha=0.7, fill="gray50", col="gray50") +
  geom_point(data = subset(plt, resection == "R2"), pch=19, cex=4.5, alpha=0.25, col="white") +
  geom_point(data = subset(plt, resection == "R2"), pch=21, cex=2.0, alpha=0.7) +
  #geom_path(data = subset(plt,survival.from.second.surgery == "<= 1yr"),arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.6 )  + 
  #geom_path(data = subset(plt,survival.from.second.surgery != "<= 1yr"),arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.25, lwd=3, col="white" )  + 
  #geom_path(data = subset(plt,survival.from.second.surgery != "<= 1yr"),arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 1 )  + 
  labs(y = "Extracellular matrix signature", x = "Time until death") +
  geom_hline(yintercept = 1.5, lty=1, lwd=3.5, col="#FFFFFFBB") +
  geom_hline(yintercept = 1.5, lty=2, col="red", lwd=0.25) +
  geom_vline(xintercept = 0, lty=1, col="black", lwd=0.35) +
  geom_vline(xintercept = -6, lty=1, col="black", lwd=0.35) +
  scale_x_reverse() +
  scale_fill_manual(values = c('high'='#009E74', 'low'='#CB75A4') ) +
  youri_gg_theme +
  labs(col = 'Survival from 2nd surgery')


ggsave("output/figures/survival_x_em_res_svvl_col.pdf",width=15,height=6)



decision <- 2.5

svvl <- plt.tight %>%
  dplyr::mutate(EM.Signature.status = as.factor(ifelse(extracellular.matrix.component.R2 > decision, "high", "low"))) %>%
  dplyr::mutate(event=1)

surv_object <- Surv(time = svvl$survivalFromSecondSurgeryDays, event=svvl$event)
fit1 <- survfit(surv_object ~  `EM.Signature.status` , data = svvl)
ggsurvplot(fit1, data = svvl, pval = TRUE, palette = c(
  'EM.Signature.status=high'=alpha('#009E74',0.7),
  'EM.Signature.status=low'=alpha('#CB75A4',0.7)
  ),xlab="Survival time from 2nd resection")

surv_pvalue(fit1)
surv_pvalue(fit1)$pval


fit.cox <- coxph(surv_object  ~  `EM.Signature.status` , data = svvl)
ggforest(fit.cox)


ggsave("output/figures/survival_x_em_res_km.pdf",width=15/2,height=6)


decision <- 2.5


svvl <- plt.tight %>%
  dplyr::mutate(EM.Signature.status = as.factor(ifelse(extracellular.matrix.component.R1 > decision, "high", "low"))) %>%
  dplyr::mutate(event=1)

library(survival)
library(survminer)
surv_object <- Surv(time = svvl$survivalDays, event=svvl$event)
fit1 <- survfit(surv_object ~  `EM.Signature.status` , data = svvl)
ggsurvplot(fit1, data = svvl, pval = TRUE)

surv_pvalue(fit1)


ggsave("output/figures/survival_x_em_res_km_R1.png",width=6*1.25,height=6)






## ggplot x mgmt

ggplot(plt.tight, aes(x = mgmtRecurrent, y=extracellular.matrix.component.R2 )) + 
  geom_violin(draw_quantiles = c(0.5), col="black") +
  geom_jitter( position=position_jitter(0.2), size=0.9)


ggplot(plt.tight, aes(x = NF1, y=extracellular.matrix.component.R2 )) + 
  geom_violin(draw_quantiles = c(0.5), col="black") +
  geom_jitter( position=position_jitter(0.2), size=0.9)



# re-newed SVVL analysis [low/high] ----


components <- gsam.rna.metadata %>%
  dplyr::select(contains(".component") | contains(".signature") | 'sid' | 'pid') %>%
  dplyr::filter(!is.na(neuron.component)) %>% 
  dplyr::rename(C1.neuron.signature = neuron.component) %>% 
  dplyr::rename(C2.oligodendrocyte.signature = oligodendrocyte.component) %>% 
  dplyr::rename(C3.endothelial.signature = endothelial.component) %>% 
  dplyr::rename(C6.ECM.signature = extracellular.matrix.component)




components.paired <- dplyr::inner_join(
  components %>% 
    dplyr::filter(grepl("1", sid)) %>%
    `colnames<-`(paste0(colnames(.), ".R1")) %>%
    dplyr::rename(pid = pid.R1)
  ,
  components %>% 
    dplyr::filter(grepl("2", sid)) %>%
    `colnames<-`(paste0(colnames(.), ".R2")) %>%
    dplyr::rename(pid = pid.R2)
  , by=c('pid'='pid')) %>%
  dplyr::left_join(
    gsam.patient.metadata %>%dplyr::select(studyID, survivalDays, survivalFromSecondSurgeryDays),
    by=c('pid'='studyID')
  ) %>% 
  dplyr::mutate(event=1) %>% 
  dplyr::mutate(C6.ECM.signature.R2.status = factor(ifelse( C6.ECM.signature.R2 > 1.4, "high", "low"), levels=c("low","high")))%>% 
  dplyr::mutate(C5.signature.R2.status = factor(ifelse( C5.signature.R2 > 0, "high", "low"), levels=c("low","high"))) %>%  # plot(sort(components.paired$C6.ECM.signature.R2)) ; abline(h=1.4)
  dplyr::mutate(C4.signature.R2.status = factor(ifelse( C4.signature.R2 > -0.78, "high", "low"), levels=c("low","high"))) %>% # plot(sort(components.paired$C4.signature.R2)) ; abline(h=-0.78)
  dplyr::mutate(C3.endothelial.signature.R2.status = factor(ifelse( C3.endothelial.signature.R2 > 1.7, "high", "low"), levels=c("low","high"))) %>%  # plot(sort(components.paired$C3.endothelial.signature.R2)) ; abline(h=1.7)
  dplyr::mutate(C2.oligodendrocyte.signature.R2.status = factor(ifelse( C2.oligodendrocyte.signature.R2 > 4.8, "high", "low"), levels=c("low","high"))) %>%  # plot(sort(components.paired$C2.oligodendrocyte.signature.R2)) ; abline(h=4.8)
  dplyr::mutate(C1.neuron.signature.R2.status = factor(ifelse( C1.neuron.signature.R2 > 1, "high", "low"), levels=c("low","high"))) # plot(sort(components.paired$C1.neuron.signature.R2)) ; abline(h=1)




surv_object <- Surv(time = components.paired$survivalDays, event=components.paired$event)
fit1 <- survfit(surv_object ~  `C6.ECM.signature.R2.status` , data = components.paired)
ggsurvplot(fit1, data = components.paired, pval = TRUE)


surv_object <- Surv(time = components.paired$survivalDays, event=components.paired$event)
fit1 <- survfit(surv_object ~  `C5.signature.R2.status` , data = components.paired)
ggsurvplot(fit1, data = components.paired, pval = TRUE)



fit.cox <- coxph(surv_object ~ `C6.ECM.signature.R2.status`  , data = components.paired)
ggforest(fit.cox)

fit.cox <- coxph(surv_object ~ `C5.signature.R2.status`  , data = components.paired)
ggforest(fit.cox)

fit.cox <- coxph(surv_object ~ `C4.signature.R2.status`  , data = components.paired)
ggforest(fit.cox)

fit.cox <- coxph(surv_object ~ `C3.endothelial.signature.R2.status`  , data = components.paired)
ggforest(fit.cox)

fit.cox <- coxph(surv_object ~ `C2.oligodendrocyte.signature.R2.status` , data = components.paired)
ggforest(fit.cox)

fit.cox <- coxph(surv_object ~ `C2.oligodendrocyte.signature.R2.status` , data = components.paired)
ggforest(fit.cox)


fit.cox <- coxph(surv_object ~
                   `C6.ECM.signature.R2.status` + 
                   `C5.signature.R2.status` +
                   `C4.signature.R2.status` +
                   `C3.endothelial.signature.R2.status` +  `C2.oligodendrocyte.signature.R2.status` + `C1.neuron.signature.R2.status` , data = components.paired)
ggforest(fit.cox)


ggsave("output/figures/R2_survival_signatures.pdf",height=5.5,width=12)




plot( -components.paired$survivalFromSecondSurgeryDays , components.paired$C5.signature.R2)


plot( components.paired$C6.ECM.signature.R2 , components.paired$C5.signature.R2)
plot( components.paired$C6.ECM.signature.R2 , components.paired$C4.signature.R2)


plot( -components.paired$survivalFromSecondSurgeryDays , components.paired$C5.signature.R2)




