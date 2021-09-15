#!/usr/bin/env R 

# load libs ----

library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(cowplot)


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


      

subtypes <- read.delim("output/tables/gsam_nmf_lda_data.txt", sep = "") 

dat <- merge(x = gsam.rna.metadata, y = gsam.patient.metadata, by.x = "pid", by.y = "studyID", all.x = TRUE)



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



## plot x neuronal & oligodendrocyte scores etc. ----



plt.tight <- gsam.patient.metadata
  
  



