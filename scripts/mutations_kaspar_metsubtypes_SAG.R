#setwd("~/projects/gsam")

#### ---- load data ---- ####
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(cowplot)

source('~/projects/gsam/scripts/R/job_gg_theme.R')
source('~/projects/gsam/scripts/R/youri_gg_theme.R')
source('~/projects/gsam/scripts/R/palette.R')
source('scripts/R/gsam_rna-seq_expression.R') # recursively calls metadata
source('scripts/R/subtype_genes.R')
source('scripts/R/glass_expression_matrix.R')
source('data/wang/msig.library.12.R')  # no license w/ code provided, can't include it in source tree
source('scripts/R/wang_glioma_intrinsic_genes.R')
       
subtypes <- read.delim("output/tables/gsam_nmf_lda_data.txt", sep = "") 

dat <- merge(x = gsam.metadata, y = gsam.patient.metadata, by.x = "pid", by.y = "studyID", all.x = TRUE)

# dat <- gsam.patient.metadata %>%
#        dplyr::select(studyID,PTEN,PIK3C2G,PIK3CA,TSC2,MTOR,PIK3R1,TP53, TP53BP1, EGFR, ERBB3, FGFR3, IGF1R,
#                       EPHA3, KIT, NF1, MAP3K1, KRAS, KMT2D, KMT2C, ARID2, SETD2, ATM, BRCA2, MSH6, BRCA1,
#                       MSH2, RB1, BUB1B, JAK2, IDH1, IDH2, APC, AXIN2, mgmtStability,tertPrimary,nMutationsPrimary,nMutationsRecurrent,
#                       cnStatusEGFRs,viiiStability,daysToSecondSurgery, survivalFromSecondSurgeryDays, gender,treatedWithTMZ,treatedWithRT, HM,nMutationsPrimary)

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

#### ---- cnv plot  ---- ####
pid_no_pair <- dat %>% 
               dplyr::group_by(pid) %>%
               dplyr::summarise(n = n()) %>%
               dplyr::filter(n == 1) %>%
               dplyr::pull(pid)

dat_cnv <- dat %>%
           dplyr::filter(! pid %in% pid_no_pair) %>%
           dplyr::group_by(pid) %>% 
           dplyr::mutate(transition=length(unique(NMF.123456.PCA.SVM.class))) %>%
           dplyr::mutate(transition = ifelse(transition == 2, T, F)) %>% 
           dplyr::mutate(CL_stable=ifelse(transition==F & NMF.123456.PCA.SVM.class == "Classical" & resection == "r2",T,F)) %>%
           dplyr::mutate(MES_stable=ifelse(transition==F & NMF.123456.PCA.SVM.class == "Mesenchymal" & resection == "r2",T,F)) %>%
           dplyr::mutate(PN_stable=ifelse(transition==F & NMF.123456.PCA.SVM.class == "Proneural" & resection == "r2",T,F)) %>%
           dplyr::mutate(to_CL=ifelse(transition==T & resection == "r2" & NMF.123456.PCA.SVM.class == "Classical",T,F)) %>%
           dplyr::mutate(to_MES=ifelse(transition==T & resection == "r2" & NMF.123456.PCA.SVM.class == "Mesenchymal",T,F)) %>%
           dplyr::mutate(to_PN=ifelse(transition==T & resection == "r2" & NMF.123456.PCA.SVM.class == "Proneural",T,F)) %>%
           dplyr::select(pid,PTEN,PIK3C2G,PIK3CA,TSC2,MTOR,PIK3R1,TP53, TP53BP1, EGFR, ERBB3, FGFR3, IGF1R,
                         EPHA3, KIT, NF1, MAP3K1, KRAS, KMT2D, KMT2C, ARID2, SETD2, ATM, BRCA2, MSH6, BRCA1,
                         MSH2, RB1, BUB1B, JAK2, IDH1, IDH2, APC, AXIN2, cnStatusCDKN2ABs, transition,CL_stable,MES_stable,PN_stable,to_CL,to_MES,to_PN) 

ordering <- dat_cnv %>%
            dplyr::select(pid,transition,to_CL,to_PN,to_MES)

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

#ggsave("output/figures/figures_SAG/subtype_CNV_pattern1.pdf",width=10 ,height=6*2)

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

