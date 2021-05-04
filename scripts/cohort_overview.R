#!/usr/bin/env R

library(ggplot2)

# load data ----

source("scripts/R/palette.R")

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


source("scripts/R/ligands.R")
source("scripts/R/subtype_genes.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

source('scripts/R/wang_glioma_intrinsic_genes.R')

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set


# G-SAM ----

# we hebben max 122 pairs?

#plt.expanded %>% dplyr::filter(pid == "GAO")
#GAO = duplicated !?

plt.ids <- gsam.rna.metadata %>%
    dplyr::filter(blacklist.pca == F) %>%
    #dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
    dplyr::select(c('pid', 'sid')) %>%
    dplyr::mutate(resection = gsub("^...(.).*$","R\\1",sid)) %>%
    dplyr::mutate(pid = as.factor(as.character(pid))) # refactor


# maak eerst losse colommen en rbind die daarna

plt.expanded <- data.frame(pid = unique(plt.ids$pid) ) %>%
  dplyr::left_join(gsam.patient.metadata %>% dplyr::select(c('studyID','treatedWithTMZ')) , by = c('pid' = 'studyID')) %>%
  dplyr::left_join(gsam.patient.metadata %>% dplyr::select(c('studyID','treatedWithRT')) , by = c('pid' = 'studyID')) %>%
  dplyr::left_join(plt.ids %>% dplyr::filter(resection == "R1") %>% dplyr::rename(R1 = sid) %>% dplyr::select(c('pid', 'R1')) , by = c('pid' = 'pid') ) %>%
  dplyr::left_join(plt.ids %>% dplyr::filter(resection == "R2") %>% dplyr::rename(R2 = sid) %>% dplyr::select(c('pid', 'R2')) , by = c('pid' = 'pid') ) %>%
  dplyr::mutate(R1.status =  ifelse(is.na(R1), "NA", "+") ) %>%
  dplyr::mutate(R2.status =  ifelse(is.na(R2), "NA", "+") ) %>%
  dplyr::mutate(R1.R2.status = case_when( is.na(R1) ~ 'R2', is.na(R2) ~ 'R1',     TRUE ~ 'both'   ) ) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::filter(sid %in% plt.ids$sid & resection == "r1") %>% dplyr::select(c('pid','NMF.123456.PCA.LDA.class')) %>% dplyr::rename(R1.rna.subtype = NMF.123456.PCA.LDA.class ),  by = c('pid' = 'pid') ) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::filter(sid %in% plt.ids$sid & resection == "r2") %>% dplyr::select(c('pid','NMF.123456.PCA.LDA.class')) %>% dplyr::rename(R2.rna.subtype = NMF.123456.PCA.LDA.class ),  by = c('pid' = 'pid') ) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::filter(sid %in% plt.ids$sid & resection == "r1") %>% dplyr::select(c('pid','pat.with.IDH')) %>% dplyr::rename(patient.idh.status = pat.with.IDH ),  by = c('pid' = 'pid') ) %>%
  dplyr::mutate(patient.idh.status = ifelse(patient.idh.status == "TRUE", "+", "-") ) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::filter(sid %in% plt.ids$sid &resection == "r1") %>% dplyr::select(c('pid','tumour.percentage.dna')) %>% dplyr::rename(R1.tumour.percentage.dna = tumour.percentage.dna )
    ,  by = c('pid' = 'pid')) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::filter(sid %in% plt.ids$sid &resection == "r2") %>% dplyr::select(c('pid','tumour.percentage.dna')) %>% dplyr::rename(R2.tumour.percentage.dna = tumour.percentage.dna )
    ,  by = c('pid' = 'pid')) %>%
  dplyr::mutate(R1.tumour.percentage.status = ifelse(R1.tumour.percentage.dna >= 15,'+','-') ) %>%
  dplyr::mutate(R2.tumour.percentage.status = ifelse(R2.tumour.percentage.dna >= 15,'+','-') ) %>%
  dplyr::mutate(R1.tumour.percentage.status.order = ifelse(R1.tumour.percentage.dna >= 15,1,2) ) %>%
  dplyr::mutate(R2.tumour.percentage.status.order = ifelse(R2.tumour.percentage.dna >= 15,1,2) ) %>%
  dplyr::mutate(treatedWithTMZ = ifelse(treatedWithTMZ == "Yes",'+','-') ) %>%
  dplyr::mutate(treatedWithRT = ifelse(treatedWithRT == "Yes",'+','-') ) %>%
  dplyr::arrange(patient.idh.status, R1.R2.status, R2.tumour.percentage.status.order,  R1.tumour.percentage.status.order, R1.rna.subtype, R2.rna.subtype) %>%
  #dplyr::arrange(patient.idh.status, R1.R2.status, R2.tumour.percentage.status.order,  R1.tumour.percentage.status.order) %>%
  dplyr::mutate(order = 1:nrow(.) )
stopifnot(nrow(plt.expanded) == 185) # ensure what gets in comes out


#plt.expanded %>% dplyr::filter(pid == "GAO")



plt <- bind_rows(
    plt.expanded %>% dplyr::select(c('pid', 'R1.status')) %>% dplyr::rename(col = R1.status) %>% dplyr::mutate(y = "RNA-Seq sample R1") ,
    plt.expanded %>% dplyr::select(c('pid', 'R2.status')) %>% dplyr::rename(col = R2.status) %>% dplyr::mutate(y = "RNA-Seq sample R2") ,
    plt.expanded %>% dplyr::select(c('pid', 'R1.rna.subtype')) %>% dplyr::rename(col = R1.rna.subtype) %>% dplyr::mutate(y = "Transcriptional subtype R1") ,
    plt.expanded %>% dplyr::select(c('pid', 'R2.rna.subtype')) %>% dplyr::rename(col = R2.rna.subtype) %>% dplyr::mutate(y = "Transcriptional subtype R2") ,
    plt.expanded %>% dplyr::select(c('pid', 'patient.idh.status')) %>% dplyr::rename(col = patient.idh.status) %>% dplyr::mutate(y = "Patient IDH status") ,
    plt.expanded %>% dplyr::select(c('pid', 'R1.tumour.percentage.status')) %>% dplyr::rename(col = R1.tumour.percentage.status) %>% dplyr::mutate(y = "Tumor cells >= 15% (WES) R1") ,
    plt.expanded %>% dplyr::select(c('pid', 'R2.tumour.percentage.status')) %>% dplyr::rename(col = R2.tumour.percentage.status) %>% dplyr::mutate(y = "Tumor cells >= 15% (WES) R2") ,
    
    plt.expanded %>% dplyr::select(c('pid', 'treatedWithTMZ')) %>% dplyr::rename(col = treatedWithTMZ) %>% dplyr::mutate(y = "Treatment: TMZ") ,
    plt.expanded %>% dplyr::select(c('pid', 'treatedWithRT')) %>% dplyr::rename(col = treatedWithRT) %>% dplyr::mutate(y = "Treatment: Radio therapy")
  ) %>%
  dplyr::mutate(col = ifelse(is.na(col),'NA',col)) %>%
  dplyr::left_join(plt.expanded %>% dplyr::select('pid','order'), by=c('pid'='pid')) 
  


ggplot(plt, aes(x = reorder(pid, order), y = y, fill=col)) +
  geom_tile(col="black") + 
  scale_fill_manual(values = c('+' = 'gray40', '-' = 'red', 'NA' = 'white', subtype_colors)) +
  job_gg_theme + 
  theme(axis.text.x = element_text(angle = 90, size = 5)) + 
  labs(y=NULL, x=NULL)


ggsave("output/figures/cohort_overview_gsam.png",width=11,height=3)



# GLASS ----


c1  <- c("GLSS-SM-R068-TP-01R-RNA-0UPMYO", "GLSS-SM-R068-R1-01R-RNA-7I5H9P",
   "GLSS-SM-R071-TP-01R-RNA-OAXGI8", "GLSS-SM-R071-R1-01R-RNA-7AZ6G2",
   "GLSS-SM-R099-R1-01R-RNA-MNTPMI", 
   "GLSS-SM-R100-TP-01R-RNA-EUI7AZ", "GLSS-SM-R100-R1-01R-RNA-46UW5U")

c2 <- gsub(".","-",colnames(read.table("data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt",sep="\t",header=T,stringsAsFactors = F) %>% dplyr::mutate("target_id"=NULL,"length"=NULL) ) ,fixed=T)
c3 <- glass.gbm.rnaseq.metadata$sid



plt.ids <- data.frame( sid = unique(sort(c(c1,c2,c3))) ) %>%
  dplyr::mutate(pid = as.factor(gsub("^([^\\-]+.[^\\-]+.[^\\-]+).*$","\\1",sid))) %>%
  dplyr::mutate(resection = as.factor(gsub("^[^\\-]+.[^\\-]+.[^\\-]+.([^\\-]+).*$","\\1",sid)))
stopifnot(duplicated(plt.ids$sid) == F)



plt.expanded <- data.frame(pid = unique(plt.ids$pid) ) %>%
  dplyr::mutate(dataset = as.factor(gsub("^([^\\-]+).+","\\1", pid) ) ) %>%
  dplyr::left_join(plt.ids %>% dplyr::filter(resection == "TP") %>% dplyr::rename(TP = sid) %>% dplyr::select(c('pid', 'TP')) , by = c('pid' = 'pid') ) %>%
  dplyr::left_join(plt.ids %>% dplyr::filter(resection == "R1") %>% dplyr::rename(R1 = sid) %>% dplyr::select(c('pid', 'R1')) , by = c('pid' = 'pid') ) %>%
  dplyr::left_join(plt.ids %>% dplyr::filter(resection == "R2") %>% dplyr::rename(R2 = sid) %>% dplyr::select(c('pid', 'R2')) , by = c('pid' = 'pid') ) %>%
  dplyr::left_join(plt.ids %>% dplyr::filter(resection == "R3") %>% dplyr::rename(R3 = sid) %>% dplyr::select(c('pid', 'R3')) , by = c('pid' = 'pid') ) %>%
  
  dplyr::mutate(TP.status =  case_when(is.na(TP) ~ 'NA', TP %in% colnames(glass.gbm.rnaseq.expression) ~ "+" , TRUE ~ "-") ) %>%
  dplyr::mutate(R1.status =  case_when(is.na(R1) ~ 'NA', R1 %in% colnames(glass.gbm.rnaseq.expression) ~ "+" , TRUE ~ "-") ) %>%
  dplyr::mutate(R2.status =  case_when(is.na(R2) ~ 'NA', R2 %in% colnames(glass.gbm.rnaseq.expression) ~ "+" , TRUE ~ "-") ) %>%
  dplyr::mutate(R3.status =  case_when(is.na(R3) ~ 'NA', R3 %in% colnames(glass.gbm.rnaseq.expression) ~ "+" , TRUE ~ "-") ) %>%
  dplyr::mutate(np = (TP.status == '+') + (R1.status == '+') + (R2.status == '+') + (R3.status == '+') ) %>%
  dplyr::mutate(IDHwt = pid %in% (glass.gbm.rnaseq.metadata %>% dplyr::filter(idh_status == "IDHwt") %>% dplyr::pull(pid)) ) %>%
  dplyr::mutate(IDHwt.order = ifelse(IDHwt, 1, 2)) %>%
  dplyr::mutate(IDHwt = ifelse(IDHwt, "+", "-")) %>%
  dplyr::mutate(Grade.II = pid %in% (glass.gbm.rnaseq.metadata %>% dplyr::filter(grade == "II") %>% dplyr::pull(pid)) ) %>%
  dplyr::mutate(Grade.III = pid %in% (glass.gbm.rnaseq.metadata %>% dplyr::filter(grade == "III") %>% dplyr::pull(pid)) ) %>%
  dplyr::mutate(Grade.IV = pid %in% (glass.gbm.rnaseq.metadata %>% dplyr::filter(grade == "IV") %>% dplyr::pull(pid)) ) %>%
  dplyr::mutate(R1.R2.Grade.IV = ifelse(Grade.IV & !Grade.III & !Grade.II,"+","-")) %>%
  dplyr::mutate(R1.R2.Grade.IV.order = ifelse(Grade.IV & !Grade.III & !Grade.II,1,2)) %>%
  dplyr::left_join(glass.metadata.all %>%
      dplyr::filter(treatment_tmz %in% c('t','f') ) %>%
      dplyr::select(c("case_barcode","treatment_tmz")) %>%
      dplyr::distinct() ,  by=c('pid'='case_barcode')) %>%
  dplyr::mutate(treatment_tmz = ifelse(treatment_tmz == "t", '+', treatment_tmz) ) %>%
  dplyr::mutate(treatment_tmz = ifelse(treatment_tmz == "f", '-', treatment_tmz) ) %>%
  dplyr::arrange(dataset, IDHwt.order, R1.R2.Grade.IV.order, -np, R3.status, R2.status, R1.status) %>%
  dplyr::mutate(order = 1:nrow(.))


# add:
# Transcriptional subtype R1
# Transcriptional subtype R2-4




plt <- bind_rows(
    plt.expanded %>% dplyr::select(c('pid', 'TP.status')) %>% dplyr::rename(col = TP.status) %>% dplyr::mutate(y = "RNA-Seq sample R1 (TP)") ,
    plt.expanded %>% dplyr::select(c('pid', 'R1.status')) %>% dplyr::rename(col = R1.status) %>% dplyr::mutate(y = "RNA-Seq sample R2") ,
    plt.expanded %>% dplyr::select(c('pid', 'R2.status')) %>% dplyr::rename(col = R2.status) %>% dplyr::mutate(y = "RNA-Seq sample R3") ,
    plt.expanded %>% dplyr::select(c('pid', 'R3.status')) %>% dplyr::rename(col = R3.status) %>% dplyr::mutate(y = "RNA-Seq sample R4") ,
    plt.expanded %>% dplyr::select(c('pid', 'dataset')) %>% dplyr::rename(col = dataset) %>% dplyr::mutate(y = "GLASS Project barcode") ,
    plt.expanded %>% dplyr::select(c('pid', 'IDHwt')) %>% dplyr::rename(col = IDHwt) %>% dplyr::mutate(y = "IDH wildtype") ,
    plt.expanded %>% dplyr::select(c('pid', 'R1.R2.Grade.IV')) %>% dplyr::rename(col = R1.R2.Grade.IV) %>% dplyr::mutate(y = "R1 (TP) & R2 Grade IV (GBM)"),
    plt.expanded %>% dplyr::select(c('pid', 'treatment_tmz')) %>% dplyr::rename(col = treatment_tmz) %>% dplyr::mutate(y = "TMZ Treatment")
  ) %>%
  dplyr::mutate(col = ifelse(is.na(col),'NA',col)) %>%
  dplyr::left_join(plt.expanded %>% dplyr::select('pid','order'), by=c('pid'='pid')) 




ggplot(plt, aes(x = reorder(pid, order), y = y, fill=col)) +
  geom_tile(col="black") + 
  scale_fill_manual(values = c('+' = 'gray40', '-' = 'red', 'NA' = 'white', subtype_colors, 'TCGA'= "#6ba6e5" , 'GLSS'=  "#eab509" )) +
  job_gg_theme + 
  theme(axis.text.x = element_text(angle = 90, size = 5)) + 
  labs(y=NULL, x=NULL)



ggsave("output/figures/cohort_overview_glass.png",width=11,height=3)



