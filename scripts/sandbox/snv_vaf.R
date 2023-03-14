#setwd("~/project/G-SAM")

library(tidyverse)

if(!exists("gsam.cnv.metadata")) {
  source('scripts/R/gsam_metadata.R')
}

available_vcf_data <- read.delim("data/DNA/GBM_allSamples_PassVariants_1736_GBM_TERT_05Jan2018_flags_VAFs.txt") %>%
                      dplyr::inner_join(gsam.cnv.metadata[,c("pid","sid","donor_ID")],.,by=c("sid"="Sample")) %>%
                      dplyr::pull(donor_ID)
available_vcf_data <- unique(available_vcf_data)

vaf <- read.delim("data/DNA/GBM_allSamples_PassVariants_1736_GBM_TERT_05Jan2018_flags_VAFs.txt") %>%
       dplyr::inner_join(gsam.cnv.metadata[,c("pid","sid","donor_ID")],.,by=c("sid"="Sample")) %>% 
       dplyr::filter(Gene == "EGFR" & Effect == "missense") %>% 
       dplyr::select(pid,donor_ID,VD,Effect) %>% 
       separate(VD, c("V1", "V2","V3","V4","missense_egfr","V6","V7"), "\\|") %>% 
       dplyr::select(pid,donor_ID,missense_egfr) %>% 
       dplyr::mutate(missense_egfr = gsub("*p\\.","",missense_egfr)) %>% 
       dplyr::group_by(donor_ID) %>% 
       dplyr::mutate(egfr_snv = paste0(missense_egfr, collapse = ";")) %>% 
       dplyr::select(-missense_egfr) %>%
       dplyr::distinct(donor_ID, .keep_all = TRUE) 

vaf$A289 <- grepl("A289V|A289D|A289T",vaf$egfr_snv)
vaf$G598 <- grepl("G598V|G598A",vaf$egfr_snv)
vaf$R108 <- grepl("R108G|R108K",vaf$egfr_snv)
vaf$resection <- ifelse(substr(vaf$donor_ID,4,4)=="1","r1","r2")

vaf  <- vaf %>% 
        dplyr::group_by(pid) %>% 
        dplyr::mutate(stable=ifelse(length(pid)==2,T,F)) %>%
        dplyr::mutate(gained=ifelse(stable==F & resection == "r2",T,F)) %>%
        dplyr::mutate(lost=ifelse(stable==F & resection == "r1",T,F)) %>%
        dplyr::mutate(trans = dplyr::case_when((stable == T) ~ "stable" , 
                      (stable == F & gained == T) ~ "gained" , 
                      (stable == F & lost == T) ~ "lost" ))