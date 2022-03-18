#!/usr/bin/env R

# load libs ----


library(ggplot2)
library(ggrepel)



# load data ----

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


source("scripts/R/gsam_metadata.R")

# prepare data ----



plt <- gsam.patient.metadata %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::filter(batch != "old" & resection == 'r1') %>%
      dplyr::filter(blacklist.pca == F) %>%
      dplyr::filter(pat.with.IDH == F) %>%
      dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
      dplyr::select(c('pid','tumour.percentage.dna')) %>%
      dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna)
    , by = c('studyID'='pid') ) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::filter(batch != "old" & resection == 'r2') %>%
      dplyr::filter(blacklist.pca == F) %>%
      dplyr::filter(pat.with.IDH == F) %>%
      dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
      dplyr::select(c('pid','tumour.percentage.dna')) %>%
      dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) ,
    by = c('studyID'='pid') ) %>%
  dplyr::filter(!is.na(tumour.percentage.dna.R1) & !is.na(tumour.percentage.dna.R2) ) %>%
  dplyr::mutate(tpc.lfc = log(tumour.percentage.dna.R2 / tumour.percentage.dna.R1))  %>%
  dplyr::mutate(extentOfResectionFirstSurgery = ifelse(grepl("Complete resection",extentOfResectionFirstSurgery) , "Complete resection" , extentOfResectionFirstSurgery)  ) %>%
  dplyr::mutate(extentOfResectionSecondSurgery = ifelse(grepl("Complete resection",extentOfResectionSecondSurgery) , "Complete resection" , extentOfResectionSecondSurgery)  ) %>%
  dplyr::mutate(extentOfResectionState = paste0('R1:',extentOfResectionFirstSurgery, '->R2:' , extentOfResectionSecondSurgery) ) %>%
  dplyr::mutate(change.extend.of.resection = extentOfResectionFirstSurgery != extentOfResectionSecondSurgery )




ggplot(plt, aes(x = reorder(studyID, tpc.lfc) , y = tpc.lfc ) ) + 
  geom_point() +
  geom_point(aes(y = -3.1, col=extentOfResectionFirstSurgery )) + 
  geom_point(aes(y = -3.3, col=extentOfResectionSecondSurgery )) + 
  geom_point(aes(y = -3.5, shape=change.extend.of.resection  )) + 
  geom_point(aes(y = -3.7, fill=tumorLocation  )  , col = rgb(1,0,1,0) , pch=21) 



ggplot(plt , aes(x = extentOfResectionFirstSurgery, y = tpc.lfc)) +
  geom_violin() +
  geom_jitter( position=position_jitter(0.2), size=0.9)
  
  
ggplot(plt , aes(x = extentOfResectionSecondSurgery, y = tpc.lfc)) +
  geom_violin() +
  geom_jitter( position=position_jitter(0.2), size=0.9)


ggplot(plt , aes(x = extentOfResectionState, y = tpc.lfc)) +
  geom_violin() +
  geom_jitter( position=position_jitter(0.2), size=0.9) + 
  job_gg_theme


ggplot(plt , aes(x = tumorLocation, y = tpc.lfc)) +
  geom_violin() +
  geom_jitter( position=position_jitter(0.2), size=0.9) + 
  job_gg_theme



