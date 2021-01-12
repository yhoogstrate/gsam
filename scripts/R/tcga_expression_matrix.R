#!/usr/bin/env R 

# ---- load libs ----

library(tidyverse)


# ---- load TCGA identifiers ----


tcga.identifiers.gbm <- read.delim('data/tcga/tables/rna-seq/gdc_sample_sheet.2019-01-28.tsv',stringsAsFactors = F) %>%
  dplyr::mutate( tcga.htseq.id = paste0('tcga.',gsub("-",".",gsub(".htseq.counts.gz","",File.Name,fixed=T))) ) %>%
  dplyr::mutate( Data.Category = NULL ) %>%
  dplyr::mutate(Data.Type = NULL) %>%
  dplyr::mutate(TCGA.pat.ID = gsub("^(....)-[^-]+-[0]*([^-]+).+$","\\1.\\2",Sample.ID)) %>%
  dplyr::filter(tcga.htseq.id != "tcga.18fbf794.a85b.4b0d.9a5c.c1c0e0ede3d8") %>% # replicate of TCGA-06-0156
  dplyr::filter(tcga.htseq.id != "tcga.f67cf68d.56b3.451d.ba21.566664213691") %>% # replicate of TCGA-06-0211
  dplyr::filter((Case.ID %in% c("TCGA-06-0211", "TCGA-06-0210", "TCGA-14-1034", "TCGA-06-0125", "TCGA-06-0190", "TCGA-19-4065" ) &
                   Sample.Type == "Recurrent Tumor") == F ) # From these, exclude the secondary - as both prim and secondary are present: "TCGA-06-0211" "TCGA-06-0210" "TCGA-14-1034" "TCGA-06-0125" "TCGA-06-0190" "TCGA-19-4065"
  #dplyr::filter(Sample.Type == "Primary Tumor") %>% Recurrent tumours were also used?


# ---- load Wang / GlioVis paper data ----


TCGA_GBM <- readRDS('data/gliovis/data/datasets/TCGA_GBM.Rds')
TCGA_GBM$sample_rnaseq <- TCGA_GBM$rseq %>%
  dplyr::mutate(Sample = gsub('\\.', '-', Sample, fixed=F )) %>%
  dplyr::select(c('EGFR', 'Sample', 'Histology', 'Recurrence')) %>%
  dplyr::mutate(Sample.RNA = case_when( # kan niet met NA values dealen
    is.na(EGFR) ~ "Not RNA-sequenced" ,
    Histology != "GBM" ~ "No GBM",
    Sample %in% tcga.identifiers.gbm$Case.ID == F ~ "no identifier",
    TRUE ~ "Match")) 
TCGA_GBM$sample_rnaseq %>% dplyr::filter(Sample.RNA %in% c('Not RNA-sequenced', 'Match') == F )



TCGA_GBM <- readRDS('data/gliovis/data/datasets/TCGA_GBM.Rds')
TCGA_GBM$sample_rnaseq <- TCGA_GBM$rseq %>%
  dplyr::mutate(Sample = gsub('\\.', '-', Sample, fixed=F )) %>%
  dplyr::select(c('EGFR', 'Sample', 'Histology', 'Recurrence')) %>%
  dplyr::mutate(Sample.RNA = case_when( # kan niet met NA values dealen
    is.na(EGFR) ~ "NA" ,
    Histology != "GBM" ~ "NA",
    Sample %in% tcga.identifiers.gbm$Case.ID == F ~ "NA",
    TRUE ~ Sample)) %>%
  dplyr::mutate(Sample.RNA = ifelse(Sample.RNA == "NA" , NA , Sample.RNA)) %>%
  dplyr::select(Sample, Sample.RNA)

# Not included:
# 3   TCGA-06-0675 Non-tumor        No GBM
# 4   TCGA-06-0678 Non-tumor        No GBM
# 5   TCGA-06-0680 Non-tumor        No GBM
# 6   TCGA-06-0681 Non-tumor        No GBM
# 7   TCGA-06-5415       GBM no identifier  not in TCGA?


# Included - but are actually recurrent samples (!) :
# 1   TCGA-06-0152       GBM no identifier  recurrent?
# 2   TCGA-06-0171       GBM no identifier  recurrent?
# 8   TCGA-14-0736       GBM no identifier  recurrent?
# 9   TCGA-14-1402       GBM no identifier  recurrent?
# 10  TCGA-19-0957       GBM no identifier  recurrent?
# 11  TCGA-19-1389       GBM no identifier  recurrent?


stopifnot( (TCGA_GBM$sample_rnaseq %>%
  dplyr::filter(!is.na(Sample.RNA)) %>%
  dplyr::pull(Sample.RNA)) %in% tcga.identifiers.gbm$Case.ID)


# ---- load actual TCGA counts and match w/ Wang ----


tcga.expression.gbm <-read.delim('data/tcga/tables/rna-seq/merged-expression_TCGA-GBM.htseq.counts.txt') %>%
  `colnames<-`(gsub("^X","tcga.",colnames(.))) %>%
  `colnames<-`(gsub("^([a-f])","tcga.\\1",colnames(.))) %>%
  `colnames<-`(gsub(".htseq.counts.gz","",colnames(.),fixed=T)) %>%
  dplyr::filter(gene.id %in% c('__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned', '__too_low_aQual') == F) %>%
  tibble::column_to_rownames('gene.id') %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('sample.id') %>%
  dplyr::left_join(tcga.identifiers.gbm %>% dplyr::select(c('tcga.htseq.id', "Case.ID" )) , by=c('sample.id'='tcga.htseq.id')) %>%
  dplyr::filter(!is.na(Case.ID)) %>%
#a = tcga.expression.gbm  %>%
  #dplyr::select('sample.id' , 'Case.ID', 'ENSG00000281884.1','ENSG00000281887.1') %>%
  dplyr::mutate(sample.id = NULL) %>%
  tibble::column_to_rownames('Case.ID') %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) %>%
  dplyr::select( TCGA_GBM$sample_rnaseq$Sample.RNA %>% purrr::discard(is.na)  ) # only take those also used in gliovis
  

stopifnot(TCGA_GBM$sample_rnaseq$Sample.RNA %>% purrr::discard(is.na) == colnames(tcga.expression.gbm))



# ---- translate ENSG to symbols ----




