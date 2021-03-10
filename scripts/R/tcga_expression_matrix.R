#!/usr/bin/env R

# load libs ----

library(tidyverse)

# three variables:
# 
# tcga.identifiers.gbm
# tcga.gbm.rnaseq.expression


# check idientifiers ----



tcga.identifiers.gbm <- read.delim('data/tcga/tables/rna-seq/gdc_sample_sheet.2019-01-28.tsv',stringsAsFactors = F) %>%
  dplyr::mutate(tcga.htseq.id = paste0('tcga.',gsub("-",".",gsub(".htseq.counts.gz","", File.Name,fixed=T)))) %>%
  dplyr::mutate(Data.Category = NULL) %>%
  dplyr::mutate(Data.Type = NULL) %>%
  dplyr::filter(Project.ID == "TCGA-GBM") %>%
  dplyr::mutate(TCGA.pat.ID = gsub("^(....)-[^-]+-[0]*([^-]+).+$","\\1.\\2", Sample.ID)) 



# obtain expression data ----


tcga.gbm.rnaseq.expression <- read.delim('data/tcga/tables/rna-seq/merged-expression_TCGA-GBM.htseq.counts.txt') %>%
  `colnames<-`(gsub("X","tcga.", colnames(.))) %>%
  `colnames<-`(gsub("^([a-f])","tcga.\\1", colnames(.))) %>%
  `colnames<-`(gsub(".htseq.counts.gz","", colnames(.), fixed=T)) %>%
  dplyr::filter(grepl("^ENSG", gene.id) ) %>%
  dplyr::mutate(gene.id =  gsub("^([A-Z0-9]+).*$","\\1", gene.id) ) %>%
  tibble::column_to_rownames('gene.id')



# pick only those w/ matching R1 & R2 ----


recurrent.cases <- tcga.identifiers.gbm %>%
  dplyr::filter(Sample.Type == "Recurrent Tumor") %>%
  dplyr::pull(Case.ID)


# from glass portal:
# "TCGA-06-0125-TP-01R-RNA-RCQ5QS" "TCGA-06-0125-R1-11R-RNA-EKH06U" "TCGA-06-0190-TP-01R-RNA-S2074S" "TCGA-06-0190-R1-01R-RNA-HWTLG4"
# "TCGA-06-0210-TP-01R-RNA-JMBENW" "TCGA-06-0210-R1-01R-RNA-WB10MK" "TCGA-06-0211-TP-01R-RNA-Y9DKHR" "TCGA-06-0211-R1-02R-RNA-3YQ1G4"
# "TCGA-14-1034-TP-01R-RNA-S6SBZH" "TCGA-14-1034-R1-01R-RNA-GSNR8T" "TCGA-DH-A669-TP-12R-RNA-XG7CTS" "TCGA-DH-A669-R1-11R-RNA-KJQZGQ"
# "TCGA-DU-5870-TP-11R-RNA-CSQWQZ" "TCGA-DU-5870-R1-12R-RNA-GIZUWN" "TCGA-DU-5872-TP-11R-RNA-HS9DNI" "TCGA-DU-5872-R1-21R-RNA-0EW2ZH"
# "TCGA-DU-6397-TP-11R-RNA-NYNT65" "TCGA-DU-6397-R1-12R-RNA-N760C1" "TCGA-DU-6404-TP-11R-RNA-A7XMXJ" "TCGA-DU-6404-R1-21R-RNA-LJZDP3"
# "TCGA-DU-6407-TP-13R-RNA-6NHIOE" "TCGA-DU-6407-R1-12R-RNA-T3UMCW" "TCGA-DU-7304-TP-12R-RNA-23D1EM" "TCGA-DU-7304-R1-12R-RNA-OMXRQS"
# "TCGA-FG-5963-TP-11R-RNA-2Y6I4P" "TCGA-FG-5963-R1-12R-RNA-1EJ1SV" "TCGA-FG-5965-TP-11R-RNA-0XBPNY" "TCGA-FG-5965-R1-11R-RNA-QSA7Y3"
# "TCGA-FG-A4MT-TP-11R-RNA-K4BQG0" "TCGA-FG-A4MT-R1-11R-RNA-9OU9UC" "TCGA-TM-A7CF-TP-11R-RNA-YLP5PU" "TCGA-TM-A7CF-R1-11R-RNA-EZWWA2"
# "TCGA-TQ-A7RK-TP-11R-RNA-LRGF1D" "TCGA-TQ-A7RK-R1-11R-RNA-PKFXB3" "TCGA-TQ-A7RV-TP-21R-RNA-6F23HF" "TCGA-TQ-A7RV-R1-11R-RNA-TIRW8I"
# "TCGA-TQ-A8XE-TP-11R-RNA-0YPLAZ" "TCGA-TQ-A8XE-R1-11R-RNA-47VZAQ"

# from TCGA
# "TCGA-06-0211" "TCGA-19-4065" "TCGA-06-0210" "TCGA-19-1389" "TCGA-14-1402" "TCGA-06-0152" "TCGA-06-0221" "TCGA-06-0171" "TCGA-19-0957"
# "TCGA-14-1034" "TCGA-06-0125" "TCGA-14-0736" "TCGA-06-0190"

all.samples.of.recurrent.cases <- tcga.identifiers.gbm %>%
  dplyr::filter( Case.ID %in% recurrent.cases ) %>%
  dplyr::pull(tcga.htseq.id)



rm(recurrent.cases)
all.samples.of.recurrent.cases[all.samples.of.recurrent.cases  %in% colnames(tcga.gbm.rnaseq.expression ) == T]


tmp <- tcga.gbm.rnaseq.expression %>%
  dplyr::filter(rowSums(.) >= (ncol(.) * 1) ) 

cond <- as.factor(rep(c("A","B"),999) [1:ncol(tmp)] )
tmp <- DESeq2::DESeqDataSetFromMatrix(tmp, S4Vectors::DataFrame(cond), ~cond)
tcga.gbm.rnaseq.expression.vst <- SummarizedExperiment::assay(DESeq2::vst(tmp, blind=T))
rm(cond, tmp)


# take only those w/ prim and sec
tcga.gbm.rnaseq.expression.vst <- tcga.gbm.rnaseq.expression.vst %>% 
  as.data.frame() %>%
  dplyr::select(all.samples.of.recurrent.cases)






