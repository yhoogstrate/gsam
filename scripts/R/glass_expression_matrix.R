#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(DESeq2)


## read counts from GLASS samples ----

# possibly some match w/ Wang dataset?


# The samples are taken from a data portal:
# www.synapse.org

# expression values are deterimined by Kallisto, per transcript


# hacking this back from GLASS git code refers to gencode v19...
# gencode.v19.chr_patch_hapl_scaff.annotation.gtf << almost no mismatches [1193]


glass.gencode.v19 <- 'data/gsam/data/GLASS_GBM_R1-R2/gencode.v19.chr_patch_hapl_scaff.annotation.translate-table.txt' %>%
  read.table(header=F, stringsAsFactors = F) %>%
  dplyr::rename(gene_symbol=V1) %>%
  dplyr::rename(gene_id=V2) %>%
  dplyr::rename(transcript_id=V3)


glass.gbm.rnaseq.expression <- 'data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt' %>%
  read.delim(stringsAsFactors = F) %>%
  dplyr::mutate(length = NULL) %>% # not needed
  dplyr::left_join(glass.gencode.v19, by=c('target_id' = 'transcript_id')) %>%
  dplyr::filter(!is.na(gene_symbol) ) %>% # 1193 transcript id's not matching gtf Ensembl 64
  dplyr::mutate(target_id = NULL) %>% # aggregate @ gene_id level
  dplyr::mutate(gene_symbol = NULL) %>%
  dplyr::group_by(gene_id) %>%
  summarise(across(everything(), list(sum))) %>%
  tibble::rownames_to_column('tmp') %>% 
  dplyr::mutate(tmp=NULL) %>%
  #dplyr::filter(gene_id %in% c("ENSG00000198804", "ENSG00000198886", "ENSG00000198938","ENSG00000198712", # exclude from normalisation?
  #                             "ENSG00000198727") == F ) %>% # extreme high counts from odd genes
  tibble::column_to_rownames('gene_id') %>%
  round() %>%
  `colnames<-`( gsub(".","-",colnames(.), fixed=T) ) %>%
  `colnames<-`( gsub("_1$","",colnames(.), fixed=F) ) %>%
  dplyr::select( # extremely strong outlier samples !!
    -c("GLSS-SM-R068-TP-01R-RNA-0UPMYO", "GLSS-SM-R068-R1-01R-RNA-7I5H9P",
       "GLSS-SM-R071-TP-01R-RNA-OAXGI8", "GLSS-SM-R071-R1-01R-RNA-7AZ6G2",
       "GLSS-SM-R099-R1-01R-RNA-MNTPMI", #"GLSS-SM-R099-TP-01R-RNA-YGXA72",
       "GLSS-SM-R100-TP-01R-RNA-EUI7AZ", "GLSS-SM-R100-R1-01R-RNA-46UW5U"
    ))



  # Using the following will pick the max() transcript per gene
  # I tried this as well in the hope that it would remove the batch
  # effect but it didn't:'
  # 
  # dplyr::arrange(gene_id) %>%
  # as.data.frame() %>%
  # mutate(rs = rowSums(dplyr::select( ., !dplyr::contains("gene_id") ))) %>%
  # tibble::as_tibble() %>%
  # dplyr::group_by(gene_id) %>%
  # top_n(n=1, wt = rs) %>%
  # dplyr::mutate(rs=NULL) %>%
  # dplyr::arrange(gene_id) %>%
  # dplyr::filter(duplicated(gene_id) == F) %>% # exclude ties
  # tibble::rownames_to_column('tmp') %>% 
  # dplyr::mutate(tmp=NULL) %>%
  # tibble::column_to_rownames('gene_id') %>%
  # round() %>%
  # `colnames<-`( gsub(".","-",colnames(.), fixed=T) ) %>%
  # `colnames<-`( gsub("_1$","",colnames(.), fixed=F) )


  


#tmp = glass.gbm.rnaseq.expression.vst[rownames(glass.gbm.rnaseq.expression.vst) == "ENSG00000198727",]
#type = as.factor(gsub("^(.).*$","\\1",names(tmp)) )
#o = order(tmp)
#plot(tmp[o], pch=19,cex=0.95,col=as.numeric(type[o]) + 1)
#rm(tmp, type, o)


# Load metadata ----


glass.gbm.rnaseq.metadata <- data.frame(sid = colnames(glass.gbm.rnaseq.expression),  stringsAsFactors = F) %>%
  dplyr::mutate(sid = gsub(".","-",sid,fixed=T)) %>% # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sid = gsub("_1$","",sid,fixed=F)) %>% # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sample_barcode = gsub("^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+)-.+$","\\1",sid)) %>%
  dplyr::mutate(pid = gsub("^(............).+$","\\1",sid)) %>%
  dplyr::arrange(pid) %>%
  dplyr::mutate(resection = as.factor(gsub("^.............(..).+$","\\1",sid))) %>% # TP is primary tumour? https://github.com/fpbarthel/GLASS
  dplyr::mutate(dataset =  as.factor(gsub("^(....).+$","\\1",sid)) ) %>%
  dplyr::left_join(
    read.table('data/gsam/data/GLASS_GBM_R1-R2/GLASS.GBM.subtypes.from.Synapse.portal.tsv', header=T,stringsAsFactors = F) %>% # https://www.synapse.org/#!Synapse:syn21441635/tables/
      dplyr::select(c('aliquot_barcode','subtype')) %>%
      dplyr::mutate(subtype = as.factor(subtype)) %>%
      dplyr::rename(GBM.transcriptional.subtype.Synapse = subtype)
    , by     = c('sid' = 'aliquot_barcode' )  )  %>%
  dplyr::mutate(sid.label = gsub("^(.)...(-..-....-).(.).*$","\\1\\2\\3",sid) ) %>%
  dplyr::mutate(dataset = gsub("^(....).*$","\\1",sid) ) %>%
  dplyr::left_join(
    read.table('data/glass_clinical_surgeries.txt',sep="\t",header=T,stringsAsFactors = F) %>%
      dplyr::mutate(ROW_ID=NULL, ROW_VERSION=NULL)
      ,
    by=c('sample_barcode'='sample_barcode')) %>%
  dplyr::mutate(condition = ifelse(resection == "TP","Primary","NotPrimary")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Primary","NotPrimary") )) %>%
  dplyr::filter(idh_status == "IDHwt") # exclude the three IDH mutants according to Synapse WGS/WES VCF files


# exclude samples from patients that had a no grade IV tumor
no.grave.iv <- glass.gbm.rnaseq.metadata %>%
  dplyr::filter(grade %in% c('II','III')) %>%
  dplyr::pull(pid)

glass.gbm.rnaseq.metadata <- glass.gbm.rnaseq.metadata %>%
  dplyr::filter(pid %in% no.grave.iv == F )


stopifnot(!is.na(glass.gbm.rnaseq.metadata$GBM.transcriptional.subtype.Synapse))





## according to TCGA [data/tcga/tables/rna-seq/gdc_sample_sheet.2019-01-28.tsv] these should be recurrent GBMs:
# TCGA-06-0125	TCGA-06-0125-02A  v
# TCGA-06-0152	TCGA-06-0152-02A  [no matching primary pair]
# TCGA-06-0171	TCGA-06-0171-02A  [no matching primary pair]
# TCGA-06-0190	TCGA-06-0190-02A  v
# TCGA-06-0210	TCGA-06-0210-02A  v
# TCGA-06-0211	TCGA-06-0211-02A  v
# TCGA-06-0221	TCGA-06-0221-02A  [no matching primary pair]
# TCGA-14-0736	TCGA-14-0736-02A  [no matching primary pair]
# TCGA-14-1034	TCGA-14-1034-02B  v
# TCGA-14-1402	TCGA-14-1402-02A  [no matching primary pair]
# TCGA-19-0957	TCGA-19-0957-02A  [no matching primary pair]
# TCGA-19-1389	TCGA-19-1389-02A  [no matching primary pair]
# TCGA-19-4065	TCGA-19-4065-02A  WEL in TCGA, WEL in GLASS data sheet, NIET in expressie tabel [glass_transcript_counts.txt] ?

# present in GLASS:
# TCGA-06-0125-R1      IDHwt Glioblastoma, IDH-wildtype
# TCGA-06-0190-R1      IDHwt Glioblastoma, IDH-wildtype
# TCGA-06-0210-R1      IDHwt Glioblastoma, IDH-wildtype
# TCGA-06-0211-R1      IDHwt Glioblastoma, IDH-wildtype
# TCGA-14-1034-R1      IDHwt Glioblastoma, IDH-wildtype
# TCGA-FG-5963-R1      IDHwt Glioblastoma, IDH-wildtype    not matching TCGA?




# select only those that are eligible for the study ----

glass.gbm.rnaseq.expression <- glass.gbm.rnaseq.expression %>%
  dplyr::select(glass.gbm.rnaseq.metadata$sid)


stopifnot(colnames(glass.gbm.rnaseq.expression) %in% glass.gbm.rnaseq.metadata$sid)



## add mutations ----

# file from Synapse portal
# tmp <- read.csv('data/variants_passgeno_20190327.csv',stringsAsFactors = F) %>%
#   dplyr::mutate(pid = gsub("^([^\\-]+.[^\\-]+.[^\\-]+).+$","\\1",aliquot_barcode)) 
# 
# #idh1:chr2:209,098,951-209,121,867
# 
# tmp.idh1 <- tmp %>%
#   dplyr::filter(chrom == "2" & start >=  209098951 & end <= 209121867 ) %>%
#   dplyr::filter(ad_alt != 0) %>%
#   dplyr::filter(pid %in% glass.gbm.rnaseq.metadata$pid) %>%
#   dplyr::filter(ssm2_pass_call == 't') %>%
#   dplyr::arrange(start, aliquot_barcode)
#   
#   # %>% dplyr::filter(start == 209113112) %>% pull(pid) %>% unique()
# 
# 
# 
# 
# 
# #idh2:chr15:90,623,937-90,648,932
# 
# tmp.idh2 <- tmp %>%
#   dplyr::filter(chrom == "15" & start >=  90623937 & end <= 90648932 ) %>%
#   dplyr::filter(ad_alt != 0) %>%
#   dplyr::filter(pid %in% glass.gbm.rnaseq.metadata$pid) %>%
#   dplyr::filter(ssm2_pass_call == 't') %>%
#   dplyr::arrange(start, aliquot_barcode)





# VST transform ----


cond <- as.factor(gsub("^(.).*$","\\1",colnames(glass.gbm.rnaseq.expression)))
glass.gbm.rnaseq.expression.vst <- DESeq2::DESeqDataSetFromMatrix(glass.gbm.rnaseq.expression, data.frame(cond), ~cond)
glass.gbm.rnaseq.expression.vst <- SummarizedExperiment::assay(DESeq2::vst(glass.gbm.rnaseq.expression.vst, blind=F))
rm(cond)




# Add TPC to metadata ----

# These estimates are just not good enough:
# 
# glass.gbm.rnaseq.metadata <- glass.gbm.rnaseq.metadata %>%
#   dplyr::mutate(pid.tmp = gsub("^([^-]+-[^-]+-[^-]+-[^-]+)-.+$","\\1",sid) ) %>%
#   dplyr::left_join(
#   read.table("output/tables/cnv/tumor.percentage.estimate_glass.txt") %>%
#     dplyr::select(c('sample','pct')) %>%
#     dplyr::mutate(pid.tmp = gsub("^([^-]+-[^-]+-[^-]+-[^-]+)-.+$","\\1",sample) ) %>%
#     dplyr::rename(tumour.percentage.dna = pct) ,
#   by = c('pid.tmp' = 'pid.tmp')) %>%
#   dplyr::mutate(pid.tmp = NULL)
# 
# stopifnot(!is.na(glass.gbm.rnaseq.metadata$tumour.percentage.dna))




# Add imputed / predict TPCs ----

glass.gbm.rnaseq.metadata <- glass.gbm.rnaseq.metadata %>%
  dplyr::mutate(tumour.percentage.dna.imputed.caret = NULL) %>%
  dplyr::mutate(tumour.percentage.dna.imputed.rf = NULL) %>%
  #dplyr::left_join(read.table("output/tables/GLASS.tumour.percentage.dna.imputed.caret.txt"), by=c('sid'='sid'))%>%
  dplyr::left_join(read.table("output/tables/GLASS.tumour.percentage.dna.imputed.rf.txt"), by=c('sid'='sid'))



# pretty poor correlation, something seems odd w/ the predictions in GLASS?
#plot(glass.gbm.rnaseq.metadata$tumour.percentage.dna , glass.gbm.rnaseq.metadata$tumour.percentage.dna.imputed.caret)

# Combine table into paired data table ----







