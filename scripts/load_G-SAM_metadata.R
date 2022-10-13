#!/usr/bin/env R

# load libs ----


library(tidyverse) # for distinct() function
library(readxl)


# patient metadata ----


# three CNV samples samples not in metadata: "AMA" "AMA" "BAO" "BAO" "FAF" "FAF"
#gsam.patient.metadata <- read.csv('data/administratie/dbGSAM_PUBLIC_VERSION.csv',stringsAsFactors=F)


gsam.patient.metadata <- read.csv('data/gsam/administratie/GSAM_combined_clinical_molecular.csv',stringsAsFactors=F) |> 
  dplyr::arrange(studyID) |> 
  dplyr::mutate(initialMGMT = NULL) |> 
  dplyr::mutate(gender = ifelse(studyID %in% c('AAT', 'AAM', 'AZH', 'HAI', 'FAG'),"Male",gender)) |>  # there's a number of samples of which the gender does not fit with the omics data - omics data determined genders are the corrected ones
  dplyr::mutate(gender = as.factor(gender)) |> 
  dplyr::mutate(survival.events = dplyr::case_when(
      status == "Deceased" ~ 1,
      status == "Censored" ~ 0,
      T ~ as.numeric(NA))) |> 
  dplyr::mutate(survival.months = survivalDays / 365.0 * 12.0) |> 
  dplyr::mutate(`X.1` = NULL)


## beva ----
### before recurrence ----

gsam.patient.metadata <- gsam.patient.metadata |> 
  dplyr::mutate(bevacizumab.before.recurrence = case_when(
    studyID %in% c('GAA','GAG','CCB','CBI','CBV','AAX') ~ "Yes",
    studyID %in% c('HAD','HAF') ~ "Trial participant",
    T ~ "No"
  ))



## PTK787 ---
### before recurrence ----


gsam.patient.metadata <- gsam.patient.metadata |> 
  dplyr::mutate(PTK787.before.recurrence = case_when(
    studyID %in% c('CBR','CBW') ~ "Yes",
    T ~ "No"
  ))



# DNA exome-seq ----


# it remains unclear what files are missing


gsam.wes.samples <- data.frame(file = Sys.glob('data/gsam/DNA/dna_data_2020/request_962/PDexport/decrypted/*/mapped_sample/*.bam')) |> 
  dplyr::mutate(wesid1 = gsub("^.+/decrypted/([^/]+)/.+$","\\1", file))  |>  #   #dplyr::mutate(wesid2 = gsub("^.+_([^\\.]+)\\..+$","\\1", file))
  dplyr::mutate(batch = as.factor(gsub("_.+$","",wesid1)))  |> 
  dplyr::mutate(wesid1 = gsub("^.+_","",wesid1))





gsam.cnv.metadata <- read.delim("data/gsam/DNA/sample codes sanger gsam.txt",stringsAsFactors=FALSE) |> 
  dplyr::mutate(pid = gsub("[1-2]$","",donor_ID))  |> 
  dplyr::select(c("donor_ID", "pid", "PD_ID", "donor_sex", "donor_age_at_diagnosis","Concentration.at.QC..ng.ul.","Volume.at.QC..ul.","Amount.at.QC..ug.")) |> 
  dplyr::full_join(
    
    read.delim('data/gsam/output/tables/cnv_copynumber-ratio.cnr_log2_all.txt',stringsAsFactors=F) |>
      dplyr::select( - c('chromosome', 'start', 'end', 'gene') ) |>
      dplyr::slice_head(n=1) |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column('cnv.table.id') |>
      dplyr::mutate(V1 = NULL) |>
      dplyr::mutate(sid = gsub("\\.b[12]$","",cnv.table.id) ) |>
      dplyr::mutate(batch =  as.factor(gsub("^[^\\.]+\\.","",cnv.table.id))) |>
      dplyr::mutate(CNVKIT.output = T)
    
    , by=c('PD_ID' = 'sid')) |>
  dplyr::left_join(gsam.patient.metadata , by=c('pid' = 'studyID')) |>
  dplyr::mutate(donor_sex = NULL) |>
  dplyr::left_join(
    read.delim('data/gsam/output/tables/dna/idh_mutations.txt', stringsAsFactors = F, header=F) |>
      `colnames<-`(c('PD_ID' , 'IDH.mutation', 'IDH.mutation.call.status', 'IDH.mutation.VAF', 'IDH.mutation.count')),
    by = c('PD_ID'='PD_ID')) |>
  dplyr::select(c('donor_ID', 'pid', 'PD_ID', 'IDH1', 'IDH.mutation', 'IDH.mutation.call.status', 'IDH.mutation.VAF', 'IDH.mutation.count','CNVKIT.output')) |>
  dplyr::mutate(tmp = ifelse(is.na(IDH.mutation), 'NA' , IDH.mutation)) |>
  dplyr::mutate(tmp = dplyr::case_when( # BAW EAF FAD and JAF are all IDH1 and IDH2 negative in Kaspars paper but apparently positive in this file. actual mutation not given.
    tmp == "NA" ~ '0',
    tmp == '-' ~ '1' , 
    TRUE ~ '2'
  ))    |>
  dplyr::group_by(pid) |>
  dplyr::mutate(pat.with.IDH = max(tmp), data = dplyr::cur_data() )  |>
  dplyr::ungroup() |>
  dplyr::mutate(tmp = NULL, data=NULL) |>
  as.data.frame() |>
  dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH == 0, NA , pat.with.IDH)) |>
  dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH == 1, F , pat.with.IDH)) |>
  dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH == 2, T , pat.with.IDH))  |>
  dplyr::mutate(resection = gsub("^...(.).*$","R\\1", donor_ID) ) |>
  dplyr::mutate(resection = ifelse(resection == "RN", NA ,resection) )




# RNA-seq metadata [full] ----

## STAR alignment statistics + patient / sample identifiers ----


gsam.rna.metadata <- read.delim("data/gsam/output/tables/gsam_featureCounts_readcounts_new.txt.summary",stringsAsFactors = F,comment="#",row.names=1) %>%
  `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1",colnames(.),fixed=F)) %>%
  dplyr::filter(rowSums(.) > 0) %>%
  t() %>%
  `colnames<-`(paste0("STAR.",colnames(.))) %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::mutate(STAR.sum = (STAR.Assigned + STAR.Unassigned_Chimera + STAR.Unassigned_Duplicate + STAR.Unassigned_MultiMapping + STAR.Unassigned_NoFeatures + STAR.Unassigned_Ambiguity)) %>%
  dplyr::mutate(pct.duplicate.STAR = STAR.Unassigned_Duplicate / STAR.sum * 100) %>%
  dplyr::mutate(pct.multimap.STAR = STAR.Unassigned_MultiMapping / STAR.sum * 100) %>%
  dplyr::mutate(pct.nofeature.STAR = STAR.Unassigned_NoFeatures / STAR.sum * 100) %>%
  dplyr::mutate(duplicate.fold.STAR = 1 / (1 - (pct.duplicate.STAR / 100)) ) %>%  # 75% duplicate means 4 fold duplication
  dplyr::mutate(sid =  gsub(".","-", sample,fixed=T) ) %>%
  dplyr::mutate(pid = as.factor( gsub("^(...).*$","\\1", sid) )) %>%
  dplyr::mutate(resection = as.factor(gsub("^...(.).*$","r\\1", sid) )) %>%
  dplyr::mutate(new.batch = grepl("new", sid)) %>%
  dplyr::mutate(name.strip = ifelse(new.batch == T , gsub("-new", "", sid, fixed=T),"")) %>%
  dplyr::mutate(old.batch = gsub("-replicate","",sid) %in% name.strip) %>%
  dplyr::mutate(name.strip = NULL) %>%
  dplyr::mutate(batch = dplyr::case_when ((old.batch == T & new.batch == F) ~ "old" , (old.batch == F & new.batch == T) ~ "new" , (old.batch == F & new.batch == F) ~ "single" )) %>%
  dplyr::mutate(batch = as.factor(batch)) %>%
  dplyr::arrange(sample) %>%
  dplyr::left_join( # add chromosomal distribution of rRNA containing alternate loci
    (
      read.delim("data/gsam/output/tables/qc/idxstats/samtools.indexstats.matrix.txt",stringsAsFactors=F,row.names=1) %>%
        `colnames<-`(gsub(".samtools.idxstats.txt","",colnames(.),fixed=T)) %>%
        t() %>%
        data.frame(stringsAsFactors = F ) %>%
        tibble::rownames_to_column("sample") %>%
        dplyr::filter(sample != "ref.len") %>%
        dplyr::mutate(X. = NULL) %>%
        dplyr::mutate(idxstats.sum = rowSums( dplyr::select(. , c(-sample) ) )) %>%
        dplyr::mutate(pct.rRNA.by.chrUn.gl000220 = chrUn_gl000220 / idxstats.sum )
    ), by=c('sample' = 'sample')) %>% 
  dplyr::mutate(resection_pair=ifelse(grepl(1,sample) == T,"matching_r1","matching_r2")) %>%
  dplyr::mutate(resection_pair=replace(resection_pair,which(batch=="old"|sample%in%c("CAO1.replicate","GAS2.replicate","FAB2","FAH2","EBP1","KAE1.new","KAE1")),NA)) %>%
  dplyr::mutate(sample = NULL)




#EBP1, FAH2 and KAE1: no pair
#FAB2: FAB2.replicate contains more vIII reads 
#CAO1.replicate, GAS2.replicate: CAO1 and GAS2 contain more vIII reads

#plot(gsam.rna.metadata$pct.rRNA.by.chrUn.gl000220 , gsam.rna.metadata$pct.multimap.STAR)
#plot(gsam.rna.metadata$pct.rRNA.by.chrUn.gl000220 , gsam.rna.metadata$pct.duplicate.STAR)
#plot(gsam.rna.metadata$pct.rRNA.by.chrUn.gl000220 , gsam.rna.metadata$pct.nofeature.STAR)
#plot(gsam.rna.metadata$pct.nofeature.STAR , gsam.rna.metadata$pct.duplicate.STAR)

# tin.py qc metrics
# tmp <- read.table('output/tables/qc/tin.py/tin.py.matrix.txt',stringsAsFactors=F,header=T)
# tmp$sid <- gsub(".bam","",tmp$Bam_file,fixed=T)
# rownames(tmp) <- tmp$sid
# tmp$Bam_file <- NULL
# gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x = "sid", by.y = "sid")


#vIII rna-seq counts
tmp <- read.table('data/gsam/output/tables/v3_extract_readcounts.txt',header=T,stringsAsFactor=F)
tmp$sample <- gsub("^.+/alignments/([^/]+)/.+$","\\1",tmp$sample)
rownames(tmp) <- tmp$sample

sel <- tmp$wt.reads.v3 + tmp$vIII.reads.v3 > 15
tmp$vIII.percentage <- NA
tmp$vIII.percentage[sel] <- tmp$vIII.reads.v3[sel] / (tmp$wt.reads.v3[sel] + tmp$vIII.reads.v3[sel]) * 100

#gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x = "sid", by.y = "sample")
gsam.rna.metadata <- gsam.rna.metadata %>%
  dplyr::left_join(tmp, by=c('sid' = 'sample'))


rm(sel, tmp)






# @TODO vIII qPCR percentage 'TODO!!!!!!
# tmp <- read.csv('data/RNA/Final_qPCR_EGFR_GSAM.csv',stringsAsFactors = F)
# tmp <- tmp[,colnames(tmp) %in% c('EGFRCt002', 'vIIICt002', 'recurrent_EGFRCt002', 'recurrent_vIIICt002','recurrent_patientID')]
# 
# tmp.1 <- tmp[,match( c( 'EGFRCt002', 'vIIICt002',  'recurrent_patientID' ) , colnames(tmp) ) ]
# tmp.1$resection <- 'r1'
# colnames(tmp.1)[colnames(tmp.1) == "EGFRCt002"] <- 'qPCR.ct.EGFR.wt'
# colnames(tmp.1)[colnames(tmp.1) == "vIIICt002"] <- 'qPCR.ct.EGFR.vIII'
# 
# tmp.2 <- tmp[,match( c( 'recurrent_EGFRCt002', 'recurrent_vIIICt002','recurrent_patientID' ) , colnames(tmp) ) ]
# tmp.2$resection <- 'r2'
# colnames(tmp.2)[colnames(tmp.2) == "recurrent_EGFRCt002"] <- 'qPCR.ct.EGFR.wt'
# colnames(tmp.2)[colnames(tmp.2) == "recurrent_vIIICt002"] <- 'qPCR.ct.EGFR.vIII'
# 
# tmp <- rbind(tmp.1, tmp.2)
# tmp$resection <- as.factor(tmp$resection)
# rm(tmp.1, tmp.2)
# 
# tmp$qPCR.percent.EGFR.vIII <- 100 - (1/(1 + 2 ^ (tmp$qPCR.ct.EGFR.wt - tmp$qPCR.ct.EGFR.vIII)) * 100)
# tmp[tmp$qPCR.ct.EGFR.vIII >= 40,]$qPCR.percent.EGFR.vIII <- 0
# 
# gsam.rna.metadata$tmp.id <- paste0(gsam.rna.metadata$pid, ".", gsam.rna.metadata$resection)
# tmp$tmp.id <- paste0(tmp$recurrent_patientID, '.', tmp$resection)
# tmp$resection <- NULL
# tmp$recurrent_patientID <- NULL
# 
# gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x="tmp.id" , by.y = "tmp.id" , all.x = T) 
# gsam.rna.metadata$tmp.id <- NULL
# x-check replictes CAO1, FAB2, GAS2 => works

# blacklist by heavy DNA contamination
blacklist.too.low.assigned <- c("AKA1", "CAC1", "AAB2", "GAS1", "KAE1", "CCZ1", "GAO2", "JAN1", "BAU2", "EAV2", "AAP1", "AZE1", "HAF1", "GAM1", "HAG1", "BAX2", "EAN1", "CBQ1", "AAD2", "HAK1", "CBG2", "BAI2", "HAE1", "CDH1", "HAI1", "KAB2", "GAE1", "BAN1", "KAC1", "KAA1", "ABA1")
blacklist.heavy.dna.contamination  <- c("CAV1", "BAT2", "EBW1", "HAE1", "BAU1", "EBO1", "GAE1", "CDH1", "KAC2", "ABA1", "KAA1", "KAC1", "BAN1", "KAB2", "GAJ2", "HAI1", "AAD2", "CBG2","AAP1")
blacklist.gender.mislabeling <- c("AAM1","AAM2","AAT1","AAT2","AZH1","AZH2","FAG1","HAI1","HAI2") # metadata of these samples can't be trusted


# some samples with rather large GC bias do not completely fall of the PCA ~ this metric may indicate bad quality but is not a very severe factory in downstream analysis
blacklist.gc.bias <- unique(c("CBI1","CBV2","CBS2","EAP2","AAS1","CBQ1","GAQ2","ECG1","EAP2","AAS1","HAF1","EAG2","HAC1","CCZ2","BAX2","GAA1","CBT1","EAV2","CBT1","HAG1","AAP1","AZG1","AAP1","KAA2","AZH2","CCZ1","AZG1","CAO1","CAO1-replicate","AZH2","ECF2","JAM1","HAF2","GAO2","CBS1","EAN1","KAE1","GAM1","CBS1","AZE1","JAN1","GAS1","EBM1","HAA2","EBN1","AZE1","CAC1","BAY2","JAB1","JAB1","BAU2","BAU2","AKA1","AAB2"))

# outlier samples found in PCA [./scripts/rna/unsupervised_expression_analysis.R], that all correlate with different qc metrics
blacklist.pca <- c("AAB2", "AAD2", "AAP1", "ABA1", "AKA1", "AZE1", "AZH2", "BAN1", "BAT2", "BAU1", "BAU2", "BAX2", "BAY2", "CAC1", "CAO1", "CAO1.replicate", "CAV1", "CBG2", "CBQ1", "CBS1", "CCZ1", "CDH1", "EAN1", "EAV2", "EBL2", "EBM1", "EBN1", "EBO1", "EBW1", "ECF1", "GAE1", "GAJ2", "GAM1", "GAO2", "GAS1", "HAA2", "HAC1", "HAE1", "HAF1", "HAF2", "HAG1", "HAI1", "HAK1", "JAB1", "JAN1", "KAA1", "KAB2", "KAC1", "KAC2", "KAE1")


gsam.rna.metadata$blacklist.too.low.assigned <- gsam.rna.metadata$sid %in% blacklist.too.low.assigned
gsam.rna.metadata$blacklist.heavy.dna.contamination <- gsam.rna.metadata$sid %in% blacklist.heavy.dna.contamination
gsam.rna.metadata$blacklist.gender.mislabeling  <- gsam.rna.metadata$sid %in% blacklist.gender.mislabeling 
gsam.rna.metadata$blacklist.gc.bias  <- gsam.rna.metadata$sid %in% blacklist.gc.bias
gsam.rna.metadata$blacklist.pca <- gsam.rna.metadata$sid %in% blacklist.pca


rm(blacklist.too.low.assigned)
rm(blacklist.heavy.dna.contamination)
rm(blacklist.pca)
rm(blacklist.gc.bias)
rm(blacklist.gender.mislabeling)

# @TODO add batches 'TO DO
# tmp <- read.table('data/administratie/plate.layout.table.txt',stringsAsFactors=F,header=T)
# gsam.rna.metadata <- merge(gsam.rna.metadata, tmp , by.x="sid" , by.y = "sid")
# rm(tmp)
# gsam.rna.metadata$plate <- as.factor(gsam.rna.metadata$plate)
# gsam.rna.metadata$storage.box <- as.factor(gsam.rna.metadata$storage.box )


# add GC offset
tmp <- read.delim("data/gsam/output/tables/qc/gc_content_rmse.txt",stringsAsFactors = F)
tmp <- tmp[order(tmp$RMSE, tmp$sample.id),]
#tmp$filename <- factor(tmp$filename , levels=tmp$filename)
# take average if multiple FQ files exist ~ manual inspection indicated barely differences across multiple FQs
tmp$filename <- NULL
tmp <- aggregate(tmp[,-1], list(tmp$sample.id), mean)
tmp$Group.1 <- gsub("repl$","replicate",tmp$Group.1)

# old line, probably failed:
#gsam.rna.metadata <- merge(gsam.rna.metadata, tmp, by.x="sid", by.y = "Group.1")
gsam.rna.metadata <- gsam.rna.metadata %>%
  dplyr::left_join(tmp, by=c('sid'='Group.1'))

rm(tmp)


## GIGA sequencing facility run statistics----


# N sheets: 6
tmp <- 'data/gsam/documents/PFrench_Summary-sheet_input-DV-qPCR.xlsx'
tmp <- read_excel(tmp,sheet=1)
tmp <- tmp[!is.na(tmp$seqID),]
tmp$seqID <- gsub("_","-",tmp$seqID,fixed=T)
tmp$`#` <- NULL
tmp$...10 <- NULL
tmp$...12 <- NULL
colnames(tmp) <- paste0('giga.', colnames(tmp))
tmp$giga.Plate <- paste0('plate',tmp$giga.Plate)
#print(tmp)

#print(dim(gsam.rna.metadata))
gsam.rna.metadata <- merge(gsam.rna.metadata, tmp, by.x = 'sid', by.y= 'giga.seqID',all=TRUE)
#print(dim(gsam.rna.metadata))

rm(tmp)


## DV200 ----


tmp <- 'data/gsam/documents/PFrench_Summary-sheet_input-DV-qPCR.xlsx'
tmp.1 <- read_excel(tmp,sheet=3)
tmp.2 <- read_excel(tmp,sheet=4)
tmp.3 <- read_excel(tmp,sheet=5)
tmp.4 <- read_excel(tmp,sheet=6)

tmp.1$`#` <- NULL
tmp.2$`#` <- NULL
tmp.3$`#` <- NULL
tmp.4$`#` <- NULL

tmp.1$`RNA quantity \r\nà atteindre\r\n(selon DV200)` <- NA
tmp.2$`ng/µL Ribo` <- NA
tmp.3$`ng/µL Ribo` <- NA
tmp.4$`ng/µL Ribo` <- NA

tmp.1 <- tmp.1[,c(1, 2, 3, 4, 5, 6, 9, 7 , 8)]
tmp.2 <- tmp.2[,c(1, 2, 3, 9, 4, 5 , 6 , 7 , 8 )]
tmp.3 <- tmp.3[,c(1, 2, 3, 9, 4, 5 , 6 , 7 , 8 )]
tmp.4 <- tmp.4[,c(1, 2, 3, 9, 4, 5 , 6 , 7 , 8 )]


stopifnot(length(unique(c(colnames(tmp.1)[1],colnames(tmp.2)[1],colnames(tmp.3)[1],colnames(tmp.4)[1]))) == 1)# [1] = seqID
stopifnot(length(unique(c(colnames(tmp.1)[2],colnames(tmp.2)[2],colnames(tmp.3)[2],colnames(tmp.4)[2]))) == 1)# [2] = isolationID
stopifnot(length(unique(c(colnames(tmp.1)[3],colnames(tmp.2)[3],colnames(tmp.3)[3],colnames(tmp.4)[3]))) == 1)
stopifnot(length(unique(c(colnames(tmp.1)[4],colnames(tmp.2)[4],colnames(tmp.3)[4],colnames(tmp.4)[4]))) == 1)
stopifnot(length(unique(c(colnames(tmp.1)[5],colnames(tmp.2)[5],colnames(tmp.3)[5],colnames(tmp.4)[5]))) == 1)
stopifnot(length(unique(c(colnames(tmp.1)[6],colnames(tmp.2)[6],colnames(tmp.3)[6],colnames(tmp.4)[6]))) == 1)
stopifnot(length(unique(c(colnames(tmp.1)[7],colnames(tmp.2)[7],colnames(tmp.3)[7],colnames(tmp.4)[7]))) == 1)
stopifnot(length(unique(c(colnames(tmp.1)[8],colnames(tmp.2)[8],colnames(tmp.3)[8],colnames(tmp.4)[8]))) == 1)
stopifnot(length(unique(c(colnames(tmp.1)[9],colnames(tmp.2)[9],colnames(tmp.3)[9],colnames(tmp.4)[9]))) == 1)

tmp <- rbind(tmp.1, tmp.2, tmp.3, tmp.4)
rm(tmp.1, tmp.2, tmp.3, tmp.4)

tmp <- tmp[!is.na(tmp$seqID),]
tmp$seqID <- gsub("_","-",tmp$seqID,fixed=T)
colnames(tmp) <- paste0('giga.', colnames(tmp))
colnames(tmp) <- gsub("[\r\n ]+",".",colnames(tmp)) # avoid white spaces

print(dim(gsam.rna.metadata))
gsam.rna.metadata <- merge(gsam.rna.metadata, tmp, by.x = 'sid', by.y= 'giga.seqID',all=TRUE)
print(dim(gsam.rna.metadata))

rm(tmp)


## Tumor Percentages (DNA) ----

tmp <- read.delim("data/gsam/output/tables/cnv/tumor.percentage.estimate.txt", sep = " ") |>
  dplyr::mutate(lfc.3p = NULL) |>
  dplyr::mutate(lfc.4p = NULL) |>
  dplyr::mutate(lfc.n = NULL) |>
  dplyr::mutate(dist = NULL) |>
  dplyr::rename(tumour.percentage.dna = pct) |>
  dplyr::filter(!duplicated(sample))

gsam.rna.metadata <- gsam.rna.metadata |>
  dplyr::mutate(tmp = gsub("^(....).*$", "\\1", sid)) |> # dummy
  dplyr::left_join(tmp, by = c("tmp" = "sample")) |>
  dplyr::mutate(tmp = NULL)

rm(tmp)



## Resection / biopsy ----


tmp <- gsam.rna.metadata |>
  dplyr::select(sid, pid, resection) |>
  dplyr::left_join(
    gsam.patient.metadata |>
      dplyr::select(studyID, extentOfResectionFirstSurgery, extentOfResectionSecondSurgery),
    by = c("pid" = "studyID")
  ) |>
  dplyr::mutate(extent = ifelse(resection == "r1", extentOfResectionFirstSurgery, extentOfResectionSecondSurgery)) |>
  dplyr::mutate(extent = case_when(
    is.na(extent) ~ as.character(NA),
    extent == "Biopsy" ~ "Biopsy",
    T ~ "Resection"
  )) |>
  dplyr::mutate(pid = NULL) |>
  dplyr::mutate(resection = NULL) |>
  dplyr::mutate(extentOfResectionFirstSurgery = NULL) |>
  dplyr::mutate(extentOfResectionSecondSurgery = NULL)

stopifnot(nrow(tmp) == 399)
stopifnot(sum(is.na(tmp$extent)) == 6)


gsam.rna.metadata <- gsam.rna.metadata %>%
  dplyr::left_join(tmp, by = c("sid" = "sid"), suffix = c("", ""))
rm(tmp)



## MGMT ----


tmp <- gsam.rna.metadata |>
  dplyr::select(.data$sid, .data$pid, .data$resection) |>
  dplyr::left_join(
    gsam.patient.metadata |>
      dplyr::select(.data$studyID, .data$mgmtPrimary, .data$mgmtRecurrent),
    by = c("pid" = "studyID")
  ) |>
  dplyr::mutate(MGMT = ifelse(.data$resection == "r1", .data$mgmtPrimary, .data$mgmtRecurrent)) |>
  dplyr::mutate(pid = NULL) |>
  dplyr::mutate(resection = NULL) |>
  dplyr::mutate(mgmtPrimary = NULL) |>
  dplyr::mutate(mgmtRecurrent = NULL) |>
  dplyr::filter(!is.na(.data$MGMT))

stopifnot(nrow(tmp) == 248)

gsam.rna.metadata <- gsam.rna.metadata %>%
  dplyr::left_join(tmp, by = c("sid" = "sid"), suffix = c("", ""))
rm(tmp)



## gliovis subtypes ----
# https://gliovis.shinyapps.io/GlioVis/


tmp <- read.csv("output/tables/gliovis/GlioVis_-_Visualization_Tools_for_Glioma_Datasets_ThreeWay.csv") |>
  dplyr::rename(gliovis.svm_call = svm_call) |>
  dplyr::rename(gliovis.knn_call = knn_call) |>
  dplyr::rename(gliovis.gsea_call = gsea_call) |>
  dplyr::rename(gliovis.equal_call = equal.call) |>
  dplyr::rename(gliovis.majority_call = majority.call) |>
  dplyr::left_join(read.csv("output/tables/gliovis/GlioVis_-_Visualization_Tools_for_Glioma_Datasets_SVM.csv") |>
    dplyr::rename(svm.Classical = Classical) |>
    dplyr::rename(svm.Mesenchymal = Mesenchymal) |>
    dplyr::rename(svm.Proneural = Proneural),
  by = c("Sample" = "Sample")
  ) |>
  dplyr::left_join(read.csv("output/tables/gliovis/GlioVis_-_Visualization_Tools_for_Glioma_Datasets_KNN.csv") |>
    dplyr::rename(knn.Classical = Classical) |>
    dplyr::rename(knn.Mesenchymal = Mesenchymal) |>
    dplyr::rename(knn.Proneural = Proneural),
  by = c("Sample" = "Sample")
  ) |>
  dplyr::left_join(read.csv("output/tables/gliovis/GlioVis_-_Visualization_Tools_for_Glioma_Datasets_ssGSEA.csv") |>
    dplyr::mutate(Sample = gsub(".", "-", X, fixed = T), X = NULL) |>
    dplyr::rename(ssGSEA.Classical.score = Classical) |>
    dplyr::rename(ssGSEA.Mesenchymal.score = Mesenchymal) |>
    dplyr::rename(ssGSEA.Proneural.score = Proneural) |>
    dplyr::rename(ssGSEA.Classical = Classical_pval) |>
    dplyr::rename(ssGSEA.Mesenchymal = Mesenchymal_pval) |>
    dplyr::rename(ssGSEA.Proneural = Proneural_pval),
  by = c("Sample" = "Sample")
  ) |>
  dplyr::rename(gliovis.ssGSEA.Proneural.score = ssGSEA.Proneural.score) |>
  dplyr::rename(gliovis.ssGSEA.Classical.score = ssGSEA.Classical.score) |>
  dplyr::rename(gliovis.ssGSEA.Mesenchymal.score = ssGSEA.Mesenchymal.score) |>
  dplyr::rename(gliovis.ssGSEA.Proneural = ssGSEA.Proneural) |>
  dplyr::rename(gliovis.ssGSEA.Classical = ssGSEA.Classical) |>
  dplyr::rename(gliovis.ssGSEA.Mesenchymal = ssGSEA.Mesenchymal)


stopifnot(tmp$gliovis.svm_call == tmp$svm.subtype.call)
stopifnot(tmp$gliovis.knn_call == tmp$knn.subtype.call)
stopifnot(tmp$gliovis.gsea_call == tmp$gsea.subtype.call)


tmp <- tmp |>
  dplyr::mutate(svm.subtype.call = NULL) |>
  dplyr::mutate(knn.subtype.call = NULL) |>
  dplyr::mutate(gsea.subtype.call = NULL)



gsam.rna.metadata <- gsam.rna.metadata |>
  dplyr::left_join(tmp, by = c("sid" = "Sample"), suffix = c("", ""))


rm(tmp)



## ssGSEA [ssgsea.GBM.classification] 2022 ----

# all data have been ran in a single batch
# see: analysis_ssGSEA.R

tmp <- read.table("output/tables/ssgsea.GBM.classification_p_result_tmp.gct.txt") |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::filter(grepl("^(GLSS|TCGA)",sid) == F) |> 
  dplyr::mutate(sid = gsub(".","-",sid, fixed=T)) |> 
  dplyr::rename(ssGSEA.2022.Proneural.enrichment_score = Proneural) |> 
  dplyr::rename(ssGSEA.2022.Classical.enrichment_score  = Classical ) |> 
  dplyr::rename(ssGSEA.2022.Mesenchymal.enrichment_score = Mesenchymal) |> 
  dplyr::rename(ssGSEA.2022.Proneural_pval = Proneural_pval) |> 
  dplyr::rename(ssGSEA.2022.Classical_pval = Classical_pval) |> 
  dplyr::rename(ssGSEA.2022.Mesenchymal_pval = Mesenchymal_pval)
stopifnot(tmp$sid %in% gsam.rna.metadata$sid) # ensure match



tmp.call <- tmp |> 
  tidyr::pivot_longer(c(ssGSEA.2022.Proneural.enrichment_score, ssGSEA.2022.Classical.enrichment_score, ssGSEA.2022.Mesenchymal.enrichment_score),
                      values_to = "enrichment_score", names_to = 'subtype.es') |> 
  dplyr::mutate(subtype.es = gsub('ssGSEA.2022.','',subtype.es)) |>  
  dplyr::mutate(subtype.es = gsub('.enrichment_score','',subtype.es)) |> 
  tidyr::pivot_longer(c(ssGSEA.2022.Proneural_pval, ssGSEA.2022.Classical_pval, ssGSEA.2022.Mesenchymal_pval),
                      values_to = "pval", names_to = 'subtype.pval') |> 
  dplyr::mutate(subtype.pval = gsub('ssGSEA.2022.','',subtype.pval)) |> 
  dplyr::mutate(subtype.pval = gsub('_pval','',subtype.pval)) |> 
  dplyr::filter(subtype.es == subtype.pval) |> 
  dplyr::rename(ssGSEA.2022.subtype = subtype.es, sutype.pval = NULL ) |> 
  dplyr::group_by(sid) |> 
  dplyr::filter(pval == min(pval)) |> 
  dplyr::ungroup() |> 
  dplyr::select(c('sid','ssGSEA.2022.subtype')) |> 
  dplyr::group_by(sid) |> 
  dplyr::summarise(ssGSEA.2022.subtype = as.factor(paste(ssGSEA.2022.subtype, collapse="|"))) |> 
  dplyr::ungroup()
stopifnot(tmp.call$sid %in% tmp$sid) # ensure match
stopifnot(tmp$sid %in% tmp.call$sid) # ensure match



gsam.rna.metadata <- gsam.rna.metadata |> 
  dplyr::left_join(tmp.call , by=c('sid'='sid'), suffix=c('','')) 

gsam.rna.metadata <- gsam.rna.metadata |> 
  dplyr::left_join(tmp , by=c('sid'='sid'), suffix=c('','')) 


rm(tmp.call, tmp)




## NMF stats 2021 ----


tmp <- read.table("data/gsam/output/tables/gsam_nmf_lda_data.txt") |>
  dplyr::mutate(sid = gsub("^GSAM-", "", .data$sid)) |>
  dplyr::rename(NMF.123456.PCA.SVM.class_2021 = .data$NMF.123456.PCA.SVM.class) |>
  dplyr::select(.data$sid, .data$`NMF.123456.PCA.SVM.class_2021`)

gsam.rna.metadata <- gsam.rna.metadata |>
  dplyr::left_join(tmp, by = c("sid" = "sid"), suffix = c("", ""))

rm(tmp)


## Gravendeel class ----


#gsam.rna.metadata <- gsam.rna.metadata %>%
#  dplyr::left_join(read.table("output/tables/gravendeel_centroid_classification_gsam.txt"), by=c('sid' = 'sid'))



## Add IDH status to RNA samples ----


gsam.rna.metadata <- gsam.rna.metadata %>%
  dplyr::mutate(tmp = gsub('-new','', sid, fixed=T) ) %>%
  dplyr::mutate(tmp = gsub('-replicate','', tmp, fixed=T) ) %>%
  dplyr::left_join(
    
    gsam.cnv.metadata %>%
      dplyr::select(c('donor_ID' , 'PD_ID', 'pat.with.IDH')) %>%
      dplyr::rename(sid = donor_ID) %>%
      dplyr::filter(!is.na(pat.with.IDH)) %>%
      dplyr::mutate(pat.with.IDH = ifelse(is.na(pat.with.IDH), 0, pat.with.IDH)) %>% # collapse those with replicates
      dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH == F, 1, pat.with.IDH)) %>%
      dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH == T, 2, pat.with.IDH)) %>%
      dplyr::group_by(sid) %>%
      dplyr::summarise(pat.with.IDH = max(pat.with.IDH), .groups = 'drop') %>%
      as.data.frame() %>%
      dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH  == 0, NA, pat.with.IDH)) %>%
      dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH  == 1, F, pat.with.IDH)) %>%
      dplyr::mutate(pat.with.IDH = ifelse(pat.with.IDH  == 2, T, pat.with.IDH))
    
    , by = c('tmp' = 'sid')) %>%
  dplyr::mutate(tmp = NULL) %>%
  dplyr::mutate( pat.with.IDH = as.logical(pat.with.IDH) )


## Add PCA signatures 2021 ----
# needed to show changes in cluster numbers 


tmp <- read.table("data/gsam/output/tables/principal_DE_cluster_components.txt") |>
  dplyr::mutate(sid = gsub("^GSAM.", "", pid), pid = NULL) |>
  dplyr::rename_with(~ paste0(.x, ".2021"), .cols = matches("compon|signatu", perl = T))

gsam.rna.metadata <- gsam.rna.metadata |>
  dplyr::left_join(tmp, by = c("sid" = "sid"), suffix = c("", ""))

rm(tmp)



## Add PCA signatures 2022 ----

# neuron_oligodendrocyte_endothelial_EM_PCA_scores


tmp <- read.table('output/tables/principal_DE_cluster_components_2022.txt')
gsam.rna.metadata <- gsam.rna.metadata |> 
  dplyr::left_join(tmp, by=c('sid'='sid'), suffix=c('','')) 

rm(tmp)



## Add GITS components 2022 ----


tmp <- readRDS("cache/analysis_GITS_space.Rds")
gsam.rna.metadata <- gsam.rna.metadata |> 
  dplyr::left_join(tmp, by=c('sid'='sid'), suffix=c('','')) 

rm(tmp)


## Add GITS travel segments ----


tmp <- readRDS('cache/analysis_GITS_space.transition.segments.Rds') |> 
  dplyr::group_by(segment) |> 
  dplyr::mutate(status = ifelse(pct_transition == min(pct_transition), "from","to")) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(segment.str = paste0(pred," [",round(pct_transition * 100,2),"%, PC1-scaled=",round(x2_transition,4),", PC2-scaled=",round(y2_transition,4),"]")) |> 
  tidyr::pivot_wider(id_cols = segment, names_from = status,
                     values_from = c(
                       pid, segment.str
                     )) |> 
  dplyr::mutate(segment = NULL, pid_to = NULL) |> 
  dplyr::rename(pid = pid_from) |> 
  dplyr::mutate(`GITS travel segments` = paste0(segment.str_from , ' -> ', segment.str_to)) |> 
  dplyr::select(pid, `GITS travel segments`) |> 
  dplyr::group_by(pid) |> 
  dplyr::summarise(`GITS travel segments` = paste0(`GITS travel segments`, collapse = " -> "))

gsam.rna.metadata <- gsam.rna.metadata |> 
  dplyr::left_join(tmp, by=c('pid'='pid'), suffix=c('','')) 

rm(tmp)



## Add EPIC scores ----


tmp <- readRDS(file ='tmp/out.epic.Rds') |>
  purrr::pluck('mRNAProportions') |> 
  dplyr::rename(`EPIC: B-cells` = `EPIC.Bcells`) |> 
  dplyr::rename(`EPIC: CAFs` = `EPIC.CAFs`) |> 
  dplyr::rename(`EPIC: CD4 T-cells` = `EPIC.CD4_Tcells`) |> 
  dplyr::rename(`EPIC: CD8 T-cells` = `EPIC.CD8_Tcells`) |> 
  dplyr::rename(`EPIC: Endothelial` = `EPIC.Endothelial`) |> 
  dplyr::rename(`EPIC: Macrophages` = `EPIC.Macrophages`) |> 
  dplyr::rename(`EPIC: NK-cells` = `EPIC.NKcells`) |> 
  dplyr::rename(`EPIC: other cells` = `EPIC.otherCells`)

gsam.rna.metadata <- gsam.rna.metadata |> 
  dplyr::left_join(tmp, by=c('sid'='sid'), suffix=c('','')) 

rm(tmp)



# 〰 © Dr. Youri Hoogstrate 〰 ----


