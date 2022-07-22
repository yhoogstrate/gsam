#!/usr/bin/env R

# load libs ----


# read counts from GLASS samples ----

# possibly some match w/ Wang dataset?


# The samples are taken from a data portal:
# www.synapse.org

# expression values are deterimined by Kallisto, per transcript


# hacking this back from GLASS git code refers to gencode v19...
# gencode.v19.chr_patch_hapl_scaff.annotation.gtf << almost no mismatches [1193]


glass.gencode.v19 <- 'data/gsam/data/GLASS_GBM_R1-R2/gencode.v19.chr_patch_hapl_scaff.annotation.translate-table.txt' |> 
  read.table(header=F, stringsAsFactors = F) |> 
  dplyr::rename(gene_symbol=V1) |> 
  dplyr::rename(gene_id=V2) |> 
  dplyr::rename(transcript_id=V3)


stopifnot(!exists('glass.gbm.rnaseq.expression')) # old


# glass.gbm.rnaseq.expression <- 'data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt' %>%
#   read.delim(stringsAsFactors = F) %>%
#   dplyr::mutate(length = NULL) %>% # not needed
#   dplyr::left_join(glass.gencode.v19, by=c('target_id' = 'transcript_id')) %>%
#   dplyr::filter(!is.na(gene_symbol) ) %>% # 1193 transcript id's not matching gtf Ensembl 64
#   dplyr::mutate(target_id = NULL) %>% # aggregate @ gene_id level
#   dplyr::mutate(gene_symbol = NULL) %>%
#   dplyr::group_by(gene_id) %>%
#   summarise(across(everything(), list(sum))) %>%
#   tibble::rownames_to_column('tmp') %>%
#   dplyr::mutate(tmp=NULL) %>%
#   #dplyr::filter(gene_id %in% c("ENSG00000198804", "ENSG00000198886", "ENSG00000198938","ENSG00000198712", # exclude from normalization?
#   #                             "ENSG00000198727") == F ) %>% # extreme high counts from odd genes
#   tibble::column_to_rownames('gene_id') %>%
#   round() %>%
#   `colnames<-`( gsub(".","-",colnames(.), fixed=T) ) %>%
#   `colnames<-`( gsub("_1$","",colnames(.), fixed=F) ) %>%
#   dplyr::select( # extremely strong outlier samples !! Turn out to be at least 3 IDH mutants
#     -c("GLSS-SM-R068-TP-01R-RNA-0UPMYO", "GLSS-SM-R068-R1-01R-RNA-7I5H9P", # IDH Mut
#        "GLSS-SM-R071-TP-01R-RNA-OAXGI8", "GLSS-SM-R071-R1-01R-RNA-7AZ6G2", # IDH Mut
#        "GLSS-SM-R099-R1-01R-RNA-MNTPMI", #"GLSS-SM-R099-TP-01R-RNA-YGXA72",
#        "GLSS-SM-R100-TP-01R-RNA-EUI7AZ", "GLSS-SM-R100-R1-01R-RNA-46UW5U"  # IDH Mut
#     ))



glass.gbm.rnaseq.expression.all.samples <- 'data/gsam/data/GLASS_GBM_R1-R2/transcript_count_matrix_all_samples.tsv' |> 
  read.delim(stringsAsFactors = F) |> 
  dplyr::mutate(length = NULL) |> # not needed
  dplyr::left_join(glass.gencode.v19, by=c('target_id' = 'transcript_id')) |>
  dplyr::filter(!is.na(gene_symbol) ) |> # 1193 transcript id's not matching gtf Ensembl 64
  dplyr::mutate(target_id = NULL) |> # aggregate @ gene_id level
  dplyr::mutate(gene_symbol = NULL) |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(across(everything(), list(sum))) |>
  tibble::rownames_to_column('tmp') |>
  dplyr::mutate(tmp=NULL) |>
  tibble::column_to_rownames('gene_id') |>
  base::round() |>
  dplyr::rename_with( ~ gsub(".", "-", .x, fixed = TRUE)) |> 
  dplyr::rename_with( ~ gsub("_1$", "", .x, fixed = FALSE))
# dplyr::select( # extremely strong outlier samples !!
#   -c("GLSS-SM-R068-TP-01R-RNA-0UPMYO", "GLSS-SM-R068-R1-01R-RNA-7I5H9P", # IDH Mut
#      "GLSS-SM-R071-TP-01R-RNA-OAXGI8", "GLSS-SM-R071-R1-01R-RNA-7AZ6G2", # IDH Mut
#      "GLSS-SM-R099-R1-01R-RNA-MNTPMI", #"GLSS-SM-R099-TP-01R-RNA-YGXA72",
#      "GLSS-SM-R100-TP-01R-RNA-EUI7AZ", "GLSS-SM-R100-R1-01R-RNA-46UW5U"  # IDH Mut
#   ))



#'@todo check with old table
#' [v] identical / missing sample names
# stopifnot(colnames(glass.gbm.rnaseq.expression) %in% colnames(glass.gbm.rnaseq.expression.all.samples))
# stopifnot(rownames(glass.gbm.rnaseq.expression) == rownames(glass.gbm.rnaseq.expression.all.samples))
#' [x] identical expression values - 21 samples have changed read counts, in approximately 150-450 genes each
# a = glass.gbm.rnaseq.expression
# b = glass.gbm.rnaseq.expression.all.samples |> dplyr::select(colnames(glass.gbm.rnaseq.expression))
# 
# for(i in 1:ncol(a)) {
#   print( paste0(colnames(a)[i]," -> ",sum(a[,i] != b[,i])))
# }

#'@todo check what is with the excluded samples in the data overview plot
#' [?] batch effect?
#' [?] suspected IDH outliers?
#' [x] 30 changed subtype ? this is worrysome






# Load metadata ----


stopifnot(!exists('glass.gbm.rnaseq.metadata'))


## batch 2021 id's ----
# those that were in the 2021 Synapse release / first revision
# labels needed to find possibly batch effects
tmp.batch.2021 <- 'data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt' |> 
  read.delim(stringsAsFactors = F,nrows=1) |> 
  dplyr::mutate(target_id=NULL, length=NULL) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::mutate(V1 = NULL) |> 
  dplyr::mutate(sid = gsub(".","-",sid,fixed=T))  |>  # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sid = gsub("_1$","",sid,fixed=F)) |>  # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(batch = paste0(gsub("^([^\\-]+).+$","\\1",sid,fixed=F),".2021"))


tmp.batch.2022 <- 'data/gsam/data/GLASS_GBM_R1-R2/transcript_count_matrix_all_samples.tsv' |> 
  read.delim(stringsAsFactors = F,nrows=1) |> 
  dplyr::mutate(target_id=NULL, length=NULL) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::mutate(V1 = NULL) |> 
  dplyr::mutate(sid = gsub(".","-",sid,fixed=T))  |>  # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sid = gsub("_1$","",sid,fixed=F)) |>  # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(batch = paste0(gsub("^([^\\-]+).+$","\\1",sid,fixed=F),".2022"))



glass.gbm.rnaseq.metadata.all.samples <- tmp.batch.2022 |> 
  dplyr::full_join(tmp.batch.2021, by=c('sid'='sid'), suffix=c('.2022','.2021')) |> 
  dplyr::mutate(sample_barcode = gsub('^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+).+$','\\1',sid))


# ensure that all the 2021 samples are in the 2022 samples
stopifnot(nrow(glass.gbm.rnaseq.metadata.all.samples |>  dplyr::filter(!is.na(batch.2021) & is.na(batch.2022))) == 0)


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(batch = ifelse(is.na(batch.2021),batch.2022, batch.2021)) |> 
  dplyr::mutate(in.batch.2021 = !is.na(batch.2021)) |> 
  dplyr::mutate(batch.2021 = NULL) |> 
  dplyr::mutate(batch.2022 = NULL) |> 
  dplyr::mutate(excluded = FALSE) |> 
  dplyr::mutate(excluded.reason = NA)


# find and exclude replicates
# glass.gbm.rnaseq.metadata.all.samples |> 
#   dplyr::arrange(sample_barcode, sid) |> 
#   dplyr::group_by(sample_barcode) |> 
#   dplyr::filter(dplyr::n() > 1)


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(excluded.reason = ifelse(sid %in% c(
    "GLSS-CU-P053-R2-01R-RNA-2HP24A",
    "GLSS-CU-P101-R1-01R-RNA-FLI06Y",
    "GLSS-SM-R101-R1-01R-RNA-7ATA59",
    "GLSS-SM-R101-TP-01R-RNA-G7G4Q5",
    "GLSS-SM-R107-R1-01R-RNA-ID07M4",
    "GLSS-SM-R110-TP-01R-RNA-RSRC7U"
  ),"replicate sample", excluded.reason  )) |> 
  dplyr::mutate(excluded = !is.na(excluded.reason))


# ensure no replicate samples exist
stopifnot(
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(excluded == F) |> 
  dplyr::filter(duplicated(sample_barcode)) |> 
  nrow() == 0
)


## subtypes ----
# subtype data from Synapse portal 7 jan 2021
tmp.subtype.2021 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/GLASS.GBM.subtypes.from.Synapse.portal.tsv', header=T,stringsAsFactors = F) |>  # https://www.synapse.org/#!Synapse:syn21441635/tables/
  dplyr::select(c('aliquot_barcode','subtype')) |> 
  dplyr::mutate(subtype = as.factor(subtype)) |> 
  dplyr::rename(GBM.transcriptional.subtype.Synapse.2021 = subtype)


# new subtype data from Synapse portal 21 jul 2022
tmp.subtype.2022 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/GLASS.GBM.subtypes.from.Synapse.portal.all_samples.tsv', header=T,stringsAsFactors = F) |>
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::filter(p_value == min(p_value)) |> 
  #dplyr::filter(enrichment_score == min(enrichment_score)) |> 
  dplyr::ungroup() |> 
  dplyr::select(c('aliquot_barcode','signature_name')) |> 
  dplyr::mutate(signature_name = as.factor(signature_name)) |> 
  dplyr::rename(GBM.transcriptional.subtype.Synapse.2022 = signature_name)



# x-check same names are used consistently throughout
stopifnot(levels(tmp.subtype.2021$GBM.transcriptional.subtype.Synapse.2021) == levels(tmp.subtype.2022$GBM.transcriptional.subtype.Synapse.2022))


### check subtype mismatches ----


dplyr::left_join(tmp.subtype.2021,tmp.subtype.2022,by=c('aliquot_barcode'='aliquot_barcode')) |> 
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021)) |> 
  dplyr::filter(as.character(GBM.transcriptional.subtype.Synapse.2021) != as.character(GBM.transcriptional.subtype.Synapse.2022)) |> 
  dplyr::mutate(pri = paste0(GBM.transcriptional.subtype.Synapse.2021 , " -to-> ", GBM.transcriptional.subtype.Synapse.2022)) |> 
  dplyr::select(aliquot_barcode, pri)
warning("significant discrepancies 2021/2022 subtype calling")



## clinical / histological ----



tmp.clinical.2021 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/glass_clinical_surgeries.txt',sep="\t",header=T,stringsAsFactors = F) |> 
  dplyr::mutate(ROW_ID=NULL, ROW_VERSION=NULL) |> 
  dplyr::mutate(sample_barcode = ifelse(sample_barcode == "", NA, sample_barcode)) |> 
  dplyr::filter(!is.na(sample_barcode)) |> 
  dplyr::left_join(tmp.batch.2021, by=c('sample_barcode'='sample_barcode')) |> 
  dplyr::mutate(treatment_alkylating_agent = dplyr::case_when(
    treatment_alkylating_agent == "t" ~ TRUE,
    treatment_alkylating_agent == "f" ~ TRUE,
    TRUE ~ NA,
  )) |> 
  dplyr::mutate(treatment_tmz = dplyr::case_when(
    treatment_tmz == "t" ~ TRUE,
    treatment_tmz == "f" ~ TRUE,
    TRUE ~ NA,
  ))

  
# check if all that have expression dat also have a clinical metadata entry:
stopifnot(tmp.batch.2021$sample_barcode.short %in% tmp.clinical.2021$sample_barcode)


# remove those that do have metadata but no expression data
tmp.clinical.2021 <- tmp.clinical.2021 |> 
    dplyr::filter(!is.na(sid))



  

# clinial 2022: syn31121219
tmp.clinical.2022 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/glass_clinical_surgeries.all_samples.txt',sep="\t",header=T,stringsAsFactors = F) |> 
  dplyr::mutate(ROW_ID=NULL, ROW_VERSION=NULL) |> 
  dplyr::mutate(sample_barcode = ifelse(sample_barcode == "", NA, sample_barcode)) |> 
  dplyr::filter(!is.na(sample_barcode)) |> 
  dplyr::left_join(tmp.batch.2022, by=c('sample_barcode'='sample_barcode')) |> 
  dplyr::mutate(treatment_alkylating_agent = dplyr::case_when(
    treatment_alkylating_agent == "1" ~ TRUE,
    treatment_alkylating_agent == "0" ~ TRUE,
    TRUE ~ NA,
  )) |> 
  dplyr::mutate(treatment_tmz = dplyr::case_when(
    treatment_tmz == "1" ~ TRUE,
    treatment_tmz == "0" ~ TRUE,
    TRUE ~ NA,
  ))


# check if all that have expression dat also have a clinical metadata entry:
stopifnot(tmp.batch.2022$sample_barcode.short %in% tmp.clinical.2022$sample_barcode)


# remove those that do have metadata but no expression data
tmp.clinical.2022 <- tmp.clinical.2022 |> 
  dplyr::filter(!is.na(sid))



# check if both tables have the same columns
stopifnot(colnames(tmp.clinical.2022) == colnames(tmp.clinical.2022))



tmp.clinical.2021 <- tmp.clinical.2021 |>
  dplyr::rename_with( ~ paste0( .x, ".2021"))

tmp.clinical.2022 <- tmp.clinical.2022 |>
  dplyr::rename_with( ~ paste0( .x, ".2022"))


tmp.clinical.merge <- 
  dplyr::full_join(
    tmp.clinical.2021,
    tmp.clinical.2022,
    by=c('sid.2021'='sid.2022')
  ) |> 
  dplyr::rename(sid = sid.2021)


rm(tmp.clinical.2021, tmp.clinical.2022) # cleanup



### check discrepancies 2021 / 2022 ----


# 3 discrepancies regarding histology
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(histology.2021 != histology.2022) |> 
#   dplyr::select(sid,
#                 case_barcode.2021,
#                 histology.2021, histology.2022, 
#                 idh_status.2021, idh_status.2022,
#                 codel_status.2021, codel_status.2022,
#                 grade.2021, grade.2022
#   )


# no strange identifier mix-ups
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(sample_barcode.2021 != sample_barcode.2022) |> 
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 sample_barcode.2021, sample_barcode.2022,
#                 histology.2021, histology.2022, 
#                 idh_status.2021, idh_status.2022,
#                 codel_status.2021, codel_status.2022,
#                 grade.2021, grade.2022
#   )


# 3 discrepancies in grading
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(grade.2021 != grade.2022) |> 
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 sample_barcode.2021, sample_barcode.2022,
#                 histology.2021, histology.2022, 
#                 idh_status.2021, idh_status.2022,
#                 codel_status.2021, codel_status.2022,
#                 grade.2021, grade.2022
#   )


# exactly identical
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(idh_status.2021 != idh_status.2022 ) |> 
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 sample_barcode.2021, sample_barcode.2022,
#                 histology.2021, histology.2022, 
#                 idh_status.2021, idh_status.2022,
#                 codel_status.2021, codel_status.2022,
#                 grade.2021, grade.2022
#   )


# exactly identical
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(codel_status.2021 != codel_status.2022 ) |> 
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 sample_barcode.2021, sample_barcode.2022,
#                 histology.2021, histology.2022, 
#                 idh_status.2021, idh_status.2022,
#                 codel_status.2021, codel_status.2022,
#                 grade.2021, grade.2022
#   )


# plenty but small discrepancies:
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(surgical_interval_mo.2021 != surgical_interval_mo.2022 ) |> 
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 surgical_interval_mo.2021 , surgical_interval_mo.2022 
#   )


# exactly identical
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(treatment_alkylating_agent.2021 != treatment_alkylating_agent.2022 ) |>
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 treatment_alkylating_agent.2021 , treatment_alkylating_agent.2022 
#   )


# exactly identical
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(treatment_tmz.2021 != treatment_tmz.2022 ) |>
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 treatment_tmz.2021 , treatment_tmz.2022
#   )


# quite striking differences in alternative chemo descriptions:
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(treatment_chemotherapy_other.2021 != treatment_chemotherapy_other.2022 ) |>
#   dplyr::select(sid,
#                 treatment_chemotherapy_other.2021 , treatment_chemotherapy_other.2022
#   )



## integrate ----


glass.gbm.rnaseq.metadata.all.samples
<- data.frame(sid = colnames(glass.gbm.rnaseq.expression.all.samples),  stringsAsFactors = F) |> 
  dplyr::mutate(sid = gsub(".","-",sid,fixed=T))  |>  # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sid = gsub("_1$","",sid,fixed=F)) |>  # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(sample_barcode = gsub("^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+)-.+$","\\1",sid)) |> 
  dplyr::mutate(pid = gsub("^(............).+$","\\1",sid)) |> 
  dplyr::arrange(pid) |> 
  dplyr::mutate(resection = as.factor(gsub("^.............(..).+$","\\1",sid)))  |>  # TP is primary tumour? https://github.com/fpbarthel/GLASS
  dplyr::mutate(dataset =  as.factor(gsub("^(....).+$","\\1",sid))) |> 
  dplyr::left_join(tmp.subtype.2021, by = c('sid' = 'aliquot_barcode' ), suffix=c('','')) |> 
  dplyr::left_join(tmp.subtype.2022, by = c('sid' = 'aliquot_barcode' ), suffix=c('','')) |> 
  dplyr::mutate(sid.label = gsub("^(.)...(-..-....-).(.).*$","\\1\\2\\3",sid) ) |> 
  dplyr::mutate(dataset = gsub("^(....).*$","\\1",sid))





  # dplyr::left_join(

  #     ,
  #   by=c('sample_barcode'='sample_barcode')) %>%
  # dplyr::mutate(condition = ifelse(resection == "TP","Primary","NotPrimary")) %>%
  # dplyr::mutate(condition = factor(condition, levels = c("Primary","NotPrimary") )) %>%
  # dplyr::filter(idh_status == "IDHwt") # exclude the three IDH mutants according to Synapse WGS/WES VCF files
  # 



# exclude samples from patients that had a no grade IV tumor
no.grade.iv <- glass.gbm.rnaseq.metadata %>%
  dplyr::filter(grade %in% c('II','III')) %>%
  dplyr::pull(pid)

glass.gbm.rnaseq.metadata <- glass.gbm.rnaseq.metadata %>%
  dplyr::filter(pid %in% no.grade.iv == F )

rm(no.grade.iv)


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







