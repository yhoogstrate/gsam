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
#   tibble::remove_rownames() |> 
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
  tibble::remove_rownames() |> 
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
  tibble::rownames_to_column('aliquot_barcode') |> 
  dplyr::mutate(V1 = NULL) |> 
  dplyr::mutate(aliquot_barcode = gsub(".","-",  aliquot_barcode,fixed=T)) |> # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(aliquot_barcode = gsub("_1$","", aliquot_barcode,fixed=F)) |> # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(batch = paste0(gsub("^([^\\-]+).+$","\\1", aliquot_barcode,fixed=F),":2021"))


tmp.batch.2022 <- 'data/gsam/data/GLASS_GBM_R1-R2/transcript_count_matrix_all_samples.tsv' |> 
  read.delim(stringsAsFactors = F,nrows=1) |> 
  dplyr::mutate(target_id=NULL, length=NULL) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('aliquot_barcode') |> 
  dplyr::mutate(V1 = NULL) |> 
  dplyr::mutate(aliquot_barcode = gsub(".",  "-", aliquot_barcode,fixed=T)) |> # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(aliquot_barcode = gsub("_1$","",  aliquot_barcode,fixed=F)) |> # by convention [https://github.com/fpbarthel/GLASS]
  dplyr::mutate(batch = paste0(gsub("^([^\\-]+).+$","\\1",aliquot_barcode,fixed=F),":2022"))



glass.gbm.rnaseq.metadata.all.samples <- tmp.batch.2022 |> 
  dplyr::full_join(tmp.batch.2021, by=c('aliquot_barcode'='aliquot_barcode'), suffix=c('.2022','.2021')) |> 
  dplyr::mutate(sample_barcode = gsub('^([^\\-]+-[^\\-]+-[^\\-]+-[^\\-]+).+$','\\1',aliquot_barcode))


rm(tmp.batch.2021, tmp.batch.2022)


#'@details [v] ensure that all the 2021 samples are in the 2022 samples
stopifnot(nrow(glass.gbm.rnaseq.metadata.all.samples |>  dplyr::filter(!is.na(batch.2021) & is.na(batch.2022))) == 0)


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(batch = ifelse(is.na(batch.2021),batch.2022, batch.2021)) |> 
  dplyr::mutate(in.batch.2021 = !is.na(batch.2021)) |> 
  dplyr::mutate(batch.2021 = NULL) |> 
  dplyr::mutate(batch.2022 = NULL)



#'@details [v] check if no strange suffixes are present
stopifnot(stringr::str_detect(glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode,"_1$", negate=T) )

#'@details [v] check if no dots are present
stopifnot(grepl(".",glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode,fixed=T) == F)

#'@details [v] check if no dots are present
stopifnot(grepl(".",glass.gbm.rnaseq.metadata.all.samples$sample_barcode,fixed=T) == F)



## subtypes ----


# subtype data from Synapse portal 7 jan 2021
tmp.subtype.2021 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/GLASS.GBM.subtypes.from.Synapse.portal.tsv', header=T,stringsAsFactors = F) |>  # https://www.synapse.org/#!Synapse:syn21441635/tables/
  dplyr::select(c('aliquot_barcode','subtype')) |> 
  dplyr::mutate(subtype = as.factor(subtype)) |> 
  dplyr::rename(GBM.transcriptional.subtype.Synapse = subtype)
stopifnot(duplicated(tmp.subtype.2021$aliquot_barcode) == FALSE)


# new subtype data from Synapse portal 21 jul 2022
tmp.subtype.2022 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/GLASS.GBM.subtypes.from.Synapse.portal.all_samples.tsv', header=T,stringsAsFactors = F) |>
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::filter(p_value == min(p_value)) |> 
  #dplyr::filter(enrichment_score == min(enrichment_score)) |> 
  dplyr::ungroup() |> 
  dplyr::select(c('aliquot_barcode','signature_name')) |> 
  dplyr::mutate(signature_name = as.factor(signature_name)) |> 
  dplyr::rename(GBM.transcriptional.subtype.Synapse = signature_name) |> 
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::summarise(GBM.transcriptional.subtype.Synapse = as.factor(paste(GBM.transcriptional.subtype.Synapse, collapse="|"))) |> 
  dplyr::ungroup()
stopifnot(nrow(tmp.subtype.2022$aliquot_barcode) == 425)
stopifnot(duplicated(tmp.subtype.2022$aliquot_barcode) == FALSE) # - some have equal prob's so should have been merged/collapsed


#'@details [-] x-check same subtype names are used consistently throughout
#stopifnot(levels(tmp.subtype.2021$GBM.transcriptional.subtype.Synapse) == levels(tmp.subtype.2022$GBM.transcriptional.subtype.Synapse))
warning('2022 subtypes contain inconfident calls with equal probabilities')


#'@details [-] check if both tables have the same columns - No - different pipelines?
# stopifnot(colnames(tmp.subtype.2021) == colnames(tmp.subtype.2022))
warning('Format of 2021 and 2022 subtype data differs')


tmp.subtype.merge <- dplyr::full_join(tmp.subtype.2021,tmp.subtype.2022,
                                      by=c('aliquot_barcode'='aliquot_barcode'),
                                      suffix = c('.2021','.2022'))


rm(tmp.subtype.2021, tmp.subtype.2022)



### check subtype mismatches ----


#'@details [v] check if all 2021 labels are in the 2022 data
stopifnot(
  (tmp.subtype.merge |> 
    dplyr::filter(is.na(GBM.transcriptional.subtype.Synapse.2022)) |> 
    nrow()) == 0)

#'@details [-] check actual discordant labels
# tmp.subtype.merge |> 
#   dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021)) |> 
#   dplyr::mutate(equal.prob = grepl("|",GBM.transcriptional.subtype.Synapse.2022,fixed=T)) |> 
#   dplyr::mutate(mismatch = stringr::str_detect(as.character(GBM.transcriptional.subtype.Synapse.2022), as.character(GBM.transcriptional.subtype.Synapse.2021)) == F) |> 
#   dplyr::mutate(discrepancy.type = dplyr::case_when(mismatch ~ "discordant call", equal.prob ~ "uncertain call", T ~ "" )) |> 
#   dplyr::filter(discrepancy.type != "") |>
#   dplyr::mutate(label.switch = paste0(GBM.transcriptional.subtype.Synapse.2021 , " (2021) -to-> ", GBM.transcriptional.subtype.Synapse.2022, " (2022)")) |> 
#   dplyr::select(aliquot_barcode, discrepancy.type, label.switch)
warning("considerable discrepancies between 2021 & 2022 subtype calling")



## clinical / histological ----


tmp.clinical.2021 <- 'data/gsam/data/GLASS_GBM_R1-R2/glass_clinical_surgeries.txt' |> 
  read.table(sep="\t",header=T,stringsAsFactors = F) |> 
  dplyr::mutate(ROW_ID=NULL) |> 
  dplyr::mutate(ROW_VERSION=NULL) |> 
  dplyr::mutate(sample_barcode = ifelse(sample_barcode == "", NA, sample_barcode)) |> 
  dplyr::filter(!is.na(sample_barcode)) |> 
  dplyr::mutate(treatment_alkylating_agent = dplyr::case_when(
    treatment_alkylating_agent == "t" ~ TRUE,
    treatment_alkylating_agent == "f" ~ TRUE,
    TRUE ~ NA )) |> 
  dplyr::mutate(treatment_tmz = dplyr::case_when(
    treatment_tmz == "t" ~ TRUE,
    treatment_tmz == "f" ~ TRUE,
    TRUE ~ NA )) |> 
  dplyr::mutate(histology = as.factor(ifelse(histology == "", NA, histology))) |> 
  dplyr::mutate(grade = as.factor(ifelse(grade == "", NA, grade))) |> 
  dplyr::mutate(idh_status = as.factor(ifelse(idh_status == "", NA, idh_status))) |> 
  dplyr::mutate(codel_status = as.factor(ifelse(codel_status == "", NA, codel_status)))


#'@details [v] check if all that have expression dat also have a clinical metadata entry:
stopifnot(glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(in.batch.2021) |> 
  #dplyr::filter(!excluded) |> # not really necessary
  dplyr::pull(sample_barcode) %in% tmp.clinical.2021$sample_barcode)



# clinial 2022: syn31121219
tmp.clinical.2022 <- read.table('data/gsam/data/GLASS_GBM_R1-R2/glass_clinical_surgeries.all_samples.txt',sep="\t",header=T,stringsAsFactors = F) |> 
  dplyr::mutate(ROW_ID=NULL) |> 
  dplyr::mutate(ROW_VERSION=NULL) |> 
  dplyr::mutate(sample_barcode = ifelse(sample_barcode == "", NA, sample_barcode)) |> 
  dplyr::filter(!is.na(sample_barcode)) |> 
  dplyr::mutate(treatment_alkylating_agent = dplyr::case_when(
    treatment_alkylating_agent == "1" ~ TRUE,
    treatment_alkylating_agent == "0" ~ TRUE,
    TRUE ~ NA )) |> 
  dplyr::mutate(treatment_tmz = dplyr::case_when(
    treatment_tmz == "1" ~ TRUE,
    treatment_tmz == "0" ~ TRUE,
    TRUE ~ NA )) |> 
  dplyr::mutate(histology = as.factor(ifelse(histology == "", NA, histology))) |> 
  dplyr::mutate(grade = as.factor(ifelse(grade == "", NA, grade))) |> 
  dplyr::mutate(idh_status = as.factor(ifelse(idh_status == "", NA, idh_status))) |> 
  dplyr::mutate(codel_status = as.factor(ifelse(codel_status == "", NA, codel_status)))



#'@details [v] ensure no strange duplicate identifiers exist
stopifnot(duplicated(tmp.clinical.2022$sample_barcode) == FALSE)


#'@details [v] check if all that have expression dat also have a clinical metadata entry:
stopifnot(glass.gbm.rnaseq.metadata.all.samples$sample_barcode %in% tmp.clinical.2022$sample_barcode)


#'@details [v] check if both tables have the same columns
stopifnot(colnames(tmp.clinical.2021) == colnames(tmp.clinical.2022))


# merge
tmp.clinical.merge <- dplyr::full_join(
    tmp.clinical.2021, tmp.clinical.2022,
    by=c('sample_barcode'='sample_barcode'),
    suffix=c('.2021','.2022')
  )
rm(tmp.clinical.2021, tmp.clinical.2022) # cleanup


#'@details [v] ensure that case_barcode's (patient identifiers) match
stopifnot(
tmp.clinical.merge |>
  dplyr::filter(!is.na(case_barcode.2021) & !is.na(case_barcode.2022)) |>
  dplyr::filter(case_barcode.2021 != case_barcode.2022) |> 
  nrow() == 0)


# redundant colums
tmp.clinical.merge <- tmp.clinical.merge |> 
  dplyr::mutate(case_barcode.2021 = NULL) |> 
  dplyr::mutate(case_barcode.2022 = NULL)


colnames(tmp.clinical.merge)



### check discrepancies 2021 / 2022 ----


#'@details [-] 3 discrepancies regarding histology
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


#'@details [v] no strange identifier mix-ups
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


#'@details [-] 3 discrepancies in grading
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


#'@details [v] IDH status exactly identical
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


#'@details [v] codel status exactly identical
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


#'@details [v] plenty but small surgical interval discrepancies - rounding errors?:
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(surgical_interval_mo.2021 != surgical_interval_mo.2022 ) |> 
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 surgical_interval_mo.2021 , surgical_interval_mo.2022 
#   )


#'@details [v] alk treatment: exactly identical
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(treatment_alkylating_agent.2021 != treatment_alkylating_agent.2022 ) |>
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 treatment_alkylating_agent.2021 , treatment_alkylating_agent.2022 
#   )


#'@details [v] TMZ: exactly identical
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(treatment_tmz.2021 != treatment_tmz.2022 ) |>
#   dplyr::select(sid,
#                 sample_barcode.2022,
#                 case_barcode.2021,
#                 treatment_tmz.2021 , treatment_tmz.2022
#   )


#'@details [-] alt chemo: quite striking differences in alternative chemo descriptions:
# tmp.clinical.merge |> 
#   dplyr::filter(!is.na(histology.2021)) |> 
#   dplyr::filter(treatment_chemotherapy_other.2021 != treatment_chemotherapy_other.2022 ) |>
#   dplyr::select(sid,
#                 treatment_chemotherapy_other.2021 , treatment_chemotherapy_other.2022
#   )
warning("descriptions of alternative chemotherapy show considerable differences - columns are excluded")

tmp.clinical.merge <- tmp.clinical.merge |> 
  dplyr::mutate(treatment_chemotherapy_other.2021 = NULL,
                treatment_chemotherapy_other.2022 = NULL)



## integrate ----


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(sid.label = gsub("^(.)...(-..-....-).(.).*$","\\1\\2\\3", aliquot_barcode)) |> 
  dplyr::mutate(case_barcode = gsub("^(............).+$","\\1", aliquot_barcode)) |> 
  dplyr::mutate(resection = factor(gsub("^.............(..).+$","\\1", aliquot_barcode), levels=c('TP','R1','R2','R3','R4'))) |>  # TP is primary tumour? https://github.com/fpbarthel/GLASS
  dplyr::mutate(is.primary = resection == "TP") |> 
  dplyr::arrange(case_barcode, resection)



#'@details [v] 425/425 subtypes in expression data
stopifnot(tmp.subtype.merge$aliquot_barcode %in% glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode)

#'@details [v] 425/425 expression data samples in subtype data
stopifnot(glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode %in% tmp.subtype.merge$aliquot_barcode)


#'@details [v] check if there is clinical data for all 425/425
stopifnot(glass.gbm.rnaseq.metadata.all.samples$sample_barcode %in% tmp.clinical.merge$sample_barcode)
# stopifnot(tmp.clinical.merge$sample_barcode %in% glass.gbm.rnaseq.metadata.all.samples$sample_barcode) there is clinical data beyond those of which expression data is available


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::left_join(tmp.subtype.merge, by = c('aliquot_barcode' = 'aliquot_barcode' ), suffix=c('','')) |> 
  dplyr::left_join(tmp.clinical.merge, by = c('sample_barcode' = 'sample_barcode' ), suffix=c('','')) |> 
  dplyr::mutate(condition = factor(ifelse(resection == "TP","Primary","NotPrimary"), levels = c("Primary","NotPrimary") ))
rm(tmp.subtype.merge, tmp.clinical.merge)

stopifnot(!is.na(glass.gbm.rnaseq.metadata.all.samples$GBM.transcriptional.subtype.Synapse.2022))



### add exclusion labels ----


# the 'excluded' labels should yet be absent
stopifnot("excluded" %in% colnames(glass.gbm.rnaseq.metadata.all.samples) == FALSE)

glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(excluded = strsplit("",":")) # quickest way to get a list of empty characters?


#### replicates ----

# manually find replicates
# find and exclude replicates
# glass.gbm.rnaseq.metadata.all.samples |>
#   dplyr::arrange(sample_barcode, aliquot_barcode) |>
#   dplyr::group_by(sample_barcode) |>
#   dplyr::filter(dplyr::n() > 1)


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(excluded = ifelse(aliquot_barcode %in% c(
      "GLSS-CU-P053-R2-01R-RNA-2HP24A",
      "GLSS-CU-P101-R1-01R-RNA-FLI06Y",
      "GLSS-SM-R101-R1-01R-RNA-7ATA59",
      "GLSS-SM-R101-TP-01R-RNA-G7G4Q5",
      "GLSS-SM-R107-R1-01R-RNA-ID07M4",
      "GLSS-SM-R110-TP-01R-RNA-RSRC7U" ), 
    lapply(excluded, function(x) { return(unique(c(x,"replicate sample"))) }), excluded))


# x-check
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(as.logical(lapply(excluded , function(x){return("replicate sample" %in% x)}))) |> 
  dplyr::select(aliquot_barcode, excluded)



#'@details [v] ensure replicates are marked as such
stopifnot(
  glass.gbm.rnaseq.metadata.all.samples |> 
    dplyr::filter(lapply(excluded , length) == 0) |> 
    dplyr::filter(duplicated(sample_barcode)) |> 
    nrow() == 0
)
stopifnot(
  glass.gbm.rnaseq.metadata.all.samples |> 
    dplyr::filter(lapply(excluded , length) > 0) |> 
    nrow() == 6
)


#### codel status ----


# use any(), to mark all samples that belong to the same patient in one go
glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::mutate(drop = dplyr::case_when(
    any(codel_status.2022 == "codel") ~ "Y",
    any(codel_status.2021 == "codel") ~ "Y",
    any(!is.na(codel_status.2021) & !is.na(codel_status.2022) & codel_status.2021 != codel_status.2022) ~ "Y", # possible 2021/2022 discrepancies
    T ~ "N"
  )) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(excluded = ifelse(drop != "N", lapply(excluded, function(x) { return(unique(c(x,"1p/19q-codel"))) }), excluded)) |> 
  dplyr::mutate(drop=NULL)


glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(as.logical(lapply(excluded , function(x){return("1p/19q-codel" %in% x)})) | codel_status.2021 == "codel" | codel_status.2022 == "codel") |> 
  dplyr::select(aliquot_barcode, codel_status.2021, codel_status.2022, excluded) |> 
  dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |> 
  as.data.frame()



#### IDH status ----


# use any(), to mark all samples that belong to the same patient in one go
glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::mutate(drop1 = dplyr::case_when( # drop if there is any IDH mutant sample that belongs to a given patient
    any(idh_status.2022 == "IDHmut") ~ "Y",
    any(idh_status.2021 == "IDHmut") ~ "Y",
    T ~ "N"
  )) |>
  dplyr::mutate(drop2 = dplyr::case_when( # drop if all samples that belong to a given patient have no IDH status or are discordant
    all(is.na(idh_status.2021) & is.na(idh_status.2022)) ~ "Y",
    any(!is.na(idh_status.2021) & !is.na(idh_status.2022) & idh_status.2021 != idh_status.2022) ~ "Y", # possible 2021/2022 discrepancies
    T ~ "N"
  )) |>
  dplyr::ungroup()


# overview
# summary(as.factor(glass.gbm.rnaseq.metadata.all.samples$drop1))
# summary(as.factor(glass.gbm.rnaseq.metadata.all.samples$drop2))


# append exclusion labels
glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(excluded = ifelse(drop1 != "N", lapply(excluded, function(x) { return(unique(c(x,"IDH-mut"))) }), excluded)) |> 
  dplyr::mutate(excluded = ifelse(drop2 != "N", lapply(excluded, function(x) { return(unique(c(x,"IDH-status N/A"))) }), excluded)) |> 
  dplyr::mutate(drop1=NULL, drop2=NULL)



# x-check / view [idh-mut]
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(idh_status.2021 == "IDHmut" | idh_status.2022 == "IDHmut" | (as.logical(lapply(excluded,function(x){return("IDH-mut" %in% x)})))) |> 
  dplyr::select(aliquot_barcode, idh_status.2021, idh_status.2022, excluded) |> 
  dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |> 
  as.data.frame()

# x-check / view [idh-mut] - nice imputed example
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(case_barcode == "GLSS-19-0279") |> 
  dplyr::select(aliquot_barcode, idh_status.2021, idh_status.2022, excluded) |> 
  dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |> 
  as.data.frame()



# x-check / view [idh N/A]
glass.gbm.rnaseq.metadata.all.samples |>
  dplyr::filter(as.logical(lapply(excluded,function(x){return("IDH-status N/A" %in% x)}))) |>
  dplyr::select(aliquot_barcode, histology.2022, who_classification.2022, idh_status.2021, idh_status.2022, excluded) |>
  dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |>
  as.data.frame()




#### grade -----


# use any(), to mark all samples that belong to the same patient in one go
glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::mutate(drop1 = dplyr::case_when( # drop if there is any grade 2/3 sample that belongs to a given patient
    any(grade.2022 %in% c("II","III")) ~ "Y1",
    any(grade.2021 %in% c("II","III")) ~ "Y2",
    T ~ "N"
  )) |>
  dplyr::mutate(drop2 = dplyr::case_when( # drop if all samples that belong to a given patient have no grading or are discordant
    all( is.na(grade.2021) &  is.na(grade.2022)) ~ "Y1",
    any(!is.na(grade.2021) & !is.na(grade.2022) & grade.2021 != grade.2022) ~ "Y2", # possible 2021/2022 discrepancies
    T ~ "N"
  )) |>
  dplyr::ungroup()


# summary(as.factor(glass.gbm.rnaseq.metadata.all.samples$drop1))
# summary(as.factor(glass.gbm.rnaseq.metadata.all.samples$drop2))


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(excluded = ifelse(drop1 != "N", lapply(excluded, function(x) { return(unique(c(x,"Grade II/III"))) }), excluded)) |> 
  dplyr::mutate(excluded = ifelse(drop2 != "N", lapply(excluded, function(x) { return(unique(c(x,"Grade N/A"))) }), excluded)) |> 
  dplyr::mutate(drop1=NULL, drop2 = NULL)


# x-check / view [grade]
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(idh_status.2021 %in% c("II","III") | idh_status.2022 %in% c("II","III") | (as.logical(lapply(excluded,function(x){return("Grade II/III" %in% x)})))) |> 
  dplyr::select(aliquot_barcode, grade.2021, grade.2022, excluded) |> 
  dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |> 
  as.data.frame()



#### histology ----


# use any(), to mark all samples that belong to the same patient in one go
glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::mutate(drop1 = dplyr::case_when( # drop if there is any non-GBM sample that belongs to a given patient
    any(histology.2022 != "Glioblastoma") ~ "Y1",
    any(histology.2021 != "Glioblastoma") ~ "Y2",
    T ~ "N"
  )) |>
  dplyr::mutate(drop2 = dplyr::case_when( # drop if all samples that belong to a given patient have no histology or are discordant
    all( is.na(histology.2021) &  is.na(histology.2022)) ~ "Y1",
    any(!is.na(histology.2021) & !is.na(histology.2022) & as.character(histology.2021) != as.character(histology.2022)) ~ "Y2", # possible 2021/2022 discrepancies
    T ~ "N"
  )) |>
  dplyr::ungroup()


#summary(as.factor(glass.gbm.rnaseq.metadata.all.samples$drop1))
# summary(as.factor(glass.gbm.rnaseq.metadata.all.samples$drop2))


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::mutate(excluded = ifelse(drop1 != "N", lapply(excluded, function(x) { return(unique(c(x,"Histology not Glioblastoma"))) }), excluded)) |> 
  dplyr::mutate(excluded = ifelse(drop2 != "N", lapply(excluded, function(x) { return(unique(c(x,"Histology N/A"))) }), excluded)) |> 
  dplyr::mutate(drop1=NULL, drop2 = NULL)



# x-check / view [histology]
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(histology.2021 != "Glioblastoma" | histology.2022 != "Glioblastoma" | (as.logical(lapply(excluded,function(x){return("Histology not Glioblastoma" %in% x)})))) |> 
  dplyr::select(aliquot_barcode, histology.2021, histology.2022, excluded) |> 
  dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |> 
  as.data.frame()



#### additional / todo ----


#'@todo make a definitive call about samples with "NOS" in the WHO classification (IDH mut not specified)



#### last: take out middle resections ----


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::mutate(excluded.tmp = (lapply(excluded, length) == 0) & !is.primary & as.numeric(resection) != max(as.numeric(resection))) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(excluded = ifelse(excluded.tmp, lapply(excluded, function(x) { return(unique(c(x,"Middle resection"))) }), excluded)) |> 
  dplyr::mutate(excluded.tmp = NULL)



# x-check / view [middle resections]
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::filter(any(as.logical(lapply(excluded,function(x){return("Middle resection" %in% x)})))) |> 
  dplyr::select(aliquot_barcode, case_barcode, resection, excluded) |> 
  as.data.frame()


glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(case_barcode == "GLSS-CU-P020") |> 
  #dplyr::filter(case_barcode == "GLSS-HF-753F") |>
  #dplyr::filter(case_barcode == "GLSS-MD-0017") |>
  dplyr::select(case_barcode, resection, excluded) |> 
  as.data.frame()


## check w/ VCF mutations ----


#'@details all that have an IDH mutation in the VCF are also marked as such

# # file from Synapse portal
# #tmp <- read.csv('data/gsam/data/GLASS_GBM_R1-R2/variants_passgeno_20190327.csv',stringsAsFactors = F) |> 
# tmp <- read.csv('data/gsam/data/GLASS_GBM_R1-R2/variants_passgeno_20220531.csv',stringsAsFactors = F) |> 
#   dplyr::mutate(aliquot_barcode = as.factor(aliquot_barcode)) |> 
#   dplyr::mutate(case_barcode = as.factor(gsub("^([^\\-]+.[^\\-]+.[^\\-]+).+$","\\1",aliquot_barcode)))
# 
# # idh1 => chr2:209,098,951-209,121,867 [R132]
# tmp.idh1 <- tmp |> 
#   dplyr::filter(chrom == "2" & start >=  209113111-1 & end <= 209113111+2+1 ) |> # hg19 R132
#   #dplyr::filter(chrom == "2" & start >=  208248387-1 & end <= 208248389+2 ) |> # hg38
#   dplyr::filter(ad_alt != 0) |> 
#   dplyr::filter(ssm2_pass_call == 't') |> 
#   dplyr::arrange(start, aliquot_barcode)
# 
# 
# # idh2 => chr15:90,623,937-90,648,932 [R140/R172]
# tmp.idh2 <- tmp %>%
#   dplyr::filter((chrom == "15" & start >=  90088701-1 & end <= 90088703+1 ) | (chrom == "15" & start >=  90088605-1 & end <= 90088607+1 ) ) %>%
#   dplyr::filter(ad_alt != 0) %>%
#   dplyr::filter(ssm2_pass_call == 't') %>%
#   dplyr::arrange(start, aliquot_barcode)
# 
# 
# glass.gbm.rnaseq.metadata.all.samples |> 
#   dplyr::filter(case_barcode %in% (tmp.idh1$case_barcode)) |> 
#   dplyr::select(aliquot_barcode,resection, excluded)  |> 
#   dplyr::mutate(excluded = as.character(lapply(excluded, paste, collapse=", "))) |> 
#   as.data.frame()
# 
# 


rm(tmp,tmp.idh1,tmp.idh2)


## final check for each variable ----

### metadata ----

colnames(glass.gbm.rnaseq.metadata.all.samples)

stopifnot(duplicated(glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode) == F)
stopifnot(glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode %in% colnames(glass.gbm.rnaseq.expression.all.samples))



glass.gbm.rnaseq.metadata.all.samples |>
  dplyr::filter(as.numeric(lapply(excluded, length)) == 0) |> 
  dplyr::pull(case_barcode) |> 
  unique() |> 
  length()

### expression data ----


# plot(sort(colSums(glass.gbm.rnaseq.expression.all.samples) / 1000000)) + abline(h=0.75)



#### old stuff (?) ----

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


glass.gbm.rnaseq.metadata.all.samples <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(lapply(excluded, length) == 0) # erase all that have reasons to be excluded


glass.gbm.rnaseq.expression.all.samples <- glass.gbm.rnaseq.expression.all.samples %>%
  dplyr::select(glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode) # match table with metadata



stopifnot(colnames(glass.gbm.rnaseq.expression.all.samples) == glass.gbm.rnaseq.metadata.all.samples$aliquot_barcode)




# VST transform ----


stopifnot(!exists('glass.gbm.rnaseq.expression.vst'))


tmp <- as.factor(paste0('c',round(runif(ncol(glass.gbm.rnaseq.expression.all.samples)))+1) )

glass.gbm.rnaseq.expression.all.samples.vst <- DESeq2::DESeqDataSetFromMatrix(glass.gbm.rnaseq.expression.all.samples, data.frame(cond = tmp), ~cond) |> 
  DESeq2::vst(blind=T) |> 
  SummarizedExperiment::assay() |> 
  as.data.frame(stringsAsFactors=F)


## QC ----

### PCA ----


plt <- glass.gbm.rnaseq.expression.all.samples.vst |> 
  t() |> 
  prcomp() |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8) |> 
  tibble::rownames_to_column('aliquot_barcode') |> 
  dplyr::left_join(glass.gbm.rnaseq.metadata.all.samples, by=c('aliquot_barcode'='aliquot_barcode'), suffix=c('','')) |> 
  dplyr::mutate(project = gsub("^([^\\-]+)-[^\\-]+-.+$","\\1",aliquot_barcode)) |>  # project
  dplyr::mutate(tissue.source = gsub("^[^\\-]+-([^\\-]+)-.+$","\\1",aliquot_barcode)) # tissue.source


# extreme batch effects
ggplot(plt, aes(x=PC1, y=PC2, label=sid.label, col=project)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  facet_grid(cols = vars(tissue.source))






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

#'@todo check w/ CD163 and TMEM144 expression

# Combine table into paired data table ----







