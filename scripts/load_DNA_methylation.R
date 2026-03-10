#!/usr/bin/env R

# 1. idats

gsam.metadata.array_samples <-  list.files(path = "/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = T) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 206)# just freeze the number to terminate on unexpected behavior
    return(.)
  })() |> 
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  dplyr::filter(!grepl("/MET2017-126-014/", channel_green)) |> # stored there for historical reasons - IDH-mutant loss study
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 87)
    return(.)
  })()  |> 
  dplyr::rename_with( ~ paste0("array_", .x)) 




gsam.metadata.array_samples <-  list.files(path = "/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/", pattern = "_(Grn|Red).idat$", recursive = T) |>
  data.frame(filename = _) |>
  dplyr::mutate(filename = paste0("/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/", filename)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 206)# just freeze the number to terminate on unexpected behavior
    return(.)
  })() |> 
  assertr::verify(file.exists(filename)) |>
  dplyr::mutate(sentrix_id = gsub("^.+/([^/]+)_(Grn|Red)\\.idat$", "\\1", filename)) |>
  dplyr::mutate(channel = gsub("^.+_(Grn|Red)\\.idat$", "\\1", filename)) |>
  tidyr::pivot_wider(id_cols = sentrix_id, names_from = channel, values_from = c(filename)) |>
  dplyr::rename(channel_green = Grn) |>
  dplyr::rename(channel_red = Red) |> 
  dplyr::filter(!grepl("/MET2017-126-014/", channel_green)) |> # stored there for historical reasons - IDH-mutant loss study
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == 87)
    return(.)
  })()  |> 
  dplyr::rename_with( ~ paste0("array_", .x)) 
## link sample names ----


tmp <- rbind(
  read.csv("/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/STS/EPIC2023-387-plate1.csv", skip=8) |> 
    dplyr::filter(!is.na(Sentrix_ID)) |> 
    dplyr::filter(grepl("^GSAM", Sample_Name)) |> 
    dplyr::mutate(array_sentrix_id = paste0(Sentrix_ID, "_", Sentrix_Position )) |> 
    dplyr::mutate(resection = paste0("R",gsub("GSAM_...(.)_.+$","\\1",Sample_Name))) |> 
    dplyr::mutate(patient_id = gsub("GSAM_(...)._.+$","\\1", Sample_Name)) |> 
    dplyr::select(array_sentrix_id, resection, patient_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 4)
      return(.)
    })()
  ,
  read.csv("/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/MET2022-350-014/MET2022-350-014_IdH.csv", skip=8) |> 
    dplyr::filter(!is.na(Sentrix_ID)) |> 
    dplyr::filter(grepl("^GSAM", Sample_Name)) |> 
    dplyr::mutate(array_sentrix_id = paste0(Sentrix_ID, "_", Sentrix_Position )) |> 
    dplyr::mutate(resection = paste0("R",gsub("GSAM_...(.)_.+$","\\1",Sample_Name))) |> 
    dplyr::mutate(patient_id = gsub("GSAM_(...)._.+$","\\1", Sample_Name)) |> 
    dplyr::select(array_sentrix_id, resection, patient_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 75)
      return(.)
    })(),
  read.csv("/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/EPIC2023-387-plate2/STS/EPIC2023-387.csv", skip=8) |> 
    dplyr::filter(!is.na(Sentrix_ID)) |> 
    dplyr::filter(grepl("^GSAM", Sample_Name)) |> 
    dplyr::mutate(array_sentrix_id = paste0(Sentrix_ID, "_", Sentrix_Position )) |> 
    dplyr::mutate(resection = paste0("R",gsub("GSAM_...(.)_.+$","\\1",Sample_Name))) |> 
    dplyr::mutate(patient_id = gsub("GSAM_(...)._.+$","\\1", Sample_Name)) |> 
    dplyr::select(array_sentrix_id, resection, patient_id) |> 
    (function(.) {
      print(dim(.))
      assertthat::assert_that(nrow(.) == 8)
      return(.)
    })()
) |> 
  dplyr::mutate(IDH = patient_id %in% c("BAW","CAV","CBG","CDF","DAB","EAF","EBD","ECB","FAD","FAL","JAB","JAD","JAF","KAC")) |> # IDH mut according to "data/gsam/output/tables/dna/idh_mutations.txt" - EAF also by MNP brain Classifier
  dplyr::mutate(resection_id = paste0(patient_id, gsub("^R","",resection))) |> 
  assertr::verify(!duplicated(resection_id)) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == nrow(gsam.metadata.array_samples))
    return(.)
  })()


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(resection_id))
rm(tmp)


gsam.metadata.array_samples <- gsam.metadata.array_samples |>
  dplyr::mutate(patient_reason_excluded = ifelse(patient_id == "AAA", "SNP based sample mismatch for this patient", NA))



tmp <- read.csv('/mnt/neuro-genomic-1-ro/gsam/administratie/GSAM_combined_clinical_molecular.csv',stringsAsFactors=F) |> 
  dplyr::rename(patient_id = studyID) |> 
  dplyr::arrange(patient_id) |> 
  dplyr::mutate(initialMGMT = NULL) |> 
  dplyr::mutate(gender = ifelse(patient_id %in% c('AAT', 'AAM', 'AZH', 'HAI', 'FAG'),"Male",gender)) |>  # there's a number of samples of which the gender does not fit with the omics data - omics data determined genders are the corrected ones
  dplyr::mutate(gender = as.factor(gender)) |> 
  dplyr::mutate(survival.events = dplyr::case_when(
    status == "Deceased" ~ 1,
    status == "Censored" ~ 0,
    T ~ as.numeric(NA))) |> 
  dplyr::rename(os.event = survival.events) |> 
  dplyr::mutate(survival.months = survivalDays / 365.0 * 12.0) |> 
  dplyr::select(patient_id, gender, 
                status, os.event,
                survivalDays, 
                progressionFreeDays,
                survivalFromSecondSurgeryDays
  ) |> 
  dplyr::filter(patient_id %in% gsam.metadata.array_samples$patient_id) |> 
  (function(.) {
    print(dim(.))
    assertthat::assert_that(nrow(.) == length(unique(gsam.metadata.array_samples$patient_id)))
    return(.)
  })() |> 
  dplyr::rename_with( ~ paste0("patient_", .x), .cols=!matches("^patient_id$",perl = T))


gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('patient_id'='patient_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(patient_gender))
rm(tmp)





# 3. MNP STP


parse_mnp_predictMGMT_csv <- function(fn, prefix) {
  
  a <- read.csv(fn,header=T) |> 
    dplyr::mutate(X=NULL, Status = NULL) |>  # also returns status if it is inconfident 
    dplyr::mutate(status = dplyr::case_when(
      Estimated < Cutoff & CI_Lower < Cutoff & CI_Upper < Cutoff ~ "unmethylated",
      Estimated > Cutoff & CI_Lower > Cutoff & CI_Upper > Cutoff ~ "methylated",
      T ~ as.character(NA)
    )) |> 
    dplyr::rename_with( ~ paste0(prefix, .x)) 
  
  return(a)
  
}




tmp.ls <- list.files(path = "/home/youri/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", pattern = "_mgmt.csv", recursive = TRUE)


tmp <- tmp.ls |> 
  data.frame(tmp_heidelberg_mgmt_report = _) |> 
  dplyr::mutate(tmp_heidelberg_mgmt_report = paste0("/home/youri/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131/", tmp_heidelberg_mgmt_report)) |>
  assertr::verify(file.exists(tmp_heidelberg_mgmt_report)) |> 
  dplyr::mutate(array_sentrix_id = gsub("^.+[^0-9]([0-9]{10,12}_[A-Z][0-9]+[A-Z][0-9]+).+$", "\\1", tmp_heidelberg_mgmt_report)) |>
  assertr::verify(!is.na(array_sentrix_id))|> 
  assertr::verify(!duplicated(array_sentrix_id)) |> 
  assertr::verify(array_sentrix_id %in% gsam.metadata.array_samples$array_sentrix_id) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(tmp = parse_mnp_predictMGMT_csv(tmp_heidelberg_mgmt_report, "array_mnp_MGMT_")) |>
  dplyr::ungroup() |> 
  tidyr::unnest(tmp) |> 
  dplyr::mutate(tmp_heidelberg_mgmt_report = NULL)



gsam.metadata.array_samples <- gsam.metadata.array_samples |> 
  dplyr::left_join(tmp, by=c('array_sentrix_id'='array_sentrix_id'), suffix=c('','')) |> 
  assertr::verify(!is.na(array_mnp_MGMT_Estimated))
rm(tmp)



# 5. SVVL primary ----

library(ggplot2)


plt <- gsam.metadata.array_samples |> 
  dplyr::filter(IDH ==F) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  
  dplyr::filter(resection == "R1") |> 
  dplyr::mutate(event = ifelse(patient_os.event, "Yes", "No"))


p1 = ggplot(plt, aes(x=reorder(resection_id, patient_survivalDays), y=patient_survivalDays, shape=event, col=array_mnp_MGMT_status)) +
  geom_point() +
  scale_shape_manual(values=c('Yes' = 19, 'No'=3)) +
  theme_bw()



# 5. SVVL recurrence ----


library(ggplot2)


plt <- gsam.metadata.array_samples |> 
  dplyr::filter(IDH ==F) |> 
  dplyr::filter(is.na(patient_reason_excluded)) |> 
  
  dplyr::filter(resection == "R2") |> 
  dplyr::mutate(event = ifelse(patient_os.event, "Yes", "No"))


p2 = ggplot(plt, aes(x=reorder(resection_id, patient_survivalFromSecondSurgeryDays, shape=patient_os.event), y=patient_survivalFromSecondSurgeryDays, col=array_mnp_MGMT_status)) +
  geom_point() +
  scale_shape_manual(values=c('Yes' = 19, 'No'=3)) +
  theme_bw()


patchwork:::`/.ggplot`(p1, p2)




