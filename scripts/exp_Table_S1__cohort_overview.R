#!/usr/bin/env R


# load libs ----


library(ggplot2)


# load data ----

source("scripts/R/palette.R")

source("scripts/load_G-SAM_metadata.R")
source("scripts/load_G-SAM_expression_data.R")


source("scripts/load_GLASS_data.R") # glass & tcga validation set


# G-SAM ----

## general stats ----

# all RNA included patients:
gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |>  # er zijn patienten met 1 goede en 1 slechte resectie
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()


# all RNA included patients w/ TMZ treatment:
all <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |>  # er zijn patienten met 1 goede en 1 slechte resectie
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::pull(pid) |> 
  unique()


gsam.patient.metadata |>
  dplyr::filter(studyID %in% all) |> 
  dplyr::filter(treatedWithTMZ == "Yes") |> 
  dim()

rm(all)


# resection biopsy stats - Primary and recurrence
gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |>  # er zijn patienten met 1 goede en 1 slechte resectie
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::filter(extent == "Biopsy") |> 
  dplyr::pull(sid) 

gsam.rna.metadata |> 
  #dplyr::filter(blacklist.pca == F) |>  # er zijn patienten met 1 goede en 1 slechte resectie
  #dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('EAV1', 'ECH2', 'ECI1', 'GAA2', 'GAQ1', 'JAK2')) |> 
  dplyr::select(sid, pid, extent, tumour.percentage.dna, blacklist.pca, pat.with.IDH)


# experimental therapies?
all <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |>  # er zijn patienten met 1 goede en 1 slechte resectie
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::pull(pid) |> 
  unique()


gsam.patient.metadata |>
  dplyr::filter(studyID %in% all) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = recode(otherTreatmentsBeforeSecondSurgery,
                                                            ' '='-',
                                                            'no'='-',
                                                            'NO'='-',
                                                            'none'='-',
                                                            'Surgery for rectum carcinoma'='-'
                                                                         )) |> 
  dplyr::pull(.data$otherTreatmentsBeforeSecondSurgery) |> 
  table()



## F] Table S1  ----
### G-SAM ----

tmp.gsam <- gsam.rna.metadata |>
  dplyr::filter(blacklist.pca == F) |> # er zijn patienten met 1 goede en 1 slechte resectie
  dplyr::filter(pat.with.IDH == F) |>
  dplyr::filter(sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
  dplyr::arrange(pid, as.character(resection)) |>
  dplyr::mutate(pid = as.factor(as.character(pid))) |>
  dplyr::left_join(gsam.patient.metadata |>
    dplyr::select(
      studyID,
      treatedWithTMZ, treatedWithRT, HM,
      mgmtPrimary, mgmtRecurrent,
      bevacizumab.before.recurrence,
      
      # svvl related
      status, survivalDays, progressionFreeDays, survivalFromSecondSurgeryDays,
      
      # treatment related
      otherTreatmentsBeforeSecondSurgery, treatmentDetailsFirstPD, treatmentDetailsFourthPD, treatmentDetailsSecondPD, treatmentDetailsThirdPD
    ),
  by = c("pid" = "studyID"), suffix = c("", "")
  ) |>
  dplyr::mutate(MGMT = gsub("ylated", "", ifelse(resection == "r1", mgmtPrimary, mgmtRecurrent))) |>
  dplyr::mutate(IDH = ifelse(pat.with.IDH, "+", "-")) |>
  dplyr::rename(TMZ = treatedWithTMZ) |>
  dplyr::rename(RT = treatedWithRT) |>
  dplyr::rename(`ssGSEA subtype` = `ssGSEA.2022.subtype`) |>
  dplyr::rename(`ssGSEA CL enrichment score` = `ssGSEA.2022.Classical.enrichment_score`) |>
  dplyr::rename(`ssGSEA CL pval` = `ssGSEA.2022.Classical_pval`) |>
  dplyr::rename(`ssGSEA MES enrichment score` = `ssGSEA.2022.Mesenchymal.enrichment_score`) |>
  dplyr::rename(`ssGSEA MES pval` = `ssGSEA.2022.Mesenchymal_pval`) |>
  dplyr::rename(`ssGSEA PN enrichment score` = `ssGSEA.2022.Proneural.enrichment_score`) |>
  dplyr::rename(`ssGSEA PN pval` = `ssGSEA.2022.Proneural_pval`) |>
  dplyr::rename(`NMF:150 meta-feature 1` = `NMF:150:1`) |>
  dplyr::rename(`NMF:150 meta-feature 2` = `NMF:150:2`) |>
  dplyr::rename(`NMF:150 meta-feature 3` = `NMF:150:3`) |>
  dplyr::rename(`NMF:150 PC1` = `NMF:150:PC1`) |>
  dplyr::rename(`NMF:150 PC2` = `NMF:150:PC2`) |>
  dplyr::rename(`NMF:150 PC1 (scaled)` = `NMF:150:PC1.n`) |>
  dplyr::rename(`NMF:150 PC2 (scaled)` = `NMF:150:PC2.n`) |>
  dplyr::rename(`NMF:150 PC1-2 (scaled) eucledian dist` = `NMF:150:PCA:eucledian.dist`) |>
  dplyr::rename(`GITS NMF:150 subtype` = `GITS.150.svm.2022.subtype`) |>
  # dplyr::mutate(GlioVis.Maj = ifelse(tumour.percentage.dna < 15, NA ,GlioVis.Maj)) |>
  dplyr::mutate(
    sig.C0.fuz = round(rna.signature.C0.fuzzy.2022, 2),
    sig.C1.col = round(rna.signature.C1.collagen.2022, 2),
    sig.C2.end = round(rna.signature.C2.endothelial.2022, 2),
    sig.C3.oli = round(rna.signature.C3.oligodendrocyte.2022, 2),
    sig.C4.neu = round(rna.signature.C4.neuron.2022, 2)
  ) |>
  dplyr::rename(
    svvl.stat = status,
    svvl.r1.days = survivalDays,
    pfs.days = progressionFreeDays,
    svvl.r2.days = survivalFromSecondSurgeryDays
  ) |>
  dplyr::mutate(angio = ifelse(bevacizumab.before.recurrence, "Bevacizumab", "No")) |>
  dplyr::mutate(purity = ifelse(is.na(tumour.percentage.dna), NA, paste0(tumour.percentage.dna, "%"))) |>
  dplyr::select(
    sid, pid, resection,
    extent, purity,
    IDH, HM, MGMT,
    TMZ, RT, angio,
    svvl.stat, svvl.r1.days, pfs.days, svvl.r2.days,

    # subtype
    `ssGSEA subtype`,
    `ssGSEA CL enrichment score`, `ssGSEA CL pval`,
    `ssGSEA MES enrichment score`, `ssGSEA MES pval`,
    `ssGSEA PN enrichment score`, `ssGSEA PN pval`,
    `NMF:150 meta-feature 1`, `NMF:150 meta-feature 2`, `NMF:150 meta-feature 3`,
    `NMF:150 PC1`, `NMF:150 PC2`,
    `NMF:150 PC1 (scaled)`, `NMF:150 PC2 (scaled)`,
    `NMF:150 PC1-2 (scaled) eucledian dist`,
    `GITS NMF:150 subtype`, `GITS travel segments`,

    # epic
    `EPIC: B-cells`,
    `EPIC: CAFs`,
    `EPIC: CD4 T-cells`,
    `EPIC: CD8 T-cells`,
    `EPIC: Endothelial`,
    `EPIC: Macrophages`,
    `EPIC: NK-cells`,
    `EPIC: other cells`,

    # signatures
    sig.C0.fuz, sig.C1.col, sig.C2.end, sig.C3.oli, sig.C4.neu,

    # therapy (extra)
    otherTreatmentsBeforeSecondSurgery
  ) |>
  dplyr::mutate(sid = gsub("-new", ".n", sid)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = tolower(otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub(":[ ]+[0-9]+\\.[0-9]+\\.[0-9]+[ \\-]+[0-9]+.[0-9]+.[0-9]+","",otherTreatmentsBeforeSecondSurgery)) |> # includes date |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("^re-surgery [0-9\\.]+: no tumor, no tissue available; (.+) at first pd \\(.+\\)","\\1",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("^yes ","",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("^\\(","",otherTreatmentsBeforeSecondSurgery)) |>  
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("\\)","",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("\\)","",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub(" after progression under tmz","",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub(" 2 cycles","",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub(" leuven","",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("lomustina","lomustine",otherTreatmentsBeforeSecondSurgery)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub(", "," & ",otherTreatmentsBeforeSecondSurgery, fixed=T)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("+","&",otherTreatmentsBeforeSecondSurgery, fixed=T)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("&"," & ",otherTreatmentsBeforeSecondSurgery, fixed=T)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("[ ]+\\&[ ]+"," & ",otherTreatmentsBeforeSecondSurgery, fixed=F)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("resection & ","",otherTreatmentsBeforeSecondSurgery, fixed=T)) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = dplyr::recode(otherTreatmentsBeforeSecondSurgery,
                                                                   'no'='',
                                                                   'none'='',
                                                                   're-irradiation'='', # standard of care
                                                                   'temozolomide'='', # standard of care
                                                                   'temsirolimius'='temsirolimus',
                                                                   'lomustine bevacizumab'='lomustine & bevacizumab'
  )) |> 
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = gsub("resection & ","",otherTreatmentsBeforeSecondSurgery, fixed=T))



head(tmp.gsam)



### GLASS ----

tmp.glass <- glass.gbm.rnaseq.metadata.all.samples |>
  dplyr::rename(`ssGSEA subtype` = `ssGSEA.2022.subtype`) |>
  dplyr::rename(`ssGSEA CL enrichment score` = `ssGSEA.2022.Classical.enrichment_score`) |>
  dplyr::rename(`ssGSEA CL pval` = `ssGSEA.2022.Classical_pval`) |>
  dplyr::rename(`ssGSEA MES enrichment score` = `ssGSEA.2022.Mesenchymal.enrichment_score`) |>
  dplyr::rename(`ssGSEA MES pval` = `ssGSEA.2022.Mesenchymal_pval`) |>
  dplyr::rename(`ssGSEA PN enrichment score` = `ssGSEA.2022.Proneural.enrichment_score`) |>
  dplyr::rename(`ssGSEA PN pval` = `ssGSEA.2022.Proneural_pval`) |>
  
  dplyr::mutate(`ssGSEA subtype` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA subtype` , '')) |> 
  dplyr::mutate(`ssGSEA CL enrichment score` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA CL enrichment score` , '')) |> 
  dplyr::mutate(`ssGSEA CL pval` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA CL pval` , '')) |> 
  dplyr::mutate(`ssGSEA MES enrichment score` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA MES enrichment score` , '')) |> 
  dplyr::mutate(`ssGSEA MES pval` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA MES pval` , '')) |> 
  dplyr::mutate(`ssGSEA PN enrichment score` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA PN enrichment score` , '')) |> 
  dplyr::mutate(`ssGSEA PN pval` = ifelse(tumour.percentage.2022 >= 15, `ssGSEA PN pval` , '')) |> 
  
  dplyr::rename(`NMF:150 meta-feature 1` = `NMF:150:1`) |>
  dplyr::rename(`NMF:150 meta-feature 2` = `NMF:150:2`) |>
  dplyr::rename(`NMF:150 meta-feature 3` = `NMF:150:3`) |>
  dplyr::rename(`NMF:150 PC1` = `NMF:150:PC1`) |>
  dplyr::rename(`NMF:150 PC2` = `NMF:150:PC2`) |>
  dplyr::rename(`NMF:150 PC1 (scaled)` = `NMF:150:PC1.n`) |>
  dplyr::rename(`NMF:150 PC2 (scaled)` = `NMF:150:PC2.n`) |>
  dplyr::rename(`NMF:150 PC1-2 (scaled) eucledian dist` = `NMF:150:PCA:eucledian.dist`) |>
  dplyr::rename(`GITS NMF:150 subtype` = `GITS.150.svm.2022.subtype`) |>
  dplyr::mutate(sig.C1.col = round(rna.signature.C1.collagen.2022, 2)) |> 
  dplyr::mutate(purity = paste0(round(tumour.percentage.2022),"%")) |> 
  dplyr::rename(purity.source = tumour.percentage.2022.source) |> 
  dplyr::select(
    aliquot_barcode,
    purity, purity.source,
    `ssGSEA subtype`,
    `ssGSEA CL enrichment score`, `ssGSEA CL pval`,
    `ssGSEA MES enrichment score`, `ssGSEA MES pval`,
    `ssGSEA PN enrichment score`, `ssGSEA PN pval`,
    `NMF:150 meta-feature 1`, `NMF:150 meta-feature 2`, `NMF:150 meta-feature 3`,
    `NMF:150 PC1`, `NMF:150 PC2`,
    `NMF:150 PC1 (scaled)`, `NMF:150 PC2 (scaled)`,
    `NMF:150 PC1-2 (scaled) eucledian dist`,
    `GITS NMF:150 subtype`,
    sig.C1.col
  )



### export ----

openxlsx::createWorkbook(
  creator = "Dr. Youri Hoogstrate",
  title = "Clinical information G-SAM study",
  subject = "Clinical information",
  category = "G-SAM study"
) -> wb
openxlsx::addWorksheet(wb, "Sheet1 - Sample info G-SAM")
openxlsx::writeDataTable(wb, sheet = "Sheet1 - Sample info G-SAM", x = tmp.gsam)
openxlsx::addWorksheet(wb, "Sheet1 - Sample info GLASS")
openxlsx::writeDataTable(wb, sheet = "Sheet1 - Sample info GLASS", x = tmp.glass)

openxlsx::saveWorkbook(wb, file = "output/tables/2022_Table_S1_clinical_information.xlsx", overwrite = T)

rm(wb, tmp.gsam, tmp.glass)


