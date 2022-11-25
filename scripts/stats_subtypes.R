#!/usr/bin/env

# load data ----


if(!exists('gsam.rna.metadata')) {
  source('scripts/load_G-SAM_metadata.R')
}

if(!exists('glass.gbm.rnaseq.metadata.all.samples')) {
  source('scripts/load_GLASS_data.R')
}



# input labels ----


## G-SAM ----

### GlioVis MJ -> ssGSEA 2022 ----

tmp <- gsam.rna.metadata  |> 
  dplyr::filter(!is.na(`gliovis.majority_call`) & !is.na(`ssGSEA.2022.subtype`)) |> 
  dplyr::select(`gliovis.majority_call` , `ssGSEA.2022.subtype`) |> 
  dplyr::mutate(textual = paste0(`gliovis.majority_call` ," => ", `ssGSEA.2022.subtype`))

stopifnot(nrow(tmp) == 287)


table(tmp$gliovis.majority_call)
table(tmp$ssGSEA.2022.subtype)


table(tmp$textual)


sum(tmp$gliovis.majority_call == "Classical")
sum(tmp$ssGSEA.2022.subtype == "Classical")

sum(tmp$gliovis.majority_call == "Mesenchymal")
sum(tmp$ssGSEA.2022.subtype == "Mesenchymal")

sum(tmp$gliovis.majority_call == "Proneural")
sum(tmp$ssGSEA.2022.subtype == "Proneural")



rm(tmp)


## GLASS ----
### Synapse 2022 x ssGSEA 2022 ----

tmp <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |> 
  dplyr::filter(!is.na(ssGSEA.Synapse.subtype.2022) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::select(ssGSEA.Synapse.subtype.2022, ssGSEA.2022.subtype) |> 
  dplyr::mutate(textual = paste0(`ssGSEA.Synapse.subtype.2022` ," => ", `ssGSEA.2022.subtype`))

stopifnot(nrow(tmp) == 216)

table(tmp$textual)


### Synapse 2021 x ssGSEA 2022  ----

tmp <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |> 
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::select(GBM.transcriptional.subtype.Synapse.2021, ssGSEA.2022.subtype) |> 
  dplyr::mutate(textual = paste0(`GBM.transcriptional.subtype.Synapse.2021` ," => ", `ssGSEA.2022.subtype`))

stopifnot(nrow(tmp) == 52)

table(tmp$textual)

sum(tmp$GBM.transcriptional.subtype.Synapse.2021 == "Classical")
sum(tmp$ssGSEA.2022.subtype == "Classical")

sum(tmp$GBM.transcriptional.subtype.Synapse.2021 == "Mesenchymal")
sum(tmp$ssGSEA.2022.subtype == "Mesenchymal")

sum(tmp$GBM.transcriptional.subtype.Synapse.2021 == "Proneural")
sum(tmp$ssGSEA.2022.subtype == "Proneural")




### Synapse 2021 x Synapse 2022 ----

# tmp <- glass.gbm.rnaseq.metadata.all.samples |> 
#   dplyr::filter(tumour.percentage.2022 >= 15) |> 
#   dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021) & !is.na(ssGSEA.Synapse.subtype.2022)) |> 
#   dplyr::select(GBM.transcriptional.subtype.Synapse.2021, ssGSEA.Synapse.subtype.2022) |> 
#   dplyr::mutate(textual = paste0(`GBM.transcriptional.subtype.Synapse.2021` ," => ", `ssGSEA.Synapse.subtype.2022`))
# 
# stopifnot(nrow(tmp) == 52)
# 
# table(tmp$textual)





# output labels ----
## G-SAM ----
# stats over verandering in finale GITS subtype labels

tmp <- gsam.rna.metadata  |> 
  dplyr::filter(!is.na(`NMF.123456.PCA.SVM.class_2021`) & !is.na(`GITS.150.svm.2022.subtype`)) |> 
  dplyr::select(`NMF.123456.PCA.SVM.class_2021` , `GITS.150.svm.2022.subtype`) |> 
  dplyr::mutate(textual = paste0(`NMF.123456.PCA.SVM.class_2021` ," => ", `GITS.150.svm.2022.subtype`))

stopifnot(nrow(tmp) == 287)

table(tmp$textual)


rm(tmp)


## GITS x ssGSEA ----


tmp.1 <- gsam.rna.metadata |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::mutate(textual = paste0(`ssGSEA.2022.subtype` ," => ", `GITS.150.svm.2022.subtype`))
table(tmp.1$textual)

(2+1+4+2+3) / (125+2+1+4+90+2+3+56)


tmp.2 <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::mutate(textual = paste0(`ssGSEA.2022.subtype` ," => ", `GITS.150.svm.2022.subtype`))
table(tmp.2$textual)


# stats ----


tmp <- gsam.rna.metadata |> 
  dplyr::filter(!is.na(ssGSEA.2022.subtype)) |> 
  
  dplyr::filter(!is.na(gliovis.majority_call)) |> 
  dplyr::filter(!is.na(NMF.123456.PCA.SVM.class_2021)) |> 
  dplyr::filter(!is.na(gliovis.gsea_call)) |> 
  dplyr::filter(!is.na(gliovis.knn_call)) |>  
  dplyr::filter(!is.na(gliovis.svm_call))
  
  


sum(tmp$NMF.123456.PCA.SVM.class_2021 != tmp$gliovis.majority_call) # 24

sum(tmp$gliovis.gsea_call != tmp$gliovis.majority_call)
sum(tmp$gliovis.knn_call != tmp$gliovis.majority_call)
sum(tmp$gliovis.svm_call != tmp$gliovis.majority_call)




sum(tmp$gliovis.majority_call != tmp$gliovis.gsea_call)
sum(tmp$gliovis.knn_call != tmp$gliovis.gsea_call)
sum(tmp$gliovis.svm_call != tmp$gliovis.gsea_call)



# expression MES in bulk ----

gids = results.out |>
  dplyr::filter(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-MES") |> 
  dplyr::pull(gid)


gids = results.out |>
  dplyr::filter(hugo_symbol == "VIM") |> 
  dplyr::pull(gid)


sids =   gsam.rna.metadata |>
  
  dplyr::filter(blacklist.pca == F) %>%
  dplyr::filter(pat.with.IDH == F) %>%
  dplyr::filter(
    sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
  ) %>%
  dplyr::filter(tumour.percentage.dna >= 15) |> 
  dplyr::filter(resection == "r1") |> 
  dplyr::pull(sid)

data = gsam.gene.expression.all |> 
  dplyr::select(sids) |> 
  tibble::rownames_to_column('gid') |> 
  dplyr::filter(.data$gid %in% gids) |> 
  tibble::column_to_rownames('gid')


data$min = as.numeric(lapply(data.frame(t(data)), min) )

data |>
       dplyr::arrange(min) |> 
       dplyr::relocate(min, .before = tidyselect::everything()) |> 
  View()



sort(as.numeric(lapply(data.frame(data), median)))


# C1 x MES-NMF ----


tmp.stats <- gsam.rna.metadata |>
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F) |> 
  dplyr::filter(tumour.percentage.dna >= 15)


cor(tmp.stats$rna.signature.C1.collagen.2022, tmp.stats$`NMF:150:3`)
#cor(tmp.stats$rna.signature.C1.collagen.2022, tmp.stats$ssGSEA.2022.Mesenchymal.enrichment_score)



tmp.stats <- gsam.rna.metadata |>
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F) |> 
  dplyr::filter(tumour.percentage.dna >= 15) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::select(pid, resection, GITS.150.svm.2022.subtype, rna.signature.C1.collagen.2022) |> 
  tidyr::pivot_wider(names_from = resection, values_from = c(GITS.150.svm.2022.subtype, rna.signature.C1.collagen.2022)) |> 
  dplyr::mutate(delta.c1 = rna.signature.C1.collagen.2022_r2 - rna.signature.C1.collagen.2022_r1) |> 
  dplyr::mutate(mes.transition = ifelse(
    GITS.150.svm.2022.subtype_r2 == "Mesenchymal" & 
      GITS.150.svm.2022.subtype_r1 != "Mesenchymal" , "Yes" ,"No"
  ))

wilcox.test(tmp.stats |> 
    dplyr::filter(mes.transition == "Yes") |>
    dplyr::pull(delta.c1) ,
  tmp.stats |> 
    dplyr::filter(mes.transition == "No") |> 
    dplyr::pull(delta.c1))


# samples

gsam.rna.metadata |>
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F) |> 
  dplyr::filter(tumour.percentage.dna >= 15) |> 
  dplyr::pull(sid) |> 
  unique() |> 
  length()



gsam.rna.metadata |>
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F) |> 
  #dplyr::filter(tumour.percentage.dna >= 15) |> 
  dplyr::pull(sid) |> 
  unique() |> 
  length()


glass.gbm.rnaseq.metadata.all.samples |> 
  #dplyr::filter(tumour.percentage.2022 >= 15) |> 
  dplyr::pull(case_barcode) |> 
  unique() |> 
  length()

glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |> 
  dplyr::group_by(case_barcode) |> 
  dplyr::filter(n() >= 2) |> 
  dplyr::ungroup( ) |> 
  dplyr::pull(case_barcode) |> 
  unique() |> 
  length()






glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022.source %in% c("CNV segment profiles","RNA-imputation")) |> 
  nrow()



glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022.source %in% c("CNV segment profiles")) |> 
  nrow()


glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022.source == "RNA-imputation") |> 
  nrow()


plot(
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022.source %in% c("CNV segment profiles")) |> 
  dplyr::pull(tumor.purity.cnv.pct.2022)
,
glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022.source %in% c("CNV segment profiles")) |> 
  dplyr::pull(tumour.percentage.dna.imputed.rf.2022.all.patients.B),
)

summary(
  lm(
    `tumour.percentage.dna.imputed.rf.2022.all.patients.B` ~ `tumor.purity.cnv.pct.2022`,
    data=glass.gbm.rnaseq.metadata.all.samples |> 
      dplyr::filter(tumour.percentage.2022.source %in% c("CNV segment profiles")) 
  )
)$r.squared



summary(
lm(
  `tumour.percentage.dna.imputed.rf.2022.all.patients.B` ~ `tumor.purity.cnv.pct.2022`,
  data=glass.gbm.rnaseq.metadata.all.samples |> 
    dplyr::filter(tumour.percentage.2022.source %in% c("CNV segment profiles")) |> 
    dplyr::mutate(`tumour.percentage.dna.imputed.rf.2022.all.patients.B` = `tumor.purity.cnv.pct.2022` +  runif(210))
)
)$r.squared




# EX / INH ----

df1 <- data.frame(hugo_symbol = c(
  "SATB2",
  "SLC17A7",
  "SV2B",
  "NRGN",
  "CHN1",
  "MLIP",
  "RALYL",
  "KIAA1211L",
  "ENC1",
  "ARPP21",
  "KCNIP4",
  "NPTX1",
  "LDB2",
  "HS3ST4",
  "GAD1",
  "SLC6A1",
  "ERBB4",
  "GAD2",
  "ARX",
  "QKI",
  "IGF1",
  "GRIK1",
  "GRIP2",
  "ADARB2",
  "MAF",
  "SPOCK3",
  "GRIP1",
  "ADRA1A",
  "ANKRD55",
  "SLC24A3",
  "DLX1",
  "POU6F2",
  "SLC35F1",
  "RORB",
  "RXFP1",
  "MGLL",
  "GRM3",
  "ADCY8",
  "HS3ST2",
  "PTPRD",
  "NTNG1",
  "THEMIS",
  "PDZRN4",
  "STAC2",
  "COL5A2",
  "TRPC3",
  "RASSF4",
  "COL6A1",
  "ANXA1",
  "KIAA1456",
  "TRHDE",
  "DYSF",
  "COL12A1",
  "PARM1",
  "LHFP",
  "FNBP1L",
  "SEMA6A",
  "CPNE7",
  "VIPR2",
  "ALDH1A1",
  "VWC2",
  "COBLL1",
  "FAM198B",
  "GRM4",
  "PALM2",
  "NEUROD6",
  "SCN9A",
  "ARHGAP12",
  "EXTL2",
  "SYNPR",
  "GRIK3",
  "DPY19L1",
  "MRAP2",
  "PLCXD3",
  "SEMA3E",
  "PHLDB2",
  "COL11A1",
  "EFHD2",
  "TLL1",
  "SMYD2",
  "TSHZ2",
  "ARHGAP29",
  "KLF3",
  "GPRIN3",
  "TSHZ3",
  "LRMP",
  "SPHKAP",
  "GPR161",
  "GREM2",
  "ANKRD29",
  "PLCXD2",
  "OCA2",
  "GOLIM4",
  "ITGA11",
  "CDK18",
  "KCNK2",
  "PLA2R1",
  "PTCHD4",
  "FGFR1",
  "LRRK1",
  "PAPPA",
  "MSRB3",
  "RP11-383H13.1",
  "ARSJ",
  "TMEFF2",
  "PDGFC",
  "IFITM10",
  "RASGEF1C",
  "CLMN",
  "CDH20",
  "OPRM1",
  "RAB3B",
  "HRK",
  "SAMD5",
  "PLAGL1",
  "MDFIC",
  "DPY19L3",
  "SLC8B1",
  "CXXC4",
  "GLDN",
  "GABRA5",
  "PAH",
  "LAMA4",
  "EMB",
  "TTC9B",
  "CPNE9",
  "GABRG1",
  "CNTNAP4",
  "ACOT11",
  "ANKRD18A",
  "ETV1",
  "SPTSSB",
  "KBTBD12",
  "TDO2",
  "HPCAL1",
  "LAMP5",
  "PDZD2",
  "CBLN2",
  "CUX2",
  "SERPINE2",
  "SLIT3",
  "CDH9",
  "STXBP6",
  "LRRTM4",
  "CCDC88C",
  "PDGFD",
  "EPHA6",
  "RGS12",
  "SLC24A4",
  "C1QL2",
  "ANO2",
  "PDE7B",
  "PCDH8",
  "SH3RF2",
  "GPR83",
  "GRIK4",
  "COL27A1",
  "ADAMTS3",
  "THSD7A",
  "FAM163A",
  "NECAB2",
  "NPTX2",
  "ITPKA",
  "GLP2R",
  "PCSK6",
  "TESPA1",
  "DGKA",
  "LONRF3",
  "KCNN3",
  "ITGA4",
  "TEK",
  "FAT4",
  "MDGA1",
  "MN1",
  "RGS20",
  "FBLN7",
  "ADAMTS8",
  "CNR1",
  "C1orf95",
  "ITGB8",
  "KCNIP2",
  "PPP4R4",
  "JUN",
  "FAM84A",
  "EXPH5",
  "ADAMTS16",
  "RAMP1",
  "CNGB1",
  "PRKCD",
  "ZNF831",
  "ATP7B",
  "ZNF516",
  "ZNF608",
  "SOWAHA",
  "TMEM132E",
  "MAS1",
  "HTR3B",
  "SLCO2A1",
  "KITLG",
  "COL5A1",
  "ART3",
  "PTCH1",
  "CPNE6",
  "ANKRD62",
  "ZFHX4",
  "GULP1",
  "LRRC2",
  "CRIP2",
  "PSRC1",
  "PCDH19",
  "PTCHD2",
  "CRLF1",
  "ARL4A",
  "IGSF3",
  "PRKCG",
  "WFS1",
  "IGFBP5",
  "KCNV1",
  "SLITRK2",
  "CUX2",
  "COBLL1",
  "RORB",
  "SLC38A11",
  "MET",
  "NTNG1",
  "CCDC168",
  "VIPR2",
  "COL5A2",
  "GPR26",
  "PDE7B",
  "GOLIM4",
  "THSD7A",
  "PLCXD3",
  "CBLN4",
  "IL1RAP",
  "VSTM2A",
  "COL6A1",
  "ANXA1",
  "TNNT2",
  "ALDH1A1",
  "C14orf132",
  "PARM1",
  "EPHA6",
  "HOPX",
  "ZBTB7C",
  "PTCHD4",
  "ZBTB18",
  "EXTL2",
  "PRDM8",
  "SLC17A6",
  "PTGFRN",
  "ARL4C",
  "PLCXD2",
  "SGK223",
  "TRPC3",
  "RAMP1",
  "FBLN7",
  "FAM84A",
  "ST8SIA4",
  "POU3F2",
  "IGSF3",
  "SCN1B",
  "NEFM",
  "SLC30A3",
  "RASSF4",
  "HSPA1B",
  "RTKN2",
  "DYNLL2",
  "SCN9A",
  "CCK",
  "PRPS2",
  "ADCY8",
  "CITED2",
  "KCNS1",
  "GRIK3",
  "LRP1B",
  "HS3ST4",
  "TLL1",
  "HS3ST2",
  "DPP10",
  "TMEFF2",
  "GRIK4",
  "ETV1",
  "SLIT3",
  "RXFP1",
  "THEMIS",
  "GALNT14",
  "PDZRN4",
  "DCHS2",
  "KLHL5",
  "ST3GAL1",
  "NETO2",
  "ROBO2",
  "GRIN3A",
  "PDE1A",
  "ASAP2",
  "CDH6",
  "ARHGEF28",
  "SLIT1",
  "SEMA3E",
  "MKX",
  "PDE1C",
  "ARSJ",
  "RASAL1",
  "DGKG",
  "KCNK2",
  "CBLN2",
  "RGS12",
  "CHST8",
  "RBM20",
  "SAMD5",
  "PTPRU",
  "ITGB8",
  "RASGEF1C",
  "TRABD2A",
  "CNTNAP3",
  "ARHGAP12",
  "CPNE7",
  "HTR1E",
  "PRKCG",
  "C1orf95",
  "SEC24D",
  "ANKRD30B",
  "PRRX1",
  "BCMO1",
  "PAPPA",
  "GNB4",
  "PELI1",
  "ANKRD18A",
  "GHR",
  "TMEM255B",
  "DYSF",
  "SEMA4A",
  "MEIS3",
  "PTCHD1",
  "XXYLT1",
  "CNNM4",
  "CNTN6",
  "PRAP1",
  "PDE9A",
  "PEPD",
  "XKR7",
  "NPY1R",
  "GLRA2",
  "CDH20",
  "SPHKAP",
  "CLMN",
  "TRPC3",
  "RAMP1",
  "COL5A2",
  "RP11-758M4.1",
  "TSHZ1",
  "PLAGL1",
  "TSHZ2",
  "SLIT1",
  "GABRG1",
  "PDE1A",
  "RP11-383H13.1",
  "CNR1",
  "SYT17",
  "NEUROD6",
  "FRMPD2",
  "CXXC4",
  "AFAP1L2",
  "PCDH20",
  "TPBG",
  "STXBP6",
  "CYP26B1",
  "SCN3B",
  "COL6A1",
  "VGF",
  "ATXN7L3B",
  "ID2",
  "PLCXD2",
  "BRPF3",
  "SCG2",
  "TMEM155",
  "DUSP26",
  "HSPA1B",
  "SEMA7A",
  "RERG",
  "KCNS1",
  "NAP1L3",
  "CHGB",
  "TSPAN7",
  "RP11-618P17.4",
  "ITM2A",
  "ANXA1",
  "SEMA5B",
  "SLC17A6",
  "KCNA1",
  "ALDH1A1",
  "TRPC3",
  "ST8SIA4",
  "TMEM35",
  "FAM198B",
  "SLC24A4",
  "FKBP5",
  "SLC25A11",
  "TCERG1L",
  "ROBO2",
  "PCED1B",
  "CPNE7",
  "COL5A2",
  "L3MBTL4",
  "COL11A1",
  "ST6GALNAC5",
  "CNGB1",
  "PRSS12",
  "OVOL2",
  "TMEFF2",
  "SCPEP1",
  "SLIT1",
  "SLC26A4",
  "ARHGAP12",
  "RASAL1",
  "BAIAP3",
  "SYNPR",
  "LRRC16B",
  "KCNN3",
  "OPRM1",
  "OCA2",
  "ANKRD18A",
  "ASAP2",
  "HTR1E",
  "OTOGL",
  "TIMP2",
  "ADAMTS9",
  "CCDC68",
  "RASGEF1C",
  "RP11-383H13.1",
  "NRBP2",
  "FNBP1L",
  "TTC40",
  "GULP1",
  "TNNT2",
  "PHLDB2",
  "GAS6",
  "PAH",
  "ZBBX",
  "CLMN",
  "ADAMTS3",
  "GPR116",
  "CBLN2",
  "LRRC2",
  "NPY1R",
  "ZNF471",
  "PRAP1",
  "XXYLT1",
  "CRYM",
  "RPH3AL",
  "SEMA4D",
  "GRM2",
  "RP11-758M4.1",
  "FAM19A5",
  "THUMPD2",
  "PLXNA1",
  "LRRIQ1",
  "LGI4",
  "ZNF737",
  "PIR",
  "HDAC7",
  "COL4A2",
  "SPHKAP",
  "TMEM196",
  "GOLIM4",
  "RXFP1",
  "KIF21B",
  "BEND4",
  "KIAA0907",
  "LHX2",
  "ZNF250",
  "MAOB",
  "GPR135",
  "SLC15A4",
  "SCN3B",
  "FAM173B",
  "PEPD",
  "GPT2",
  "TECTA",
  "PIK3C2G",
  "COL12A1",
  "CLDN16",
  "TEX30",
  "ALKBH1",
  "TXNRD3",
  "SLC19A1",
  "NSUN6",
  "BIVM",
  "GFM2",
  "KLF6",
  "PLS3",
  "ADCYAP1R1",
  "GPR63",
  "PKIA",
  "SYT17",
  "B3GNTL1",
  "BAALC",
  "ARSJ",
  "EFNB2",
  "GGCT",
  "CPNE6",
  "EPHA10",
  "ATP10D",
  "FXC1",
  "TRIM66",
  "GREM2",
  "RP11-637O19.3",
  "PACRGL",
  "FAT2",
  "CDC7",
  "PBX3",
  "DPY19L1",
  "SWAP70",
  "MYO5C",
  "OBSCN",
  "ZNF331",
  "SMYD2",
  "DGKG",
  "LZTS1",
  "CDH9",
  "GRIK3",
  "ROBO2",
  "TESPA1",
  "ITGB8",
  "ASIC2",
  "SEMA5B",
  "OLFML2B",
  "PDZRN4",
  "RGS12",
  "CA10",
  "SLC4A4",
  "ITGA8",
  "DCHS2",
  "SMAD3",
  "SLIT3",
  "HCRTR2",
  "RYR3",
  "MOXD1",
  "SULF1",
  "GFRA1",
  "CBLN2",
  "PDE9A",
  "PBX3",
  "SLC24A4",
  "ZNF462",
  "SPOCK3",
  "MDGA1",
  "PTPN13",
  "ADAMTS3",
  "SERPINE2",
  "GLDC",
  "RORB",
  "EPHA3",
  "FAM19A2",
  "IL1RAPL2",
  "LRRK1",
  "PDE1C",
  "PARM1",
  "COL11A1",
  "TLL1",
  "TDO2",
  "FSTL4",
  "SLC30A3",
  "TRABD2A",
  "CTNNA2",
  "PTCHD4",
  "COBLL1",
  "FSTL1",
  "PLCXD3",
  "KCNIP4",
  "RNF182",
  "TSHZ2",
  "ST6GALNAC5",
  "CCDC68",
  "LHFP",
  "RHBDL3",
  "ARHGAP29",
  "DNAJC5G",
  "ARSJ",
  "PTGFRN",
  "VWC2",
  "TCERG1L",
  "CDH20",
  "KIAA1456",
  "LHX2",
  "PCED1B",
  "PDGFRB",
  "PHLDB2",
  "FNDC5",
  "NTNG1",
  "CMTM8",
  "KIRREL3",
  "CDK18",
  "PTPRK",
  "KCNIP4",
  "PDZRN3",
  "NWD2",
  "TMEM132D",
  "CCK",
  "SLC24A4",
  "TLE3",
  "MGAT5B",
  "ZNF804B",
  "PTPRO",
  "ITGB8",
  "THEMIS",
  "CTNNA2",
  "RTN4RL1",
  "TLE1",
  "FRMD6",
  "NPTXR",
  "CBLN2",
  "NTNG2",
  "MOXD1",
  "PPP4R4",
  "RTN4RL2",
  "C1QL3",
  "CACNA1G",
  "BACE2",
  "POSTN",
  "COL27A1",
  "NTNG1",
  "TRIM54",
  "SLC30A3",
  "AMOT",
  "ITGA7",
  "DNAJC5G",
  "LHX2",
  "CHST8",
  "GPR153",
  "GALNT18",
  "MAN1C1",
  "SMYD1",
  "ANKRD30B",
  "SLC26A4",
  "RGS12",
  "GULP1",
  "SH2D4A",
  "PTPRQ",
  "SNED1",
  "GALNT14",
  "HEY1",
  "FSTL4",
  "KIRREL3",
  "GRASP",
  "PDE4C",
  "UNC5B",
  "PLA2G4A",
  "TRAF5",
  "PRLR",
  "GPR63",
  "C14orf23",
  "ADAMTS3",
  "SCN1B",
  "ANXA5",
  "MUM1L1",
  "CD63",
  "FMNL3",
  "RAVER2",
  "TNS3",
  "ZNF697",
  "ADAMTS8",
  "SLC20A1",
  "CPNE9",
  "COL5A3",
  "DGKA",
  "NPTX2",
  "SGK223",
  "TMEM196",
  "SEMA3E",
  "HS3ST4",
  "KIAA1456",
  "ASIC2",
  "MCC",
  "GRM4",
  "RYR3",
  "KLHL5",
  "C3orf70",
  "GHR",
  "MAP3K8",
  "TRPM3",
  "ITPR2",
  "SULF1",
  "VWA2",
  "NXPH3",
  "KIAA1211",
  "SPHKAP",
  "IL1RAP",
  "ROBO2",
  "PTCHD4",
  "KCNH5",
  "VANGL2",
  "RAB3B",
  "GRIK3",
  "ATP1B2",
  "MAPKAPK2",
  "COL11A1",
  "ITGA8",
  "RGMA",
  "ACOT11",
  "PDE4B",
  "RNF220",
  "AMOTL1",
  "LHFP",
  "KCNK2",
  "RBM20",
  "NT5DC3",
  "HIP1R",
  "SYT6",
  "CYP26B1",
  "ANKRD62",
  "NXPH2",
  "IMPA2",
  "MGAT4C",
  "SLC4A4",
  "HCRTR2",
  "GRIN3A",
  "ARSJ",
  "HPCAL1",
  "FRMD3",
  "HTR2C",
  "THEMIS2",
  "KIAA1644",
  "PDE9A",
  "SEL1L3",
  "SPTB",
  "FNBP1L",
  "LGR6",
  "ABO",
  "RHOJ",
  "GRIK3",
  "CTNNA2",
  "PDGFC",
  "ST6GALNAC5",
  "SLIT3",
  "TMEM108",
  "PTPRD",
  "GRIK4",
  "ARHGEF28",
  "LRP1B",
  "SLC30A3",
  "GRM3",
  "KIRREL3",
  "CNR1",
  "RGS4",
  "NKAIN2",
  "PRRX1",
  "OPCML",
  "SYT17",
  "PDE4B",
  "GULP1",
  "PCDH20",
  "LRRTM4",
  "PVRL1",
  "TRIM54",
  "TNS3",
  "PPP4R4",
  "TSHZ3",
  "TPO",
  "NPTX2",
  "RXRA",
  "OPRM1",
  "ABCG1",
  "SEL1L3",
  "SLC35A1",
  "ACTN2",
  "CHST8",
  "LRRC2",
  "DNAJC5G",
  "C4orf32",
  "ROBO2",
  "TMEM44",
  "CPNE8",
  "KBTBD2",
  "FSTL4",
  "ABHD13",
  "C1QL3",
  "SLITRK4",
  "TCERG1L",
  "ARMCX3",
  "IL17REL",
  "RAVER2",
  "WSCD2",
  "ADCY8",
  "KITLG",
  "FEM1C",
  "ANP32A",
  "SH2D1B",
  "CPNE9",
  "SLCO4C1",
  "PEX3",
  "ETV1",
  "HEY1",
  "TESPA1",
  "DYNLL2",
  "KCNG2",
  "CAMK1G",
  "SMYD2",
  "CNOT11",
  "AMER3",
  "TXNRD3",
  "WDR1",
  "PDE4D",
  "RP11-383H13.1",
  "PAPPA",
  "CHGA",
  "BAIAP3",
  "CHMP4B",
  "COL6A1",
  "RGS12",
  "POSTN",
  "SMYD1",
  "NTNG2",
  "SYNPR",
  "NR4A2",
  "MCTP2",
  "STK32B",
  "RNF152",
  "GFRA1",
  "SPOCK3",
  "KCNIP4",
  "OLFML2B",
  "PDLIM5",
  "B3GAT2",
  "HS3ST4",
  "GAS2L3",
  "LSAMP",
  "NTNG1",
  "RERG",
  "VWC2",
  "TGFA",
  "ATP10A",
  "ST8SIA2",
  "SCN7A",
  "PTGER3",
  "PRLR",
  "ADRA2A",
  "ITPKB",
  "TRPC7",
  "CDK18",
  "TMEM200A",
  "PTPRQ",
  "GDPD5",
  "RIN2",
  "RASSF3",
  "MAMDC2",
  "FBLN5",
  "ERCC4",
  "PLD5",
  "HGF",
  "CDH20",
  "LAMA4",
  "DGKD",
  "CNNM4",
  "PDLIM7",
  "ITGA11",
  "GLUL",
  "NHS",
  "CDC14A",
  "ABCC3",
  "SCN1B",
  "CSPG4",
  "NPNT",
  "TMEM132D",
  "PLA2G4A",
  "STEAP4",
  "TMEM196",
  "PDIA5",
  "ASPHD2",
  "CALB2",
  "RGN",
  "ELK3",
  "KCNG3",
  "SEPT10",
  "SPTBN5",
  "LRIG1",
  "B2M",
  "ZNF684",
  "C10orf11",
  "GNB4",
  "WNT2B",
  "SYT10",
  "DPYS",
  "TLL1",
  "HAAO",
  "GCNT4",
  "TPMT",
  "PLAGL1",
  "PPFIBP2",
  "FMNL3",
  "SPTB",
  "UBR7",
  "JUP",
  "RGS9",
  "GPR153",
  "FAP",
  "TRPC3",
  "NPY1R",
  "CYSTM1",
  "VSTM2A"
)) |> 
  dplyr::filter(!duplicated(hugo_symbol)) |> 
  dplyr::mutate(status="exhibitory")


df2 <- data.frame(hugo_symbol = c(
  "GRIK3",
  "OPRD1",
  "SPARCL1",
  "TSHZ3",
  "KIAA1211L",
  "DPP10",
  "SLC9A9",
  "PDE1A",
  "COL19A1",
  "LHX6",
  "ZBBX",
  "STXBP6",
  "SLC36A1",
  "RGS4",
  "CNTNAP3",
  "ST8SIA4",
  "TAC1",
  "RGS5",
  "BCL11A",
  "PDGFC",
  "BRINP3",
  "PDZRN4",
  "KIAA0226L",
  "PVALB",
  "MAN1A1",
  "SULF1",
  "SLC30A3",
  "NEFM",
  "C3orf58",
  "SAMD5",
  "TMEM132C",
  "SLIT2",
  "ELL2",
  "PRDM1",
  "SCN9A",
  "ADAMTS20",
  "ANKRD18A",
  "TMEM155",
  "SST",
  "HGF",
  "TRPC4",
  "EPHB6",
  "ST6GALNAC5",
  "KCNA1",
  "ADARB2",
  "NFIB",
  "RGS12",
  "NFIX",
  "PROX1",
  "EGFR",
  "KIT",
  "NR2F2",
  "AP1S2",
  "CNTNAP4",
  "DOCK10",
  "CHRNA2",
  "CNR1",
  "PKP2",
  "ARPP21",
  "RERG",
  "NECAB2",
  "SV2C",
  "GALNTL6",
  "VIP",
  "LAMP5",
  "SALL1",
  "CCK",
  "CDHR3",
  "FBXL7",
  "ASIC4",
  "LINGO2",
  "PDGFD",
  "IGFBP5",
  "SEZ6L",
  "SHISA8",
  "SYNPR",
  "NR2E1",
  "TP53I11",
  "ATP11C",
  "PLD5",
  "HPCAL1",
  "KCNJ2",
  "NKAIN3",
  "IGFBP3",
  "CHRNA7",
  "L3MBTL4",
  "TRIM36",
  "CHRDL1",
  "ADRA1B",
  "KIRREL",
  "CRIM1",
  "RGS8",
  "CAMKV",
  "TRPC3",
  "UNC5B",
  "NTNG1",
  "CA3",
  "TSHZ2",
  "ERBB4",
  "GOLIM4",
  "PLCXD3",
  "SLIT2",
  "ZMAT4",
  "DPP10",
  "GAD2",
  "LYPD6",
  "LYPD6B",
  "FNBP1L",
  "RGS5",
  "SLC26A4",
  "PVALB",
  "SLIT1",
  "RPH3AL",
  "TAC1",
  "NPPC",
  "GABRD",
  "SULF1",
  "RNF144B",
  "LRRC38",
  "LDB2",
  "C8orf4",
  "MYO1B",
  "SYNPR",
  "GRIK1",
  "FSTL5",
  "KCNIP1",
  "ELFN1",
  "CHRM3",
  "MAN1A1",
  "LHFPL3",
  "FAT4",
  "GULP1",
  "TRPC6",
  "GAD1",
  "PRNP",
  "CLU",
  "TMEM66",
  "ROBO2",
  "TSPAN7",
  "GRIK1",
  "ARPP21",
  "FIGN",
  "GPRC6A",
  "SKAP1",
  "MET",
  "RORB",
  "SGCZ",
  "NCAM2",
  "CBLN4",
  "PDE4B",
  "PENK",
  "STK32B",
  "DSP",
  "KCNIP4",
  "SLC24A4",
  "DOCK10",
  "SLC45A2",
  "KIRREL3",
  "CALB1",
  "CNTN6",
  "HPGD",
  "NABP1",
  "MKX",
  "AOX1",
  "CPVL",
  "MAML2",
  "ST6GALNAC5",
  "ZMAT4",
  "SAMD5",
  "PAWR",
  "ERBB4",
  "L3MBTL4",
  "GALNTL6",
  "ARPP21",
  "PDE4B",
  "SEZ6L",
  "CDH10",
  "SV2C",
  "FGF13",
  "LAMP5",
  "PTCHD4",
  "CACNA2D1",
  "SPHKAP",
  "GOLIM4",
  "MYO16",
  "FBXL7",
  "TRPC3",
  "EYA4",
  "SGCZ",
  "MPPED1",
  "GRIN2A",
  "UNC5D",
  "ST8SIA4",
  "KIT",
  "SERTM1",
  "HRH1",
  "PDGFD",
  "CBLN4",
  "GRM5",
  "MGAT4C",
  "ARHGAP31",
  "TMEM132D",
  "PMEPA1",
  "OLFM3",
  "IGFBP5",
  "COL11A1",
  "CHRNA2",
  "OSBPL3",
  "GPR149",
  "LHFP",
  "PTGDS",
  "RARB",
  "NPR3",
  "RYR3",
  "ZMAT4",
  "SLC7A11",
  "CXCL12",
  "COL9A1",
  "SLC22A3",
  "CBLN1",
  "EXPH5",
  "KCNK1",
  "CSMD1",
  "DUSP10",
  "SLC12A7",
  "ASB13",
  "EDNRA",
  "IL1RAP",
  "SYT2",
  "GULP1",
  "RNF144B",
  "AHRR",
  "IGSF21",
  "KIT",
  "SHISA8",
  "CYP27A1",
  "ARAP2",
  "PDE9A",
  "KCNK2",
  "PAPPA",
  "FOSL2",
  "ACTN2",
  "TSHZ2",
  "LYPD6",
  "VIP",
  "TNS1",
  "FAM126A",
  "GRM5",
  "MN1",
  "MAML3",
  "PDZRN3",
  "FSTL5",
  "ST6GALNAC5",
  "CHSY3",
  "CACNA1I",
  "MAN1A1",
  "NETO2",
  "BCL11A",
  "EPHA7",
  "EPHB6",
  "B3GALTL",
  "PTCHD2",
  "SEMA5B",
  "CACNG3",
  "PRKG2",
  "OCA2",
  "WIF1",
  "SOX5",
  "MPPED1",
  "RGMB",
  "FOXO1",
  "BAIAP3",
  "SLITRK4",
  "ST6GAL2",
  "GABRA5",
  "THBS1",
  "ZNF385B",
  "ZNF462",
  "PDZRN4",
  "NDNF",
  "MEIS3",
  "STMN1",
  "GLCE",
  "LMO4",
  "RELN",
  "GABRG1",
  "SPHKAP",
  "DDR2",
  "DPY19L1",
  "KCNK9",
  "CNN3",
  "FAM19A4",
  "FRAS1",
  "MGAT4C",
  "ALCAM",
  "CDH9",
  "COBLL1",
  "PLCXD3",
  "KCNQ5",
  "CDH8",
  "CDHR3",
  "PRR16",
  "PCDH10",
  "ANGPT1",
  "SLC7A11",
  "KCND3",
  "CNR1",
  "HEG1",
  "SCN7A",
  "NETO1",
  "SYT17",
  "SP8",
  "GREM2",
  "SOCS2",
  "TRPC4",
  "ABI3BP",
  "PLS3",
  "FAM126A",
  "CDH11",
  "CNTFR",
  "ADRA2A",
  "SLC6A16",
  "CDKN2D",
  "RGMB",
  "RP11-58C22.1",
  "EDNRA",
  "RXRG",
  "TMPO",
  "EPHA7",
  "ADRA1D",
  "ANKFN1",
  "CHL1",
  "PRKAA1",
  "NECAB2",
  "NDN",
  "KBTBD11",
  "H2AFZ",
  "JPH1",
  "NUPL2",
  "ITGB8",
  "RASSF2",
  "C5orf30",
  "GRPEL1",
  "LSS",
  "SDC2",
  "INSIG1",
  "CHRNA7",
  "SPARCL1",
  "PRKCG",
  "RFC3",
  "ARAP2",
  "IGSF21",
  "KCNIP4",
  "KCNJ2",
  "EXPH5",
  "PREX1",
  "PTGFR",
  "SLC24A4",
  "SLC12A7",
  "CDK6",
  "SCTR",
  "PCDH19",
  "SYT2",
  "SORL1",
  "ANXA11",
  "VWDE",
  "TOX3",
  "PDE9A",
  "RP11-664D7.4",
  "NXPH2",
  "PROS1",
  "KIRREL3",
  "CNTN6",
  "HTR2C",
  "ENC1",
  "C16orf95",
  "HPSE",
  "WFS1",
  "ESYT3",
  "ARHGEF26",
  "SYT10",
  "BTBD3",
  "CACNG3",
  "NR2E1",
  "SLC22A3",
  "CD24",
  "RGS16",
  "ISM1",
  "KIRREL",
  "NEXN",
  "IER5",
  "KCNN3",
  "PTN",
  "PLXNB1",
  "VSTM2B",
  "LRIG1",
  "EIF3I",
  "PRR16",
  "FRAS1",
  "FREM1",
  "RELN",
  "POU6F2",
  "BMP6",
  "CNR1",
  "FIGN",
  "PLCXD3",
  "CRH",
  "GDA",
  "FAT4",
  "BMP3",
  "NDNF",
  "IL1RAP",
  "SERPINE2",
  "SYT13",
  "PARM1",
  "ITGA8",
  "LHX6",
  "NOS1",
  "FREM2",
  "GPC6",
  "SYNDIG1",
  "CA8",
  "CPNE8",
  "HCRTR2",
  "CA1",
  "TMEM255A",
  "ST6GALNAC5",
  "VAV3",
  "RND3",
  "TAC1",
  "GPR126",
  "NRP2",
  "CA13",
  "ROR1",
  "OPRM1",
  "TXK",
  "CALB1",
  "PDGFD",
  "ZBTB7C",
  "FBN2",
  "RET",
  "GULP1",
  "ZNF804B",
  "LHX2",
  "CA3",
  "NTNG1",
  "JADE3",
  "SLC7A3",
  "ARHGAP6",
  "KCTD12",
  "MARCH3",
  "FAM179A",
  "NWD2",
  "CNTN6",
  "GFRA2"
)) |> 
  dplyr::filter(!duplicated(hugo_symbol)) |> 
  dplyr::mutate(status="inhibitory")

df <- rbind(df1,df2) |> 
  pivot_wider(names_from = status, values_from=status) |> 
  dplyr::filter(is.na(exhibitory) | is.na(inhibitory)) |> 
  dplyr::mutate(status = ifelse(is.na(exhibitory),"inh","exh")) |> 
  dplyr::mutate(exhibitory = NULL, inhibitory=NULL) |> 
  dplyr::mutate(C4 = ifelse(hugo_symbol %in% (results.out |> dplyr::filter(C4.2022) |> dplyr::pull(hugo_symbol)),"in C4","not DE over time"))


plot(as.factor(df$status),as.factor(df$C4))



exp <- gsam.rna.metadata |> dplyr::select( sid,
                                    rna.signature.C1.collagen.2022, 
                                    rna.signature.C4.neuron.2022)
saveRDS(exp, file="/tmp/gsam-signatures.Rds")

write.csv(exp, file="/tmp/gsam-signatures.csv")

