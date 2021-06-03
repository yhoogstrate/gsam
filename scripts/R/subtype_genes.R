#!/usr/bin/env R 

# load expression data for matching ----

if(!exists('gencode.31')) {
  source('scripts/R/gsam_rna-seq_expression.R')
}


# TCGA wang/verhaak 2017 subtypes ----


subtype.classical.tt2 <- {{}}
subtype.classical.tt2['1'] <- 'PTPRA'
subtype.classical.tt2['2'] <- 'ELOVL2'
subtype.classical.tt2['3'] <- 'MLC1'
subtype.classical.tt2['4'] <- 'SOX9'
subtype.classical.tt2['5'] <- 'ARNTL'
subtype.classical.tt2['6'] <- 'DENND2A'
subtype.classical.tt2['7'] <- 'BBS1'
subtype.classical.tt2['8'] <- 'ABLIM1'
subtype.classical.tt2['9'] <- 'PAX6'
subtype.classical.tt2['10'] <- 'ZHX3'
subtype.classical.tt2['11'] <- 'USP8'
subtype.classical.tt2['12'] <- 'PLCG1'
subtype.classical.tt2['13'] <- 'CDH4'
subtype.classical.tt2['14'] <- 'RASGRP1'
subtype.classical.tt2['15'] <- 'ACSBG1'
subtype.classical.tt2['16'] <- 'CST3'
subtype.classical.tt2['17'] <- 'BCKDHB'
subtype.classical.tt2['18'] <- 'LHFPL6' # 'LHFP' << in the origial paper
subtype.classical.tt2['19'] <- 'VAV3'
subtype.classical.tt2['20'] <- 'ACSL3'
subtype.classical.tt2['21'] <- 'EYA2'
subtype.classical.tt2['22'] <- 'SEPT11'
subtype.classical.tt2['23'] <- 'SLC4A4'
subtype.classical.tt2['24'] <- 'SLC20A2'
subtype.classical.tt2['25'] <- 'DGLUCY' # 'C14orf159' << in the origial paper
subtype.classical.tt2['26'] <- 'CTNND1'
subtype.classical.tt2['27'] <- 'ZFHX4'
subtype.classical.tt2['28'] <- 'SPRY2'
subtype.classical.tt2['29'] <- 'ZNF45'
subtype.classical.tt2['30'] <- 'NCOA1'
subtype.classical.tt2['31'] <- 'PLCE1'
subtype.classical.tt2['32'] <- 'DTNA'
subtype.classical.tt2['33'] <- 'POLRMT'
subtype.classical.tt2['34'] <- 'SALL1'
subtype.classical.tt2['35'] <- 'TYK2'
subtype.classical.tt2['36'] <- 'TJP1'
subtype.classical.tt2['37'] <- 'MEOX2'
subtype.classical.tt2['38'] <- 'FGFR3'
subtype.classical.tt2['39'] <- 'STXBP3'
subtype.classical.tt2['40'] <- 'GRIK1'
subtype.classical.tt2['41'] <- 'GATM'
subtype.classical.tt2['42'] <- 'UPF1'
subtype.classical.tt2['43'] <- 'NPEPL1'
subtype.classical.tt2['44'] <- 'EFCAB14' # 'KIAA0494' << in the origial paper
subtype.classical.tt2['45'] <- 'RBCK1'
subtype.classical.tt2['46'] <- 'PHKB'
subtype.classical.tt2['47'] <- 'SLC3A2'
subtype.classical.tt2['48'] <- 'PPARGC1A'
subtype.classical.tt2['49'] <- 'PNPLA6'
subtype.classical.tt2['50'] <- 'MYO5C'


subtype.proneural.tt2 <- {{}}
subtype.proneural.tt2['1'] <- 'JPT1' # 'HN1' << in the origial paper
subtype.proneural.tt2['2'] <- 'RAB33A'
subtype.proneural.tt2['3'] <- 'HDAC2'
subtype.proneural.tt2['4'] <- 'MYT1'
subtype.proneural.tt2['5'] <- 'MTSS1'
subtype.proneural.tt2['6'] <- 'HOXD3'
subtype.proneural.tt2['7'] <- 'GPR17'
subtype.proneural.tt2['8'] <- 'PTTG1'
subtype.proneural.tt2['9'] <- 'KLRC3'
subtype.proneural.tt2['10'] <- 'PLAAT1' # 'HRASLS' << in the origial paper
subtype.proneural.tt2['11'] <- 'TCP1'
subtype.proneural.tt2['12'] <- 'NPPA'
subtype.proneural.tt2['13'] <- 'PFDN2'
subtype.proneural.tt2['14'] <- 'CA10'
subtype.proneural.tt2['15'] <- 'EPHB1'
subtype.proneural.tt2['16'] <- 'UGT8'
subtype.proneural.tt2['17'] <- 'PAK5' # 'PAK7' << in the origial paper
subtype.proneural.tt2['18'] <- 'SLC1A1'
subtype.proneural.tt2['19'] <- 'NARF'
subtype.proneural.tt2['20'] <- 'DCTN3'
subtype.proneural.tt2['21'] <- 'SMPD3'
subtype.proneural.tt2['22'] <- 'ZNF804A'
subtype.proneural.tt2['23'] <- 'RASL11B'
subtype.proneural.tt2['24'] <- 'MYB'
subtype.proneural.tt2['25'] <- 'PDGFRA'
subtype.proneural.tt2['26'] <- 'ERBB3'
subtype.proneural.tt2['27'] <- 'CLGN'
subtype.proneural.tt2['28'] <- 'SOX10'
subtype.proneural.tt2['29'] <- 'BCL11A'
subtype.proneural.tt2['30'] <- 'NMU'
subtype.proneural.tt2['31'] <- 'ZFP69B' # 'ZNF643' << in the origial paper
subtype.proneural.tt2['32'] <- 'CDKN1C'
subtype.proneural.tt2['33'] <- 'JPH3'
subtype.proneural.tt2['34'] <- 'PCDHA9'
subtype.proneural.tt2['35'] <- 'IL1RAPL1'
subtype.proneural.tt2['36'] <- 'MAST1'
subtype.proneural.tt2['37'] <- 'VIPR2'
subtype.proneural.tt2['38'] <- 'SIM2'
subtype.proneural.tt2['39'] <- 'BAMBI'
subtype.proneural.tt2['40'] <- 'PKMYT1'
subtype.proneural.tt2['41'] <- 'PLCB4'
subtype.proneural.tt2['42'] <- 'SLC17A6'
subtype.proneural.tt2['43'] <- 'KLRK1'
subtype.proneural.tt2['44'] <- 'CENPJ'
subtype.proneural.tt2['45'] <- 'NHLH1'
subtype.proneural.tt2['46'] <- 'GABRB3'
subtype.proneural.tt2['47'] <- 'KLRC4'
subtype.proneural.tt2['48'] <- 'KCNK3'
subtype.proneural.tt2['49'] <- 'GRID2'
subtype.proneural.tt2['50'] <- 'DACH1'


subtype.mesenchymal.tt2 <- {{}}
subtype.mesenchymal.tt2['1'] <- 'ARPC1B'
subtype.mesenchymal.tt2['2'] <- 'S100A11'
subtype.mesenchymal.tt2['3'] <- 'CTSC'
subtype.mesenchymal.tt2['4'] <- 'GLIPR1'
subtype.mesenchymal.tt2['5'] <- 'NNMT'
subtype.mesenchymal.tt2['6'] <- 'VDR'
subtype.mesenchymal.tt2['7'] <- 'RGS2'
subtype.mesenchymal.tt2['8'] <- 'CTSB'
subtype.mesenchymal.tt2['9'] <- 'TGFBI'
subtype.mesenchymal.tt2['10'] <- 'PLAUR'
subtype.mesenchymal.tt2['11'] <- 'LY96'
subtype.mesenchymal.tt2['12'] <- 'BCL3'
subtype.mesenchymal.tt2['13'] <- 'TNFAIP8'
subtype.mesenchymal.tt2['14'] <- 'IER3'
subtype.mesenchymal.tt2['15'] <- 'PRSS23'
subtype.mesenchymal.tt2['16'] <- 'IL7R'
subtype.mesenchymal.tt2['17'] <- 'RAB27A'
subtype.mesenchymal.tt2['18'] <- 'RUNX1'
subtype.mesenchymal.tt2['19'] <- 'P4HA2'
subtype.mesenchymal.tt2['20'] <- 'CYP1B1'
subtype.mesenchymal.tt2['21'] <- 'BACE2'
subtype.mesenchymal.tt2['22'] <- 'ACPP'
subtype.mesenchymal.tt2['23'] <- 'FTL'
subtype.mesenchymal.tt2['24'] <- 'SLPI'
subtype.mesenchymal.tt2['25'] <- 'RAC2'
subtype.mesenchymal.tt2['26'] <- 'RARRES1'
subtype.mesenchymal.tt2['27'] <- 'SYNGR2'
subtype.mesenchymal.tt2['28'] <- 'THBS1'
subtype.mesenchymal.tt2['29'] <- 'IL6'
subtype.mesenchymal.tt2['30'] <- 'CAV1'
subtype.mesenchymal.tt2['31'] <- 'PI3'
subtype.mesenchymal.tt2['32'] <- 'CDCP1'
subtype.mesenchymal.tt2['33'] <- 'ITGB1'
subtype.mesenchymal.tt2['34'] <- 'LOX'
subtype.mesenchymal.tt2['35'] <- 'CD72'
subtype.mesenchymal.tt2['36'] <- 'COL1A2'
subtype.mesenchymal.tt2['37'] <- 'ANPEP'
subtype.mesenchymal.tt2['38'] <- 'MMP7'
subtype.mesenchymal.tt2['39'] <- 'SPAG4'
subtype.mesenchymal.tt2['40'] <- 'BNC2'
subtype.mesenchymal.tt2['41'] <- 'NDRG1'
subtype.mesenchymal.tt2['42'] <- 'CNN2'
subtype.mesenchymal.tt2['43'] <- 'LUM'
subtype.mesenchymal.tt2['44'] <- 'PTGS2'
subtype.mesenchymal.tt2['45'] <- 'COL3A1'
subtype.mesenchymal.tt2['46'] <- 'COL5A1'
subtype.mesenchymal.tt2['47'] <- 'SDC1'
subtype.mesenchymal.tt2['48'] <- 'COL1A1'
subtype.mesenchymal.tt2['49'] <- 'GPRC5A'
subtype.mesenchymal.tt2['50'] <- 'COL15A1'


## ensure match in gencode 31 ----


stopifnot(length(subtype.classical.tt2[subtype.classical.tt2 %in% gencode.31$GENE == F]) == 0)
stopifnot(length(subtype.proneural.tt2[subtype.proneural.tt2 %in% gencode.31$GENE == F]) == 0)
stopifnot(length(subtype.mesenchymal.tt2[subtype.mesenchymal.tt2 %in% gencode.31$GENE == F]) == 0)



subtype.mesenchymal <-data.frame(symbol=subtype.mesenchymal.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)


subtype.proneural <-data.frame(symbol=subtype.proneural.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)


subtype.classical <-data.frame(symbol=subtype.classical.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)



stopifnot(subtype.mesenchymal$id %in% rownames(gsam.rnaseq.expression))
stopifnot(subtype.proneural$id %in% rownames(gsam.rnaseq.expression))
stopifnot(subtype.classical$id %in% rownames(gsam.rnaseq.expression))


# GLiTS Redux set ----


# https://www.nature.com/articles/s41374-020-0437-0/figures/1?proof=t


classical.glits.redux.tt2 <- {{}}
classical.glits.redux.tt2['1'] <- 'EGFR'
classical.glits.redux.tt2['2'] <- 'NOS2'
classical.glits.redux.tt2['3'] <- 'GAS1'
classical.glits.redux.tt2['4'] <- 'KCNF1'
classical.glits.redux.tt2['5'] <- 'ELOVL2'
classical.glits.redux.tt2['6'] <- 'VAV3'
classical.glits.redux.tt2['7'] <- 'SOCS2'
classical.glits.redux.tt2['8'] <- 'MEOX2'
classical.glits.redux.tt2['9'] <- 'RGS6'
classical.glits.redux.tt2['10'] <- 'HS3ST3B1'
classical.glits.redux.tt2['11'] <- 'PDGFA'
classical.glits.redux.tt2['12'] <- 'GRIK1'
classical.glits.redux.tt2['13'] <- 'EYA2'
classical.glits.redux.tt2['14'] <- 'CAMK2B'
classical.glits.redux.tt2['15'] <- 'SLC4A4'
classical.glits.redux.tt2['16'] <- 'CDH4'
classical.glits.redux.tt2['17'] <- 'WSCD1'
classical.glits.redux.tt2['18'] <- 'ACSBG1'
classical.glits.redux.tt2['19'] <- 'MLC1'
classical.glits.redux.tt2['20'] <- 'NES'
classical.glits.redux.tt2['21'] <- 'KLHL4'
classical.glits.redux.tt2['22'] <- 'ZFHX4'
classical.glits.redux.tt2['23'] <- 'KLHDC8A'
classical.glits.redux.tt2['24'] <- 'SPRY2'
classical.glits.redux.tt2['25'] <- 'DENND2A'
classical.glits.redux.tt2['26'] <- 'MEIS1'
classical.glits.redux.tt2['27'] <- 'SOX9'
classical.glits.redux.tt2['28'] <- 'ARAP3'

mesenchymal.glits.redux.tt2 <- {{}}
mesenchymal.glits.redux.tt2['1'] <- 'KLRC4'
mesenchymal.glits.redux.tt2['2'] <- 'GRID2'
mesenchymal.glits.redux.tt2['3'] <- 'UGT8'
mesenchymal.glits.redux.tt2['4'] <- 'EPHB1'
mesenchymal.glits.redux.tt2['5'] <- 'MYT1'
mesenchymal.glits.redux.tt2['6'] <- 'DCX'
mesenchymal.glits.redux.tt2['7'] <- 'KLRC3'
mesenchymal.glits.redux.tt2['8'] <- 'TMSB15A'
mesenchymal.glits.redux.tt2['9'] <- 'CNTN1'
mesenchymal.glits.redux.tt2['10'] <- 'ERBB3'
mesenchymal.glits.redux.tt2['11'] <- 'PLAAT1' # formerly known as 'HRASLS' [https://www.genecards.org/cgi-bin/carddisp.pl?gene=PLAAT1]
mesenchymal.glits.redux.tt2['12'] <- 'DNM3'
mesenchymal.glits.redux.tt2['13'] <- 'SCN3A'
mesenchymal.glits.redux.tt2['14'] <- 'PAK3'
mesenchymal.glits.redux.tt2['15'] <- 'GNG4'
mesenchymal.glits.redux.tt2['16'] <- 'AMOTL2'
mesenchymal.glits.redux.tt2['17'] <- 'CHD7'
mesenchymal.glits.redux.tt2['18'] <- 'CRMP1'
mesenchymal.glits.redux.tt2['19'] <- 'SOX4'
mesenchymal.glits.redux.tt2['20'] <- 'SATB1'

proneural.glits.redux.tt2 <- {{}}
proneural.glits.redux.tt2['1`'] <- 'LOX'
proneural.glits.redux.tt2['2'] <- 'SERPINE1'
proneural.glits.redux.tt2['3'] <- 'COL1A1'
proneural.glits.redux.tt2['4'] <- 'CRYBG1' # formerly known as 'AIM1' [https://www.genecards.org/cgi-bin/carddisp.pl?gene=CRYBG1]
proneural.glits.redux.tt2['5'] <- 'FCGR2B'
proneural.glits.redux.tt2['6'] <- 'COL1A2'
proneural.glits.redux.tt2['7'] <- 'FHL2'
proneural.glits.redux.tt2['8'] <- 'IGFBP6'
proneural.glits.redux.tt2['9'] <- 'TGFBI'
proneural.glits.redux.tt2['10'] <- 'PLAU'
proneural.glits.redux.tt2['11'] <- 'P4HA2'
proneural.glits.redux.tt2['12'] <- 'CSTA'
proneural.glits.redux.tt2['13'] <- 'DSE'
proneural.glits.redux.tt2['14'] <- 'LAMB1'
proneural.glits.redux.tt2['15'] <- 'SRPX2'
proneural.glits.redux.tt2['16'] <- 'TIMP1'
proneural.glits.redux.tt2['17'] <- 'S100A11'
proneural.glits.redux.tt2['18'] <- 'PLAUR'
proneural.glits.redux.tt2['19'] <- 'FPR3'
proneural.glits.redux.tt2['20'] <- 'CLEC2B'
proneural.glits.redux.tt2['21'] <- 'DCBLD2'
proneural.glits.redux.tt2['22'] <- 'MYOF'
proneural.glits.redux.tt2['23'] <- 'MAFB'
proneural.glits.redux.tt2['24'] <- 'ARPC1B'
proneural.glits.redux.tt2['25'] <- 'RAB27A'
proneural.glits.redux.tt2['26'] <- 'LY75'
proneural.glits.redux.tt2['27'] <- 'PLBD1'
proneural.glits.redux.tt2['28'] <- 'LY96'


## ensure match in gencode 31 ----




stopifnot(classical.glits.redux.tt2 %in% gencode.31$GENE) # classical.glits.redux.tt2[classical.glits.redux.tt2 %in% gencode.31$GENE == F]
stopifnot(mesenchymal.glits.redux.tt2 %in% gencode.31$GENE) # mesenchymal.glits.redux.tt2[mesenchymal.glits.redux.tt2 %in% gencode.31$GENE == F]
stopifnot(proneural.glits.redux.tt2 %in% gencode.31$GENE) # proneural.glits.redux.tt2[proneural.glits.redux.tt2 %in% gencode.31$GENE == F]




mesenchymal.glits.redux <-data.frame(symbol=mesenchymal.glits.redux.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)


proneural.glits.redux <-data.frame(symbol=proneural.glits.redux.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)


classical.glits.redux <-data.frame(symbol=classical.glits.redux.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)



stopifnot(mesenchymal.glits.redux$id %in% rownames(gsam.rnaseq.expression))
stopifnot(proneural.glits.redux$id %in% rownames(gsam.rnaseq.expression))
stopifnot(classical.glits.redux$id %in% rownames(gsam.rnaseq.expression))



