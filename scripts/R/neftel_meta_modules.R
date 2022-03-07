#!/usr/bin/env R

# load expression data for matching ----


if(!exists('gencode.31')) {
  source('scripts/R/gsam_rna-seq_expression.R')
}



# load meta modules ----

# https://www.sciencedirect.com/science/article/pii/S0092867419306877?via%3Dihub#mmc2
neftel.meta.modules <- read_xlsx('data/1-s2.0-S0092867419306877-mmc2.xlsx', skip=4) %>%
  as.data.frame()


neftel.meta.modules.MES1.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`MES1`)) %>%
  dplyr::pull(`MES1`) %>%
  gsub("^C8orf4$","TCIM",.) %>%
  purrr::discard(. %in% c("S100A16","ENO2","TUBA1C")) # non-mutex
stopifnot(neftel.meta.modules.MES1.tt2 %in% gencode.31$GENE) # neftel.meta.modules.MES1.tt2[neftel.meta.modules.MES1.tt2 %in% gencode.31$GENE == F]

neftel.meta.modules.MES2.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`MES2`)) %>%
  dplyr::pull(`MES2`) %>%
  gsub("^ERO1L$","ERO1A",.) %>%
  purrr::discard(. %in% c("S100A16","ENO2","TUBA1C")) # non-mutex
stopifnot(neftel.meta.modules.MES2.tt2 %in% gencode.31$GENE) # neftel.meta.modules.MES2.tt2[neftel.meta.modules.MES2.tt2 %in% gencode.31$GENE == F]
 
neftel.meta.modules.AC.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`AC`)) %>%
  dplyr::pull(`AC`) %>%
  gsub("^PPAP2B$","PLPP3",.)
stopifnot(neftel.meta.modules.AC.tt2 %in% gencode.31$GENE) # neftel.meta.modules.AC.tt2[neftel.meta.modules.AC.tt2 %in% gencode.31$GENE == F]

neftel.meta.modules.OPC.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`OPC`)) %>%
  dplyr::pull(`OPC`) %>%
  gsub("^LPPR1$","PLPPR1",.) %>%
  gsub("^HRASLS$","PLAAT1",.)
stopifnot(neftel.meta.modules.OPC.tt2 %in% gencode.31$GENE) # neftel.meta.modules.OPC.tt2[neftel.meta.modules.OPC.tt2 %in% gencode.31$GENE == F]

neftel.meta.modules.NPC1.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`NPC1`)) %>%
  dplyr::pull(`NPC1`) %>%
  purrr::discard(. %in% c("CD24")) %>% # no entry of CD24 in gencode 31
  gsub("^HN1$","JPT1",.) %>%
  gsub("^GPR56$","ADGRG1",.)
stopifnot(neftel.meta.modules.NPC1.tt2 %in% gencode.31$GENE) # neftel.meta.modules.NPC1.tt2[neftel.meta.modules.NPC1.tt2 %in% gencode.31$GENE == F]

neftel.meta.modules.NPC2.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`NPC2`)) %>%
  dplyr::pull(`NPC2`) %>%
  purrr::discard(. %in% c("CD24")) %>% # no entry of CD24 in gencode 31
  gsub("^HMP19$","NSG2",.) %>%
  gsub("^HN1$","JPT1",.) %>%
  gsub("^LOC150568$","LINC01102",.) 
stopifnot(neftel.meta.modules.NPC2.tt2 %in% gencode.31$GENE) # neftel.meta.modules.NPC2.tt2[neftel.meta.modules.NPC2.tt2 %in% gencode.31$GENE == F]

neftel.meta.modules.G1.S.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`G1/S`)) %>%
  dplyr::pull(`G1/S`) %>%
  gsub("^KIAA0101$","PCLAF",.) %>%
  gsub("^MLF1IP$","CENPU",.) 
stopifnot(neftel.meta.modules.G1.S.tt2 %in% gencode.31$GENE) # neftel.meta.modules.G1.S.tt2[neftel.meta.modules.G1.S.tt2 %in% gencode.31$GENE == F]

neftel.meta.modules.G2.M.tt2 <- neftel.meta.modules %>%
  dplyr::filter(!is.na(`G2/M`)) %>%
  dplyr::pull(`G2/M`)
stopifnot(neftel.meta.modules.G2.M.tt2 %in% gencode.31$GENE) # neftel.meta.modules.G2.M.tt2[neftel.meta.modules.G2.M.tt2 %in% gencode.31$GENE == F]



rm(neftel.meta.modules)



# modules as subentry of gencode ----



neftel.meta.modules.MES1 <-data.frame(symbol=neftel.meta.modules.MES1.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.MES2 <-data.frame(symbol=neftel.meta.modules.MES2.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.AC <-data.frame(symbol=neftel.meta.modules.AC.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.OPC <-data.frame(symbol=neftel.meta.modules.OPC.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.NPC1 <-data.frame(symbol=neftel.meta.modules.NPC1.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.NPC2 <-data.frame(symbol=neftel.meta.modules.NPC2.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.G1.S <-data.frame(symbol=neftel.meta.modules.G1.S.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)

neftel.meta.modules.G2.M <-data.frame(symbol=neftel.meta.modules.G2.M.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)








