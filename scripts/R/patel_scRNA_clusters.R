#!/usr/bin/env R


# scRNA-seq gene clusters from: https://science.sciencemag.org/content/344/6190/1396


# load expression data for matching ----


if(!exists('gencode.31')) {
  source('scripts/R/gsam_rna-seq_expression.R')
}


# Cell Cycle ----


patel.cell.cycle.tt2 <- c(
  "TOP2A",
  "RRM2",
  "CDK1",
  "PBK",
  "BIRC5",
  "AURKB",
  "KIF15",
  "KIF4A",
  "UBE2T",
  "NUSAP1",
  "CENPF",
  "CENPU", # described as 'MLF1IP' in the manuscript
  "PCLAF", # described as 'KIAA0101' in the manuscript
  "CCNB2",
  "FANCI",
  "CDC6",
  "BUB1",
  "TPX2",
  "DTL",
  "FANCD2",
  "MAD2L1",
  "NCAPG2",
  "CLSPN",
  "CCNB1",
  "XRCC2",
  "KIF20A",
  "CCNE2",
  "DHFR",
  "RFC3",
  "DSN1",
  "HMGB2",
  "SPAG5",
  "POLQ",
  "NCAPD2",
  "RACGAP1",
  "ATAD2",
  "TIMELESS",
  "FEN1",
  "KNTC1",
  "CDC7"
)

stopifnot(patel.cell.cycle.tt2 %in% gencode.31$GENE) # patel.cell.cycle.tt2[patel.cell.cycle.tt2 %in% gencode.31$GENE == F]




# Hypoxia ----

# https://science.sciencemag.org/content/sci/344/6190/1396/F2.large.jpg


patel.hypoxia.tt2 <- c(
  "VEGFA",
  "ADM",
  "IGFBP5",
  "AKAP12",
  "HILPDA",
  "NDRG1",
  'ERO1A', # described as 'ERO1L' in the manuscript
  "NRN1",
  "ZNF395",
  "SLC2A3",
  "PGK1",
  "ALDOA",
  "MT1X",
  "TREM1",
  "UBC",
  "DHRS3",
  "TCAF2", # described as 'FAM115C' in the manuscript
  "EPAS1",
  "ENO1",
  "RPS27",
  "LDHA",
  "TPI1",
  "MT2A",
  "LGALS3",
  "SLC2A1",
  "PFKP",
  "SLC39A14",
  "TMSB10",
  "GBE1",
  "RPS18",
  "SMIM3", # described as 'C5orf62' in the manuscript
  "RPL10",
  "IGFBP3",
  "SDC4",
  # "FTL", appears in both hypoxia & complete/imm resp
  "BHLHE40",
  "EEF1A1",
  "ACTG1",
  "ARRDC3",
  "MIF"
)

stopifnot(patel.hypoxia.tt2 %in% gencode.31$GENE) # patel.hypoxia.tt2[patel.hypoxia.tt2 %in% gencode.31$GENE == F]


# complement/immune response ----

# http://www.sciencemag.org/content/344/6190/1396/suppl/DC1
# Fig S10

patel.complement.immune.response.tt2 <- c(
  "CHI3L2",
  "C3",
  "SERPINA3",
  "SOD2",
  "MGST1",
  "A2M",
  "CHI3L1",
  "MAOB",
  "GBP1",
  "SAA1",
  "GBP2",
  "SKAP2",
  "C1S",
  "ID3",
  "CADPS",
  "ACTN1",
  "GLRX",
  "TMBIM1",
  "TAGLN",
  "CD44",
  "GLIPR1",
  "SERPING1",
  "NNMT",
  "MYBPC1",
  "HSPB8",
  "MYLK",
  "PIRT",
  "GASK1B", # described as 'FAM198B' in the manuscript
  "NUPR1",
  "PGAM2",
  "HLA-A",
  "C1R",
  "S1PR1",
  "MYL6",
  "GFAP",
  "MIR4709",
  "FN1",
  "FAM107A",
  "HLA-DRA"
  # ,"FTL" appears in both hypoxia & complete/imm resp
)


stopifnot(patel.complement.immune.response.tt2 %in% gencode.31$GENE) # patel.complement.immune.response.tt2[patel.complement.immune.response.tt2 %in% gencode.31$GENE == F]


# check mutex ----


# require mutual exclusivity
stopifnot(patel.cell.cycle.tt2 %in% c(patel.complement.immune.response.tt2, patel.hypoxia.tt2) == F )
stopifnot(patel.complement.immune.response.tt2 %in% c(patel.cell.cycle.tt2, patel.hypoxia.tt2) == F )
stopifnot(patel.hypoxia.tt2 %in% c(patel.complement.immune.response.tt2, patel.cell.cycle.tt2) == F )

# "FTL" is in both patel.complement.immune.response.tt2 and patel.hypoxia.tt2



# obtain gids etc. from gencode ----


patel.complement.immune.response <-data.frame(symbol=patel.complement.immune.response.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)


patel.hypoxia <- data.frame(symbol=patel.hypoxia.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)


patel.cell.cycle <- data.frame(symbol=patel.cell.cycle.tt2) %>%
  dplyr::left_join(gencode.31, by=c('symbol' = 'GENE')) %>%
  dplyr::mutate(id = paste0(ENSG , "|" , symbol , "|" , V1 , ":" , V4, "-" , V5 ,"(", V7 , ")")) # "ENSG00000141736.13_3|ERBB2|chr17:37844167-37886679(+)







