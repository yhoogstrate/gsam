#!/usr/bin/env R

# settings ----


options(warnPartialMatchDollar = TRUE) # https://stackoverflow.com/questions/32854683/data-frames-in-r-name-autocompletion


# load libs ----


library(tidyverse)
library(DESeq2)


# load data ----


source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set + metedata



# prepare data ----

## GSAM ----




gsam.metadata.all <- gsam.rna.metadata %>%
  dplyr::filter(blacklist.pca == F) %>%
  dplyr::filter(pat.with.IDH == F) %>%
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
  dplyr::filter(tumour.percentage.dna >= 15) %>%
  dplyr::mutate(tpc = 1 - (tumour.percentage.dna / 100))


gsam.metadata.all.paired <- gsam.metadata.all %>%
  dplyr::filter(pid %in% 
                  (gsam.metadata.all %>%
                     dplyr::group_by(pid) %>%
                     dplyr::tally() %>%
                     dplyr::filter(n == 2) %>% 
                     dplyr::ungroup() %>%
                     dplyr::filter(!duplicated(pid)) %>%
                     dplyr::pull(pid))
  ) %>%
  dplyr::mutate(pid = as.factor(as.character(pid))) # re-factor?





gsam.gene.expression.all <- gsam.rnaseq.expression %>%
  dplyr::select(gsam.metadata.all$sid)
stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)





## GLASS ----



glass.metadata.all  <- glass.gbm.rnaseq.metadata


glass.gene.expression.all <- glass.gbm.rnaseq.expression %>%
  dplyr::select(glass.metadata.all$sid)


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)




## Combined ----


combined.gene.expression <- dplyr::inner_join(
  gsam.gene.expression.all %>%
    `colnames<-`(paste0("GSAM.", colnames(.))) %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::mutate(gid = gsub('^(ENSG[0-9]+).+$','\\1', gid)) ,
  glass.gene.expression.all %>%
   tibble::rownames_to_column('gid'),
  by=c('gid'='gid') ) %>%
  tibble::column_to_rownames('gid') %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  limma::removeBatchEffect(as.factor(gsub("^(..).*$","\\1",colnames(.)) == "GS")) %>% # remove batch effect :)
  as.data.frame()


# neuron correlated genes ----



neuron.genes <- c("PTPRN","ICAM5","VWA7","PANX2","STRC","NPAS4","EGR4","OTOF","L1CAM","DOC2B","STXBP6","REPS2","DYNC1I1","AMPH",
                  "PWWP3B","KCNQ5","ZNF98","KCNK12","BSPRY","YPEL4","PDZD7","ZNF365","SLC6A15","RBM11","CDKL2","KCTD4","RASGRF2",
                  "IGFL2","TRIM54","TNNT2","ESYT3","RBP4","FSTL5","EPHB6","CAMKK1","SNCG","MATK","HS3ST2","FABP3","CBLN2","ATP8A2",
                  "SULT4A1","ATP2B3","TMEM130","SSTR1","CNNM1","CABP1","SHISAL1","FGF13","CBLN1","NEFH","KCNS1","RYR2","NEFM",
                  "ANKRD30B","KHDRBS2","SYCE1","SLC35F3","CDH18","CDH12","TRHDE","TCERG1L","FAM153CP","FAM153B","FAM153A","SOHLH1",
                  "MTUS2","DLGAP2","FRMPD4","CCK","HOOK1","MAP7D2","GALNTL6","CDH9","NEGR1","LRFN5","GRM1","ACVR1C","PPP4R4","DRD1",
                  "OCA2","CNGB1","RTN4RL1","GLP2R","MEPE","VSNL1","GABRB2","KRT222","HTR2A","SLITRK4","RXFP1","VSTM2L","SLC6A7",
                  "ISLR2","CYP4X1","COL26A1","CABLES1","SYNGR3","ARHGDIG","PSD","NPM2","SGSM1","OLFM1","SRRM4","NEUROD2","UBE2QL1",
                  "SVOP","CELF4","PCSK2","ATP1A3","SCN3B","CHGB","JPH3","DLGAP3","CKMT1B","KCNH3","DOC2A","KCNN1","UNC5A","SHANK1",
                  "TMEM151B","VWA5B2","STX1B","FBXL16","ADAM11","TMEM63C","PDE2A","TMEM246","PRKCZ","SNAP91","SH3GL2","KCTD16",
                  "NMNAT2","STXBP5L","TAC3","RIMBP2","CRHR2","CUX2","CLEC2L","SLC12A5","LRTM2","IQSEC3","WSCD2","HTR5A","GPR83",
                  "NRGN","KCNS2","SLC30A3","DNAJC5G","SYT7","GABRD","TGFBR3L","KCNC2","GABRG2","GABRA1","HCN1","CREG2","CACNG3",
                  "C1QL3","SLC6A17","SLC17A7","GABRG3","GABRA5","KCNV1","GPR26","GAD2","GPR22","CAMK1G","PDYN","SERTM1","HPCA",
                  "CPLX2","SYT5","SLC8A2","CHRM1","RAB3A","GRIN1","CHGA","SNCB","FAM163B","SYN1","FSTL4","GPR61","SLC32A1","CPLX1",
                  "SPTB","GLS2","SYT13","SNAP25","SYT1","SV2C","SYNPR","OLFM3","RBFOX1","UNC13C","RBFOX3","SYT2","PHYHIP","KCNA4",
                  "HPCAL4","MPPED1","EMX1","PACSIN1","CALY","CACNA1B","TMEM132D","SV2B","CCKBR","RASAL1","DDN","C4orf50","HRH3",
                  "CACNA1I","PHF24","MFSD4A","CAMK2A","SST","PRKCG","TBR1","SLC4A10","ARHGAP44","AJAP1","KCNT1","CHD5","GALNT9",
                  "GRM4","NGEF","RTN4R","SCN2B","NAP1L2","DMTN","BRINP1","GRIN2A","GDA","SNCA","PRKCB","SERPINI1","NECAB1","KCNK1",
                  "AK5","GABRA2","PPP2R2C","CPNE6","SOWAHA","ADARB2","NPY","KIRREL3","PNMA8B","FAIM2","DLG2","KIAA0319","SLC39A12",
                  "ETNPPL","HPSE2")

#results.out %>% dplyr::filter(in.gsam == T & in.glass == T) %>% 
#         dplyr::filter(hugo_symbol %in% neuron.genes) %>% dplyr::pull('ensembl_id')

neuron.genes.ens <- c("ENSG00000187730","ENSG00000067606","ENSG00000196581","ENSG00000116254","ENSG00000121769","ENSG00000121905","ENSG00000116544","ENSG00000116983",
                      "ENSG00000186377","ENSG00000134709","ENSG00000172260","ENSG00000154027","ENSG00000118733","ENSG00000156097","ENSG00000197106","ENSG00000157064",
                      "ENSG00000118194","ENSG00000143858","ENSG00000174514","ENSG00000008118","ENSG00000135750","ENSG00000183780","ENSG00000198626","ENSG00000163032",
                      "ENSG00000115155","ENSG00000115194","ENSG00000163793","ENSG00000138100","ENSG00000184261","ENSG00000135638","ENSG00000135625","ENSG00000175874",
                      "ENSG00000123612","ENSG00000136535","ENSG00000144290","ENSG00000054356","ENSG00000066248","ENSG00000187094","ENSG00000163630","ENSG00000145087",
                      "ENSG00000158220","ENSG00000163536","ENSG00000145198","ENSG00000157005","ENSG00000168993","ENSG00000181215","ENSG00000074211","ENSG00000151834",
                      "ENSG00000138769","ENSG00000152595","ENSG00000145335","ENSG00000164089","ENSG00000171509","ENSG00000168843","ENSG00000174473","ENSG00000215218",
                      "ENSG00000145526","ENSG00000154162","ENSG00000113100","ENSG00000164588","ENSG00000122012","ENSG00000113319","ENSG00000198944","ENSG00000053108",
                      "ENSG00000183775","ENSG00000011083","ENSG00000070808","ENSG00000145864","ENSG00000022355","ENSG00000113327","ENSG00000184845","ENSG00000145920",
                      "ENSG00000182230","ENSG00000074317","ENSG00000113763","ENSG00000170074","ENSG00000204677","ENSG00000137261","ENSG00000204396","ENSG00000124493",
                      "ENSG00000124507","ENSG00000178233","ENSG00000112232","ENSG00000185760","ENSG00000065609","ENSG00000152822","ENSG00000122585","ENSG00000106113",
                      "ENSG00000078053","ENSG00000158560","ENSG00000166448","ENSG00000160963","ENSG00000172209","ENSG00000236279","ENSG00000106123","ENSG00000157219",
                      "ENSG00000198010","ENSG00000158806","ENSG00000158856","ENSG00000168490","ENSG00000104722","ENSG00000123119","ENSG00000156486","ENSG00000164794",
                      "ENSG00000107295","ENSG00000122733","ENSG00000119125","ENSG00000165152","ENSG00000119411","ENSG00000078725","ENSG00000196990","ENSG00000130558",
                      "ENSG00000165643","ENSG00000107147","ENSG00000176884","ENSG00000148408","ENSG00000185736","ENSG00000165985","ENSG00000148482","ENSG00000136750",
                      "ENSG00000138311","ENSG00000173267","ENSG00000138207","ENSG00000172987","ENSG00000119946","ENSG00000186862","ENSG00000059915","ENSG00000154478",
                      "ENSG00000176769","ENSG00000130643","ENSG00000171772","ENSG00000110148","ENSG00000182255","ENSG00000019505","ENSG00000166793","ENSG00000011347",
                      "ENSG00000168539","ENSG00000174576","ENSG00000186642","ENSG00000150672","ENSG00000123901","ENSG00000149575","ENSG00000166257","ENSG00000154146",
                      "ENSG00000149571","ENSG00000120645","ENSG00000166159","ENSG00000181418","ENSG00000135519","ENSG00000135472","ENSG00000135423","ENSG00000166863",
                      "ENSG00000072657","ENSG00000166006","ENSG00000067715","ENSG00000072041","ENSG00000075035","ENSG00000166111","ENSG00000111249","ENSG00000111344",
                      "ENSG00000139767","ENSG00000157782","ENSG00000151952","ENSG00000060709","ENSG00000182870","ENSG00000132932","ENSG00000132938","ENSG00000180440",
                      "ENSG00000180332","ENSG00000102468","ENSG00000100884","ENSG00000168952","ENSG00000139874","ENSG00000165379","ENSG00000070182","ENSG00000165548",
                      "ENSG00000100604","ENSG00000119698","ENSG00000186297","ENSG00000182256","ENSG00000104044","ENSG00000237289","ENSG00000242866","ENSG00000137766",
                      "ENSG00000167178","ENSG00000185518","ENSG00000242173","ENSG00000127585","ENSG00000127561","ENSG00000078328","ENSG00000183454","ENSG00000122254",
                      "ENSG00000166501","ENSG00000006116","ENSG00000149927","ENSG00000099365","ENSG00000102924","ENSG00000070729","ENSG00000154118","ENSG00000272636",
                      "ENSG00000185924","ENSG00000004660","ENSG00000065325","ENSG00000006740","ENSG00000171532","ENSG00000213424","ENSG00000073670","ENSG00000167281",
                      "ENSG00000180777","ENSG00000134508","ENSG00000101489","ENSG00000141668","ENSG00000007264","ENSG00000260001","ENSG00000105376","ENSG00000105642",
                      "ENSG00000105649","ENSG00000197360","ENSG00000105409","ENSG00000204866","ENSG00000204851","ENSG00000118160","ENSG00000104888","ENSG00000161681",
                      "ENSG00000126583","ENSG00000129990","ENSG00000101327","ENSG00000089199","ENSG00000132639","ENSG00000125851","ENSG00000132821","ENSG00000101438",
                      "ENSG00000124134","ENSG00000124140","ENSG00000101180","ENSG00000185272","ENSG00000040608","ENSG00000167037","ENSG00000100285","ENSG00000100346",
                      "ENSG00000186732","ENSG00000130540","ENSG00000138944","ENSG00000073150","ENSG00000169933","ENSG00000169891","ENSG00000184368","ENSG00000008056",
                      "ENSG00000186462","ENSG00000157502","ENSG00000129682","ENSG00000179542","ENSG00000067842","ENSG00000198910")




combined.gene.expression.neuron <- combined.gene.expression %>%
  tibble::rownames_to_column('ens') %>%
  dplyr::filter(ens %in% neuron.genes.ens) %>%
  tibble::column_to_rownames('ens')

pca.neuron <- prcomp(t(combined.gene.expression.neuron))
#plot(pca.neuron)
#library(ggbiplot)
#ggbiplot(pca.neuron)

# correlates nicely to the individual genes
#plot(as.numeric(combined.gene.expression.neuron[4,]), as.numeric(pca.neuron$x %>% as.data.frame %>% dplyr::pull('PC1') )  )


pca.neuron$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pid') %>%
  dplyr::select(c('pid','PC1')) %>%
  dplyr::rename(neuron.component = PC1)



# oligodendrocyte correlated genes ----


oligodendrocyte.genes <- c("RASGRF1","FBXO2","PPP1R16B","TPPP","SEC14L5","TMEM151A","LGI3","TMCC2","HHATL","RAB11FIP4","PDIA2",
                           "HCN2","KCNJ9","DNAJC6","TUBB4A","ADCY5","CSDC2","AC118754.1","PLIN4","HSPA2","PI16","PTGDS","CDK18",
                           "FA2H","AATK","NKX6-2","MAG","PLCH2","FAM131C","PPP1R14A","CHADL","TMEM88B","ABCA2","PLPP2","CERCAM",
                           "BOK","ACP7","CYS1","ANXA3","TNFSF9","AVPI1","MYH11","ADH1B","TUBA4A","CORO6","IL12RB2","TESPA1",
                           "MPP7","RSPO3","KCNJ12","OPN4","MKX","FRAS1","CPNE7")

results.out %>% dplyr::filter(in.gsam == T & in.glass == T) %>% 
  dplyr::filter(hugo_symbol %in% oligodendrocyte.genes) %>% dplyr::pull('ensembl_id')

oligodendrocyte.genes.ens <- c("ENSG00000205116","ENSG00000149527","ENSG00000116661","ENSG00000185519","ENSG00000116675",
                               "ENSG00000081985","ENSG00000162728","ENSG00000133069","ENSG00000117266","ENSG00000205795",
                               "ENSG00000127824","ENSG00000176720","ENSG00000010282","ENSG00000173175","ENSG00000138759",
                               "ENSG00000138772","ENSG00000196616","ENSG00000171368","ENSG00000164530","ENSG00000146374",
                               "ENSG00000168481","ENSG00000167123","ENSG00000107317","ENSG00000107331","ENSG00000150051",
                               "ENSG00000150054","ENSG00000122375","ENSG00000119986","ENSG00000148826","ENSG00000179292",
                               "ENSG00000135426","ENSG00000126803","ENSG00000058335","ENSG00000185615","ENSG00000103184",
                               "ENSG00000133392","ENSG00000103089","ENSG00000178773","ENSG00000183018","ENSG00000184185",
                               "ENSG00000167549","ENSG00000131242","ENSG00000181409","ENSG00000141934","ENSG00000099822",
                               "ENSG00000167676","ENSG00000104833","ENSG00000125657","ENSG00000105695","ENSG00000167641",
                               "ENSG00000183760","ENSG00000101445","ENSG00000100399","ENSG00000172346")




combined.gene.expression.oligodendrocyte <- combined.gene.expression %>%
  tibble::rownames_to_column('ens') %>%
  dplyr::filter(ens %in% oligodendrocyte.genes.ens) %>%
  tibble::column_to_rownames('ens')

pca.oligodendrocyte <- prcomp(t(combined.gene.expression.oligodendrocyte))
plot(pca.oligodendrocyte)
#library(ggbiplot)
#ggbiplot(pca.oligodendrocyte)

# correlates nicely to the individual genes
#plot(as.numeric(combined.gene.expression.oligodendrocyte[4,]), as.numeric(pca.oligodendrocyte$x %>% as.data.frame %>% dplyr::pull('PC1') )  )


pca.oligodendrocyte$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pid') %>%
  dplyr::select(c('pid','PC1')) %>%
  dplyr::rename(oligodendrocyte.component = PC1)



# endothelial correlated genes ----






endothelial.genes <- c("VWF","TIE1","HIGD1B","MMRN1","CYSLTR2","MMP25","FLT4","BCL6B","GRAP","LAMC3","DPEP1","PXDNL","ANGPT2",
                       "PALD1","ADGRD1","GBP6","SLC52A3","CLDN5","VWA2","ABCB1","THSD7B","SPINK8","FOXQ1","ZIC3","NODAL")

results.out %>% dplyr::filter(in.gsam == T & in.glass == T) %>% 
  dplyr::filter(hugo_symbol %in% endothelial.genes) %>% dplyr::pull('ensembl_id')

endothelial.genes.ens <- c("ENSG00000066056","ENSG00000183347","ENSG00000144229","ENSG00000229453","ENSG00000138722",
                           "ENSG00000037280","ENSG00000164379","ENSG00000085563","ENSG00000091879","ENSG00000147485",
                           "ENSG00000050555","ENSG00000156574","ENSG00000107719","ENSG00000165816","ENSG00000110799",
                           "ENSG00000111452","ENSG00000152207","ENSG00000008516","ENSG00000015413","ENSG00000161940",
                           "ENSG00000154016","ENSG00000131097","ENSG00000101276","ENSG00000184113","ENSG00000156925")




combined.gene.expression.endothelial <- combined.gene.expression %>%
  tibble::rownames_to_column('ens') %>%
  dplyr::filter(ens %in% endothelial.genes.ens) %>%
  tibble::column_to_rownames('ens')

pca.endothelial <- prcomp(t(combined.gene.expression.endothelial))
plot(pca.endothelial)
#library(ggbiplot)
#ggbiplot(pca.endothelial)

# correlates nicely to the individual genes
#plot(as.numeric(combined.gene.expression.endothelial[4,]), as.numeric(pca.endothelial$x %>% as.data.frame %>% dplyr::pull('PC1') )  )


pca.endothelial$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pid') %>%
  dplyr::select(c('pid','PC1')) %>%
  dplyr::rename(endothelial.component = PC1)



# extra-cellular matrix correlated genes ----


extracellular.matrix.genes <- c("PRF1","ARHGAP9","FCMR","LXN","KCNE3","NR5A2","FPR2","CCL13","MMP7","CALCR","LRG1","SAA2","PI3","LIF","HSPA6","CRABP2","CILP2","DPT","FGF7","COL10A1","FBN1","GLT8D2","IRX3","MFAP5","MFAP4","COL8A2","FNDC1","MMP11","MFAP2","COL1A2","COL1A1","COL5A1","ADAMTS2","TPSB2","KRT8","OMD","OGN","MME","MLPH","MRC1L1","PTGFR","TWIST2","C5orf46","TNNT3","ASS1","PERP","KLHDC7B","CCL8")

results.out %>% dplyr::filter(in.gsam == T & in.glass == T) %>% 
  dplyr::filter(hugo_symbol %in% extracellular.matrix.genes) %>% dplyr::pull('ensembl_id')

extracellular.matrix.genes.ens <- c("ENSG00000117122","ENSG00000171812","ENSG00000122420","ENSG00000143320","ENSG00000173110","ENSG00000143196","ENSG00000116833","ENSG00000162894","ENSG00000115648","ENSG00000233608","ENSG00000196549","ENSG00000079257","ENSG00000178776","ENSG00000087116","ENSG00000123500","ENSG00000112378","ENSG00000164694","ENSG00000004948","ENSG00000164692","ENSG00000106809","ENSG00000127083","ENSG00000130707","ENSG00000130635","ENSG00000183748","ENSG00000180644","ENSG00000130595","ENSG00000134339","ENSG00000175538","ENSG00000137673","ENSG00000197614","ENSG00000170421","ENSG00000123329","ENSG00000120820","ENSG00000166147","ENSG00000140285","ENSG00000197253","ENSG00000177508","ENSG00000166482","ENSG00000108700","ENSG00000181374","ENSG00000108821","ENSG00000171236","ENSG00000160161","ENSG00000171049","ENSG00000124102","ENSG00000099953","ENSG00000128342","ENSG00000130487")




combined.gene.expression.extracellular.matrix <- combined.gene.expression %>%
  tibble::rownames_to_column('ens') %>%
  dplyr::filter(ens %in% extracellular.matrix.genes.ens) %>%
  tibble::column_to_rownames('ens')

pca.extracellular.matrix <- prcomp(t(combined.gene.expression.extracellular.matrix))
plot(pca.extracellular.matrix)
#library(ggbiplot)
#ggbiplot(pca.extracellular.matrix)

# correlates nicely to the individual genes
#plot(as.numeric(combined.gene.expression.extracellular.matrix[4,]), as.numeric(pca.extracellular.matrix$x %>% as.data.frame %>% dplyr::pull('PC1') )  )


pca.extracellular.matrix$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pid') %>%
  dplyr::select(c('pid','PC1')) %>%
  dplyr::rename(extracellular.matrix.component = PC1)



# combine ----


principal_de_components <- pca.neuron$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pid') %>%
  dplyr::select(c('pid','PC1')) %>%
  dplyr::rename(neuron.component = PC1) %>%
  dplyr::left_join(
    
    pca.oligodendrocyte$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column('pid') %>%
      dplyr::select(c('pid','PC1')) %>%
      dplyr::rename(oligodendrocyte.component = PC1)
    
    , by=c('pid'='pid')
  ) %>% dplyr::left_join(
    
    pca.endothelial$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column('pid') %>%
      dplyr::select(c('pid','PC1')) %>%
      dplyr::rename(endothelial.component = PC1),
    
    by=c('pid'='pid')
    
  ) %>% dplyr::left_join(
    
    pca.extracellular.matrix$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column('pid') %>%
      dplyr::select(c('pid','PC1')) %>%
      dplyr::rename(extracellular.matrix.component = PC1),
    
    by=c('pid'='pid')
  )


write.table(principal_de_components , "output/tables/principal_DE_cluster_components.txt")
  
  



