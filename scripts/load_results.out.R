#!/usr/bin/env R

# load generic file ----

results.out <- readRDS(file = 'tmp/results.out.Rds') |> 
  
  dplyr::mutate(`p.value.glass.cor.tpc` = NULL) |>        # obsolete: 2021 rev
  dplyr::mutate(`estimate.glass.cor.tpc` = NULL) |>       # obsolete: 2021 rev
  dplyr::mutate(`statistic.glass.cor.tpc` = NULL) |>      # obsolete: 2021 rev
  
  dplyr::mutate(`baseMean.glass.res` = NULL) |>           # obsolete: 2021 rev
  dplyr::mutate(`baseMean.glass.tpc.res` = NULL) |>       # obsolete: 2021 rev
  
  dplyr::mutate(`lfcSE.glass.res` = NULL) |>              # obsolete: 2021 rev
  dplyr::mutate(`lfcSE.glass.tpc.res` = NULL) |>          # obsolete: 2021 rev
  
  dplyr::mutate(`log2FoldChange.glass.res` = NULL) |>     # obsolete: 2021 rev
  dplyr::mutate(`log2FoldChange.glass.tpc.res` = NULL) |> # obsolete: 2021 rev
  
  dplyr::mutate(`padj.glass.res` = NULL) |>               # obsolete: 2021 rev
  dplyr::mutate(`padj.glass.tpc.res` = NULL) |>           # obsolete: 2021 rev
  
  dplyr::mutate(`pvalue.glass.res` = NULL) |>             # obsolete: 2021 rev
  dplyr::mutate(`pvalue.glass.tpc.res` = NULL) |>         # obsolete: 2021 rev
  
  dplyr::mutate(`stat.glass.res` = NULL) |>               # obsolete: 2021 rev
  dplyr::mutate(`stat.glass.tpc.res` = NULL)              # obsolete: 2021 rev
  




# DGE GLASS 2022 ----


tmp <- readRDS("tmp/analysis_DGE_GLASS-2022.tpc.Rds") 
#|> 
#  dplyr::rename_with( ~ gsub(".res", ".tpc.res", .x, fixed = TRUE))


results.out <- results.out |> 
  dplyr::left_join(
    tmp,
    by=c('ensembl_id'='ensembl_id'), suffix=c('','')
  )

rm(tmp)


# cor GLASS 2022 ----

tmp <- readRDS('tmp/analysis_cor_purity_expression_GLASS-2022.Rds')

results.out <- results.out |> 
  dplyr::left_join(
    tmp,
    by=c('ensembl_id'='ensembl_id'), suffix=c('','')
  )

rm(tmp)

# go.0031012 COLLAGEN pathway ----

go.0031012 <- c("MMP25","ANOS1","DCN","SEMA3B","SERPINB1","MYOC","TIMP2","VCAN","TNC","USH2A","COL9A2",
                "ADAMTS6","LTBP1","ELN","LAMC3","COL23A1","LAMA3","ENTPD2","CBLN4","ITIH4","ITIH1","F7","
                LAMC2","COL11A1","WNT8A","GPC1","CCN5","CDON","NTN1","COL17A1","FGFR2","PKM","TGFBR3",
                "FGF10","SPP2","ADAM11","NTN4","LMAN1","DLG1","GPC4","FBLN1","COL5A3","COL4A4","IMPG2",
                "COL19A1","EPYC","COL16A1","FCN1","OVGP1","WNT11","ACHE","ADAMTS2","MMP2","NID2","ATRN",
                "LTBP4","ICAM1","P3H2","GLG1","LAMB4","LAMB1","APOH","ANGPT2","CCDC80","CMA1","COL9A3","TGFB2",
                "CRISP3","AMELY","IGFALS","HNRNPM","MMP11","LGALS1","TFIP11","TIMP3","PDGFB","CHADL","CTSG","COCH",
                "PAPLN","MMP9","BMP7","CTSZ","COL20A1","ANGPT4","LAMA1","MXRA5","SRPX","F9","TIMP1","SRPX2","CBLN1",
                "MMP15","CRISPLD2","ZP2","CTSH","SFRP1","CCN4","IL7","FGL1","CLC","APLP1","TGFB1","COMP","TFPI2",
                "WNT2","PTPRZ1","PCOLCE","SERPINE1","PLOD3","RARRES2","AEBP1","MEGF9","OGN","ASPN","ECM2","AMBP",
                "ATRNL1","CXCL12","SPOCK2","KAZALD1","WNT3","LGALS3BP","COL1A1","VTN","SOD3","CTSC","TECTA","HPX",
                "APOA4","APOC3","MDK","VWF","WNT5B","MGP","COL12A1","COL9A1","GPLD1","SMOC2","IMPG1","VEGFA","CCN6",
                "LAMA4","ERBIN","LOX","SPARC","THBS4","FGF1","KNG1","HRG","WNT5A","COL7A1","EFEMP1","FN1","WNT6",
                "TNR","ANGPTL1","PRG4","NID1","ZP4","MFAP2","P3H1","F3","SERPINC1","MUC5B","CTSD","MMP8","APOA1",
                "CCN2","LTBP2","TGFB3","TNN","TGFBI","CLU","A1BG","FMOD","PLG","SERAC1","ANXA11","MMP19","COL10A1",
                "PI3","SLPI","MATN4","NDP","F13A1","SERPINB6","COL21A1","AMELX","FLRT3","MMP24","CFP","PZP","OMD",
                "FGL2","LRRC17","LOXL1","APOE","NCAN","PXDN","GDF15","COL5A1","LAMA5","FIBCD1","ANGPTL6","F12",
                "LGALS3","PODNL1","RTBDN","MATN3","EMILIN2","SERPINF1","ENAM","ITGB4","MATN2","BCAN","HAPLN2",
                "APCS","ANGPTL3","CHI3L1","POSTN","LOXL2","ADAMDEC1","CD180","WNT2B","COL4A2","ADAMTS8","ANXA1",
                "CTSL","ADAM19","AGT","LAMC1","SERPINE2","ADAMTS7","CHAD","ANGPTL2","CCN3","TINAG","LRRC32",
                "SULF1","MMP7","MMP20","MMP27","MMP13","THBS1","EMILIN1","ANXA7","ADAMTS14","CILP","SEMA7A",
                "MMRN1","FRAS1","FBN2","COL2A1","LRIG3","INHBE","PTPRQ","LUM","KERA","FBLN5","MFAP1","ADAMTS17",
                "PCSK6","LMAN1L","HAPLN3","MFGE8","TGFB1I1","ADAMTS18","CDH13","GFOD2","COL6A1","COL6A2","NTN5",
                "ADAMTS10","FBN3","FCN3","HSPG2","CCN1","TINAGL1","DPT","HDGF","HMCN1","ECM1","ADAMTSL4",
                "S100A8","S100A7","FLG","LEFTY2","FBLN7","TMEFF2","LRIG1","COL8A1","MUC4","AHSG","SFRP2",
                "ADAMTS16","HAPLN1","POMZP3","ZAN","CASK","GPC3","HMCN2","HSD17B12","SERPING1","SERPINH1",
                "NCAM1","ZP1","MMP3","FREM2","ADAMTS12","ITIH2","DST","SPARCL1","DSPP","DMP1","MEPE",
                "ABI3BP","ANGPT1","MMP21","ADAMTS1","ADAMTS5","WNT7A","MMP16","ADAMTS3","ADAMTSL3","PRG3",
                "TIMP4","MMP14","ACAN","COLEC12","ADAMTS4","LAD1","C1QC","SPON2","CSTB","ADAMTS13","FCN2",
                "AZGP1","COL26A1","NTN3","MATN1","SDC3","ALPL","WNT4","NTNG1","OLFML2B","LRRTM1","CTSS","S100A9",
                "COL6A3","IGFBP7","IHH","FBLN2","ADAMTS9","PTX3","CRELD1","PF4","CPA3","TGM4","CLEC3B","ANXA5",
                "EDIL3","EGFLAM","TLR3","BMPER","RELL2","SHH","COL1A2","CTSB","TNFRSF11B","SBSPON","CTHRC1","FREM1",
                "MBL2","UCMA","VWA2","OTOGL","HTRA1","ADAMTS15","FBN1","SERPINB8","MFAP4","HSP90B1","SERPINB12","MMP10",
                "NAV2","ELFN2","GREM1","MMP26","SERPINF2","KRT1","ANGPTL4","SOST","LTBP3","SCARA3","VASN","TNXB","COL3A1",
                "NPNT","LGALS9","COL4A3","THBS3","COL22A1","SDC2","LINGO1","MUC17","LRRN2","SERPINB9","CDH2","COL24A1",
                "FGG","FGA","FGB","LGALS4","COL8A2","ANGPTL7","LAMB2","LRRC15","TPSAB1","EFEMP2","COL6A5","EGFL7",
                "HPSE2","HPSE","LRRN3","ADAMTS20","MMRN2","C1QA","NDNF","DAG1","CSPG4","PHOSPHO1","CTSF","PODN",
                "LINGO2","FGFBP3","CD248","ZG16","NPPA","A2M","LRRN1","LRRTM4","CLEC14A","CD151","CPN2","CALR",
                "GPC5","VWA1","LRRC3B","SSC5D","F2","ADIPOQ","OTOL1","BGN","MXRA7","SNORC","ANXA2","COL18A1",
                "FREM3","GPC6","CCBE1","DGCR6","EMILIN3","EFNA5","MUC6","FLRT2","GP1BA","OLFML2A","ADAMTSL5",
                "RTN4RL1","THBS2","PRG2","C17orf58","RTN4RL2","EMID1","ANGPTL5","BCAM","COL4A1","HAPLN4","AMTN",
                "THSD4","DMBT1","COL14A1","EYS","COL4A5","AGRN","OTOG","PLSCR1","ZP3","IFNA2","SERPINA5","COL25A1",
                "VWC2","OPTC","PRELP","NYX","RELN","MMP23B","SPOCK3","SERPINA3","S100A4","NTNG2","PRTN3","LAMA2",
                "MMP1","COL27A1","LAMB3","ANXA4","ANXA6","SERPINA1","TPSB2","COL13A1","ELANE","COL4A6","MFAP5",
                "PSAP","S100A10","ADAMTSL2","HRNR","S100A6","MMP17","SMOC1","LRRTM3","EGFL6","LRIG2","L1CAM",
                "TGM2","COL11A2","COL5A2","COL15A1","C6orf15","LRRC3C","PRSS1","VIT","DEFA1","COL6A6","COLQ","LINGO4","GPC2","ANG","COL28A1","MUC5AC","LINGO3","ELFN1","ORM2","ORM1","DEFA1B","TMEFF1","MARCOL","OC90","LRRC24","TRIL","GH1","MMP12","SPON1","ANXA8","RBP3","GDF10","MMP28","PRSS2")


results.out <- results.out |> 
  dplyr::mutate(EM.GO.0031012 = hugo_symbol %in% go.0031012)


rm(go.0031012)



# RNA-seq rCor clusters  ----

# C6.2021 <- c('CRABP2', 'CILP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
#              'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
#              "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
#              "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
#              "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")



tmp.1 <- readRDS('cache/h.2022.Rds')
tmp.2 <- data.frame(gid = tmp.1$labels[rev(tmp.1$order)]) |> 
  dplyr::mutate(C0.start = which(gid == "ENSG00000119714.11_5|GPR68|chr14:91698876-91720269(-)")) |> 
  dplyr::mutate(C0.end = which(gid == "ENSG00000105697.9_4|HAMP|chr19:35771619-35776046(+)")) |> 
  
  dplyr::mutate(C1.start = which(gid == "ENSG00000130487.8_2|KLHDC7B|chr22:50984632-50989452(+)")) |> 
  dplyr::mutate(C1.end = which(gid == "ENSG00000122420.10_5|PTGFR|chr1:78769568-79006386(+)")) |>  
  
  dplyr::mutate(C2.start = which(gid == "ENSG00000128917.8_4|DLL4|chr15:41221538-41231271(+)")) |> 
  dplyr::mutate(C2.end = which(gid == "ENSG00000085563.14_5|ABCB1|chr7:87132949-87342639(-)")) |> 
  
  dplyr::mutate(C3.start = which(gid == "ENSG00000132702.13_3|HAPLN2|chr1:156589123-156595517(+)")) |> 
  dplyr::mutate(C3.end = which(gid == "ENSG00000099822.3_3|HCN2|chr19:589881-617159(+)")) |> 

  # SEPT12 does not fit with c3 or c4
  
  dplyr::mutate(C4.start = which(gid == "ENSG00000185272.14_3|RBM11|chr21:15588451-15600693(+)")) |> 
  dplyr::mutate(C4.end = which(gid == "ENSG00000178773.15_6|CPNE7|chr16:89642166-89663654(+)")) |> 
  
  
  dplyr::mutate(C4s1.start = which(gid == "ENSG00000185272.14_3|RBM11|chr21:15588451-15600693(+)")) |> 
  dplyr::mutate(C4s1.end = which(gid == "ENSG00000054356.14_5|PTPRN|chr2:220154345-220174370(-)")) |> 

  dplyr::mutate(C4s2.start = which(gid == "ENSG00000133392.18_5|MYH11|chr16:15796992-15950885(-)")) |> 
  dplyr::mutate(C4s2.end = which(gid == "ENSG00000178773.15_6|CPNE7|chr16:89642166-89663654(+)")) |> 
  
  dplyr::mutate(cluster.2022 = case_when(
    dplyr::row_number() >= C0.start & dplyr::row_number() <= C0.end ~ "C0",
    dplyr::row_number() >= C1.start & dplyr::row_number() <= C1.end ~ "C1",
    dplyr::row_number() >= C2.start & dplyr::row_number() <= C2.end ~ "C2",
    dplyr::row_number() >= C3.start & dplyr::row_number() <= C3.end ~ "C3",
    dplyr::row_number() >= C4.start & dplyr::row_number() <= C4.end ~ "C4",
    dplyr::row_number() >= C4.start & dplyr::row_number() <= C4.end ~ "C4",
    T ~ "-"
  )) |> 
  
  dplyr::mutate(C0.start = NULL, C0.end = NULL,
                C1.start = NULL, C1.end = NULL,
                C2.start = NULL, C2.end = NULL,
                C3.start = NULL, C3.end = NULL,
                C4.start = NULL, C4.end = NULL) |> 
  
  dplyr::mutate(cluster.2022 = factor(cluster.2022, levels = c("C0", "C1", "C2", "C3", "C4", "-"))) |> 
  
  dplyr::mutate(C0.2022 = cluster.2022 == "C0") |> 
  dplyr::mutate(C1.2022 = cluster.2022 == "C1") |> 
  dplyr::mutate(C2.2022 = cluster.2022 == "C2") |> 
  dplyr::mutate(C3.2022 = cluster.2022 == "C3") |> 
  dplyr::mutate(C4.2022 = cluster.2022 == "C4") |> 
  
  dplyr::mutate(C4s1.2022 = dplyr::row_number() >= C4s1.start & dplyr::row_number() <= C4s1.end ) |> 
  dplyr::mutate(C4s2.2022 = dplyr::row_number() >= C4s2.start & dplyr::row_number() <= C4s2.end ) |> 
  
  dplyr::mutate(C4s1.start = NULL, C4s1.end = NULL, 
                C4s2.start = NULL, C4s2.end = NULL)

rm(tmp.1)


results.out <- results.out |> 
  dplyr::left_join(tmp.2, by=c('gid'='gid'), suffix=c('','') )


rm(tmp.2)






