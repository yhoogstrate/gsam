#!/usr/bin/env R

# settings ----


options(warnPartialMatchDollar = TRUE) # https://stackoverflow.com/questions/32854683/data-frames-in-r-name-autocompletion


# load libs ----


library(tidyverse)
library(DESeq2)


library(ggplot2)
library(ggrepel)

#library(pheatmap)
library(fgsea)
#library(limma)

library(EnhancedVolcano)
library(patchwork)

library(ggsignif)
library(xlsx)


# load data ----

#source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
#ensembl_genes <- get_ensembl_hsapiens_gene_ids()

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


source("scripts/R/ligands.R")
source("scripts/R/subtype_genes.R") # Verhaak/Wang/TCGA + GliTS Redux
source("scripts/R/patel_scRNA_clusters.R")
source("scripts/R/neftel_meta_modules.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

source('scripts/R/wang_glioma_intrinsic_genes.R')

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set + metedata

source("scripts/R/cor_cor_plot.R")


#m <- c("ward.D", "ward.D2", "single", "complete", "average" , "mcquitty" , "median" , "centroid")
m <- "use correlation of correlation instead?!"


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




gsam.metadata.r1 <- gsam.metadata.all %>%
  dplyr::filter(resection == "r1")

gsam.metadata.r2 <- gsam.metadata.all %>%
  dplyr::filter(resection == "r2")


# uncomment out the minimum TPC
wilcox.test(gsam.metadata.r1 %>% dplyr::pull(tumour.percentage.dna),
            gsam.metadata.r2 %>% dplyr::pull(tumour.percentage.dna), alternative = "two.sided")




gsam.gene.expression.all <- gsam.rnaseq.expression %>%
  dplyr::select(gsam.metadata.all$sid)
stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)


gsam.gene.expression.all.paired <- gsam.gene.expression.all %>%
  dplyr::select(gsam.metadata.all.paired$sid)
stopifnot(colnames(gsam.gene.expression.all.paired) == gsam.metadata.all.paired$sid)



# gsam.bfg.expression.all <- gsam.rnaseq.expression %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   dplyr::filter( # bona fide glioma genes (BFGs)
#     (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
#       (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
#   ) %>%
#   tibble::column_to_rownames('gid') %>%
#   dplyr::select(gsam.metadata.all$sid)
# 
# stopifnot(colnames(gsam.bfg.expression.all) == gsam.metadata.all$sid)


# gsam.gene.expression.all.vst
# TODO gsam.gene.expression.R1.vst
# TODO gsam.gene.expression.R2.vst




gsam.gene.expression.all.vst <- gsam.gene.expression.all %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay()




## GLASS ----



glass.metadata.all  <- glass.gbm.rnaseq.metadata


glass.gene.expression.all <- glass.gbm.rnaseq.expression %>%
  dplyr::select(glass.metadata.all$sid)


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)



glass.metadata.tp <- glass.metadata.all %>%
  dplyr::filter(resection == "TP")

glass.metadata.tr <- glass.metadata.all %>%
  dplyr::filter(resection != "TP")


# uncomment out the minimum TPC
wilcox.test(glass.metadata.tp %>% dplyr::pull(tumour.percentage.dna.imputed.rf) ,
            glass.metadata.tr %>% dplyr::pull(tumour.percentage.dna.imputed.rf) , alternative="two.sided")


# ggplot(glass.metadata.all %>%
#          dplyr::mutate(resection.status = ifelse(resection == "TP", "primary" , "recurrent")),
#        aes(x = resection.status, y = tumour.percentage.dna.imputed.rf)) +
#   geom_violin() +
#   #geom_segment(aes(x=0.55, y=m1.gsam, xend=1.45, yend=m1.gsam),lty=1,lwd=0.15, col="gray40") +
#   #geom_segment(aes(x=1.55, y=m2.gsam, xend=2.45, yend=m2.gsam),lty=1,lwd=0.15, col="gray40") +
#   geom_jitter( position=position_jitter(0.2), size=0.9) +
#   ylim(25, 80) +
#   labs(x = NULL, col = NULL, y = "WES estimate tumour cell percentage" ) +
#   job_gg_theme



## per-gene results table ----

# how can these appear:
#31        TGFA      ENSG00000163235.16_5|TGFA|chr2:70674416-70781325(-) oligodendrocyte
#32        TGFA                                                     <NA> oligodendrocyte
# 2 distinct ENS id's


results.out <- dplyr::full_join(
    gsam.gene.expression.all %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::select(gid) %>%
    dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
    dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
    dplyr::mutate(in.gsam = T) %>%
    dplyr::mutate(chr = as.factor(gsub("^.+(chr[^:]+):.+$","\\1",gid)))
    ,
    glass.gene.expression.all %>%
      tibble::rownames_to_column('ensembl_id') %>%
      dplyr::select(ensembl_id) %>% 
      dplyr::mutate(in.glass = T)
    , by = c('ensembl_id'='ensembl_id')) %>%
  dplyr::left_join(
    glass.gencode.v19 %>%
      dplyr::arrange( transcript_id) %>%
      dplyr::select(c('gene_symbol','gene_id')) %>%
      dplyr::filter(!duplicated(gene_id)) %>%
      dplyr::rename(ensembl_id = gene_id) %>%
      dplyr::rename(hugo_symbol.gencode = gene_symbol) ,
    by = c ('ensembl_id'='ensembl_id')) %>%
  dplyr::mutate(hugo_symbol = ifelse(is.na(hugo_symbol) , hugo_symbol.gencode , hugo_symbol )) %>%
  dplyr::mutate(hugo_symbol.gencode = NULL) %>%
  dplyr::mutate(is.bfg = (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) | (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol))
    


stopifnot(sum(duplicated(results.out$ensembl_id)) == 0)



### add pharmacogenomics genes ----

# look for transporter genes / influx / efflux genes within a volcano like plot [hedgehog signalling?]
# PharmGKB: ABCB1, OPRM1, INHBA, ITGBL1, FCGR2B

r1 <- read.delim('data/pharmgkb/relationships.tsv',sep="\t",stringsAsFactors = F) %>%
  dplyr::filter((Entity1_type == "Gene" & Entity2_type == "Chemical")) %>%
  dplyr::filter(Association == "not associated") %>%
  dplyr::pull(Entity1_name) %>%
  unique()

r2<- read.delim('data/pharmgkb/relationships.tsv',sep="\t",stringsAsFactors = F) %>%
  dplyr::filter((Entity2_type == "Gene" & Entity1_type == "Chemical")) %>%
  dplyr::filter(Association == "not associated") %>%
  dplyr::pull(Entity2_name) %>%
  unique()


results.out <- results.out %>%
  dplyr::mutate(pharma.relation = hugo_symbol %in% c(r1, r2))

rm(r1, r2)








### McKenzie cell type labels ----

dim(results.out)

#### top_all_enrich

# n = 5809
# tmp <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_all_enrich',skip=2) %>%
#   dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
#   dplyr::arrange(desc(grand_mean)) %>%
#   dplyr::filter(gene %in% results.out$hugo_symbol ) %>%
#   dplyr::mutate(Celltype = case_when(
#     Celltype == "ast" ~ "astrocyte" ,
#     Celltype == "end" ~ "endothelial",
#     Celltype == "mic" ~ "microglia", 
#     Celltype == "neu" ~ "neuron",
#     Celltype == "oli" ~ "oligodendrocyte",
#     Celltype == "opc" ~ "OPC")) %>%
#   dplyr::group_by(Celltype) %>%
#   dplyr::slice_head(n=75) %>%
#   dplyr::filter(gene %in% unique(.[["gene"]][duplicated(.[["gene"]])]) == F) %>% # only those that match 1 cell type
#   dplyr::rename(McKenzie_celltype_top_all_enrich = Celltype) %>%
#   dplyr::mutate(grand_mean = NULL)
# dim(tmp)
# 
# results.out$McKenzie_celltype_top_all_enrich = NULL
# results.out <- results.out %>%
#   dplyr::left_join(tmp, by=c('hugo_symbol'='gene'))
# rm(tmp)
# dim(results.out)
# 


#### top_human_specificity ----


# monocyten/macrofagen uit bone marrow: CD163"
# T-cells: "ITGA5", "ITGB1", "MSN", "FAS", "FLNA", "CD44", "RUNX1", "RUNX2"
# TAMs: "NRP1", "SPP1", "LYN", "LIMS1", "C5AR1", "PLAUR", "CEBPB"
# 'key immune markers': "CD3", "CD68"
# angiogenesis: 'FLT1', 'MMP14', 'ENG', 'SERPINE1'

# n = 3943
tmp <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_specificity') %>%
  dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
  dplyr::arrange(desc(grand_mean)) %>%
  dplyr::filter(gene %in% results.out$hugo_symbol ) %>%
  dplyr::mutate(Celltype = case_when(
    Celltype == "ast" ~ "astrocyte" ,
    Celltype == "end" ~ "endothelial",
    Celltype == "mic" ~ "microglia/TAM", 
    Celltype == "neu" ~ "neuron",
    Celltype == "oli" ~ "oligodendrocyte",
    Celltype == "opc" ~ "OPC")) %>%
  dplyr::group_by(Celltype) %>%
  dplyr::slice_head(n=200) %>%
  dplyr::filter(gene %in% unique(.[["gene"]][duplicated(.[["gene"]])]) == F) %>% # only those that match 1 cell type
  dplyr::ungroup() %>%
  dplyr::add_row(gene = "CD163",Celltype = "microglia/TAM") %>%
  dplyr::add_row(gene = "RBFOX3",Celltype = "neuron") %>%
  dplyr::add_row(gene = "FLT1",Celltype = "endothelial") %>%
  dplyr::mutate(show.marker = ifelse(gene %in% c(
    "TBX3","ERG", "FLT1","RERGL","VCAM1" , # endothelial
    "CD74", "CD14", "CD84", "CD53", "CD163", # TAM/Microglia
    "OPALIN","PLP1","CNDP1","CNTN2","TMEM144", # oligodendrocyte
    "PAX6","ELOVL2","FGFR3","SOX9", "IL33","GLI3", "SLC4A4", # astrocyte
    
    "BCAS1","ERBB3","GFAP","GABRA1","GABRB2","RBFOX3","TGFA","VIP"
    
    ) , T , F) ) %>%
  dplyr::rename(McKenzie_celltype_top_human_specificity = Celltype) %>%
  dplyr::mutate(grand_mean = NULL)
#dim(tmp)

results.out$McKenzie_celltype_top_human_specificity <- NULL
results.out$show.marker <- NULL
results.out <- results.out %>%
  dplyr::left_join(tmp, by=c('hugo_symbol'='gene')) %>%
  dplyr::mutate(show.marker.chr7 = ifelse(hugo_symbol %in% c("PDGFA","GNA12","ETV1","CBX3","CREB5","TRIL","GLI3","COA1","EGFR","PTPN12","CDK6","DOCK4","PTPRZ1","TRIM24"), T , F) )
rm(tmp)
dim(results.out)


# top_human_expression = non informative



#### top_human_enrich

# # n = 4380
# tmp <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_enrich') %>%
#   dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
#   dplyr::arrange(desc(grand_mean)) %>%
#   dplyr::filter(gene %in% results.out$hugo_symbol ) %>%
#   dplyr::mutate(Celltype = case_when(
#     Celltype == "ast" ~ "astrocyte" ,
#     Celltype == "end" ~ "endothelial",
#     Celltype == "mic" ~ "microglia",
#     Celltype == "neu" ~ "neuron",
#     Celltype == "oli" ~ "oligodendrocyte",
#     Celltype == "opc" ~ "OPC")) %>%
#   dplyr::group_by(Celltype) %>%
#   dplyr::slice_head(n=75) %>%
#   dplyr::filter(gene %in% unique(.[["gene"]][duplicated(.[["gene"]])]) == F) %>% # only those that match 1 cell type
#   dplyr::rename(McKenzie_celltype_top_human_enrich = Celltype) %>%
#   dplyr::mutate(grand_mean = NULL)
# #dim(tmp)
# 
# results.out$McKenzie_celltype_top_human_enrich = NULL
# results.out <- results.out %>%
#   dplyr::left_join(tmp, by=c('hugo_symbol'='gene'))
# rm(tmp)
# dim(results.out)


### TCGA / Wang / Verhaak 2017 subtype signature ----


results.out <- results.out %>%
  dplyr::mutate(TCGA.subtype.marker = NULL) %>% # remove if already exists
  dplyr::left_join(rbind(
    data.frame(
      TCGA.subtype.marker = "TCGA-CL",
      ensembl_id = gsub("^([A-Z0-9]+).+$", "\\1", subtype.classical$ENSG)
    ),
    data.frame(
      TCGA.subtype.marker = "TCGA-PN",
      ensembl_id = gsub("^([A-Z0-9]+).+$", "\\1", subtype.proneural$ENSG)
    ),
    data.frame(
      TCGA.subtype.marker = "TCGA-MES",
      ensembl_id = gsub("^([A-Z0-9]+).+$", "\\1", subtype.mesenchymal$ENSG)
    )
  ) ,
  by = c('ensembl_id' = 'ensembl_id'))



### GliTS redux signature ----



results.out <- results.out %>%
  dplyr::mutate(GliTS.reduxsubtype.marker = NULL) %>% # remove if already exists
  dplyr::left_join(rbind(
    data.frame(
      GliTS.reduxsubtype.marker = "GliTS-CL",
      ensembl_id = gsub("^([A-Z0-9]+).+$", "\\1", classical.glits.redux$ENSG)
    ),
    data.frame(
      GliTS.reduxsubtype.marker = "GliTS-PN",
      ensembl_id = gsub("^([A-Z0-9]+).+$", "\\1", proneural.glits.redux$ENSG)
    ),
    data.frame(
      GliTS.reduxsubtype.marker = "GliTS-MES",
      ensembl_id = gsub("^([A-Z0-9]+).+$", "\\1", mesenchymal.glits.redux$ENSG)
    )
  ) ,
  by = c('ensembl_id' = 'ensembl_id'))



### Patel scRNA-seq clusters ----



results.out <- results.out %>%
  dplyr::mutate(patel.scRNAseq.cluster = NULL) %>% # reset
  dplyr::mutate(patel.scRNAseq.cluster = case_when(
    hugo_symbol %in% patel.cell.cycle.tt2 ~ "Cell cycle",
    hugo_symbol %in% patel.complement.immune.response.tt2 ~ "Complete/Immune response",
    hugo_symbol %in% patel.hypoxia.tt2 ~ "Hypoxia",
    T ~ "")) %>%
  dplyr::mutate(patel.scRNAseq.cluster = ifelse(patel.scRNAseq.cluster == "", NA, patel.scRNAseq.cluster))


results.out %>% dplyr::filter(!is.na(patel.scRNAseq.cluster) )





### Neftel scRNA-seq 6 meta clusters ----

# too many overlaps in the clusters, also AC/ODS, no mutex state


results.out <- results.out %>%
  dplyr::mutate(neftel.meta.module.MES1 = hugo_symbol %in% neftel.meta.modules.MES1.tt2)  %>%
  dplyr::mutate(neftel.meta.module.MES2 = hugo_symbol %in% neftel.meta.modules.MES2.tt2)  %>%
  dplyr::mutate(neftel.meta.module.AC = hugo_symbol %in% neftel.meta.modules.AC.tt2)  %>%
  dplyr::mutate(neftel.meta.module.OPC = hugo_symbol %in% neftel.meta.modules.OPC.tt2)  %>%
  dplyr::mutate(neftel.meta.module.NPC1 = hugo_symbol %in% neftel.meta.modules.NPC1.tt2)  %>%
  dplyr::mutate(neftel.meta.module.NPC2 = hugo_symbol %in% neftel.meta.modules.NPC2.tt2)  %>%
  dplyr::mutate(neftel.meta.module.G1.S = hugo_symbol %in% neftel.meta.modules.G1.S.tt2)  %>%
  dplyr::mutate(neftel.meta.module.G2.M = hugo_symbol %in% neftel.meta.modules.G2.M.tt2) 



### GO:0005201 extracellular matrix structural constituent ----


results.out <- results.out %>%
  dplyr::mutate(EM.struct.constituent = ensembl_id %in%
  (read.csv('data/GO0005201_extracellular_matrix_structural_constituent.csv') %>%
      dplyr::pull(converted_alias))
  )



## plotting subset type ----

# voeg toe:
# - present in non-malignant cells [x] << gebruik ofwel levi's data ofwel een van die twee studies
# - present in malignent cells [x]

# McKenzie markers:





results.out <- results.out %>% 
  dplyr::mutate(primary.marker.genes = '') %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("RBFOX3", "GABRB2","GABRA1","GABRG2") , "neuron", primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("CACHD1","AHCYL1","GPR37L1","BMPR1B", "GFAP") , "astrocyte", primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("PLP1","OPALIN", "TMEM144","MOG") , "oligodendrocyte", primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("TIE1","RGS5","NOSTRIN","FLT1") , "endothelial", primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("CD163","CD14","C1QA","THEMIS2") , "microglia/TAM", primary.marker.genes))

stopifnot(
  results.out %>%
    dplyr::filter(primary.marker.genes != "") %>%
    dplyr::select(c(hugo_symbol, primary.marker.genes, `McKenzie_celltype_top_human_specificity` , show.marker)) %>%
    dplyr::pull(primary.marker.genes)
  ==
  results.out %>%
    dplyr::filter(primary.marker.genes != "") %>%
    dplyr::select(c(hugo_symbol, primary.marker.genes, `McKenzie_celltype_top_human_specificity` , show.marker)) %>%
    dplyr::pull(McKenzie_celltype_top_human_specificity)
  )



# custom markers

results.out <- results.out %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("CREB5","TRIM24","ETV1","COA1","SETD5","SMAD5") , "chr7 gained (tumor)", primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("CD45") , "non-malignant", primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in% c("CD2", "CD3D", "CD3E","CD8A") , "TIL / T-cell", primary.marker.genes)) # Van Levi



# subtype markers

results.out %>%
  dplyr::filter(!is.na(TCGA.subtype.marker)) %>%
  dplyr::select(c( gid, McKenzie_celltype_top_human_specificity, TCGA.subtype.marker)) %>%
  dplyr::filter(TCGA.subtype.marker=='TCGA-CL')



results.out <- results.out %>%
  dplyr::mutate(primary.marker.genes = NA) %>% 
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in%  c("COL1A1", "LAMB1", "RUNX1", "S100A11"), paste0(primary.marker.genes, ",MES subtype"), primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in%  c("ELOVL2", "VAV3", "SOCS2", "MLC1"), paste0(primary.marker.genes, ",CL subtype"), primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = ifelse(hugo_symbol %in%  c("PLAAT1", "DNM3", "ERBB3","SOX10"), paste0(primary.marker.genes, ",PN subtype"), primary.marker.genes)) %>%
  dplyr::mutate(primary.marker.genes = gsub("^[,]+","", primary.marker.genes) ) 



# aesthetics

results.out <- results.out %>%
  dplyr::mutate(primary.marker.genes = ifelse(primary.marker.genes == "", NA, primary.marker.genes))



# final view?

results.out %>% dplyr::filter(primary.marker.genes != "" & !is.na(primary.marker.genes))



# 
# #"PDGFA", "PDGFRA", "OLIG1", "OLIG2", "OLIG3", # TCGA/PN?
# 
# 
# ###, "MME", "ERG", "FCER2", "EPCAM", "EREG" << !!
# ,"EGFR"
# ,"EREG","AREG", "BTC","HBEGF","NGF","TGFA","EGF","EPGN",
# 
# #"ARHGAP28","RHGEF26","BVES-AS1","BVES","CACNA2D1","CDH4","CNGA3","COL11A1","ELOVL2","ETV4","EVA1A","FGFR3","GNAI1","LFNG","LHFPL6","POPDC3","PPARGC1A","RFX4","RNF180","ROBO2","SLC24A3","SOCS2-AS1","SOCS2","SOX9","SULF1","TACR1","TAP1","VAV3"
# 
# "TOP2A", "CDK1", "DTL", "CCNB1", "XRCC2", "CCNE2", "DSN1", "TIMELESS", # Cell Cycle genes Patel/Bernstein
# "VEGFA", "ADM", "TREM1", "ENO1", "LDHA", "NRN1", "UBC", "GBE1", "MIF" # Hypoxia genes Patel Bernstein



# corr TPC [GSAM] ----


tmp <- gsam.gene.expression.all.vst %>%
  as.data.frame %>%
  dplyr::select(everything(),-contains("-new")) %>%
  as.matrix

tmp.tpc <- gsam.metadata.all %>%
  dplyr::filter(sid %in% colnames(tmp)) %>%
  dplyr::pull(tumour.percentage.dna, name=sid)


stopifnot(colnames(tmp) == names(tmp.tpc))



gsam.gene.expression.all.cor.estimate <- data.frame(apply(tmp, 1, function (x)  cor.test(x, tmp.tpc) %>% purrr::pluck('estimate') ), stringsAsFactors = F) %>%
  `colnames<-`("estimate") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor.statistic <- data.frame(apply(tmp, 1, function (x)  cor.test(x, tmp.tpc) %>% purrr::pluck('statistic') ), stringsAsFactors = F) %>%
  `colnames<-`("statistic") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor.p.value <- data.frame(apply(tmp, 1, function (x)  cor.test(x, tmp.tpc) %>% purrr::pluck('p.value') ), stringsAsFactors = F) %>%
  `colnames<-`("p.value") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor <- gsam.gene.expression.all.cor.estimate %>%
  dplyr::left_join(gsam.gene.expression.all.cor.statistic , by=c('gid' = 'gid') ) %>%
  dplyr::left_join(gsam.gene.expression.all.cor.p.value , by=c('gid' = 'gid') ) %>%
  `colnames<-`(paste0(colnames(.), ".gsam.cor.tpc")) %>%
  dplyr::rename(gid = gid.gsam.cor.tpc)

rm(tmp, tmp.tpc, gsam.gene.expression.all.cor.estimate , gsam.gene.expression.all.cor.statistic , gsam.gene.expression.all.cor.p.value)



results.out <- results.out %>%
  dplyr::left_join(gsam.gene.expression.all.cor , by = c('gid' = 'gid')) 

stopifnot("statistic.gsam.cor.tpc" %in% colnames(gsam.gene.expression.all.cor))
stopifnot("statistic.gsam.cor.tpc" %in% colnames(results.out))




# 
# 
# plot(gsam.gene.expression.all.cor$statistic,abs(gsam.gene.expression.all.cor$estimate) + runif(nrow(gsam.gene.expression.all.cor), 0, 0.2) , pch=19,cex=0.05)
# 
# 
# 
# stopifnot(gsam.res.tpc.all$gid == gsam.gene.expression.all.cor$gid)
# plot(gsam.res.tpc.all$log2FoldChange , gsam.gene.expression.all.cor$statistic , pch=19,cex=0.05)
# 
# 
# 
# a = gsam.gene.expression.all.cor %>%
#   dplyr::arrange(-p.value) %>%
#   dplyr::top_n(500, -p.value) %>%
#   dplyr::pull(hugo_symbol)
# 
# b = gsam.gene.expression.all.cor %>%
#   dplyr::arrange(-p.value) %>%
#   dplyr::top_n(500, -p.value) %>%
#   dplyr::pull(ensembl_id)



# corr TPC [GLASS] ----




#glass.gbm.rnaseq.expression.vst.old <- glass.gbm.rnaseq.expression.vst
#glass.metadata.all.old <- glass.metadata.all



# sample with '99' pct tumor in wrong determination
#glass.gbm.rnaseq.expression.vst <- glass.gbm.rnaseq.expression.vst %>%
#  as.data.frame %>%
#  dplyr::select(-c('GLSS-HF-2919-TP'))
#glass.metadata.all <- glass.metadata.all %>%
#  dplyr::filter(sid != 'GLSS-HF-2919-TP' )


#glass.gbm.rnaseq.expression.vst <- glass.gbm.rnaseq.expression.vst %>%
#  as.data.frame %>%
#  dplyr::select(contains('GLSS'))
#glass.metadata.all <- glass.metadata.all %>%
#  dplyr::filter(grepl("GLSS",sid))



stopifnot(colnames(glass.gbm.rnaseq.expression.vst) == glass.metadata.all$sid)
glass.gene.expression.all.cor.estimate <- data.frame(apply(glass.gbm.rnaseq.expression.vst, 1, function (x)  cor.test(x, glass.metadata.all %>% dplyr::pull(tumour.percentage.dna.imputed.rf)) %>% purrr::pluck('estimate') ), stringsAsFactors = F) %>%
  `colnames<-`("estimate") %>%
  tibble::rownames_to_column('gid')

glass.gene.expression.all.cor.statistic <- data.frame(apply(glass.gbm.rnaseq.expression.vst, 1, function (x)  cor.test(x, glass.metadata.all %>% dplyr::pull(tumour.percentage.dna.imputed.rf)) %>% purrr::pluck('statistic') ), stringsAsFactors = F) %>%
  `colnames<-`("statistic") %>%
  tibble::rownames_to_column('gid')

glass.gene.expression.all.cor.p.value <- data.frame(apply(glass.gbm.rnaseq.expression.vst, 1, function (x)  cor.test(x, glass.metadata.all %>% dplyr::pull(tumour.percentage.dna.imputed.rf)) %>% purrr::pluck('p.value') ), stringsAsFactors = F) %>%
  `colnames<-`("p.value") %>%
  tibble::rownames_to_column('gid')

glass.gene.expression.all.cor <- glass.gene.expression.all.cor.estimate %>%
  dplyr::left_join(glass.gene.expression.all.cor.statistic , by=c('gid' = 'gid') ) %>%
  dplyr::left_join(glass.gene.expression.all.cor.p.value , by=c('gid' = 'gid') ) %>%
  `colnames<-`(paste0(colnames(.), ".glass.cor.tpc")) %>%
  dplyr::rename(gid = gid.glass.cor.tpc)

rm(glass.gene.expression.all.cor.estimate , glass.gene.expression.all.cor.statistic , glass.gene.expression.all.cor.p.value)



results.out <- results.out %>%
  dplyr::left_join(glass.gene.expression.all.cor , by = c('ensembl_id' = 'gid')) 

stopifnot("statistic.gsam.cor.tpc" %in% colnames(gsam.gene.expression.all.cor))
stopifnot("statistic.gsam.cor.tpc" %in% colnames(results.out))




#plot(results.out$stat.gsam.res , results.out$stat.glass.res)
a <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc) &!is.na(statistic.glass.cor.tpc))
cor(a$statistic.gsam.cor.tpc , a$statistic.glass.cor.tpc)
#plot(a$statistic.gsam.cor.tpc , a$statistic.glass.cor.tpc,pch=19,cex=0.01)


# all samples: 0.493
# without the 99% sample: 0.4999
# without the TCGA samples: 0.507





# plot(results.out$statistic.gsam.cor.tpc, results.out$statistic.glass.cor.tpc, pch=19, cex=0.4, ylim = c(-10,10), col="gray")
# abline(h=0, col="red")
# abline(v=0, col="red")
# a <- results.out %>% dplyr::filter(statistic.gsam.cor.tpc < -13)
# text(a$statistic.gsam.cor.tpc, a$statistic.glass.cor.tpc, a$hugo_symbol, pos=4 ,col="blue", cex=0.8)
# a <- results.out %>% dplyr::filter(statistic.gsam.cor.tpc > 9.75)
# a <- results.out %>% dplyr::filter(hugo_symbol %in% c('EGFR','SEC61G','ITGB8',
#   'LFNG','PDGFA','DENND2A','CDKN2C','ZNF558','PRPF31','SNX5','NCOA3','POU3F2','NFIA','SOX9','BTC','SOCS2',
#   'EREG','AREG','EGF','TGFA'
#   ))
# text(a$statistic.gsam.cor.tpc, a$statistic.glass.cor.tpc, a$hugo_symbol, pos=2 ,col="blue", cex=0.8)






# DE unpaired [G-SAM] ----


gsam.gene.res.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                              dplyr::filter(rowSums(.) > ncol(.) * 3)
                                            , gsam.metadata.all, ~resection ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>%
  as.data.frame() %>%
  #dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  `colnames<-`(paste0(colnames(.),".gsam.res")) %>%
  dplyr::rename(gid = gid.gsam.res)

results.out <- results.out %>%
  dplyr::mutate(baseMean.gsam.res = NULL) %>%
  dplyr::mutate(log2FoldChange.gsam.res = NULL) %>%
  dplyr::mutate(lfcSE.gsam.res = NULL) %>%
  dplyr::mutate(stat.gsam.res = NULL) %>%
  dplyr::mutate(pvalue.gsam.res = NULL) %>%
  dplyr::mutate(padj.gsam.res = NULL) %>%
  dplyr::mutate(significant.gsam.res = NULL) %>%
  dplyr::left_join(gsam.gene.res.res , by = c('gid' = 'gid'))



# DE unpaired ~ TPC [G-SAM] ----


gsam.gene.res.tpc.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                                  dplyr::filter(rowSums(.) > ncol(.) * 3),
                                                gsam.metadata.all, ~tpc + resection ) %>% # + resection; corrected for tpc
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  #dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj, log2FoldChange) %>%
  tibble::rownames_to_column('gid') %>%
  `colnames<-`(paste0(colnames(.),".gsam.tpc.res")) %>%
  dplyr::rename(gid = gid.gsam.tpc.res)

results.out <- results.out %>%
  dplyr::mutate(baseMean.gsam.tpc.res = NULL) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = NULL) %>%
  dplyr::mutate(lfcSE.gsam.tpc.res = NULL) %>%
  dplyr::mutate(stat.gsam.tpc.res = NULL) %>%
  dplyr::mutate(pvalue.gsam.tpc.res = NULL) %>%
  dplyr::mutate(padj.gsam.tpc.res = NULL) %>%
  dplyr::mutate(significant.gsam.tpc.res = NULL) %>%
  dplyr::left_join(gsam.gene.res.tpc.res , by = c('gid' = 'gid'))



# DE unpaired all [GLASS] ----


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)

glass.gene.res.res <- DESeqDataSetFromMatrix(glass.gene.expression.all %>%
                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
                                             , glass.metadata.all, ~condition ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  #dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5) %>%
  dplyr::arrange(padj, log2FoldChange) %>%
  tibble::rownames_to_column('ensembl_id') %>% 
  `colnames<-`(paste0(colnames(.),".glass.res")) %>%
  dplyr::rename(ensembl_id = ensembl_id.glass.res)

colnames(glass.gene.res.res)


results.out <- results.out %>%
  dplyr::mutate(baseMean.glass.res = NULL) %>%
  dplyr::mutate(log2FoldChange.glass.res = NULL) %>%
  dplyr::mutate(lfcSE.glass.res = NULL) %>%
  dplyr::mutate(stat.glass.res = NULL) %>%
  dplyr::mutate(pvalue.glass.res = NULL) %>%
  dplyr::mutate(padj.glass.res = NULL) %>%
  dplyr::left_join(glass.gene.res.res, by = c('ensembl_id' = 'ensembl_id'))


glass.gene.expression.all.vst <- DESeqDataSetFromMatrix(glass.gene.expression.all %>%
                                                         dplyr::filter(rowSums(.) > ncol(.) * 3)
                                                       , glass.metadata.all, ~condition ) %>% # + resection
  assay() %>%
  vst(blind=T)




# DE unpaired ~ TPC [GLASS : RNA imputed TPC] ----



stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)

glass.gene.res.tpc.res <- DESeqDataSetFromMatrix(glass.gene.expression.all %>%
                                                   dplyr::filter(rowSums(.) > ncol(.) * 3),
                                                 glass.metadata.all %>%
                                                   dplyr::mutate(tpc = 1 - (tumour.percentage.dna.imputed.rf / 100)),
                                                 ~tpc + condition ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  #dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5) %>%
  dplyr::arrange(padj, log2FoldChange) %>%
  tibble::rownames_to_column('ensembl_id') %>% 
  `colnames<-`(paste0(colnames(.),".glass.tpc.res")) %>%
  dplyr::rename(ensembl_id = ensembl_id.glass.tpc.res)


results.out <- results.out %>%
  dplyr::mutate(baseMean.glass.tpc.res = NULL) %>%
  dplyr::mutate(log2FoldChange.glass.tpc.res = NULL) %>%
  dplyr::mutate(lfcSE.glass.tpc.res = NULL) %>%
  dplyr::mutate(stat.glass.tpc.res = NULL) %>%
  dplyr::mutate(pvalue.glass.tpc.res = NULL) %>%
  dplyr::mutate(padj.glass.tpc.res = NULL) %>%
  dplyr::left_join(glass.gene.res.tpc.res, by = c('ensembl_id' = 'ensembl_id'))




















# ::::::::::::: ----
# Export (!) ----

# uncomment when updating code above
# saveRDS(results.out, file="tmp/results.out.Rds")

# ::::::::::::: ----


# __________ ----

# :::::::::::::::::::: ----
# Import (quick) ----


results.out <- readRDS(file = 'tmp/results.out.Rds')


# stats ----

# n significant G-SAM


stat <- results.out %>%
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") )


stat %>%
  dplyr::filter(padj.gsam.tpc.res <= 0.01 & abs(log2FoldChange.gsam.tpc.res) > 0.5 & direction.gsam.tpc.res == direction.glass.tpc.res) %>%
  dim

stat %>%
  dplyr::filter(padj.glass.tpc.res <= 0.01 & abs(log2FoldChange.glass.tpc.res) > 0.5 & direction.gsam.tpc.res == direction.glass.tpc.res) %>%
  dim



# :::::::::::::::::::: ----


# hypeR enrichment ---- 
# 
# # ensembl to entrez?
# # data(examplePathways)
# # data(exampleRanks)
# # # https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# 
# 
# ## https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# 
# genesets <- list()
# 
# y <- results.out %>% 
#   dplyr::filter(!is.na(stat.gsam.tpc.res)) %>%
#   dplyr::filter(!is.na(hugo_symbol)) %>%
#   dplyr::filter(lfcSE.gsam.tpc.res < 0.3) %>%
#   #dplyr::filter(baseMean >= 10) %>% 
#   dplyr::distinct() %>% 
#   dplyr::group_by(hugo_symbol) %>% 
#   dplyr::summarize(stat.gsam.tpc.res=mean(stat.gsam.tpc.res))
# 
# y.ordered <- (y %>% dplyr::select(hugo_symbol, stat.gsam.tpc.res) %>% dplyr::arrange(stat.gsam.tpc.res)) %>% tibble::deframe()
# y.ordered.abs <- (y %>% dplyr::select(hugo_symbol, stat.gsam.tpc.res) %>% dplyr::arrange(abs(stat.gsam.tpc.res))) %>% tibble::deframe()
# 
# 
# for(a in names(genesets)) {
#   for(b in names(genesets[[a]])) {
#     if(grepl("blood vessel",b,ignore.case = T) |   grepl("adasd vascul",b,ignore.case = T) ) {
#       print(a)
#       print(b)
#       print("")
#     }
#   }
# }
# 
# 
# grepl("asd","AsD",ignore.case = T)


## Hallmark genes ----

genesets$HALLMARK     <- hypeR::msigdb_gsets("Homo sapiens", "H", "", clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$HALLMARK, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)
tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)

z$plots[1]
z$plots[2]
z$plots[3]


z$plots[[which(rownames(z$as.data.frame()) == "Coagulation")]]


## Positional ----

# TODO: chr19 GLASS?

genesets$POSITIONAL   <- hypeR::msigdb_gsets("Homo sapiens", "C1", "", clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$POSITIONAL, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)

tibble::as_tibble(z$as.data.frame())

tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)


z$plots[[which(rownames(z$as.data.frame()) == "Chr19q13")]]#genesets$POSITIONAL$Chr19q13
z$plots[[which(rownames(z$as.data.frame()) == "Chr19p13")]]
z$plots[[which(rownames(z$as.data.frame()) == "Chr1p36")]]#genesets$POSITIONAL$Chr1p36
z$plots[[which(rownames(z$as.data.frame()) == "Chr7q22")]]#genesets$POSITIONAL$Chr7q22



z$plots[[which(rownames(z$as.data.frame()) == "Chr7p21")]] # ETV1, TWIST1, MEOX2?
#z$plots[[which(rownames(z$as.data.frame()) == "Chrxq13")]]


## CGP: chemical and genetic perturbations ----

genesets$CGP          <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CGP', clean = T)$genesets

#z <- hypeR::hypeR(genesets = genesets$CGP, signature = y.ordered, background = nrow(y), test = 'hypergeometric', absolute = F, quiet = T, plotting = T)
z <- hypeR::hypeR(genesets = genesets$CGP, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)


tibble::as_tibble(z$as.data.frame())


z$plots[[which(rownames(z$as.data.frame()) == "Gobert Oligodendrocyte Differentiation Up")]]
z$plots[[which(rownames(z$as.data.frame()) == "Kobayashi Egfr Signaling 24hr Dn")]]#genesets$CGP$`Kobayashi Egfr Signaling 24hr Dn`



tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)



# Pca2 - sterk
z$plots[[which(rownames(z$as.data.frame()) == "Nakayama Soft Tissue Tumors Pca2 Up")]] 
z$plots[[which(rownames(z$as.data.frame()) == "Creighton Endocrine Therapy Resistance 2")]]




## [x] Curated canonical: Biocarta ----
## TOO small gene sets?! ~20 genes each
# genesets$BIOCARTA     <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:BIOCARTA', clean = T)$genesets
# 
# z <- hypeR::hypeR(genesets = genesets$BIOCARTA, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
# hyp_dots(z)
# 
# 
# 
# tibble::as_tibble(z$as.data.frame())
# 
# z$plots[[which(rownames(z$as.data.frame()) == "Rb Pathway")]] #genesets$BIOCARTA$`Rb Pathway`
# 
# 
# tibble::as_tibble(z$as.data.frame()) %>%
#   filter(score < 0)


## Curated canonical: Reactome ----

genesets$REACTOME     <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:REACTOME', clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$REACTOME, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = F)
hyp_dots(z)

tibble::as_tibble(z$as.data.frame())

tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)

## Curated canonical: Kegg ----


genesets$KEGG         <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:KEGG', clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$KEGG, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = F)
hyp_dots(z)


tibble::as_tibble(z$as.data.frame())

tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)



## Curated canonical: Wiki pathways ----


genesets$WIKIPATHWAYS <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:WIKIPATHWAYS', clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$WIKIPATHWAYS, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = F)
hyp_dots(z)


tibble::as_tibble(z$as.data.frame())

tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)



# compare w/:
# "Erbb Signaling Pathway"
# "Egfegfr Signaling Pathway"
# "Egfr Tyrosine Kinase Inhibitor Resistance"
# "Erbb Signaling Pathway"

# for(n in names(genesets$WIKIPATHWAYS)) {
#   p = genesets$WIKIPATHWAYS[[n]]
#   if("EGFR" %in% p) {
#     print(n)
#   }
# }


## TFT: transcription factor targets ----

genesets$TFT <- hypeR::msigdb_gsets(species='Homo sapiens', category='C3',  subcategory='TFT:GTRD', clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$TFT, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)

# https://pubmed.ncbi.nlm.nih.gov/18852215/
# Hsd17b8
# Cebpz: CCAAT


tibble::as_tibble(z$as.data.frame())

tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)


z$plots[[which(rownames(z$as.data.frame()) == "Nfkbia Target Genes")]]#genesets$TFT$`Nfkbia Target Genes`



## Computed: Cancer modules ----


genesets$COMP <- hypeR::msigdb_gsets(species='Homo sapiens', category='C4',  subcategory='CM', clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$COMP, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)

# module 54: Cell cycle (expression cluster) - https://www.gsea-msigdb.org/gsea/msigdb/cards/MODULE_54.html - http://robotics.stanford.edu/~erans/cancer/modules/module_54


tibble::as_tibble(z$as.data.frame())

tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)


z$plots[[which(rownames(z$as.data.frame()) == "54")]]




## GO:BP (Biological Process) [g:profiler?] ----

genesets$GO_BP        <- hypeR::msigdb_gsets(species='Homo sapiens', category='C5', subcategory='GO:BP', clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$GO_BP, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)




tibble::as_tibble(z$as.data.frame())


z$plots[[which(rownames(z$as.data.frame()) == "Dna Replication")]]
z$plots[[which(rownames(z$as.data.frame()) == "Dna Dependent Dna Replication")]]#genesets$GO_BP$`Dna Dependent Dna Replication`


tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)




## oncogenic signature gene sets ----

# heel grote db? heel traag(!)

# nice that up and down are separated

genesets$ONCOSIG      <- hypeR::msigdb_gsets(species='Homo sapiens', category='C6',   clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$ONCOSIG, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)
#hyp_emap(z)

# CSR = fibroblast? - https://www.gsea-msigdb.org/gsea/msigdb/cards/CSR_LATE_UP.V1_UP.html
# PRC2 = fibroblast? - https://www.gsea-msigdb.org/gsea/msigdb/cards/PRC2_EED_UP.V1_DN.html


# doen niks: !!
# EGFR_UP.V1_DN
# EGFR_UP.V1_UP
# ERBB2_UP.V1_DN
# ERBB2_UP.V1_UP
# PTEN_DN.V1_DN
# PTEN_DN.V1_UP
# PTEN_DN.V2_DN
# PTEN_DN.V2_UP
# GLI1_UP.V1_DN
# GLI1_UP.V1_UP
# z$plots[[1]]
# z$plots[[2]]
# 
# z$plots[[3]]
# z$plots[[6]]
# 
# z$plots[[9]]
# 
# 
# z$plots[[which(rownames(z$as.data.frame()) == "Gli1 Up.v1 Dn")]]


tibble::as_tibble(z$as.data.frame())

z$plots[[which(rownames(z$as.data.frame()) == "Csr Late Up.v1 Up")]]#genesets$ONCOSIG$`Csr Late Up.v1 Up`
z$plots[[which(rownames(z$as.data.frame()) == "Prc2 Eed Up.v1 Dn")]]#genesets$ONCOSIG$`Prc2 Eed Up.v1 Dn`
z$plots[[which(rownames(z$as.data.frame()) == "Rb P130 Dn.v1 Up")]]#genesets$ONCOSIG$`Rb P130 Dn.v1 Up`



tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)


z$plots[[which(rownames(z$as.data.frame()) == "Prc2 Ezh2 Up.v1 Dn")]]#genesets$ONCOSIG$`Prc2 Ezh2 Up.v1 Dn`
z$plots[[which(rownames(z$as.data.frame()) == "Pten Dn.v2 Dn")]]#genesets$ONCOSIG$`Pten Dn.v2 Dn`



## immunologic signature gene sets  ----

genesets$IMMUSIG      <- hypeR::msigdb_gsets(species='Homo sapiens', category='C7',   clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$IMMUSIG, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)




tibble::as_tibble(z$as.data.frame())

# TODO:
#z$plots[[which(rownames(z$as.data.frame()) == "Day6 Vs Day10 Traf6ko Eff Cd8 Tcell Up")]]#genesets$IMMUSIG$`Day6 Vs Day10 Traf6ko Eff Cd8 Tcell Up`
#z$plots[[which(rownames(z$as.data.frame()) == "Gse14415 Natural Treg Vs Tconv Dn")]]#genesets$IMMUSIG$`Gse14415 Natural Treg Vs Tconv Dn`



tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)

# TODO:
#z$plots[[which(rownames(z$as.data.frame()) == "Gse45365 Wt Vs Ifnar Ko Bcell Dn")]]



## cell type signature gene sets ----

genesets$CELLTYPE      <- hypeR::msigdb_gsets(species='Homo sapiens', category='C8',   clean = T)$genesets

z <- hypeR::hypeR(genesets = genesets$CELLTYPE, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
hyp_dots(z)
hyp_emap(z)


tibble::as_tibble(z$as.data.frame())


tibble::as_tibble(z$as.data.frame()) %>%
  filter(score < 0)



###





# z <- hypeR::hypeR(genesets = genesets$BIOCARTA, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = F)
# zz <- tibble::as_tibble(z$as.data.frame())
# 

z <- hypeR::hypeR(genesets = genesets$REACTOME, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = F)
zz <- tibble::as_tibble(z$as.data.frame())
hyp_emap(z)
hyp_dots(z)



z <- hypeR::hypeR(genesets = genesets$KEGG, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
zz <- tibble::as_tibble(z$as.data.frame())


# hyp_show(z)
# hyp_dots(z)
# hyp_emap(z)
# hyp_hmap(z)

z$plots[[1]]
z$plots[[2]]
z$plots[[4]]



z <- hypeR::hypeR(genesets = genesets$GO_BP, signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = F)
zz <- tibble::as_tibble(z$as.data.frame())
hyp_dots(z)



z <- hypeR::hypeR(genesets = genesets$GO_BP[c("Chemical Synaptic Transmission Postsynaptic", names(genesets$GO_BP)[1:25])]
                    , signature = y.ordered, background = nrow(y), test = 'kstest', absolute = F, quiet = T, plotting = T)
zz <- tibble::as_tibble(z$as.data.frame())
hyp_dots(z)
z$plots[[which(zz$label == "Chemical Synaptic Transmission Postsynaptic")]]

z$plots[[3]]




# < : : : deprecated code : : : > -----

# [x] DE unpaired ~ TPC as numeric [GSAM] ----

# quite sketchy test - sensitive to outliers and many NA corrected padj's
# 
# stopifnot(colnames(gsam.gene.expression.all) == gsam.metadata.all$sid)
# 
# gsam.res.tpc.all <- DESeqDataSetFromMatrix(gsam.gene.expression.all, gsam.metadata.all, ~tpc ) %>% # + tpc
#   DESeq(parallel=T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   `colnames<-`(paste0(colnames(.),".gsam.tpc"))  %>%
#   dplyr::rename(gid = gid.gsam.tpc)
#   
# #sum(is.na ( gsam.res.tpc.all$padj.tpc ))
# #sum(is.na ( gsam.res.tpc.all$p.value.tpc ))
# 
# results.out <- results.out %>%
#   dplyr::left_join(gsam.res.tpc.all , by = c('gid' = 'gid')) 
# 


# 
# gsam.res.tpc.all %>%
#   dplyr::top_n(500, padj.tpc) %>% dplyr::pull(hugo_symbol.tpc)




# [x] zelfde maar dan voor alleen R1 ----

# 
# 
# expression <- gsam.rnaseq.expression %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   #  dplyr::filter( # bona fide glioma genes (BFGs)
#   #    (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
#   #      (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
#   #  ) %>%
#   tibble::column_to_rownames('gid') %>%
#   dplyr::select(metadata$sid)
# 
# stopifnot(colnames(expression) == metadata$sid)
# 
# 
# res.tpc.r1 <- DESeqDataSetFromMatrix(expression, metadata, ~tpc ) %>% # + tpc
#   DESeq(parallel=T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   #dplyr::arrange(padj) %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".r1"))
# 
# 
# 
# [x] zelfde maar dan voor alleen R2 ----
# 
# rm(metadata, expression)
# 
# metadata <- gsam.rna.metadata %>%
#   dplyr::filter(blacklist.pca == F) %>%
#   dplyr::filter(pat.with.IDH == F) %>%
#   dplyr::filter(resection == 'r2') %>%
#   dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
#   dplyr::filter(tumour.percentage.dna >= 25) %>%
#   dplyr::mutate(tpc = 1 - (tumour.percentage.dna / 100))
# 
# 
# expression <- gsam.rnaseq.expression %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   #  dplyr::filter( # bona fide glioma genes (BFGs)
#   #    (ensembl_id  %in% wang.glioma.intrinsic.genes$ENSG.short) |
#   #      (hugo_symbol %in% wang.glioma.intrinsic.genes$Gene_Symbol)
#   #  ) %>%
#   tibble::column_to_rownames('gid') %>%
#   dplyr::select(metadata$sid)
# 
# stopifnot(colnames(expression) == metadata$sid)
# 
# 
# res.tpc.r2 <- DESeqDataSetFromMatrix(expression, metadata, ~tpc ) %>% # + tpc
#   DESeq(parallel=T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   #dplyr::arrange(padj) %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".r2"))
# 
# 
# 
# plt <- res.tpc.r1 %>%
#   dplyr::left_join( res.tpc.r2 , by = c('gid.r1'='gid.r2'))
# 
# ggplot(plt, aes(x =log2FoldChange.r1 , y = log2FoldChange.r2, label =  hugo_symbol.r2 )) +
#   geom_text_repel(data = subset(plt, abs(log2FoldChange.r1 - log2FoldChange.r2) > 8.5 ) , size=2) +
#   geom_point(cex=0.5) 
# 



# [x] DE unpaired BFGs ----
# 
# # TODO exclude IDH mutants [check]
# # TODO x-check replicates/duplicates? [check] << taken out using 'gsam.rnaseq.expression'
# # TODO exclude low tumour percentage ? << not if BFGs are explitily used & double test?
# 
# 
# gsam.bfg.res.tpc <- DESeqDataSetFromMatrix(gsam.bfg.expression.all, gsam.metadata.all, ~tpc ) %>% # + tpc
#   DESeq(parallel = T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::arrange(padj) %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".tpc"))
# 
# 
# gsam.bfg.res.res <- DESeqDataSetFromMatrix(gsam.bfg.expression.all, gsam.metadata.all, ~resection ) %>% # + resection
#   DESeq(parallel = T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   dplyr::arrange(padj) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) )%>%
#   `colnames<-`(paste0(colnames(.),".res"))
# 
# 
# gsam.bfg.res.tpc.res <- DESeqDataSetFromMatrix(gsam.bfg.expression.all, gsam.metadata.all, ~tpc + resection ) %>% # + resection; corrected for tpc
#   DESeq(parallel = T) %>%
#   results() %>% 
#   as.data.frame() %>%
#   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
#   dplyr::arrange(padj) %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
#   dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
#   `colnames<-`(paste0(colnames(.),".tpc.res"))
# 
# 
# 
# gsam.bfg.res.combined <- gsam.bfg.res.tpc %>%
#   dplyr::full_join(gsam.bfg.res.res, by=c('gid.tpc' = 'gid.res')) %>%
#   dplyr::full_join(gsam.bfg.res.tpc.res, by=c('gid.tpc' = 'gid.tpc.res')) %>%
#   dplyr::full_join(gsam.gene.expression.all.cor, by=c('gid.tpc' = 'gid')) %>%
#   dplyr::filter(!is.na(padj.res))%>%
#   dplyr::filter(!is.na(padj.tpc))%>%
#   dplyr::filter(!is.na(padj.tpc.res))
# 
# 
# 
# 
# p1 <- ggplot(gsam.bfg.res.combined, aes(x=log2FoldChange.res ,
#                                         y=statistic,
#                                         col=significant.res,
#                                         label=hugo_symbol.tpc.res ) ) + 
#   geom_point(data=subset(gsam.bfg.res.combined, significant.res == F ),pch=19,cex=0.05) +
#   geom_point(data=subset(gsam.bfg.res.combined, significant.res == T ),pch=19,cex=0.5) +
#   scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
#   geom_text_repel(data = subset(gsam.bfg.res.combined, significant.res == T & abs(log2FoldChange.res) > 1), size = 2.4 )  +
#   labs(x = "log2FC R1 vs. R2 (unpaired)",
#        y="Correlation t-statistic with tumour percentage",
#        col="Difference significant (R1 ~ R2)"
#   ) +
#   youri_gg_theme
# 
# p2 <- ggplot(gsam.bfg.res.combined, aes(x=log2FoldChange.tpc.res ,
#                                         y= statistic,
#                                         col=significant.tpc.res,
#                                         label=hugo_symbol.tpc.res ) ) + 
#   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == F ),pch=19,cex=0.05) +
#   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == T ),pch=19,cex=0.5) +
#   scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
#   geom_text_repel(data = subset(gsam.bfg.res.combined, significant.tpc.res == T & abs(log2FoldChange.tpc.res) > 1), size = 2.4 )  +
#   labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage",
#        col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
#   youri_gg_theme
# 
# p1 + p2
# 
# # 
# # 
# # p3 <- ggplot(gsam.bfg.res.combined, aes(x=log2FoldChange.tpc.res, y= log2FoldChange.tpc, col=significant.tpc.res, label=hugo_symbol.tpc.res ) ) + 
# #   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == F ),pch=19,cex=0.05) +
# #   geom_point(data=subset(gsam.bfg.res.combined, significant.tpc.res == T ),pch=19,cex=0.3) +
# #   scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
# #   #geom_text_repel(data = subset(gsam.bfg.res.combined, significant.tpc.res == T & abs(log2FoldChange.tpc.res) > 1), size = 2.4 )  +
# #   youri_gg_theme
# # 


# [x] DE paired all [G-SAM] ----

# https://support.bioconductor.org/p/59481/
# fitType='local' or 'mean' 

if(!file.exists('tmp/gsam.gene.res.res.paired.Rds')) {
  
  stopifnot(F)
  
  
  # gsam.gene.res.res.paired <- DESeqDataSetFromMatrix(gsam.gene.expression.all.paired %>%
  #                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
  #                                             , gsam.metadata.all.paired , ~pid + resection ) %>% # + resection
  #   DESeq(parallel = T, fitType="local") %>%
  #   results() %>% 
  #   as.data.frame() %>%
  #   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 )  %>%
  #   dplyr::arrange(padj) %>% 
  #   tibble::rownames_to_column('gid') %>%
  #   `colnames<-`(paste0(colnames(.),".res.paired")) %>%
  #   dplyr::rename(gid = gid.res.paired)
  # 
  # saveRDS(gsam.gene.res.res.paired, file = "tmp/gsam.gene.res.res.paired.Rds")
  
  
  
  stopifnot(gsam.metadata.all.paired$sid == colnames(gsam.gene.expression.all.paired))
  # gsam.gene.res.tpc.res.paired <- DESeqDataSetFromMatrix(gsam.gene.expression.all.paired %>%
  #                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
  #                                             , gsam.metadata.all.paired , ~tpc + pid + resection ) %>% # + resection
  #   DESeq(parallel = T, fitType="local") %>%
  #   results() %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  #   dplyr::arrange(padj) %>%
  #   tibble::rownames_to_column('gid') %>%
  #   `colnames<-`(paste0(colnames(.),".tpc.res.paired")) %>%
  #   dplyr::rename(gid = gid.tpc.res.paired)
  # 
  # saveRDS(gsam.gene.res.tpc.res.paired, file = "tmp/gsam.gene.res.tpc.res.paired.Rds")
  
  stopifnot(gsam.metadata.all$sid == colnames(gsam.gene.expression.all))
  test.a <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , gsam.metadata.all , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.a")) %>%
    dplyr::rename(gid = gid.test.a)
  
  test.a.ntpc <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                          dplyr::filter(rowSums(.) > ncol(.) * 3)
                                        , gsam.metadata.all , ~ resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.a.ntpc")) %>%
    dplyr::rename(gid = gid.test.a.ntpc)
  
  
  test.b <- DESeqDataSetFromMatrix(gsam.gene.expression.all.paired %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , gsam.metadata.all.paired , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.b")) %>%
    dplyr::rename(gid = gid.test.b)
  
  
  m <- gsam.metadata.all %>%
    dplyr::filter() %>% 
    dplyr::filter(sid %in% c('FAF2', 'AIA2', 'HAF2-new') == F )  
  e <- gsam.gene.expression.all %>% dplyr::select(m$sid)
  stopifnot(m$sid == colnames(e))
  test.c <- DESeqDataSetFromMatrix(e %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , m , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.c")) %>%
    dplyr::rename(gid = gid.test.c)
  
  
  m <- gsam.metadata.all %>%
    dplyr::filter() %>% 
    dplyr::filter(sid %in% c('FAF2', 'AIA2', 'HAF2-new', "AOA2","BAE2","AAX1","EBB2") == F )  
  e <- gsam.gene.expression.all %>% dplyr::select(m$sid)
  stopifnot(m$sid == colnames(e))
  test.d <- DESeqDataSetFromMatrix(e %>%
                                     dplyr::filter(rowSums(.) > ncol(.) * 3)
                                   , m , ~tpc + resection ) %>% # + resection
    DESeq(parallel = T, fitType="local") %>%
    results() %>%
    as.data.frame() %>%
    dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
    dplyr::arrange(padj) %>%
    tibble::rownames_to_column('gid') %>%
    `colnames<-`(paste0(colnames(.),".test.d")) %>%
    dplyr::rename(gid = gid.test.d)
  
  
  
  
  
  subset(gsam.gene.res.tpc.res.paired, gid %in% a)$log2FoldChange.tpc.res.paired %>% length()
  subset(test.a, gid %in% a)$log2FoldChange.test.a %>% length()
  subset(test.a.ntpc, gid %in% a)$log2FoldChange.test.a.ntpc %>% length()
  subset(test.b, gid %in% a)$log2FoldChange.test.b %>% length()
  subset(test.c, gid %in% a)$log2FoldChange.test.c %>% length()
  subset(test.d, gid %in% a)$log2FoldChange.test.d %>% length()
  #%>% length()
  
  plt <- gsam.gene.res.tpc.res.paired %>%
    dplyr::left_join(test.a , by  = c('gid'='gid')) %>%
    dplyr::mutate(col = ifelse(lfcSE.test.a  > 0.3 , "outlier", "regular") )
  ggplot(plt , aes(x=stat.tpc.res.paired,
                   y=stat.test.a , 
                   col=col )) +
    geom_point(data = subset(plt, col == "regular"), pch=19,cex=0.5) +
    geom_point(data = subset(plt, col == "outlier"), pch=19,cex=0.85)
  
  
  aa = gsam.gene.res.tpc.res%>% dplyr::filter( padj.tpc.res < 0.2 ) %>% pull(lfcSE.tpc.res)
  ab = gsam.gene.res.tpc.res %>% dplyr::filter( gid %in% a ) %>% pull(lfcSE.tpc.res)
  plot(density(aa))
  abline(v=ab)
  abline(v=0.3, col="red")
  
  
  # a = 
  # [1] "ENSG00000162571.13_5|TTLL10|chr1:1109264-1133315(+)"       "ENSG00000137975.8_3|CLCA2|chr1:86889854-86922236(+)"      
  # [3] "ENSG00000143196.5_4|DPT|chr1:168664706-168698444(-)"       "ENSG00000159173.19_4|TNNI1|chr1:201372896-201398994(-)"   
  # [5] "ENSG00000122180.5_3|MYOG|chr1:203052257-203055140(-)"      "ENSG00000168530.16_4|MYL1|chr2:211154874-211179898(-)"    
  # [7] "ENSG00000287059.1_1|AC090004.2|chr3:14080147-14088681(-)"  "ENSG00000185290.4_4|NUPR2|chr7:56182374-56184110(-)"      
  # [9] "ENSG00000285670.1_2|AC006970.3|chr7:56282672-56288440(+)"  "ENSG00000147573.17_5|TRIM55|chr8:67039131-67087720(+)"    
  # [11] "ENSG00000215182.6|MUC5AC|chr11:1151580-1222364(+)"         "ENSG00000129152.4_3|MYOD1|chr11:17741118-17743683(+)"     
  # [13] "ENSG00000230657.6_3|PRB4|chr12:11460017-11463369(-)"       "ENSG00000251655.6_3|PRB1|chr12:11504757-11548500(-)"      
  # [15] "ENSG00000121335.12_7|PRB2|chr12:11544474-11653975(-)"      "ENSG00000279134.1_5|AC090643.1|chr12:58488697-58491796(-)"
  # [17] "ENSG00000283361.2_7|CFAP97D2|chr13:114920166-114988559(+)" "ENSG00000226777.7_5|FAM30A|chr14:106383838-106398502(+)"  
  # [19] "ENSG00000260496.3_6|AC009041.1|chr16:1041151-1050926(-)"   "ENSG00000262152.7_5|LINC00514|chr16:3038257-3052017(+)"   
  # [21] "ENSG00000260034.1_6|LCMT1-AS2|chr16:25151898-25160353(-)"  "ENSG00000133020.4_3|MYH8|chr17:10293639-10325267(-)"      
  # [23] "ENSG00000128422.17_6|KRT17|chr17:39775694-39781094(-)"     "ENSG00000175894.18_7|TSPEAR|chr21:45917776-46131487(-)"   
  # [25] "ENSG00000187268.12_5|FAM9C|chrX:13053736-13062801(-)"     
  
  
} else {
  
  gsam.gene.res.res.paired <- readRDS("tmp/gsam.gene.res.res.paired.Rds")
  gsam.gene.res.tpc.res.paired <- readRDS("tmp/gsam.gene.res.tpc.res.paired.Rds")
  
}


results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res.paired, by = c('gid' = 'gid')) %>%
  dplyr::left_join(gsam.gene.res.tpc.res.paired, by = c('gid' = 'gid'))






plt <- results.out  %>%
  dplyr::mutate(show.label = !is.na(log2FoldChange.tpc.res.paired) & 
                  !is.na(log2FoldChange.tpc.res) & 
                  abs(log2FoldChange.tpc.res.paired - log2FoldChange.tpc.res) > 2
  )


#%>%
#  dplyr::filter(!is.na(log2FoldChange.res) & !is.na(statistic.cor.tpc)) %>%
#dplyr::mutate(is.limited.res = as.character(log2FoldChange.res > 3)) %>% # change pch to something that is limited
#dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3 , log2FoldChange.res)) %>%
#dplyr::mutate(is.limited.tpc.res = as.character(log2FoldChange.tpc.res > 3)) %>% # change pch to something that is limited
#dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3 , log2FoldChange.tpc.res))


# MYL1 & OPALIN & MYOD1 ??? 


p1 <- ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc , label=hugo_symbol
                      #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                      #shape = is.limited.res ,
                      #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  #geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  #scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  #scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  geom_text_repel(data =subset(plt, show.label == T),col="red") +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage") +
  xlim(-2.5,4)


p2 <- ggplot(plt, aes(x = log2FoldChange.tpc.res.paired , y = statistic.cor.tpc , label=hugo_symbol
                      #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                      #shape = is.limited.res ,
                      #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  #geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  #scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  #scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  geom_text_repel(data =subset(plt, show.label == T),col='red') +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage") +
  xlim(-2.5,4)


p1 + p2 




results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res.paired, by = c('gid' = 'gid')) %>%
  dplyr::left_join(gsam.gene.res.tpc.res.paired, by = c('gid' = 'gid'))
plt <- results.out  
p.corr <- ggplot(plt, aes(x = log2FoldChange.tpc.res.paired , y = log2FoldChange.tpc.res , label=hugo_symbol
                          #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                          #shape = is.limited.res ,
                          #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(pch=19,cex=0.05) +
  geom_point(pch=19,cex=0.25, data =subset(plt, gid %in% a),col="red") +
  #geom_smooth(data = subset(plt, padj.tpc.res > 0.05),method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  #scale_shape_manual(values = c('TRUE'=4, 'FALSE' = 19)    ) +
  #scale_size_manual(values = c('TRUE'=0.75, 'FALSE' = 0.05)    ) +
  geom_text_repel(data =subset(plt, gid %in% a),col="red") +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage") + 
  xlim(-2,4) +
  ylim(-2,4) 

p.uncorr + p.corr
subset(plt, gid %in% a) $log2FoldChange.tpc.res
subset(plt, gid %in% a) $log2FoldChange.tpc.res.paired

# .paired
# > subset(plt, gid %in% a) $log2FoldChange.tpc.res
# [1]  2.616752  2.705314  2.839108        NA        NA        NA  4.413780 -2.573126  2.801220        NA  5.968508        NA  2.413724  8.257744  8.114305
# [16]  3.549533  3.056017  3.897042  4.550233 -2.336730  2.838127        NA  3.528260  2.586812        NA





lm(  log2FoldChange.tpc.res.paired ~ log2FoldChange.tpc.res  , data=plt)
lm(  log2FoldChange.tpc.res ~ log2FoldChange.tpc.res.paired   , data=plt)

c <- plt %>%
  dplyr::filter(!is.na(log2FoldChange.tpc.res) & !is.na(log2FoldChange.tpc.res.paired))
#& abs(log2FoldChange.tpc.res) < 2 & abs(log2FoldChange.tpc.res.paired) < 2 )
cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired  )
cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired, method = "spearman"  )

cor(  c$stat.tpc.res , c$stat.tpc.res.paired  )
cor(  c$stat.tpc.res , c$stat.tpc.res.paired, method = "spearman"  )



cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired  )
cor(  c$stat.tpc.res , c$stat.tpc.res.paired  )
c <- c %>% dplyr::filter(!is.na(padj.tpc.res) & !is.na (padj.tpc.res.paired  ) )
cor(  c$padj.tpc.res , c$padj.tpc.res.paired , method="" )
cor(  c$log2FoldChange.tpc.res , c$log2FoldChange.tpc.res.paired, method = "spearman"  )
#cor(  c$stat.tpc.res , c$stat.tpc.res.paired, method = "spearman"  )





#ggplot(plt, aes(x = stat.tpc.res.paired , y = stat.tpc.res , label=hugo_symbol
ggplot(plt, aes(x = stat.tpc.res.paired , y = stat.tpc.res , label=hugo_symbol
                #ggplot(plt, aes(x = log2FoldChange.res.paired , y = statistic.cor.tpc, 
                #shape = is.limited.res ,
                #size = is.limited.res 
) ) +
  #p1 <- ggplot(plt, aes(x = stat.res , y = statistic.cor.tpc  ) ) +
  geom_point(data=subset(plt, show.label == F), pch=19,cex=0.05) +
  geom_point(data=subset(plt, show.label == T), pch=19,cex=0.25,col="red") +
  geom_text_repel(data =subset(plt, show.label == T),col="red") +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (paired)",
       y="Correlation t-statistic with tumour percentage")
#+
#xlim(-3,10)+
#ylim(-3,10)



a = plt %>% dplyr::filter(show.label == T ) %>% pull(gid)


b = as.data.frame(gsam.gene.expression.all.vst) %>% 
  dplyr::filter(rownames(.) %in% a)

c = prcomp(t(b))
screeplot(c)


d = as.data.frame(c$x) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(in.paired.data = gid  %in% colnames(gsam.gene.expression.all.paired) ) %>%
  dplyr::mutate(res = gsub("^...(.).*$","R\\1",gid))


ggplot(d, aes(x=PC1, y=PC3, label=gid, col=in.paired.data) ) + 
  geom_point() + 
  geom_text_repel(data = subset(d, in.paired.data == F))



#AOA2
#BAE2
#AAX1
#EBB2

b = data.frame(
  gen = as.numeric(gsam.gene.expression.all.vst %>% as.data.frame() %>% dplyr::filter(rownames(.) == a[9]) ),
  sid=colnames(gsam.gene.expression.all.vst) ) %>%
  dplyr::mutate(in.paired.data = sid  %in% colnames(gsam.gene.expression.all.paired) ) %>%
  dplyr::mutate(res = gsub("^...(.).*$","R\\1",sid)) %>%
  dplyr::left_join(gsam.metadata.all %>% dplyr::select(c('sid', 'tpc')) , by = c('sid' = 'sid')) %>%
  dplyr::mutate(pid = as.factor(gsub("^(...).*$","\\1",sid)))


# AOA2, BAE2
ggplot(b, aes(x = tpc, y = gen, col=sid %in% c("AOA2","BAE2"), label=sid, shape = in.paired.data, group = pid)) +
  #geom_line(lwd=0.25,size=0.4,col="gray60") +
  geom_point(size=2)  +
  geom_text_repel(data = subset(b, gen > 5))



b1 <- b %>% dplyr::filter(in.paired.data) %>% dplyr::filter(res == 'R1')
b2 <- b %>% dplyr::filter(in.paired.data) %>% dplyr::filter(res == 'R2')
colnames(b1) <- paste0(colnames(b1), ".R1")
colnames(b2) <- paste0(colnames(b2), ".R2")
c = dplyr::left_join(b1, b2, by=c('pid.R1' = 'pid.R2')) %>%
  dplyr::mutate(gen.avg = (gen.R2 + gen.R1) / 2) %>%
  dplyr::mutate(dist = gen.R2 - gen.R1) %>%
  dplyr::mutate(dist.tpc = tpc.R2 - tpc.R1) %>%
  dplyr::arrange(dist) %>%
  dplyr::mutate(col = ifelse(sid.R1 %in% c('FAF1', 'AIA1','HAF1-new') , 1 , 2) )

## strong diff FAF2 / AIA2 / HAF2-new ?

plot(c$gen.R1 - c$gen.avg, y = 1:nrow(c), pch=19,col="red")
points(c$gen.R2 - c$gen.avg, y = 1:nrow(c), pch=19,col="blue")


plot(c$dist.tpc , c$dist, col=c$col, pch=19) # deze zou het meeste moeten zeggen
d <- b %>% dplyr::filter(in.paired.data == F)
points(d$tpc - median(c(c$tpc.R1 , c$tpc.R2)), d$gen - median(c(c$gen.R1 , c$gen.R2)) ,col=4, pch=19)
abline(h=0, col="gray30",lwd=0.3)
#plot((c$tpc.R1 + c$tpc.R2) / 2 , c$dist)
#plot(c$tpc.R1 , c$dist)
#plot(c$tpc.R2 , c$dist)


plot(b$tpc, b$gen, col=ifelse(b$res == "R1",2,4), pch=19)


# < : / : deprecated code : : : > -----

# plots ----




## 1 :: gene-set corr plot ----

# to find whether the scRNA-seq clusters are (largely) identical to verhaak etc.


tmp <- rbind(
  results.out %>%
    dplyr::filter(!is.na(GliTS.reduxsubtype.marker) & GliTS.reduxsubtype.marker == "GliTS-CL") %>%
    dplyr::mutate(label = paste0("Subtype: ",GliTS.reduxsubtype.marker)) %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(!is.na(GliTS.reduxsubtype.marker) & GliTS.reduxsubtype.marker == "GliTS-MES") %>%
    dplyr::mutate(label = paste0("Subtype: ",GliTS.reduxsubtype.marker)) %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(!is.na(GliTS.reduxsubtype.marker) & GliTS.reduxsubtype.marker == "GliTS-PN") %>%
    dplyr::mutate(label = paste0("Subtype: ",GliTS.reduxsubtype.marker)) %>%
    head(n=14)
  ,
  
  results.out %>%
    dplyr::filter(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Complete/Immune response") %>%
    dplyr::mutate(label = "Patel: Compl/Imm.res") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Cell cycle") %>%
    dplyr::mutate(label = paste0("Patel: ",patel.scRNAseq.cluster))%>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Hypoxia") %>%
    dplyr::mutate(label = paste0("Patel: ",patel.scRNAseq.cluster))%>%
    head(n=14)
  ,
  
  results.out %>%
    dplyr::filter(neftel.meta.module.MES1 == T) %>%
    dplyr::mutate(label = "Neftel: MES1") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.MES2 == T) %>%
    dplyr::mutate(label = "Neftel: MES2") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.AC == T) %>%
    dplyr::mutate(label = "Neftel: AC") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.OPC == T) %>%
    dplyr::mutate(label = "Neftel: OPC") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.NPC1 == T) %>%
    dplyr::mutate(label = "Neftel: NPC1") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.NPC2 == T) %>%
    dplyr::mutate(label = "Neftel: NPC2") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.G1.S == T) %>%
    dplyr::mutate(label = "Neftel: G1.S") %>%
    head(n=14)
  ,
  results.out %>%
    dplyr::filter(results.out$neftel.meta.module.G2.M == T) %>%
    dplyr::mutate(label = "Neftel: G2.M") %>%
    head(n=14),
  
  results.out %>%
    dplyr::filter(grepl("astrocyte",primary.marker.genes)) %>%
    dplyr::mutate(label = "McKenzy: ast") %>%
    head(n=14),
  results.out %>%
    dplyr::filter(grepl("neuron",primary.marker.genes)) %>%
    dplyr::mutate(label = "McKenzy: neu") %>%
    head(n=14),
  results.out %>%
    dplyr::filter(grepl("oligodendrocyte",primary.marker.genes)) %>%
    dplyr::mutate(label = "McKenzy: oli") %>%
    head(n=14),
  results.out %>%
    dplyr::filter(grepl("microglia",primary.marker.genes)) %>%
    dplyr::mutate(label = "McKenzy: tam") %>%
    head(n=14)
    
) %>%
  dplyr::filter(!is.na(gid)) %>%
  dplyr::select(c(gid, label, hugo_symbol)) %>%
  dplyr::left_join(
    gsam.gene.expression.all.vst %>%
      as.data.frame %>%
      tibble::rownames_to_column('gid'),
    by=c('gid'='gid')
  ) %>%
  dplyr::mutate(
    label = factor(label, levels = 
                   c(  "Patel: Cell cycle" ,   "Neftel: G1.S" , "Neftel: G2.M", 
                       "Subtype: GliTS-MES",   "Neftel: MES1" ,  "Neftel: MES2"  , "Patel: Hypoxia" , "Patel: Compl/Imm.res",
                       "Neftel: AC" , "Subtype: GliTS-CL" ,   "Neftel: OPC" ,
                       "Neftel: NPC1" ,  "Neftel: NPC2","Subtype: GliTS-PN" )
                     )) %>%
  dplyr::arrange(label) %>%
  dplyr::mutate(label=as.character(label))



#plt <- tmp %>% dplyr::select(-c(gid,label,hugo_symbol)) %>% as.matrix %>% t() %>% cor() %>% `rownames<-`(tmp %>%   dplyr::mutate(label = ifelse(duplicated(label), "", label)) %>% dplyr::pull(label))
#corrplot::corrplot(plt,tl.cex=0.7)
 




plt <- tmp %>%
  dplyr::filter(!duplicated(hugo_symbol)) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::select(-c(label, gid)) %>%
  as.matrix %>% 
  t() %>%
  cor()

#h <- hclust(as.dist(1 - plt))
h <- hclust( as.dist(1 - cor(plt)) ) # Geniale manier om te clusteren!!!
o <- h$labels[h$order]

plt <- plt %>%
  as.data.frame %>%
  dplyr::select(o) %>%
  t() %>%
  as.data.frame %>%
  dplyr::select(o) %>%
  t() %>%
  as.matrix

corrplot::corrplot(plt)
  




plt.expanded <- data.frame()
for(i in 1:nrow(plt)) {
  x <- rownames(plt)[i]
  
  for(j in 1:ncol(plt)) {
    y <- colnames(plt)[j]
    
    plt.expanded <- dplyr::bind_rows(plt.expanded, c(x = x, x.order = i,
                                                     y = y, y.order = nrow(plt) - j,
                                                     order = (i - 1) * nrow(plt) + j,
                                                     value = plt[i,j],
                                                     type = 'cor'))
  }
}

plt.expanded <- plt.expanded %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::mutate(x.order = as.numeric(x.order)) %>%
  dplyr::mutate(y.order = as.numeric(y.order)) %>%
  dplyr::mutate(order = as.numeric(order))


# add metadata
tmp.meta <- results.out %>% dplyr::filter(!is.na(GliTS.reduxsubtype.marker) & GliTS.reduxsubtype.marker == "GliTS-CL") %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 13) %>%
  dplyr::mutate(y = "GliTS-CL")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(!is.na(GliTS.reduxsubtype.marker) & GliTS.reduxsubtype.marker == "GliTS-PN") %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 6) %>%
  dplyr::mutate(y = "GliTS-PN")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(!is.na(GliTS.reduxsubtype.marker) & GliTS.reduxsubtype.marker == "GliTS-MES") %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 5) %>%
  dplyr::mutate(y = "GliTS-MES")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Complete/Immune response") %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 4) %>%
  dplyr::mutate(y = "Patel: Comp/Imm.res")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Cell cycle") %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 10) %>%
  dplyr::mutate(y = "Patel: Cell Cycle")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Hypoxia") %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 1) %>%
  dplyr::mutate(y = "Patel: Hypoxia")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.MES1 == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 3) %>%
  dplyr::mutate(y = "Neftel: MES1")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.MES2 == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 2) %>%
  dplyr::mutate(y = "Neftel: MES2")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.AC == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 14) %>%
  dplyr::mutate(y = "Neftel: AC")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.OPC == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 7) %>%
  dplyr::mutate(y = "Neftel: OPC")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)


# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.NPC1 == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 9) %>%
  dplyr::mutate(y = "Neftel: NPC1")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.NPC2 == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 8) %>%
  dplyr::mutate(y = "Neftel: NPC2")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.G1.S == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 12) %>%
  dplyr::mutate(y = "Neftel: G1+S")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(neftel.meta.module.G2.M == T) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 11) %>%
  dplyr::mutate(y = "Neftel: G2+M")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(grepl("astrocyte",primary.marker.genes)) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 15) %>%
  dplyr::mutate(y = "McKenzy: ast")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(grepl("neuron",primary.marker.genes)) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 16) %>%
  dplyr::mutate(y = "McKenzy: neu")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(grepl("oligodendrocyte",primary.marker.genes)) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 17) %>%
  dplyr::mutate(y = "McKenzy: oli")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)

# add metadata
tmp.meta <- results.out %>% dplyr::filter(grepl("microglia",primary.marker.genes)) %>% dplyr::pull(hugo_symbol)
tmp.meta <- plt.expanded %>%
  dplyr::filter(x %in% tmp.meta & x == y) %>%
  dplyr::mutate(y.order = nrow(plt ) + 18) %>%
  dplyr::mutate(y = "McKenzy: tam")  %>%
  dplyr::mutate(type = y) %>% 
  dplyr::mutate(value=0)
dim(tmp.meta)
plt.expanded <- dplyr::bind_rows(plt.expanded, tmp.meta)



plt.expanded <- plt.expanded %>%
  dplyr::mutate(mckenzy = grepl("McKenzy", type))




col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                   "#4393C3", "#2166AC", "#053061"))



p1 <- ggplot(plt.expanded %>% dplyr::filter(type != 'cor'), aes(x=reorder(x, x.order), y=reorder(y, y.order),  fill=type)) +
  geom_tile(col="black") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "", col=NULL, fill=NULL) +
  guides(fill=FALSE) +
  theme(axis.text = element_text(),
        axis.line = element_line(colour="black"),
        text = element_text(size=15)) +
  labs(y = "") + 
  facet_grid(rows = vars(mckenzy), scales = "free", space = "free") +
  theme(strip.text.y = element_blank())


p2 <- ggplot(plt.expanded %>% dplyr::filter(type == 'cor'), aes(x=reorder(x, x.order), y=reorder(y, y.order),  col=type, fill=value)) +
  geom_tile(col = "gray") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = "") +
  scale_fill_gradientn(
    colours = col2(200),
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
    limits = c(-1,1)
  ) +
  theme(legend.position = 'bottom',  axis.text = element_text(),
        axis.line = element_line(colour="black"),
        text = element_text(size=5)) +
  labs(y = "")

p1 / p2 +  plot_layout(heights = c(1,  3))


#ggsave("/tmp/cor_clusters.pdf", width = 11, height = 11)
ggsave("/tmp/cor_clusters.png", width = 11 , height = 11 )





#scale_y_discrete(labels = NULL, breaks = NULL) +









## 2 :: [E:S1] :: TPC violin ----





plt <- rbind(
  gsam.metadata.all %>%
    dplyr::mutate(resection = ifelse(resection == 'r1', "Primary Res.", "Recurrent Res.")) %>%
    dplyr::select(c(sid, resection, tumour.percentage.dna)) %>%
    dplyr::mutate(dataset = "(A) G-SAM")
  ,
  glass.metadata.all %>%
    dplyr::mutate(resection = ifelse(resection == 'TP', "Primary Res.", "Recurrent Res.")) %>%
    dplyr::mutate(tumour.percentage.dna = tumour.percentage.dna.imputed.rf) %>%
    dplyr::select(c(sid, resection, tumour.percentage.dna)) %>%
    dplyr::mutate(dataset = "(B) GLASS")) %>%
  dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
  dplyr::left_join(
    gsam.patient.metadata %>% dplyr::select(studyID, HM),
    by = c('pid' = 'studyID')
  ) 

# HM does not clearly change, seems actually more stable with small N - exclude from fig..

ggplot(plt, aes(x = resection, y = tumour.percentage.dna, fill=resection)) +
  facet_grid(cols=vars(dataset), scales = "free", space="free_y") + 
  geom_violin(draw_quantiles = c(0.6), col="black", alpha=0.2) +
  geom_jitter( position=position_jitter(0.2), size=2.5, pch=21, col="black") +
  ylim(0, 100) +
  labs(x = NULL, col = NA, y = "Tumor cell percentage" ) +
  geom_signif(
    comparisons = list(c("Primary Res." ,  "Recurrent Res.")),
    test="wilcox.test",
    col="black"
  ) + 
  scale_fill_manual(values =  resection_colors ) +
  job_gg_theme




# @TODO include Wilcox pvalues
ggsave('output/figures/figure-tcp_estimate.pdf', width = 5.7 * 1.4, height = 4.3 * 1.1)





## 3 :: TopCor genes ----


sel <- rbind(
  plt %>%
    dplyr::arrange(-statistic.gsam.cor.tpc) %>%
    head(100),
  plt %>%
    dplyr::arrange(statistic.gsam.cor.tpc) %>%
    head(100)
  )


# sel <- plt %>% 
#   dplyr::filter(hugo_symbol %in% c(
#     "ACTN4", "AEBP1", "AHR", "AIRE", "ALX1", "ALX3", "ALX4", "ANHX", "APP", "AR", "ARGFX", "ARNT", "ARNT2", "ARNTL", "ARNTL2", "ARX", "ASCL1", "ASCL2", "ASCL3", "ASCL4", "ASCL5", "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6", "ATF6B", "ATF7", "ATOH1", "ATOH7", "ATOH8", "ATXN3L", "BACH1", "BACH2", "BARHL1", "BARHL2", "BARX1", "BARX2", "BATF", "BATF2", "BATF3", "BBX", "BCL11A", "BCL11B", "BCL6", "BCL6B", "BCOR", "BHLHA15", "BHLHA9", "BHLHE22", "BHLHE23", "BHLHE40", "BHLHE41", "BMI1", "BORCS8-MEF2B", "BPTF", "BSX", "CALCOCO1", "CARF", "CASZ1", "CC2D1A", "CC2D1B", "CCAR1", "CDC5L", "CDK9", "CDX1", "CDX2", "CDX4", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CHD2", "CHD7", "CIART", "CIC", "CLOCK", "CPHXL", "CREB1", "CREB3", "CREB3L1", "CREB3L2", "CREB3L3", "CREB3L4", "CREB5", "CREBRF", "CREM", "CRX", "CRY1", "CTCF", "CTCFL", "CUX1", "CUX2", "DACH1", "DACH2", "DBP", "DDIT3", "DEAF1", "DHX36", "DHX9", "DLX1", "DLX2", "DLX3", "DLX4", "DLX5", "DLX6", "DMBX1", "DMRT1", "DMRT2", "DMRT3", "DMRTA1", "DMRTA2", "DMRTB1", "DMRTC2", "DMTF1", "DNMT3A", "DPRX", "DRGX", "DUX4", "DUXA", "DUXB", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8", "E4F1", "EBF1", "EBF2", "EBF3", "EBF4", "EGR1", "EGR2", "EGR3", "EGR4", "EHF", "EHMT2", "ELF1", "ELF3", "ELF4", "ELK1", "ELK3", "ELK4", "EMX1", "EMX2", "EN1", "EN2", "ENO1", "EOMES", "EPAS1", "ERG", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "ESX1", "ETS1", "ETS2", "ETV1", "ETV2", "ETV3", "ETV4", "ETV5", "ETV6", "ETV7", "EVX1", "EVX2", "EZH2", "FERD3L", "FEZF1", "FEZF2", "FIGLA", "FLI1", "FOS", "FOSB", "FOSL1", "FOSL2", "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXB2", "FOXC1", "FOXC2", "FOXD1", "FOXD2", "FOXD3", "FOXD4", "FOXD4L1", "FOXD4L3", "FOXD4L4", "FOXD4L5", "FOXD4L6", "FOXE1", "FOXE3", "FOXF1", "FOXF2", "FOXI1", "FOXI2", "FOXI3", "FOXJ1", "FOXJ2", "FOXJ3", "FOXK1", "FOXK2", "FOXL1", "FOXL2", "FOXL3", "FOXM1", "FOXO1", "FOXO3", "FOXO4", "FOXO6", "FOXP1", "FOXP2", "FOXP3", "FOXP4", "FOXQ1", "FOXS1", "GABPA", "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6", "GBX1", "GBX2", "GCM1",
#     "GCM2", "GFI1", "GFI1B", "GLI1", "GLI2", "GLI3", "GLI4", "GLIS1", "GLIS2", "GLIS3", "GMEB1", "GMEB2", "GRHL1", "GRHL2", "GRHL3", "GSC", "GSC2", "GSX1", "GTF2A1", "GTF2B", "GTF2IRD1", "GZF1", "H2AZ1", "H3-3A", "H3-3B", "HAND1", "HAND2", "HBP1", "HDAC1", "HDAC4", "HDAC5", "HDAC6", "HDX", "HELT", "HES1", "HES2", "HES3", "HES4", "HES5", "HES6", "HES7", "HESX1", "HEY1", "HEY2", "HEYL", "HHEX", "HIC1", "HIC2", "HIF1A", "HIF3A", "HINFP", "HIVEP1", "HIVEP2", "HIVEP3", "HLF", "HMGA2", "HMX1", "HMX2", "HMX3", "HNF1A", "HNF1B", "HNF4A", "HNF4G", "HNRNPU", "HOXA10", "HOXA11", "HOXA13", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXB1", "HOXB13", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9", "HOXC10", "HOXC11", "HOXC13", "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXD1", "HOXD10", "HOXD11", "HOXD13", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HSF1", "HSF2", "HSF4", "HSF5", "HSFX1", "HSFX2", "HSFX3", "HSFX4", "HSFY1", "HSFY2", "IFI16", "IKZF1", "IKZF2", "IKZF3", "IKZF4", "IKZF5", "initial_alias,", "INSM1", "INSM2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "IRX1", "IRX2", "IRX3", "IRX4", "IRX5", "IRX6", "ISL1", "ISX", "JDP2", "JUN", "JUNB",
#     "JUND", "KAT2B", "KCNIP3", "KDM2B", "KDM6A", "KDM6B", "KLF1", "KLF10", "KLF11", "KLF12", "KLF13", "KLF14", "KLF15", "KLF16", "KLF17", "KLF18", "KLF2", "KLF3", "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9", "LEF1", "LEUTX", "LHX1", "LHX2", "LHX3", "LHX4", "LHX5", "LHX6", "LHX8", "LHX9", "LITAF", "LMX1A", "LMX1B", "LRRFIP1", "LYL1", "MACROH2A1", "MACROH2A2", "MAF", "MAFA", "MAFB", "MAFF", "MAFG", "MAFK", "MAX", "MAZ", "MECOM", "MED1", "MED12", "MED8", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "MEIS1", "MEIS2", "MEIS3", "MEOX1", "MEOX2", "MESP1", "MESP2", "MGA", "MITF", "MIXL1", "MKX", "MLX", "MLXIP", "MLXIPL", "MNT", "MSC", "MSGN1", "MSX1", "MSX2", "MTA1", "MTF1", "MTF2", "MUC1", "MXD1", "MXD3", "MXD4", "MXI1", "MYB", "MYBBP1A", "MYBL1", "MYBL2", "MYC", "MYCL", "MYCN", "MYF5", "MYF6", "MYNN", "MYOD1", "MYOG", "MYPOP", "MYT1", "MYT1L", "MZF1", "NACC1", "NACC2", "nan", "nan", "nan", "nan", "nan", "nan", "nan", "NANOG", "NANOGP8", "NCOA1", "NCOA2", "NDN", "NEUROD1", "NEUROD2", "NEUROD4", "NEUROD6", "NEUROG1", "NEUROG2", "NEUROG3", "NFAT5", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "NFE2", "NFE2L1", "NFE2L2", "NFE2L3", "NFIA", "NFIB", "NFIC", "NFIL3", "NFIX", "NFKB1", "NFKB2", "NFX1", "NFXL1", "NFYA", "NFYB", "NFYC", "NHLH1", "NHLH2", "NKRF", "NKX1-1", "NKX1-2", "NKX2-1", "NKX2-2", "NKX2-3", "NKX2-4", "NKX2-5", "NKX2-6", "NKX2-8", "NKX3-1", "NKX3-2", "NKX6-1", "NKX6-2", "NKX6-3", "NLRC5", "NME1", "NOBOX", "NOTO", "NPAS1", "NPAS2", "NPAS3", "NPAS4", "NR1D1", "NR1D2", "NR1H2", "NR1H3", "NR1H4", "NR1I2", "NR1I3", "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1", "NRF1", "NRIP1", "NRL", "NSD1", "OLIG1", "OLIG2", "OLIG3", "ONECUT1", "ONECUT2", "ONECUT3", "OSR1", "OSR2", "OTP", "OTX1", "OTX2", "OVOL1", "OVOL2", "OVOL3", "PATZ1", "PAX1", "PAX2",
#     "PAX3", "PAX4", "PAX5", "PAX6", "PAX7", "PAX8", "PAX9", "PBX1", "PBX2", "PBX3", "PBX4", "PDX1", "PEG3", "PER1", "PGBD1", "PGR", "PHOX2A", "PHOX2B", "PITX1", "PITX2", "PITX3", "PKNOX1", "PKNOX2", "PLAG1", "PLAGL1", "POU1F1", "POU2F1", "POU2F2", "POU2F3", "POU3F1", "POU3F2", "POU3F3", "POU3F4", "POU4F1", "POU4F2", "POU4F3", "POU5F1", "POU5F1B", "POU5F2", "POU6F1", "POU6F2", "PPARA", "PPARD", "PPARG", "PRDM1", "PRDM14", "PRDM15", "PRDM16", "PRDM2", "PRDM4", "PRDM5", "PRMT5", "PROP1", "PROX1", "PROX2", "PRRX1", "PRRX2", "PTF1A", "PURA", "PURB", "PURG", "RAD23B", "RARA", "RARB", "RARG", "RAX", "RAX2", "RB1", "RBAK", "RBBP4", "RBL1", "RBL2", "RBMX", "RBPJ", "RBPJL", "REL", "RELA", "RELB", "REST", "RFX1", "RFX2", "RFX3", "RFX4", "RFX5", "RFX6", "RFX7", "RFX8", "RFXANK", "RFXAP", "RHOXF1", "RHOXF2", "RHOXF2B", "RORA", "RORB", "RORC", "RPS3", "RREB1", "RUNX1", "RUNX2", "RUNX3", "RUVBL2", "RXRA", "RXRB", "RXRG", "SAFB", "SALL1", "SALL2", "SALL3", "SALL4", "SARS1", "SATB1", "SATB2", "SCRT1", "SCRT2", "SCX", "SIM1", "SIM2", "SIRT1", "SIX1", "SIX2", "SIX3", "SIX4", "SIX5", "SIX6", "SKI", "SKIL", "SKOR1", "SKOR2", "SLC2A4RG", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD9", "SMYD3", "SNAI1", "SNAI2", "SNAI3", "SNAPC3", "SOHLH1", "SOHLH2", "SOX1", "SOX10", "SOX11", "SOX12", "SOX13", "SOX14", "SOX15", "SOX17", "SOX18", "SOX2", "SOX21", "SOX3", "SOX30", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8", "SOX9",
#     "SP1", "SP2", "SP3", "SP4", "SP5", "SP6", "SP7", "SP8", "SP9", "SPI1", "SPIB", "SPIC", "SREBF1", "SREBF2", "SRF", "SRY", "ST18", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "STK16", "STOX1", "SUB1", "TAF1", "TAL1", "TAL2", "TARDBP", "TBP", "TBPL1", "TBR1", "TBX1", "TBX10", "TBX15", "TBX18", "TBX19", "TBX2", "TBX20", "TBX21", "TBX22", "TBX3", "TBX4", "TBX5", "TBX6", "TBXT", "TCF12", "TCF15", "TCF21", "TCF23", "TCF24", "TCF3", "TCF4", "TCF7", "TCF7L1", "TCF7L2", "TCFL5", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "TEF", "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E", "TFAP4", "TFCP2", "TFCP2L1", "TFDP1", "TFDP2", "TFE3", "TFEB", "TFEC", "TGIF1", "TGIF2", "TGIF2LX", "TGIF2LY", "THAP1", "THAP11", "THRA", "THRAP3", "THRB", "TLX1", "TLX2", "TLX3", "TOP1", "TP53", "TP63", "TP73", "TPRX1", "TRIM24", "TRPS1", "TWIST1", "TWIST2", "UBP1", "USF1", "USF2", "USP3", "UTY", "VAX1", "VAX2", "VDR", "VENTX", "VEZF1", "VSX2", "WBP2", "WIZ", "WT1", "XBP1", "YAP1", "YBX3", "YY1", "YY2", "ZBED1", "ZBTB1", "ZBTB10", "ZBTB11", "ZBTB12", "ZBTB14", "ZBTB16", "ZBTB17", "ZBTB18", "ZBTB2", "ZBTB20", "ZBTB22", "ZBTB24", "ZBTB25", "ZBTB26", "ZBTB3", "ZBTB32", "ZBTB33", "ZBTB34", "ZBTB37", "ZBTB38", "ZBTB39",
#     "ZBTB4", "ZBTB40", "ZBTB41", "ZBTB42", "ZBTB43", "ZBTB45", "ZBTB46", "ZBTB47", "ZBTB48", "ZBTB49", "ZBTB5", "ZBTB6", "ZBTB7A", "ZBTB7B", "ZBTB7C", "ZBTB8A", "ZBTB8B", "ZBTB9", "ZC3H8", "ZEB1", "ZEB2", "ZFAT", "ZFHX2", "ZFHX3", "ZFHX4", "ZFP1", "ZFP14", "ZFP2", "ZFP28", "ZFP3", "ZFP30", "ZFP37", "ZFP41", "ZFP42", "ZFP57", "ZFP62", "ZFP69", "ZFP69B", "ZFP82", "ZFP90", "ZFP92", "ZGPAT", "ZHX3", "ZIC1", "ZIC2", "ZIC3", "ZIC4", "ZIC5", "ZIK1", "ZIM2", "ZIM3", "ZKSCAN1", "ZKSCAN2", "ZKSCAN3", "ZKSCAN4", "ZKSCAN5", "ZKSCAN7", "ZKSCAN8", "ZNF10", "ZNF100", "ZNF101", "ZNF107", "ZNF112", "ZNF114", "ZNF117", "ZNF12", "ZNF121", "ZNF124", "ZNF131", "ZNF132", "ZNF133", "ZNF134", "ZNF135", "ZNF136", "ZNF138", "ZNF14", "ZNF140", "ZNF141", "ZNF142", "ZNF143", "ZNF146", "ZNF148", "ZNF154", "ZNF155", "ZNF157", "ZNF16", "ZNF160", "ZNF165", "ZNF169", "ZNF17", "ZNF174", "ZNF175", "ZNF177", "ZNF18", "ZNF180", "ZNF181", "ZNF182", "ZNF184", "ZNF189", "ZNF19", "ZNF195", "ZNF197", "ZNF2", "ZNF20", "ZNF202", "ZNF205", "ZNF208", "ZNF211", "ZNF212", "ZNF213", "ZNF214", "ZNF215", "ZNF217", "ZNF219", "ZNF22", "ZNF221", "ZNF222", "ZNF223", "ZNF224", "ZNF225", "ZNF226", "ZNF227", "ZNF229", "ZNF23", "ZNF230", "ZNF232",
#     "ZNF233", "ZNF234", "ZNF235", "ZNF236", "ZNF239", "ZNF24", "ZNF248", "ZNF25", "ZNF250", "ZNF251", "ZNF253", "ZNF254", "ZNF256", "ZNF257", "ZNF26", "ZNF260", "ZNF263", "ZNF264", "ZNF266", "ZNF267", "ZNF268", "ZNF273", "ZNF274", "ZNF275", "ZNF277", "ZNF28", "ZNF280A", "ZNF280B", "ZNF280C", "ZNF280D", "ZNF281", "ZNF282", "ZNF283", "ZNF284", "ZNF285", "ZNF286A", "ZNF287", "ZNF296", "ZNF3", "ZNF30", "ZNF300", "ZNF302", "ZNF304", "ZNF311", "ZNF317", "ZNF319", "ZNF32", "ZNF320", "ZNF322", "ZNF324", "ZNF324B", "ZNF329", "ZNF331", "ZNF333", "ZNF334", "ZNF335", "ZNF337", "ZNF33A", "ZNF33B", "ZNF34", "ZNF341", "ZNF343", "ZNF345", "ZNF347", "ZNF35", "ZNF350", "ZNF354A", "ZNF354B", "ZNF354C", "ZNF358", "ZNF362", "ZNF366", "ZNF367", "ZNF37A", "ZNF382", "ZNF383", "ZNF384", "ZNF391", "ZNF394", "ZNF395", "ZNF396", "ZNF397", "ZNF398", "ZNF404", "ZNF408", "ZNF41", "ZNF415", "ZNF416", "ZNF417", "ZNF418", "ZNF419", "ZNF420", "ZNF423", "ZNF425", "ZNF426", "ZNF429", "ZNF43", "ZNF430", "ZNF431", "ZNF432", "ZNF433", "ZNF436", "ZNF438", "ZNF439", "ZNF44", "ZNF440", "ZNF441", "ZNF442", "ZNF443", "ZNF444", "ZNF445", "ZNF446", "ZNF449", "ZNF45", "ZNF454", "ZNF460", "ZNF461", "ZNF467", "ZNF468", "ZNF470", "ZNF471", "ZNF473", "ZNF479", "ZNF48", "ZNF480", "ZNF483", "ZNF484", "ZNF485", "ZNF486", "ZNF490", "ZNF491", "ZNF492", "ZNF493", "ZNF496", "ZNF497", "ZNF500", "ZNF501", "ZNF502", "ZNF506", "ZNF510",
#     "ZNF514", "ZNF516", "ZNF517", "ZNF518A", "ZNF518B", "ZNF519", "ZNF521", "ZNF524", "ZNF526", "ZNF527", "ZNF528", "ZNF529", "ZNF530", "ZNF534", "ZNF536", "ZNF540", "ZNF543", "ZNF544", "ZNF546", "ZNF547", "ZNF548", "ZNF549", "ZNF550", "ZNF551", "ZNF552", "ZNF554", "ZNF555", "ZNF556", "ZNF557", "ZNF558", "ZNF559", "ZNF560", "ZNF561", "ZNF562", "ZNF563", "ZNF564", "ZNF565", "ZNF566", "ZNF567", "ZNF568", "ZNF569", "ZNF57", "ZNF570", "ZNF571", "ZNF572", "ZNF573", "ZNF574", "ZNF575", "ZNF577", "ZNF578", "ZNF580", "ZNF581", "ZNF582", "ZNF583", "ZNF584", "ZNF585A", "ZNF585B", "ZNF586", "ZNF587", "ZNF587B", "ZNF589", "ZNF594", "ZNF595", "ZNF596", "ZNF597", "ZNF599", "ZNF600", "ZNF605", "ZNF606", "ZNF607", "ZNF610", "ZNF611", "ZNF613", "ZNF614", "ZNF615", "ZNF616", "ZNF619", "ZNF620", "ZNF621", "ZNF623", "ZNF624", "ZNF625", "ZNF626", "ZNF627", "ZNF628", "ZNF629", "ZNF630", "ZNF639", "ZNF641", "ZNF644", "ZNF646", "ZNF648", "ZNF649", "ZNF652", "ZNF655", "ZNF658", "ZNF66", "ZNF660", "ZNF662", "ZNF664", "ZNF665", "ZNF667", "ZNF668", "ZNF669", "ZNF670", "ZNF671", "ZNF674", "ZNF675", "ZNF676", "ZNF677", "ZNF678", "ZNF679", "ZNF680", "ZNF681", "ZNF682", "ZNF683", "ZNF684", "ZNF688", "ZNF689", "ZNF69",
#     "ZNF691", "ZNF692", "ZNF695", "ZNF696", "ZNF697", "ZNF699", "ZNF7", "ZNF70", "ZNF700", "ZNF701", "ZNF704", "ZNF705A", "ZNF705B", "ZNF705D", "ZNF705E", "ZNF705G", "ZNF707", "ZNF708", "ZNF709", "ZNF71", "ZNF710", "ZNF713", "ZNF716", "ZNF718", "ZNF721", "ZNF723", "ZNF724", "ZNF726", "ZNF727", "ZNF728", "ZNF729", "ZNF730", "ZNF732", "ZNF735", "ZNF736", "ZNF737", "ZNF74", "ZNF740", "ZNF746", "ZNF749", "ZNF750", "ZNF75A", "ZNF75D", "ZNF76", "ZNF761", "ZNF763", "ZNF764", "ZNF765", "ZNF768", "ZNF77", "ZNF770", "ZNF771", "ZNF772", "ZNF773", "ZNF774", "ZNF775", "ZNF776", "ZNF777", "ZNF778", "ZNF780A", "ZNF780B", "ZNF782", "ZNF783", "ZNF784", "ZNF785", "ZNF786", "ZNF789", "ZNF79", "ZNF790", "ZNF791", "ZNF792", "ZNF793", "ZNF799", "ZNF8", "ZNF80", "ZNF805", "ZNF808", "ZNF81", "ZNF813", "ZNF814", "ZNF816", "ZNF821", "ZNF823", "ZNF829", "ZNF83", "ZNF835", "ZNF836", "ZNF837", "ZNF84", "ZNF841", "ZNF844", "ZNF845", "ZNF846", "ZNF85", "ZNF850", "ZNF852", "ZNF853", "ZNF860", "ZNF865", "ZNF875", "ZNF878", "ZNF879", "ZNF880", "ZNF888", "ZNF891", "ZNF90", "ZNF91", "ZNF92", "ZNF93", "ZNF98", "ZNF99", "ZSCAN1", "ZSCAN10", "ZSCAN12", "ZSCAN16", "ZSCAN18", "ZSCAN2", "ZSCAN20", "ZSCAN21", "ZSCAN22", "ZSCAN23", "ZSCAN25", "ZSCAN26", "ZSCAN29", "ZSCAN30", "ZSCAN31", "ZSCAN32", "ZSCAN4", "ZSCAN5A", "ZSCAN5B", "ZSCAN5C", "ZSCAN9"
#   ))


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(log2FoldChange.gsam.res > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2, 2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(col.chr7 = ifelse(gid %in% sel$gid  ,"at chr7",  "not at chr7")) 


ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = col.chr7 , size = col.chr7, col = col.chr7  ,fill =  col.chr7  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(data = subset(plt, col.chr7 != "at chr7")) +
  geom_point(data = subset(plt, col.chr7 == "at chr7")) +
  scale_shape_manual(values = c('truncated'= 4, 'not at chr7' = 19 , 'at chr7'= 21 )    ) +
  scale_size_manual(values = c( 'not at chr7' = 0.1 , 'at chr7'= 0.8 )    ) +
  scale_color_manual(values = c('truncated'= 'black', 'not at chr7' = rgb(0,0,0,0.35), 'at chr7'= rgb(0,0,0,0.5)  )    ) +
  scale_fill_manual(values = c('truncated'= 'black', 'not at chr7' = rgb(0,0,0,0.35), 'at chr7'= rgb(1,0,0,1)  )    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") + 
  xlim(-2,2) + ylim(-16.5,16.5)





## 4 :: [E:1] :: FC res x corr TPC G-SAM  [diag; chr7 ; chr10 ] ----


results.out %>%
  dplyr::filter(hugo_symbol %in% c('SMAD5','SETD5')) %>%
  dplyr::select(c('hugo_symbol', 'estimate.gsam.cor.tpc', 'estimate.glass.cor.tpc'))


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(log2FoldChange.gsam.res > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(col.chr7 = case_when(is.limited.gsam.res == T ~ "truncated" , chr == "chr7" ~ "at chr7", chr != "chr7" ~ "not at chr7")) %>%
  dplyr::mutate(col.chr10 = case_when(is.limited.gsam.res == T ~ "truncated" , chr == "chr10" ~ "at chr10", chr != "chr10" ~ "not at chr10")) 


p1 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = is.limited.gsam.res ,
                      size = is.limited.gsam.res  ,
                      col = is.limited.gsam.res   ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.gsam.res > 0.05),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.35))    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)

p2 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = col.chr7,
                      size = col.chr7  ,
                      col = col.chr7  ,
                      fill =  col.chr7  
) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(data = subset(plt, chr != "chr7")) +
  geom_point(data = subset(plt, chr == "chr7", stroke=0.001)) +
  scale_shape_manual(values = c('truncated'= 4, 'not at chr7' = 19 , 'at chr7'= 21 )    ) +
  scale_size_manual(values = c('truncated'= 0.85, 'not at chr7' = 0.1 , 'at chr7'= 0.8 )    ) +
  scale_color_manual(values = c('truncated'= 'black', 'not at chr7' = rgb(0,0,0,0.35), 'at chr7'= rgb(0,0,0,0.5)  )    ) +
  scale_fill_manual(values = c('truncated'= 'black', 'not at chr7' = rgb(0,0,0,0.35), 'at chr7'= rgb(1,0,0,1)  )    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)

p3 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y = statistic.gsam.cor.tpc, 
                      shape = col.chr10,
                      size = col.chr10  ,
                      col = col.chr10  ,
                      fill =  col.chr10  
) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(data = subset(plt, chr != "chr10")) +
  geom_point(data = subset(plt, chr == "chr10", stroke=0.001)) +
  scale_shape_manual(values = c('truncated'= 4, 'not at chr10' = 19 , 'at chr10'= 21 )    ) +
  scale_size_manual(values = c('truncated'= 0.85, 'not at chr10' = 0.1 , 'at chr10'= 0.8 )    ) +
  scale_color_manual(values = c('truncated'= 'black', 'not at chr10' = rgb(0,0,0,0.35), 'at chr10'= rgb(0,0,0,0.5)  )    ) +
  scale_fill_manual(values = c('truncated'= 'black', 'not at chr10' = rgb(0,0,0,0.35), 'at chr10'= rgb(1,0,0,1)  )    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)

p1 / p2 / p3


ggsave('output/figures/figure-DGE_log2FoldChange.res_x_statistic.cor.tpc.png', width = 5.7 * 1.4 , height = 4 * 1.4 * 3)




## 2.5 :: G-SAM.res & GLASS.res ----


### a :: GLASS res x G-SAM TPC corr ----


plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>% # TPC correlation G-SAM
  dplyr::filter(!is.na(statistic.glass.cor.tpc)) %>% # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(log2FoldChange.glass.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 2.5, 2.5 , log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res < -2.5, -2.5 , log2FoldChange.glass.tpc.res)) 



p1 <- ggplot(plt, aes(x = log2FoldChange.glass.res ,
                      y =  statistic.glass.cor.tpc,
                      shape = is.limited.glass.res ,
                      size = is.limited.glass.res  ,
                      col = is.limited.glass.res  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.glass.res > 0.05 &  is.limited.glass.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)

p2 <- ggplot(plt, aes(x = log2FoldChange.glass.tpc.res ,
                      y =  statistic.glass.cor.tpc,
                      shape = is.limited.glass.tpc.res ,
                      size = is.limited.glass.tpc.res  ,
                      col = is.limited.glass.tpc.res  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.glass.tpc.res > 0.05 &
                              is.limited.glass.tpc.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)


p1 + p2 





# ---


p3 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                      y =  statistic.gsam.cor.tpc,
                      shape = is.limited.gsam.res ,
                      size = is.limited.gsam.res  ,
                      col = is.limited.gsam.res  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.gsam.res > 0.05 &  is.limited.gsam.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)

p4 <- ggplot(plt, aes(x = log2FoldChange.gsam.tpc.res ,
                      y =  statistic.gsam.cor.tpc,
                      shape = is.limited.gsam.tpc.res ,
                      size = is.limited.gsam.tpc.res  ,
                      col = is.limited.gsam.tpc.res  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.05 &  is.limited.gsam.tpc.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5) +
  ylim(-16.5,16.5)

(p3 + p4) / (p1 + p2)



ggsave("output/figures/figure-DE_with_tcp_correction_gsam.png",width = 5.7 * 1.4  * 2, height = 4 * 1.4 * 2 )

## Figure 3D-G ----


plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>% # TPC correlation G-SAM
  dplyr::filter(!is.na(statistic.glass.cor.tpc)) %>% # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(log2FoldChange.glass.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 2.5, 2.5 , log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res < -2.5, -2.5 , log2FoldChange.glass.tpc.res)) 

plt.expanded <- rbind(
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.res) %>%
    dplyr::rename(padj = padj.gsam.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(padj = padj.gsam.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.tpc.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = T) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.res) %>%
    dplyr::rename(padj = padj.glass.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.tpc.res) %>%
    dplyr::rename(padj = padj.glass.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.tpc.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = T) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj)
  ) %>%
  dplyr::mutate(DESeq2.tcp.corrected = ifelse(DESeq2.tcp.corrected == T, "Tumor cell-% corrected" , "Tumor cell-% uncorrected")) %>%
  dplyr::mutate(DESeq2.tcp.corrected = factor(DESeq2.tcp.corrected, levels=
                                                c("Tumor cell-% uncorrected" , "Tumor cell-% corrected"))) %>%
  dplyr::mutate(is.limited = is.limited == 'TRUE') 


ggplot(plt.expanded, aes(x = log2FoldChange ,
                y =  cor.stat,
                shape = is.limited ,
                col = is.limited ,
                #size = is.limited  ,
                fill = dataset  ) ) +
  facet_grid(rows = vars(dataset), cols=vars(DESeq2.tcp.corrected), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2) +
  geom_smooth(data = subset(plt.expanded, padj > 0.05 &  is.limited == FALSE),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="#ff2929", show.legend=F ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e6935622', 'GLASS' = '#69a4d522') ) + 
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5)


ggsave("output/figures/DESeq2_x_tumor_cell_percentage.pdf", width = 5.7 * 2, height = 4 * 2)


## Figure 4AB ----


# show CD14, CD74, CD163,
# CD33, CD84, FCER1G, GPR34, TMEM119



plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>% # TPC correlation G-SAM
  dplyr::filter(!is.na(statistic.glass.cor.tpc)) %>% # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.mg.or.tam = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == "microglia/TAM") 




plt.expanded <- rbind(
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.res) %>%
    dplyr::rename(padj = padj.gsam.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam, hugo_symbol)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.res) %>%
    dplyr::rename(padj = padj.glass.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam, hugo_symbol)
) %>%
  dplyr::mutate(is.limited = is.limited == 'TRUE') %>%
  dplyr::mutate(col = case_when(
    is.mg.or.tam == F ~ "no-mg-or-tam",
    is.mg.or.tam == T ~ dataset
  )) %>%
  dplyr::mutate(show.label = hugo_symbol %in% c("CD33", "CD84", "FCER1G", "GPR34", "CD14","CD74", "CD163"))


ggplot(plt.expanded, aes(x = log2FoldChange ,
                         y =  cor.stat,
                         shape = is.limited ,
                         col =  is.limited ,
                         label = hugo_symbol,
                         #size = is.limited  ,
                         fill = col  ) ) +
  facet_grid(cols = vars(dataset), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.mg.or.tam` == F)) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.mg.or.tam` == T)) +
  geom_smooth(data = subset(plt.expanded, padj > 0.05 &  is.limited == FALSE),
              aes(group=1),
              method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="#ff2929", show.legend=F ) +
  geom_text_repel(data=subset(plt.expanded, show.label  & log2FoldChange > 0), size=2.5 , 
                  segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, 
                  direction = "y", hjust = "left", col="black", nudge_y = -4.8 ) +
  geom_text_repel(data=subset(plt.expanded, show.label  & log2FoldChange < 0), size=2.5 ,
                  segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, 
                  direction = "y", hjust = "right", col="black", nudge_y = -4.8 ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e69356', 'GLASS' = '#69a4d5', 'no-mg-or-tam'='white') ) + 
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5)



ggsave('output/figures/geiser_plot_mckenzy_immune_tam.pdf', height=5.75 * 1.1,width=12 * 1.1)

## Figure 4AB restyled double ----


# show CD14, CD74, CD163,
# CD33, CD84, FCER1G, GPR34, TMEM119



plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>% # TPC correlation G-SAM
  dplyr::filter(!is.na(statistic.glass.cor.tpc)) %>% # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(log2FoldChange.glass.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 2.5, 2.5 , log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res < -2.5, -2.5 , log2FoldChange.glass.tpc.res)) %>% 
  dplyr::mutate(is.mg.or.tam = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == "microglia/TAM")



# pivot longer is complicated because both the log2fc & cor.stat need to be pivotted

plt.expanded <- rbind(
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.res) %>%
    dplyr::rename(padj = padj.gsam.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(padj = padj.gsam.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.tpc.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = T) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.res) %>%
    dplyr::rename(padj = padj.glass.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.tpc.res) %>%
    dplyr::rename(padj = padj.glass.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.tpc.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = T) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
) %>%
  dplyr::mutate(DESeq2.tcp.corrected = ifelse(DESeq2.tcp.corrected == T, "Tumor cell-% corrected" , "Tumor cell-% uncorrected")) %>%
  dplyr::mutate(DESeq2.tcp.corrected = factor(DESeq2.tcp.corrected, levels=
                                                c("Tumor cell-% uncorrected" , "Tumor cell-% corrected"))) %>%
  dplyr::mutate(is.limited = is.limited == 'TRUE') %>% 
  dplyr::mutate(col = case_when(
    is.mg.or.tam == F ~ "no-mg-or-tam",
    is.mg.or.tam == T ~ dataset
  )) %>%
  dplyr::mutate(show.label = hugo_symbol %in% c("CD33", "CD84", "FCER1G", "GPR34", "CD14","CD74", "CD163")) %>% 
  dplyr::rename(`Purity corrected` = DESeq2.tcp.corrected) 




ggplot(plt.expanded, aes(x = log2FoldChange ,
                         y =  cor.stat,
                         shape = is.limited ,
                         col =  is.limited ,
                         label = hugo_symbol,
                         #size = is.limited  ,
                         fill = col  ) ) +
  facet_grid(rows = vars(dataset), cols=vars(`Purity corrected`), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.mg.or.tam` == F)) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.mg.or.tam` == T)) +
  geom_smooth(data = subset(plt.expanded, padj > 0.05 &  is.limited == FALSE),
              aes(group=1),
              method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="#ff2929", show.legend=F ) +
  geom_text_repel(data=subset(plt.expanded, show.label  & log2FoldChange > 0), size=2.5 , 
                  segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, 
                  direction = "y", hjust = "left", col="black", nudge_y = -4.8 ) +
  geom_text_repel(data=subset(plt.expanded, show.label  & log2FoldChange < 0), size=2.5 ,
                  segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, 
                  direction = "y", hjust = "right", col="black", nudge_y = -4.8 ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e69356', 'GLASS' = '#69a4d5', 'no-mg-or-tam'='white') ) + 
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5)



#ggsave('output/figures/geiser_plot_mckenzy_immune_tam_double.pdf', height=5.75 * 1.1,width=12 * 1.1)
ggsave('output/figures/geiser_plot_mckenzy_immune_tam_double.pdf', height=5.75 * 1.1 * 1.1,width=12 * 1.1)





## Figure 6 [EM] ----



plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>% # TPC correlation G-SAM
  dplyr::filter(!is.na(statistic.glass.cor.tpc)) %>% # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(is.em = !is.na(EM.struct.constituent) & EM.struct.constituent == T)



plt.expanded <- rbind(
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.res) %>%
    dplyr::rename(padj = padj.gsam.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.em)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.res) %>%
    dplyr::rename(padj = padj.glass.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.em)
) %>%
  dplyr::mutate(is.limited = is.limited == 'TRUE') %>%
  dplyr::mutate(col = case_when(
    is.em == F ~ "no-EM",
    is.em == T ~ dataset
  ))


ggplot(plt.expanded, aes(x = log2FoldChange ,
                         y =  cor.stat,
                         shape = is.limited ,
                         col =  is.limited ,
                         #size = is.limited  ,
                         fill = col  ) ) +
  facet_grid(cols = vars(dataset), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.em` == F)) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.em` == T)) +
  geom_smooth(data = subset(plt.expanded, padj > 0.05 &  is.limited == FALSE),
              aes(group=1),
              method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="#ff2929", show.legend=F ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e69356', 'GLASS' = '#69a4d5', 'no-EM'='white') ) + 
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5)


ggsave('output/figures/geiser_plot_mckenzy_immune_tam.pdf', height=5.75 * 1.1,width=12 * 1.1)



### C6 & pericytes ----

C6 <- c('CRABP2', 'CILP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
        'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
        "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
        "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
        "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")



plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) %>% # TPC correlation G-SAM
  dplyr::filter(!is.na(statistic.glass.cor.tpc)) %>% # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) %>%
  dplyr::filter(!is.na(log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(log2FoldChange.glass.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 2.5, 2.5 , log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res < -2.5, -2.5 , log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(is.C6 = hugo_symbol %in% C6) %>% 
  dplyr::mutate(pericyte.marker = hugo_symbol %in% c('PDGFRB','CD248','RGS5','HEYL','CFH'))



plt.expanded <- rbind(
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(padj = padj.gsam.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.tpc.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.C6, pericyte.marker)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.tpc.res) %>%
    dplyr::rename(padj = padj.glass.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.glass.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.glass.tpc.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.C6, pericyte.marker)
) %>%
  dplyr::mutate(is.limited = is.limited == 'TRUE') %>%
  dplyr::mutate(col = case_when(
    pericyte.marker == T ~ "pericyte marker",
    is.C6 == T ~ "C6",
    T ~ "other"
  ))


ggplot(plt.expanded, aes(x = log2FoldChange ,
                         y =  cor.stat,
                         shape = is.limited ,
                         col =  is.limited ,
                         #size = is.limited  ,
                         label=hugo_symbol,
                         fill = col  ) ) +
  facet_grid(cols = vars(dataset), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded, `col` == "other")) +
  geom_point(size=2.2, data = subset(plt.expanded, `col` == "C6")) +
  geom_point(size=2.2, data = subset(plt.expanded, `col` == "pericyte marker")) +
  #geom_smooth(data = subset(plt.expanded, padj > 0.05 &  is.limited == FALSE),
  #            aes(group=1),
  #            method="lm",
  #            se = FALSE,  formula=y ~ x, orientation="y", col="#ff2929", show.legend=F ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = c('C6' = '#e69356', 'pericyte marker' = '#69a4d5', 'other'='white') ) + 
  geom_text_repel(data=subset(plt.expanded, col == "C6"), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, direction = "y", hjust = "left", col="gray30" ) +
  geom_text_repel(data=subset(plt.expanded, col == "pericyte marker"), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, direction = "y", hjust = "right",col="black" ) +
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.75,2.75)

ggsave('output/figures/geiser_plot_C6_pericyte.pdf', height=5.75 * 1.1,width=12 * 1.1)
ggsave('output/figures/geiser_plot_C6_pericyte.png', height=5.75 * 1.1,width=12 * 1.1)



### b :: GLASS res x G-SAM res ----


plt <- results.out %>%
  dplyr::filter(!is.na(stat.glass.res) & !is.na(stat.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(stat.gsam.res) > 8)) %>% # change pch to something that is limited
  dplyr::mutate(stat.gsam.res = ifelse(stat.gsam.res > 8, 8 , stat.gsam.res)) %>%
  dplyr::mutate(stat.gsam.res = ifelse(stat.gsam.res < -8, -8 , stat.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(stat.glass.res) > 8)) %>% # change pch to something that is limited
  dplyr::mutate(stat.glass.res = ifelse(stat.glass.res > 8, 8 , stat.glass.res))  %>%
  dplyr::mutate(stat.glass.res = ifelse(stat.glass.res < -8, -8 , stat.glass.res))  %>%
  dplyr::mutate(is.limited = is.limited.gsam.res == "TRUE" | is.limited.glass.res == "TRUE")



p2 <- ggplot(plt, aes(x = stat.gsam.res ,
                y =  stat.glass.res ,
                shape = is.limited ,
                size = is.limited  ,
                col = is.limited  ) ) +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point() +
  geom_smooth(data = subset(plt, is.limited == F), method="lm", se = FALSE,  formula=y ~ x, col="red" , size=0.8) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
  scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
  scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
  youri_gg_theme + 
  labs(x = "DESeq2 stat G-SAM",
       y="DESeq2 stat GLASS",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") + 
  xlim(-8,8) +
  ylim(-8,8)



cor(plt$stat.gsam.res, plt$stat.glass.res) # r=0.63


p1 / p2



ggsave('output/figures/paper_dge_glass_x_gsam_uncorercted.png', width = 5.7 * 1.4 , height = 4 * 1.4 * 2)



# ggplot(plt, aes(x = stat.gsam.res ,
#                 y =  stat.glass.res ,
#                 shape = is.limited ,
#                 size = is.limited  ,
#                 col = is.limited  ) ) +
#   geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
#   geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
#   geom_point() +
#   #geom_smooth(data = subset(plt, is.limited == F), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
#   scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19)   )  +
#   scale_size_manual(values = c('TRUE'= 0.85, 'FALSE' = 0.1)    ) +
#   scale_color_manual(values = c('TRUE'= 'black', 'FALSE' = rgb(0,0,0,0.65))    ) +
#   youri_gg_theme + 
#   labs(x = "log2FC R1 vs. R2 G-SAM",
#        y="log2FC R1 vs. R2 GLASS",
#        shape = "Truncated at x-axis",
#        size="Truncated at x-axis",
#        col="Truncated at x-axis") +
#   xlim(-2.5,2.5) +
#   ylim(-2.5, 2.5)



## 2.6 :: GLASS & ch19 ----

# wees zeker dat beide plots exact dezelfde genen bevatten


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res) & !is.na(log2FoldChange.gsam.res) & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res))  %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res))  %>%
  dplyr::mutate(is.chr19 = as.character(chr == "chr19")) %>%
  dplyr::mutate(size.type.glass = case_when(
    is.limited.glass.res == "TRUE" ~ "limited",
    is.chr19 == "TRUE" ~ "chr19",
    TRUE ~ "not chr19")) %>%
  dplyr::mutate(size.type.gsam = case_when(
    is.limited.gsam.res == "TRUE" ~ "limited",
    is.chr19 == "TRUE" ~ "chr19",
    TRUE ~ "not chr19"))



### a :: GLASS ----



plt$is.limited.gsam.res %>% as.factor %>% summary

plt$is.limited %>% as.factor %>% summary
plt$size.type %>% as.factor %>% summary


p1 <- ggplot(plt, aes(x = log2FoldChange.glass.res ,
                      y =  statistic.gsam.cor.tpc,
                      shape = is.limited.glass.res ,
                      size = size.type.glass  ,
                      col = is.chr19))   +
  geom_point(data=subset(plt, is.chr19 == 'FALSE')) +
  geom_point(data=subset(plt, is.chr19 == 'TRUE')) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19))  +
  scale_size_manual(values = c('chr19'= 0.65,
                               'not chr19' = 0.1,
                               'limited' = 0.86)) +
  scale_color_manual(values = c('FALSE'= 'black', 'TRUE' = rgb(1,0,0,0.85)))  +
  youri_gg_theme + 
  labs(x = "log2FC GLASS",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis")  +
  xlim(-2.5,2.5)



### b :: GSAM ----


p2 <- ggplot(plt, aes(x = log2FoldChange.gsam.res ,
                y =  statistic.gsam.cor.tpc,
                shape = is.limited.gsam.res ,
                size = size.type.gsam  ,
                col = is.chr19))   +
  geom_point(data=subset(plt, is.chr19 == 'FALSE')) +
  geom_point(data=subset(plt, is.chr19 == 'TRUE')) +
  scale_shape_manual(values = c('TRUE'= 4, 'FALSE' = 19))  +
  scale_size_manual(values = c('chr19'= 0.65,
                               'not chr19' = 0.1,
                               'limited' = 0.86)) +
  scale_color_manual(values = c('FALSE'= 'black', 'TRUE' = rgb(1,0,0,0.85)))  +
  youri_gg_theme + 
  labs(x = "log2FC G-SAM",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis")  +
  xlim(-2.5,2.5)



p1 + p2




ggsave('output/figures/paper_dge_glass_x_gsam_chr19.png', width = 5.7 * 1.1 * 2 , height = 4 * 1.6 )




## 2.7 :: G-SAM FC corrected  ----



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(label.chr7 = case_when(chr == "chr7" ~ "chr7", TRUE ~ "remaining" )) %>%
  dplyr::mutate(label.chr10 = case_when(chr == "chr10" ~ "chr10", TRUE ~ "remaining" ))


plt <- rbind(
  plt %>%
    dplyr::mutate(lfc = log2FoldChange.gsam.res ) %>%
    dplyr::mutate(cor = statistic.gsam.cor.tpc) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DGE.tumor.cell.percentage.correction =  'Tumor cell percentage: Not Corrected')
  ,
  plt %>%
    dplyr::mutate(lfc = log2FoldChange.gsam.tpc.res ) %>%
    dplyr::mutate(cor = statistic.gsam.cor.tpc) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DGE.tumor.cell.percentage.correction =  'Tumor cell percentage: Corrected')
  ,
  plt %>%
    dplyr::mutate(lfc = log2FoldChange.glass.res ) %>%
    dplyr::mutate(cor = statistic.glass.cor.tpc) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DGE.tumor.cell.percentage.correction =  'Tumor cell percentage: Not Corrected')
  ,
  plt %>%
    dplyr::mutate(lfc = log2FoldChange.glass.tpc.res ) %>%
    dplyr::mutate(cor = statistic.glass.cor.tpc) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DGE.tumor.cell.percentage.correction =  'Tumor cell percentage: Corrected')
) %>% 
  dplyr::mutate(DGE.tumor.cell.percentage.correction = factor(DGE.tumor.cell.percentage.correction, levels = c('Tumor cell percentage: Not Corrected', 'Tumor cell percentage: Corrected')))


### _ :: gg + facet all ----

ggplot(plt, aes(x=lfc, y=cor ) ) + 
  geom_point(cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.05 &  is.limited.gsam.tpc.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_color_manual(values = c('remaining'='black','chr7'='red') ) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       #,col="Difference significant (R1 ~ R2)"
  ) +
  facet_grid(cols = vars(DGE.tumor.cell.percentage.correction), rows=vars(dataset), scales = "free", space = "free") +
  youri_gg_theme + 
  xlim(-2.5,2.5)



ggsave('output/figures/paper_dge_gsam_corrected_uncorrected.png', width = 5.7 * 1.1 * 2 , height = 4 * 1.6 )



## 2.8.5 Geizerplot ----



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2.5, 2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2.5, -2.5 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate( show.label =  
                   #log2FoldChange.glass.res < -1.5 | 
                  hugo_symbol %in% 
                  #c("TOP2A", "CDK1", "DTL", "CCNB1", "XRCC2", "CCNE2", "DSN1", "TIMELESS") # Cell Cycle genes Patel/Bernstein
                  #c("VEGFA", "ADM", "TREM1", "ENO1", "LDHA", "NRN1", "UBC", "GBE1", "MIF")# Hypoxia genes Patel Bernstein
                  c("PTPRZ1", "CCND1", "SCG2", "MAP2", "IDH1", "ENO2", "PLTP", "NOON", "METTL7B", "PIK3R3", "TSPAN3")# Stem cell signature Patel?
                  )


p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.res ,
                y=statistic.gsam.cor.tpc ,
                label = hugo_symbol , 
                col=show.label)) + 
  geom_point(data=subset(plt, show.label == F), size = 0.85, col = 'gray30') +
  geom_point(data=subset(plt, show.label == T), size = 0.85) +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res > 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, direction = "y", hjust = "left" ) +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res < 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, direction = "y", hjust = "right" ) +
  scale_size_manual(values = c('not significant'=0.85, 'significant'=0.85, 'GABA'=1.75, 'neuronal'=1.75 ,'TRUE'=2)) +
  scale_color_manual(values = c('not significant'= rgb(0,0,0,0.15),'significant'= rgb(0,0,0,0.15),'TRUE'='red','GABA'='red','neuronal'=rgb(0.2,0.2,1.0))) +
  labs(x = "log2FC R1 vs. R2 G-SAM (tumor cell-% corrected)",
       y="Correlation t-statistic with tumour percentage GSAM"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + xlim(-3.5, 3.5)


p2 <- ggplot(plt, aes(x=log2FoldChange.gsam.res ,
                y=statistic.glass.cor.tpc ,
                label = hugo_symbol )) + 
  geom_point(size = 0.85, col = 'gray30') +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res > 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, direction = "y", hjust = "left" ) +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res < 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, direction = "y", hjust = "right" ) +
  scale_size_manual(values = c('not significant'=0.85, 'significant'=0.85, 'GABA'=1.75, 'neuronal'=1.75 )) +
  scale_color_manual(values = c('not significant'= rgb(0,0,0,0.15),'significant'= rgb(0,0,0,0.15),'GABA'='red','neuronal'=rgb(0.2,0.2,1.0))) +
  labs(x = "log2FC R1 vs. R2 G-SAM (tumor cell-% corrected)",
       y="Correlation t-statistic with tumour percentage GLASS"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + xlim(-3.5, 3.5)


p3 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc ,
                      label = hugo_symbol )) + 
  geom_point(size = 0.85, col = 'gray30') +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res > 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, direction = "y", hjust = "left" ) +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res < 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, direction = "y", hjust = "right" ) +
  scale_size_manual(values = c('not significant'=0.65, 'significant'=0.85, 'GABA'=1.75, 'neuronal'=1.75 )) +
  scale_color_manual(values = c('not significant'= rgb(0,0,0,0.15),'significant'= rgb(0,0,0,0.15),'GABA'='red','neuronal'=rgb(0.2,0.2,1.0))) +
  labs(x = "log2FC R1 vs. R2 GLASS (tumor cell-% corrected)",
       y="Correlation t-statistic with tumour percentage G-SAM"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  geom_vline(xintercept = 0, col="red") +
  youri_gg_theme + xlim(-3.5, 3.5)


p4 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.glass.cor.tpc ,
                      label = hugo_symbol )) + 
  geom_point(size = 0.85, col = 'gray30') +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res > 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, direction = "y", hjust = "left" ) +
  geom_text_repel(data=subset(plt, show.label  & log2FoldChange.glass.res < 0), size=2.5 , segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, direction = "y", hjust = "right" ) +
  scale_size_manual(values = c('not significant'=0.65, 'significant'=0.85, 'GABA'=1.75, 'neuronal'=1.75 )) +
  scale_color_manual(values = c('not significant'= rgb(0,0,0,0.15),'significant'= rgb(0,0,0,0.15),'GABA'='red','neuronal'=rgb(0.2,0.2,1.0))) +
  labs(x = "log2FC R1 vs. R2 GLASS (tumor cell-% corrected)",
       y="Correlation t-statistic with tumour percentage GLASS"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  geom_vline(xintercept = 0, col="red") +
  youri_gg_theme + xlim(-3.5, 3.5)


#p1

(p1 + p2 ) / (p3 + p4)




## 2.9 Corrected LFC + GAB(R)A labels ----


plt <- results.out %>%
  #dplyr::filter(is.na(padj.glass.res) | (!is.na(padj.glass.res) & padj.glass.res  < 0.01) ) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  #dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res <= 0.01 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) %>%
  dplyr::mutate(show.label.neuronal.label = hugo_symbol %in% c("RBFOX3", "SATB2", "SLC17A7", "RORB", "GAD1", "GAD2", "SST", "LHX6","ADARB2","VIP")) %>%
  dplyr::mutate(show.label.neuronal = hugo_symbol %in% c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3","GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST","VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5","NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK","ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A","RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87","ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2","DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2","CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4","ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2","CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1","GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2","RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")) %>%
  dplyr::mutate(show.label.gaba = hugo_symbol %in% c("GABBR2","GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6","GABRB1","GABRB2","GABRB3","GABRD","GABRE","GABRG1","GABRG2","GABRG3","GABRP","GABRQ","GABRR1","GABRR2","GABRR3")) %>%
  dplyr::mutate(show.label = as.factor ( case_when(
    show.label.gaba == T ~ 'GABA',
    show.label.neuronal == T & show.label.gaba == F ~ 'neuronal',
    significant == T ~ 'significant',
    T ~ 'not significant')))

p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.res ,
                y=statistic.gsam.cor.tpc ,
                col=show.label ,
                size = show.label,
                label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == 'not significant')) +
  geom_point(data=subset(plt, show.label == 'significant')) +
  geom_point(data=subset(plt, show.label == 'neuronal')) +
  geom_point(data=subset(plt, show.label == 'GABA')) +
  scale_size_manual(values = c('not significant'=0.85, 'significant'=0.85, 'GABA'=1.75, 'neuronal'=1.75 )) +
  scale_color_manual(values = c('not significant'= rgb(0,0,0,0.15),'significant'= rgb(0,0,0,0.15),'GABA'='red','neuronal'=rgb(0.2,0.2,1.0))) +
  geom_text_repel(data=subset(plt, (show.label == 'GABA' | show.label.neuronal.label ) & log2FoldChange.gsam.tpc.res > 0), size=2.5 , segment.size = 0.25, segment.linetype = 1,
                  nudge_x = 3.1, direction = "y", hjust = "left" ) +
  geom_text_repel(data=subset(plt, (show.label == 'GABA' | show.label.neuronal.label ) & log2FoldChange.gsam.tpc.res < 0), size=2.5 , segment.size = 0.25, segment.linetype = 1,
                  nudge_x = -3.1, direction = "y", hjust = "right" ) +
  labs(x = "log2FC R1 vs. R2 G-SAM (tumor cell-% corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + xlim(-3.5, 3.5)





# orientation
plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res) & !is.na(log2FoldChange.glass.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 3, 3 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -3, -3 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(significant = padj.glass.res < 0.01 & lfcSE.glass.res < 0.78 & abs(log2FoldChange.glass.res) > 0.5 ) %>%
  #dplyr::mutate(significant = padj.glass.res < 0.01 ) %>%
  dplyr::mutate(show.label.neuronal.label = hugo_symbol %in% c("RBFOX3", "SATB2", "SLC17A7", "RORB", "GAD1", "GAD2", "SST", "LHX6","ADARB2","VIP")) %>%
  dplyr::mutate(show.label.neuronal = hugo_symbol %in% c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3","GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST","VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5","NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK","ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A","RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87","ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2","DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2","CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4","ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2","CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1","GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2","RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")) %>%
  dplyr::mutate(show.label.gaba = hugo_symbol %in% c("GABBR2","GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6","GABRB1","GABRB2","GABRB3","GABRD","GABRE","GABRG1","GABRG2","GABRG3","GABRP","GABRQ","GABRR1","GABRR2","GABRR3")) %>%
  dplyr::mutate(show.label = as.factor ( case_when(
    show.label.gaba == T ~ 'GABA',
    show.label.neuronal == T & show.label.gaba == F ~ 'neuronal',
    significant == T ~ 'significant',
    T ~ 'not significant')))


p2 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc ,
                      col=show.label,
                      size=show.label,
                      label = hugo_symbol)) + 
  geom_point(data=subset(plt, show.label == 'not significant')) +
  geom_point(data=subset(plt, show.label == 'significant')) +
  geom_point(data=subset(plt, show.label == 'neuronal')) +
  geom_point(data=subset(plt, show.label == 'GABA')) +
  scale_size_manual(values = c('not significant'=0.85, 'significant'=0.85, 'GABA'=1.75, 'neuronal'=1.75 )) +
  scale_color_manual(values = c('not significant'= rgb(0,0,0,0.15),'significant'= rgb(0,0,0,0.15),'GABA'='red','neuronal'=rgb(0.2,0.2,1.0))) +
  geom_text_repel(data=subset(plt, (show.label == 'GABA' | show.label.neuronal.label ) == T & log2FoldChange.glass.res > 0), size=2.5 ,
                  nudge_x = 3.1, direction = "y", hjust = "left" , segment.size = 0.25, segment.linetype = 1
                  ) +
  geom_text_repel(data=subset(plt, (show.label == 'GABA' | show.label.neuronal.label ) == T & log2FoldChange.glass.res < 0), size=2.5 ,
                  nudge_x = -3.1, direction = "y", hjust = "right" , segment.size = 0.25 , segment.linetype = 1) +
  labs(x = "log2FC R1 vs. R2 GLASS",
       y="Correlation t-statistic with tumour percentage (in G-SAM!)"
       ,col="Difference significant (R1 ~ R2)" ) +
  youri_gg_theme + 
  xlim(-3.5, 3.5)



p1 + p2



ggsave("output/figures/paper_dge_GABA-genes.png",height=5.7 * 1.1,width=4 * 1.6 * 2)



## 2.10 Corrected LFC [TNNT] ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c(#'CDK4','MDM2','GLI1','GLIS1',
    "TNNI3K","TNN","TNNT2","TNNI1","TNNC1","TNNI2","TNNT3","TNNT1","TNNC2"
  ))

p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res ,
                      y=statistic.gsam.cor.tpc ,
                      col=significant,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label == T),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,
                  nudge_x = 2.9, direction = "y", hjust = "left") + #, lwd=0.5
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,
                  nudge_x = -2.9, direction = "y", hjust = "right")+ #, lwd=0.5
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)




plt <- results.out %>%
  #dplyr::filter(is.na(padj.glass.res) | (!is.na(padj.glass.res) & padj.glass.res  < 0.01) ) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.glass.res < 0.01 & lfcSE.glass.res < 0.6 & abs(log2FoldChange.glass.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c(
    "TNNI3K","TNN","TNNT2","TNNI1","TNNC1","TNNI2","TNNT3","TNNT1","TNNC2"
  ))


p2 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc  ,
                      col=significant,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label == T),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,
                  nudge_x = 3.1, direction = "y", hjust = "left" #, lwd=0.5
  ) +
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,
                  nudge_x = -3.1, direction = "y", hjust = "right" #, lwd=0.5
  ) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 GLASS",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)



p1 + p2



#ggsave("/tmp/gabra.png",height=10 * 1.3,width=4.5 * 1.3)




## 2.** test purposes for other gene sets ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c(#'CDK4','MDM2','GLI1','GLIS1',
    # TIL
    #"ACAP1","ARHGAP15","ARHGAP25","ARHGAP30","ARHGAP9","ARHGAP9","C1orf38","CCL5","CCR2","CCR7","CD2","CD247","CD27","CD28","CD38","CD3D","CD3E","CD3G","CD40","CD48","CD52","CD52","CD53","CD6","CD79A","CD86","CD8A","CFH ","CFHR1","CLEC2D","CORO1A","CSF2RB","CST7","CYBB","DOCK11","DOCK2","DOK2","EVI2B","F5","FAM65B","FCRL3","FCRL5","FYB","GIMAP4","GIMAP5","GIMAP5","GIMAP6","GIMAP7","GLYR1 ","GPR171","GPR18","GPSM3","GVINP1","GZMK","HCLS1","HCST","ICOS","IFFO1","IGHA1 ","IGHA2 ","IGHD ","IGHG1 ","IGHG1 ","IGHG3 ","IGHG4 ","IGHM","IGHM","IGHM ","IGHM ","IGHM ","IGHV4-31 ","IGKV1-5","IGKV4-1","IGLL5 ","IGLV2-11","IKZF1","IKZF1","IL10RA","IL16","IL2RB","IL2RG","IL7R","IL7R","INPP5D","IRF4","ITGA4","ITK","ITM2C","KLHL6","KLRD1","KLRK1","LAT ","LAX1","LCK","LCK","LCP2","LILRB1","LOC100133862","LOC100133862","LOC100133862","LPXN","LST1","LST1","LY9","LY9","MAL","MFNG","MGC29506","MPEG1","MPEG1","MS4A1","MS4A1","MS4A6A","NCF4","NCKAP1L","NLRC3","P2RY8","PAG1","PAG1","PARVG","PAX5","PIK3CD","PLAC8","PLEK","PLEK","PRKCB","PRKCQ","PTPRC","PTPRCAP","PVRIG","SASH3","SELL","SELPLG","SEPT6","SH2D1A","SIRPG","SIT1","SLAMF1","SLAMF6","SPNS1","STAT4","STK10","STK10","TARP","TARP ","TARP ","TARP ","TBC1D10C","TBX21","TCL1A","TIGIT","TNFRSF4","TRAC","TRAC","TRAC ","TRAF3IP3","TRAJ17 ","TRAT1","TRAV20","TRBC1","TRBC1","TRBC1 ","TRBC2","TRGC2","TRGC2","TRGC2","VAMP5","XCL1","XCL1 ","XCL2"
    #"ACAP1","ARHGAP15","ARHGAP25","C1orf38","CCR2","CCR7","CD2","CD247","CD27","CD28","CD38","CD3D","CD3E","CD3G","CD40","CD48","CD52","CD52","CD53","CD6","CD86","CD8A","CFH","CFHR1","CLEC2D","CORO1A","CSF2RB","CST7","CYBB","DOCK2","DOK2","EVI2B","F5","FAM65B","GIMAP4","GIMAP5","GIMAP5","GIMAP6","GLYR1","SEPT6","GPR171","GPR18","GPSM3","GVINP1","GZMK","HCLS1","ICOS","IFFO1","IGHA1","IGHA2","IGHD","IGHG1","IGHG3","IGHG4","IGHM","IGHV4-31","LOC100133862","IGHG1","IGHM","LOC100133862","IGHM","IGHM","IGHM","LOC100133862","IGKV1-5","IGKV4-1","IGLL5","IGLV2-11","IL10RA","IL16","IL2RB","IL2RG","IL7R","INPP5D","IRF4","ITGA4","ITK","ITM2C","KLRD1","KLRK1","LAT","SPNS1","LAX1","LCK","LCK","LCP2","LILRB1","LPXN","LST1","LST1","LY9","MAL","MFNG","MGC29506","MS4A1","NCF4","NCKAP1L","PAX5","PIK3CD","PLAC8","PLEK","PLEK","PRKCB","PRKCQ","PTPRC","PTPRCAP","PVRIG","SASH3","SELL","SELPLG","SH2D1A","SIRPG","SIT1","SLAMF1","STAT4","STK10","STK10","TARP","TARP","TRGC2","TARP","TRGC2","TARP","TRGC2","TBX21","TCL1A","TNFRSF4","TRAC","TRAC","TRAC","TRAJ17","TRAV20","TRAF3IP3","TRAT1","TRBC1","TRBC1","TRBC1","TRBC2","VAMP5","XCL1","XCL1","XCL2"
    #"IL2RG","CXCR6","CD3D","CD2","ITGAL","TAGAP","CIITA","HLA-DRA","PTPRC","CXCL9","CCL5","NKG7","GZMA","PRF1","CCR5","CD3E","GZMK","IFNG","HLA-E","GZMB","PDCD1","SLAMF6","CXCL13","CXCL10","IDO1","LAG3","STAT1","CXCL11"
    
    #"PDGFRA", "PDGFRB","CSPG4", "ACTA2"
    "GPR34", "MME"
  ))

p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.res ,
                      y=statistic.gsam.cor.tpc ,
                      col=significant,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label == T),col="black",cex=0.65) +
  # geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,
  #                 nudge_x = 2.9, direction = "y", hjust = "left") + #, lwd=0.5
  # geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,
  #                 nudge_x = -2.9, direction = "y", hjust = "right")+ #, lwd=0.5
  scale_color_manual(values = c('TRUE'='gray60','FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)




plt <- results.out %>%
  #dplyr::filter(is.na(padj.glass.res) | (!is.na(padj.glass.res) & padj.glass.res  < 0.01) ) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.glass.res < 0.01 & lfcSE.glass.res < 0.6 & abs(log2FoldChange.glass.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c(
    # TIL
    #"ACAP1","ARHGAP15","ARHGAP25","ARHGAP30","ARHGAP9","ARHGAP9","C1orf38","CCL5","CCR2","CCR7","CD2","CD247","CD27","CD28","CD38","CD3D","CD3E","CD3G","CD40","CD48","CD52","CD52","CD53","CD6","CD79A","CD86","CD8A","CFH ","CFHR1","CLEC2D","CORO1A","CSF2RB","CST7","CYBB","DOCK11","DOCK2","DOK2","EVI2B","F5","FAM65B","FCRL3","FCRL5","FYB","GIMAP4","GIMAP5","GIMAP5","GIMAP6","GIMAP7","GLYR1 ","GPR171","GPR18","GPSM3","GVINP1","GZMK","HCLS1","HCST","ICOS","IFFO1","IGHA1 ","IGHA2 ","IGHD ","IGHG1 ","IGHG1 ","IGHG3 ","IGHG4 ","IGHM","IGHM","IGHM ","IGHM ","IGHM ","IGHV4-31 ","IGKV1-5","IGKV4-1","IGLL5 ","IGLV2-11","IKZF1","IKZF1","IL10RA","IL16","IL2RB","IL2RG","IL7R","IL7R","INPP5D","IRF4","ITGA4","ITK","ITM2C","KLHL6","KLRD1","KLRK1","LAT ","LAX1","LCK","LCK","LCP2","LILRB1","LOC100133862","LOC100133862","LOC100133862","LPXN","LST1","LST1","LY9","LY9","MAL","MFNG","MGC29506","MPEG1","MPEG1","MS4A1","MS4A1","MS4A6A","NCF4","NCKAP1L","NLRC3","P2RY8","PAG1","PAG1","PARVG","PAX5","PIK3CD","PLAC8","PLEK","PLEK","PRKCB","PRKCQ","PTPRC","PTPRCAP","PVRIG","SASH3","SELL","SELPLG","SEPT6","SH2D1A","SIRPG","SIT1","SLAMF1","SLAMF6","SPNS1","STAT4","STK10","STK10","TARP","TARP ","TARP ","TARP ","TBC1D10C","TBX21","TCL1A","TIGIT","TNFRSF4","TRAC","TRAC","TRAC ","TRAF3IP3","TRAJ17 ","TRAT1","TRAV20","TRBC1","TRBC1","TRBC1 ","TRBC2","TRGC2","TRGC2","TRGC2","VAMP5","XCL1","XCL1 ","XCL2"
    #"ACAP1","ARHGAP15","ARHGAP25","C1orf38","CCR2","CCR7","CD2","CD247","CD27","CD28","CD38","CD3D","CD3E","CD3G","CD40","CD48","CD52","CD52","CD53","CD6","CD86","CD8A","CFH","CFHR1","CLEC2D","CORO1A","CSF2RB","CST7","CYBB","DOCK2","DOK2","EVI2B","F5","FAM65B","GIMAP4","GIMAP5","GIMAP5","GIMAP6","GLYR1","SEPT6","GPR171","GPR18","GPSM3","GVINP1","GZMK","HCLS1","ICOS","IFFO1","IGHA1","IGHA2","IGHD","IGHG1","IGHG3","IGHG4","IGHM","IGHV4-31","LOC100133862","IGHG1","IGHM","LOC100133862","IGHM","IGHM","IGHM","LOC100133862","IGKV1-5","IGKV4-1","IGLL5","IGLV2-11","IL10RA","IL16","IL2RB","IL2RG","IL7R","INPP5D","IRF4","ITGA4","ITK","ITM2C","KLRD1","KLRK1","LAT","SPNS1","LAX1","LCK","LCK","LCP2","LILRB1","LPXN","LST1","LST1","LY9","MAL","MFNG","MGC29506","MS4A1","NCF4","NCKAP1L","PAX5","PIK3CD","PLAC8","PLEK","PLEK","PRKCB","PRKCQ","PTPRC","PTPRCAP","PVRIG","SASH3","SELL","SELPLG","SH2D1A","SIRPG","SIT1","SLAMF1","STAT4","STK10","STK10","TARP","TARP","TRGC2","TARP","TRGC2","TARP","TRGC2","TBX21","TCL1A","TNFRSF4","TRAC","TRAC","TRAC","TRAJ17","TRAV20","TRAF3IP3","TRAT1","TRBC1","TRBC1","TRBC1","TRBC2","VAMP5","XCL1","XCL1","XCL2"
    #"IL2RG","CXCR6","CD3D","CD2","ITGAL","TAGAP","CIITA","HLA-DRA","PTPRC","CXCL9","CCL5","NKG7","GZMA","PRF1","CCR5","CD3E","GZMK","IFNG","HLA-E","GZMB","PDCD1","SLAMF6","CXCL13","CXCL10","IDO1","LAG3","STAT1","CXCL11"
    
    # Bcell
    #"IGG","IGG1","CD27","CD38","CD78","CD138","CD319","IL-6","CD138","IGA","IGG","IGE","IGA1","IGG1","IGE1","CD20","CD27","CD40","CD80","PDCD1LG2"
    
    # MG
    #"TMEM119", "CD33",
    
    "GPR34", "MME"
    #Fcer1g, Gpr34, Adora3, C1qb, C3, Ly86, P2ry13, Tbxas1, and Tlr7
  ))


p2 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc  ,
                      col=significant,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label == T),col="black",cex=0.65) +
  # geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,
  #                 nudge_x = 3.1, direction = "y", hjust = "left" #, lwd=0.5
  # ) +
  # geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,
  #                 nudge_x = -3.1, direction = "y", hjust = "right" #, lwd=0.5
  # ) +
  scale_color_manual(values = c('TRUE'='gray60','FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 GLASS",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)



p1 + p2





## 2.11 Corrected LFC + vascular/angio ----



# signi in both corrected and uncorrected, but lfcSE outlier in corrected
# "TTLL10"     
# "MYO3A"     
# "AC090791.1" "AC090124.2" "AC090643.1" "AC009041.1" "LINC00514"  "KRT9"       "CKM"        "PI3"   
#'EGFR','MHMT','CD4', "CXCL12", "BLNK", "DDB2","RBP1", "PLXNB1",
#"SOX2", "NANOG",
#"CDKN2A", "CDKN2B", "APEX1", "NF1", "TP53", "CD40",
#"GSTM1", "SOCS2", "BTC", "FGFR3", 
#"OCT4", "NOS1",
#"POU5F1"

# presynapse 
synapse <- read.csv('/tmp/gProfiler_hsapiens_4-2-2021_4-39-17 PM.csv')



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c(#'CDK4','MDM2','GLI1','GLIS1',
    #"CPEB1",
    #"GABARAP","GABBR2","GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6","GABRB1","GABRB2","GABRB3","GABRD","GABRE","GABRG1","GABRG2","GABRG3","GABRP","GABRQ","GABRR1","GABRR2","GABRR3"
    #"GRIN1","GRIN2A","GRIN2B","GRIN2C","GRIN2D","GRIN3A","GRIN3B","GRM1","GRM2","GRM3","GRM4","GRM5","GRM6","GRM7"
    #"HTR1A","HTR1B","HTR1D","HTR1E","HTR1F","HTR2A","HTR2B","HTR2C","HTR3A","HTR3B","HTR3C","HTR3D","HTR3E","HTR4","HTR5A","HTR6","HTR7","HTT"
    #"PCDH15","PCDH17","PCDH8","PCDHB10","PCDHB11","PCDHB13","PCDHB14","PCDHB16","PCDHB2","PCDHB3","PCDHB4","PCDHB5","PCDHB6","PCDHB9"
    #"SLC12A4","SLC12A5","SLC12A6","SLC12A7","SLC16A1","SLC16A3","SLC17A5","SLC17A6","SLC17A7","SLC17A8","SLC18A1","SLC18A2","SLC18A3","SLC1A1","SLC1A2","SLC1A3","SLC1A4","SLC1A6","SLC1A7","SLC22A1","SLC22A2","SLC22A3","SLC29A1","SLC29A2","SLC29A4","SLC2A1","SLC2A4","SLC2A8","SLC30A1","SLC30A3","SLC32A1","SLC3A2","SLC40A1","SLC4A10","SLC4A7","SLC4A8","SLC5A7","SLC6A1","SLC6A11","SLC6A17","SLC6A2","SLC6A3","SLC6A4","SLC6A5","SLC6A6","SLC6A9","SLC8A1","SLC8A2","SLC8A3","SLC9A6","SLC9B2"
    #"ANKRD1","ANKRD18A","ANKRD18B","ANKRD29","ANKRD30B"
    #"SYTL1","SYT6","SYT11","SYT2","SYT14","SYTL3","SYT8","SYT9","SYT13","SYT7","SYT12","SYTL2","SYT10","SYT1","SYT16","SYT17","SYT4","SYT3","SYT5","SYTL5","SYTL4"
    # TNNT1:HP:0000707	Abnormality of the nervous system
    #genesets$CELLTYPE$"Fan Embryonic Ctx Microglia 1"
    #genesets$CELLTYPE$`Fan Embryonic Ctx Big Groups Microglia`
    #genesets$CELLTYPE$`Fan Embryonic Ctx Big Groups Brain Immune`
    #genesets$CELLTYPE$`Fan Embryonic Ctx Opc` # Oligodendrocyte progenitor cells
    #genesets$CGP$`Nakayama Soft Tissue Tumors Pca2 U`
    #genesets$CGP$`Kobayashi Egfr Signaling 24hr Dn`
    #genesets$COMP$`54` # RB/EGFR overlap?
    #genesets$ONCOSIG$`Csr Late Up.v1 Up` # rb cell cycle overlap?
    #genesets$ONCOSIG$`Prc2 Ezh2 Up.v1 Dn`
    #genesets$IMMUSIG$`Day6 Vs Day10 Traf6ko Eff Cd8 Tcell Up`
    #genesets$IMMUSIG$`Gse14415 Natural Treg Vs Tconv Dn`
    
    genesets$GO_BP$`Vasculature Development`
    #genesets$GO_BP$`Blood Vessel Morphogenesis` 
    
  ))


genesets$GO_BP$`Vasculature Development` %>% length()
genesets$GO_BP$`Blood Vessel Development` %>% length()
genesets$GO_BP$`Blood Vessel Morphogenesis` %>% length()

ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res ,
                y=statistic.gsam.cor.tpc ,
                col=significant,
                label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  #geom_point(data=subset(plt, significant == T),cex=0.45) +
  #geom_point(data=subset(plt, significant != T),cex=0.35) +
  #geom_point(data=subset(plt, significant == T),cex=0.45) +
  #geom_point(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < -0.25),col="red",cex=0.65) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,
  #                nudge_x = 3.1, direction = "y", hjust = "left") + #, lwd=0.5
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,
  #                nudge_x = 3.1, direction = "y", hjust = "left")+ #, lwd=0.5
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)




plt <- results.out %>%
  #dplyr::filter(is.na(padj.glass.res) | (!is.na(padj.glass.res) & padj.glass.res  < 0.01) ) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.glass.res < 0.01 & lfcSE.glass.res < 0.6 & abs(log2FoldChange.glass.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) 
  dplyr::mutate(show.label = hugo_symbol %in% c(
    "TNNI3K","TNN","TNNT2","TNNI1","TNNC1","TNNI2","TNNT3","TNNT1","TNNC2"
  ))


p2 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc  ,
                      col=significant,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label == T),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,
                  nudge_x = 3.1, direction = "y", hjust = "left" #, lwd=0.5
  ) +
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,
                  nudge_x = -3.1, direction = "y", hjust = "right" #, lwd=0.5
  ) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 GLASS",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)



p1 + p2



ggsave("/tmp/gabra.png",height=10 * 1.3,width=4.5 * 1.3)




## 2.12 Corrected LFC + oncogenes ----



plt <- results.out %>%
  dplyr::mutate(log2FoldChange.gsam.res = NULL) %>% # force using corrected 
  dplyr::mutate(log2FoldChange.glass.res = NULL) %>% # force using corrected 
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%

  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(significant = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 & #padj.glass.tpc.res < 0.05
                  abs(log2FoldChange.glass.tpc.res) > 0.5 & 
                  direction.gsam.tpc.res == direction.glass.tpc.res #lfcSE.gsam.tpc.res < 0.3 &
  ) %>%
  dplyr::mutate(show.label.gains = hugo_symbol %in% c(
    "AKT1","AKT3","EGFR", "PDGFRA","MET", "PIK3C2B", "MDM2","MDM4", "CDK4",
    "CDK6","SOX2","FGFR3","MYCN","MYC","CCND1","CCND2","BMI1" # gains
  )) %>%
  dplyr::mutate(show.label.other = hugo_symbol %in% c( 
    "TP53", "PIK3CA","PIK3R1","NF1","SPTA1","GABRA6","ABCC6","CXCL12",
    "LTBP4", "TGFB1","PREX1","MSH6", "MSH2", "MLH1","VEGFA", "PTPN11",
    "STAT3","DCC", "MGMT", "PTCH1","VEGFA","GLI1"
  )) %>%
  dplyr::mutate(show.label.losses = hugo_symbol %in% c(
    "CDKN2A","CDKN2B","PTEN","RB1","NF1","QKI","CDKN2C","TP53","MTAP","ELAVL2","PARK2"
  ))


p1a <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.other == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.other == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.other == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.other == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Other genes of interest in GBM") +
  youri_gg_theme + xlim(-3, 3)
p1b <- ggplot(plt, aes(x=log2FoldChange.glass.tpc.res , y=statistic.glass.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.other == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.other == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.other == T & log2FoldChange.glass.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.other == T & log2FoldChange.glass.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Other genes of interest in GBM") +
  youri_gg_theme + xlim(-3, 3)


p2a <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.gains == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.gains == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.gains == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.gains == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Commonly gained genes in GBM") +
  youri_gg_theme + xlim(-3, 3)
p2b <- ggplot(plt, aes(x=log2FoldChange.glass.tpc.res , y=statistic.glass.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.gains == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.gains == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.gains == T & log2FoldChange.glass.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.gains == T & log2FoldChange.glass.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Commonly gained genes in GBM") +
  youri_gg_theme + xlim(-3, 3)



p3a <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.losses == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.losses == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.losses == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.losses == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Commonly lost genes in GBM") +
  youri_gg_theme + xlim(-3, 3)
p3b <- ggplot(plt, aes(x=log2FoldChange.glass.tpc.res , y=statistic.glass.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.losses == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.losses == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.losses == T & log2FoldChange.glass.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.losses == T & log2FoldChange.glass.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Commonly lost genes in GBM") +
  youri_gg_theme + xlim(-3, 3)


(p2a + p1a + p3a) / (p2b + p1b + p3b) 


ggsave("output/figures/figure-DGE_glioblastoma_relevant_genes.png", width=15,height=10)


## Figure S2 CDE-FGH? oncogenes ----


plt <- results.out %>%
  dplyr::mutate(log2FoldChange.gsam.res = NULL) %>% # force using corrected 
  dplyr::mutate(log2FoldChange.glass.res = NULL) %>% # force using corrected 
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(significant = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 & #padj.glass.tpc.res < 0.05
                  abs(log2FoldChange.glass.tpc.res) > 0.5 & 
                  direction.gsam.tpc.res == direction.glass.tpc.res #lfcSE.gsam.tpc.res < 0.3 &
  ) %>%
  dplyr::mutate(show.label.gains = hugo_symbol %in% c(
    "AKT1","AKT3","EGFR", "PDGFRA","MET", "PIK3C2B", "MDM2","MDM4", "CDK4",
    "CDK6","SOX2","FGFR3","MYCN","MYC","CCND1","CCND2","BMI1" # gains
  )) %>%
  dplyr::mutate(show.label.other = hugo_symbol %in% c( 
    "TP53", "PIK3CA","PIK3R1","NF1","SPTA1","GABRA6","ABCC6","CXCL12",
    "LTBP4", "TGFB1","PREX1","MSH6", "MSH2", "MLH1","VEGFA", "PTPN11",
    "STAT3","DCC", "MGMT", "PTCH1","VEGFA","GLI1",
    "NODAL",
    "PTPN11", # mut [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5627776/]
    "LTBP4", "MSH6", "PRDM2", "IGF1R"
  )) %>%
  dplyr::mutate(show.label.losses = hugo_symbol %in% c(
    "CDKN2A","CDKN2B","PTEN","RB1","NF1","QKI","CDKN2C","TP53","MTAP","ELAVL2","PARK2"
  ))


#PTPN11, LTBP4, MSH6, PRDM2, IGF1R



plt.expanded <- rbind(
  plt %>% 
    dplyr::rename(cor.tpc = statistic.gsam.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(show.label = show.label.gains) %>% 
    dplyr::mutate(marker.genes = "Common gains") %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.gsam.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(show.label = show.label.other) %>% 
    dplyr::mutate(marker.genes = "Other") %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.gsam.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(show.label = show.label.losses) %>% 
    dplyr::mutate(marker.genes = "Common losses") %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.glass.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.tpc.res) %>%
    dplyr::rename(show.label = show.label.gains) %>% 
    dplyr::mutate(marker.genes = "Common gains") %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.glass.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.tpc.res) %>%
    dplyr::rename(show.label = show.label.other) %>% 
    dplyr::mutate(marker.genes = "Other") %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.glass.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.glass.tpc.res) %>%
    dplyr::rename(show.label = show.label.losses) %>% 
    dplyr::mutate(marker.genes = "Common losses") %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
) %>%
  dplyr::mutate(status = case_when(show.label == F & significant == F ~ "other",
                                   show.label == F & significant == T ~ "significant",
                                   show.label == T & significant == F ~ "show label",
                                   show.label == T & significant == T ~ "label & show significant")) %>%
  dplyr::filter(!is.na(cor.tpc) & !is.na(log2FoldChange)) %>%
  dplyr::mutate(marker.genes = factor(marker.genes, levels = c("Common gains", "Other", "Common losses"))) %>%
  dplyr::mutate(limited = abs(log2FoldChange) > 2.5) %>%
  dplyr::mutate(log2FoldChange = ifelse(limited & log2FoldChange < 0, -2.5, log2FoldChange)) %>%
  dplyr::mutate(log2FoldChange = ifelse(limited & log2FoldChange > 0, 2.5, log2FoldChange)) %>% 
  dplyr::mutate(limited = as.character(limited))



ggplot(plt.expanded, aes(x = log2FoldChange, y = cor.tpc, shape=limited,
                         fill = status, col = status, label=hugo_symbol)) +
  facet_grid(cols = vars(marker.genes), rows = vars(dataset), scales = "free") +
  geom_point(data = subset(plt.expanded, status == "other"), size=1.8, col="#00000044") +
  geom_point(data = subset(plt.expanded, status == "significant"), size=2, col="black") +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  geom_text_repel(data=subset(plt.expanded, grepl("label",status) & log2FoldChange > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt.expanded, grepl("label",status) & log2FoldChange < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  geom_point(data = subset(plt.expanded, status == "show label"), size=2.75, col="black") +
  geom_point(data = subset(plt.expanded, status == "label & show significant"), size=2.75, col="black") +
  scale_shape_manual(values = c('TRUE' = 23, 'FALSE' = 21)) +
  scale_fill_manual(values = c("other" = "white", "significant" = "#eab509DD", "show label" = "#6ba6e5DD", "label & show significant" = "red")) + 
  scale_color_manual(values = c("other" = "white", "significant" = "#eab509DD", "show label" = "#6ba6e5DD", "label & show significant" = "red")) + 
  youri_gg_theme +
  xlim(-3, 3) +
  labs(x="log2FC R1. vs R2 (unpaired; tumour-% corrected)", y="Correlation t-statistic with tumour-%")


ggsave("output/figures/paper_dge_gains_losses_other.png", width=10,height=7.5)




w## Figure S2 AB chr7, chr10, PTEN + EGFR ----



plt <- results.out %>%
  dplyr::mutate(mark.chr7 = chr == "chr7") %>%
  dplyr::mutate(mark.chr10 = chr == "chr10") %>%
  dplyr::mutate(show.label.chr7 = grepl("EGFR",hugo_symbol)) %>%
  dplyr::mutate(show.label.chr10 = grepl("PTEN|MGMT",hugo_symbol)) # TET1 really interesting

# plt <- plt %>% 
#   dplyr::mutate(show.label.chr10 = chr=="chr10" & statistic.gsam.cor.tpc > 2.5)


plt.expanded <- rbind(
  plt %>% 
    dplyr::mutate(status = case_when(show.label.chr7 == T ~ "chr7 label", mark.chr7 == T ~ "chr7", TRUE ~ "other")) %>% 
    dplyr::mutate(panel = "chr7 genes") %>%
    dplyr::select(hugo_symbol, statistic.gsam.cor.tpc, log2FoldChange.gsam.res, status, panel, chr)
  ,
  plt %>% 
    dplyr::mutate(status = case_when(show.label.chr10 == T ~ "chr10 label", mark.chr10 == T ~ "chr10", TRUE ~ "other")) %>% 
    dplyr::mutate(panel = "chr10 genes") %>%
    dplyr::select(hugo_symbol, statistic.gsam.cor.tpc, log2FoldChange.gsam.res, status, panel, chr)
) %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc) & !is.na(log2FoldChange.gsam.res)) %>%
  dplyr::mutate(limited = abs(log2FoldChange.gsam.res) > 2.5) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(limited & log2FoldChange.gsam.res < 0, -2.5, log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(limited & log2FoldChange.gsam.res > 0, 2.5, log2FoldChange.gsam.res)) %>%
  dplyr::mutate(panel = factor(panel, levels = c('chr7 genes', 'chr10 genes') ))



ggplot(plt.expanded, aes(x = log2FoldChange.gsam.res, y = statistic.gsam.cor.tpc, shape=limited,
                         fill = status, col = status, label=hugo_symbol)) +
  facet_grid(cols = vars(panel), scales = "free") +
  geom_point(data = subset(plt.expanded, status == "other"), size=1.8, col="#00000044") +
  geom_point(data = subset(plt.expanded, status %in% c("chr7", "chr10")), size=2.5, col="black") +
  geom_point(data = subset(plt.expanded, grepl("label", status)), size=2.75, fill="red", col="black") +
  geom_text_repel(data=subset(plt.expanded, status == "chr7 label"),  size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  geom_text_repel(data=subset(plt.expanded, status == "chr10 label"),  size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  scale_shape_manual(values = c('TRUE' = 23, 'FALSE' = 21)) +
  scale_fill_manual(values=c('chr7'='#6ba6e5AA','chr10'='#eab509AA','other'='white' , "chr10 label"='red',"chr7 label"='red')) + 
  scale_color_manual(values=c('chr7'='#6ba6e5AA','chr10'='#eab509AA','other'='white' , "chr10 label"='red',"chr7 label"='red')) + 
  youri_gg_theme +
  xlim(-2.75, 2.75)


ggsave("output/figures/paper_dge_chr7_chr10_geiser_uncorrected.png",width=10,height=5)



## 2.13 Uncorrected LFC + cell types ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.glass.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2, 2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2, -2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 3, 3 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -3, -3 , log2FoldChange.glass.res)) %>%
  
  # dplyr::mutate(show.label.neurons = hugo_symbol %in% c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3","GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST","VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5","NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK","ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A","RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87","ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2","DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2","CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4","ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2","CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1","GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2","RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")) %>%
  # dplyr::mutate(show.label.endothelial = hugo_symbol %in% c("APOLD1","FLT1","RGS5","PTPRB","TM4SF1","ABCB1","ITM2A","SDPR","SLCO1A2","FN1","EMCN","ESAM","NOSTRIN","CD34","SLC38A5","CYYR1","PODXL","CDH5","VWF","MECOM","CD93","ABCG2","TEK","PALMD","ERG","CLDN5","PECAM1","KDR","ITGA1","ICAM2","ATP10A","ANXA3","CA4","MYCT1","GIMAP6","ANXA1","PTRF","KIAA1462","EBF1","HMCN1","ENG","IGFBP7","ARHGAP29","ANXA2","OCLN","HIGD1B","SLC2A1","GNG11","SLC19A3","EPAS1","TBX3","SRGN","SOX7","SLC16A4","CAV1","CLIC5","VIM","HEG1","CCDC141","C10ORF10","EDN1","ROBO4","TMEM204","PROM1","IFITM1","LEF1","COBLL1","WWTR1","HBB","ETS1","SLC39A8","COL4A1","OSMR","ADCY4","TIE1","EDN3","THBD","BSG","AHNAK","MYO1B","IL1R1","CXCL12","CLEC14A","GATA2","SGPP2","SHE","PLTP","SPARC","ACVRL1","MMRN2","NID1","TNFSF10","FOXC1","UACA","CGNL1","MFSD2A","NET1","ABCC9","FLI1","C1ORF54")) %>%
  # dplyr::mutate(show.label.microglia = hugo_symbol %in% c("CCL4","CCL3","CSF1R","CX3CR1","P2RY12","C1QB","RGS1","GPR183","GPR34","CTSS","LAPTM5","CD53","IL1A","C3AR1","PLEK","FCGR2A","CD83","ITGAM","P2RY13","CD86","TREM2","TYROBP","FCER1G","NCKAP1L","SELPLG","SLC2A5","CD14","C1QC","C1QA","MPEG1","HAVCR2","PTAFR","LY86","AIF1","ALOX5AP","LPCAT2","SLA","PTPRC","FCGR1A","CCL2","BLNK","IL10RA","BCL2A1","C5AR1","RHOH","CD84","CSF3R","TLR7","TLR2","HPGDS","LCP1","CD300A","FYB","MRC1","FAM105A","IRF8","LCP2","RGS10","CD74","PTPN6","TBXAS1","LYZ","DOCK2","TMEM119","NLRP3","ARHGDIB","CCRL2","IKZF1","ARHGAP25","DOCK8","HEXB","THEMIS2","SAMSN1","HK2","PLD4","APBB1IP","ITGB2","RUNX1","SLCO2B1","TLR1","FGD2","HCLS1","GPR84","OLFML3","MAFB","PIK3CG","SIGLEC7","IL1B","PIK3R5","IL6R","CXCL16","CLEC4A","PTGS1","SUSD3","LYN","VAV1","SLC11A1","RBM47","SYK","C10ORF128")[1:50]) %>%
  # dplyr::mutate(show.label.oligodendrocytes = hugo_symbol %in% c("PLP1","MOBP","CLDN11","MBP","UGT8","ERMN","MOG","MAG","OPALIN","CNP","MAL","GPR37","TF","MYRF","GJB1","ASPA","ENPP2","BCAS1","LPAR1","FA2H","ENPP6","APOD","CNTN2","CRYAB","KLK6","ERBB3","ANLN","SEPT4","PLEKHB1","TMEFF2","ST18","PTGDS","PEX5L","SLAIN1","QDPR","PLLP","TMEM125","HHIP","LGI3","TUBB4A","PLEKHH1","S1PR5","MAP6D1","GSN","EVI2A","EDIL3","CMTM5","GJC3","CA14","NFASC","TPPP","TMEM88B","TRIM59","CDH19","APLP1","NIPAL4","ADAMTS4","STMN4","S100B","CA2","PRR18","OLIG1","FOLH1","NINJ2","NDRG1","SLC24A2","SGK2","GALNT6","KCNA1","SH3TC2","TTLL7","SH3GL3","DOCK5","SCD","FEZ1","SLC44A1","RHOU","PPP1R16B","TSPAN2","C10ORF90","TNFAIP6","NKAIN2","MOB3B","PRKCQ","PPP1R14A","PLA2G16","DBNDD2","CDK18","PCDH9","ANO4","AGPAT4","OMG","FGFR2","TMEM63A","GLTP","CCP110","PLEKHG3","RAB33A","PSAT1","ZNF536")) %>%
   dplyr::mutate(show.label.oligodendrocyte.progenitor.cells = hugo_symbol %in% c("PDGFRA","TNR","PCDH15","SHC4","VCAN","LHFPL3","NEU4","GPR17","PTPRZ1","OLIG1","MMP16","DSCAM","C8ORF46","SEMA5A","MATN4","UGT8","GRIA3","CNTN1","BCAS1","SULF2","LUZP2","GJC3","NXPH1","APOD","MEGF11","LRRTM3","BRINP3","GALNT13","GRIA4","MYT1","SUSD5","LRRN1","SOX10","PRKCQ","SOX6","ITGB8","TMEM255A","GFRA1","RLBP1","PNLIP","XYLT1","GPSM2","TMEM255B","SEZ6L","STK32A","C14ORF37","LPPR5","SEMA3D","CSPG4","CSMD3","TMEM132B","SCRG1","KCNH8","CACNG4","UGDH","DPP6","BCAT1","PLLP","ERBB3","RNF43","S100B","SORCS1","OLIG2","CHRNA4","KCNJ16","PPAPDC1A","CSMD1","OPCML","PRKG2","COBL","FIGN","ACAN","TGFA","NLGN1","SLC6A13","EMID1","CHST6","TMEM100","GAL3ST1","EDIL3","KCNJ10","SLITRK3","SNTG1","CSPG5","ERBB4","SLC35F1","B3GAT2","C1QL1","SERINC5","CKAP2","LRRTM4","DPYD","SLITRK1","NCALD","CALCRL","SPP1","ZNF488","ADAM12","SULF1","HAS2")) %>% # grote overlap met tumor cel / gbm subtype genen?
  # dplyr::mutate(show.label.astrocytes = hugo_symbol %in% c("AQP4","GJA1","GJB6","SLC4A4","SLC1A2","F3","BMPR1B","FGFR3","SLC39A12","CLDN10","DIO2","ALDOC","ALDH1L1","SLC1A3","CLU","ATP13A4","SLCO1C1","SLC14A1","CHRDL1","GPR37L1","ACSBG1","ATP1A2","SLC25A18","EDNRB","PPAP2B","GFAP","SOX9","SDC4","PPP1R3C","NCAN","MLC1","GLI3","SLC7A11","ACSL6","RFX4","ID4","AGT","SFXN5","GABRG1","PAX6","RORB","GRM3","PTPRZ1","PSD2","SLC6A11","ATP1B2","NTSR2","S1PR1","SLC15A2","ELOVL2","TRIL","SCARA3","MGST1","KIAA1161","FAM107A","BCAN","SPARCL1","NWD1","NTRK2","SLC7A10","SCG3","ACOT11","KCNN3","MFGE8","RANBP3L","GPC5","EZR","ADHFE1","GABRB1","TMEM47","PAMR1","CPE","FABP7","LIX1","SLC13A5","IL33","SLC7A2","EGFR","PREX2","NDRG2","DTNA","ABCD2","HEPACAM","RGS20","ARHGEF26","GPAM","CHI3L1","ADCYAP1R1","GDPD2","SLC1A4","POU3F2","ETNPPL","MEGF10","MT3","TTYH1","PRODH","PLCD4","DDAH1","LGR4","HTRA1")) %>%
  
  dplyr::mutate(show.label.astrocytes = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == 'astrocyte') %>%
  dplyr::mutate(show.label.endothelial = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == 'endothelial') %>%
  dplyr::mutate(show.label.microglia = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == 'microglia/TAM') %>%
  dplyr::mutate(show.label.neurons = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == 'neuron') %>%
  dplyr::mutate(show.label.oligodendrocytes = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == 'oligodendrocyte')



plt <- plt %>% dplyr::mutate(show.label = show.label.neurons)
p1a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Neuron marker genes") +
  youri_gg_theme + xlim(-2.5, 2.5)
p1b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Neuron marker genes") +
  youri_gg_theme + xlim(-3.5, 3.5)



plt <- plt %>% dplyr::mutate(show.label = show.label.endothelial)
p2a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Endothelial marker genes") +
  youri_gg_theme + xlim(-2.5 ,2.5)
p2b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Endothelial marker genes") +
  youri_gg_theme + xlim(-3.5, 3.5)



plt <- plt %>% dplyr::mutate(show.label = show.label.microglia)
p3a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="TAM/Microglia marker genes") +
  youri_gg_theme + xlim(-2.5, 2.5)
p3b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="TAM/Microglia marker genes") +
  youri_gg_theme + xlim(-3.5, 3.5)



plt <- plt %>% dplyr::mutate(show.label = show.label.oligodendrocytes)
p4a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Oligodendrocyte marker genes") +
  youri_gg_theme + xlim(-2.5, 2.5)
p4b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Oligodendrocyte marker genes") +
  youri_gg_theme + xlim(-3.5, 3.5)


plt <- plt %>% dplyr::mutate(show.label = show.label.astrocytes)
p5a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Astrocyte marker genes") +
  youri_gg_theme + xlim(-2.5, 2.5)
p5b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Astrocyte marker genes") +
  youri_gg_theme + xlim(-3.5, 3.5)




plt <- plt %>% dplyr::mutate(show.label = chr == "chr7" & baseMean.gsam.res > 300 & baseMean.glass.res > 300)
p6a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker.chr7 == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker.chr7 == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 

  #geom_label_repel(data=subset(plt, statistic.gsam.cor.tpc > 7 &  chr == "chr7" & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  #geom_label_repel(data=subset(plt, statistic.gsam.cor.tpc > 7 &  chr == "chr7" & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="chr7 genes (chr7 amplification marker; baseMean > 300)") +
  youri_gg_theme + xlim(-2.5, 2.5)
p6b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & show.marker.chr7 == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 , nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T & show.marker.chr7 == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="chr7 genes (chr7 amplification marker; baseMean > 300)") +
  youri_gg_theme + xlim(-3.5, 3.5)




# plt <- plt %>% dplyr::mutate(show.label = show.label.oligodendrocyte.progenitor.cells)
# p6a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
#   geom_point(data=subset(plt, show.label == F),cex=0.35) +
#   geom_smooth(data = subset(plt, padj.glass.res > 0.05 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col="blue" , size=0.4) +
#   geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
#   scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
#   labs(x = "log2FC R1 vs. R2 [G-SAM]",
#        y="Corr t-stat tumour-%"
#        ,col="Oligod. Prog. marker genes") +
#   youri_gg_theme + xlim(-2, 2)
# p6b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
#   geom_point(data=subset(plt, show.label == F),cex=0.35) +
#   geom_smooth(data = subset(plt, padj.glass.res > 0.05 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col="blue" , size=0.4) +
#   geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
#   scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
#   labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
#        y="Corr t-stat tumour-%"
#        ,col="Oligod. Prog. marker genes") +
#   youri_gg_theme + xlim(-3, 3)


(p1a + p1b) /
(p2a + p2b) / 
(p3a + p3b) / 
(p4a + p4b) / 
(p5a + p5b) / 
(p6a + p6b)

ggsave("output/figures/paper_dge_cell-type_genes.png",width=2*6, height=6*6)


## 2.14 Corrplot marker genes ----


### A :: each patient separately ----


labels <- results.out %>% 
  dplyr::filter(!is.na(primary.marker.genes) | hugo_symbol %in% c(
    # "SLC35E3", "DACT2", "TLX2", "SEPT12",
    # "TLR10", "SCIN", "HAMP", "CXCR2",
    # "GLI1", "PRB2", "MAFA", "HAPLN1"
    #"BIRC5", "CD44", "DANCR", "EZH2", "HIF1A", "ID1", "ID2", "IGFBP2", "ITGA6", "MECAM", "MET",
    #"MYC", "NOS2", "OLIG2", "PDGFRA", "PDPN", "PI3", "POSTN", "PROM1", "TGFBR2", "TNFAIP3",
    #"MFAP2","COL12A1"
  )) %>%
  dplyr::filter(!is.na(gid)) %>%
  dplyr::filter(primary.marker.genes %in% c('CL subtype', 'MES subtype', 'PN subtype') == F ) %>%
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical)


plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL)


plt.tpc <- data.frame(sid = colnames(plt)) %>%
  dplyr::left_join(gsam.metadata.all %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
  dplyr::rename('tumor-% DNA' = tumour.percentage.dna) %>%
  tibble::column_to_rownames('sid') %>%
  t()


plt <- rbind(plt, plt.tpc)


labels <- labels %>% 
  `rownames<-`(NULL) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL)  %>%
  dplyr::mutate_all(is.na)


cor_cor_plot(plt, labels %>%
               dplyr::select(
                 `chr7 gained (tumor)`,
                 `tumor-% DNA`,
                 `astrocyte`,
                 `endothelial`,
                 `oligodendrocyte`,
                 `neuron`,
                 `TIL / T-cell`,
                 `microglia/TAM`
               )
             )



ggsave("output/figures/paper_dge_corrplot_expression_gene_per_patient.png",width = 1200 * 2, height = 900 * 2 ,units="px" )



## gli1/MGMT plots ----

# 
# 
# plt <- gsam.gene.expression.all.vst %>%
#   as.data.frame %>%
#   tibble::rownames_to_column('gid') %>%
#   dplyr::filter(grepl('GLI1|MGMT',gid)) %>%
#   tibble::column_to_rownames('gid') %>%
#   t() %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column('sid') %>%
#   dplyr::rename(GLI1.vst = `ENSG00000111087.10_4|GLI1|chr12:57853568-57866051(+)`) %>%
#   dplyr::rename(MGMT.vst = `ENSG00000170430.10_4|MGMT|chr10:131265454-131569247(+)`) %>%
#   dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
#   dplyr::mutate(resection = gsub("^...(.).*$","R\\1",sid)) %>%
#   dplyr::left_join(
#     gsam.patient.metadata %>%
#       dplyr::select(studyID, treatedWithTMZ, nCyclesTMZ,daysToLastCycleTMZ),
#     by = c('pid'='studyID')
#   ) %>%
#   dplyr::mutate(order =  match(1:nrow(.), order(resection, treatedWithTMZ, GLI1.vst)))
# 
# 
# ggplot(plt, aes(x = reorder(sid, order), y = GLI1.vst, col=treatedWithTMZ)) +
#   geom_point()
# 
# ggplot(plt, aes(x = resection, y = GLI1.vst, col=treatedWithTMZ, group=pid)) +
#   facet_grid(cols = vars(treatedWithTMZ), scales = "free") +
#   geom_point() +
#   geom_line()
# 
# ggplot(plt, aes(x = nCyclesTMZ, y = GLI1.vst, col=treatedWithTMZ, group=pid)) +
#   facet_grid(cols = vars(resection), scales = "free") +
#   geom_point() +
#   geom_line()
# 
# 
# ggplot(plt, aes(x = nCyclesTMZ, y = MGMT.vst, col=treatedWithTMZ, group=pid)) +
#   facet_grid(cols = vars(resection), scales = "free") +
#   geom_point() +
#   geom_line()
# 
# 
# ggplot(plt, aes(x = daysToLastCycleTMZ, y = GLI1.vst, col=treatedWithTMZ, group=pid)) +
#   facet_grid(cols = vars(resection), scales = "free") +
#   geom_point() +
#   geom_line()
# 
# ggplot(plt, aes(x = daysToLastCycleTMZ, y = MGMT.vst, col=treatedWithTMZ, group=pid)) +
#   facet_grid(cols = vars(resection), scales = "free") +
#   geom_point() +
#   geom_line()



### Figure 7|S7 LFC per gepaard sample ----


Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2", Pericyte)) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>%
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>%
  dplyr::mutate(`pericyte` = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>%
  dplyr::select(c("gid","hugo_symbol",
                  "oligodendrocyte",
                  'PN subtype',
                  "neuron",
                  'CL subtype',
                  "astrocyte",
                  "tumor-% DNA",
                  "chr7 gained (tumor)",
                  'pericyte',
                  'Extracellular Matrix',
                  'MES subtype',
                  "microglia/TAM",
                  "TIL / T-cell",
                  "endothelial",
  )) %>% 
  dplyr::rename(`NMF:V[1] (MES)` = `MES subtype`) %>% 
  dplyr::rename(`NMF:V[2] (CL)` = `CL subtype`) %>% 
  dplyr::rename(`NMF:V[3] (PN)` = `PN subtype`) %>% 
  dplyr::rename(`OD` = `oligodendrocyte`) %>% 
  dplyr::rename(`NE` = `neuron`) %>% 
  dplyr::rename(`AC` = `astrocyte`) %>% 
  dplyr::rename(`EN` = `endothelial`) %>%
  dplyr::rename(`PE` = `pericyte`) %>% 
  dplyr::rename(`TIL/TC` = `TIL / T-cell`) %>% 
  dplyr::rename(`TAM/MG` = `microglia/TAM`) %>% 
  dplyr::rename(`ECM` = `Extracellular Matrix`) %>% 
  dplyr::rename(`T% (DNA)` = `tumor-% DNA`) %>% 
  dplyr::rename(`chr7 (T)` = `chr7 gained (tumor)`)





plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(

    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



# Add NMF 1 - 3
# Add ECM signature
# Add EPIC macrophage score [not necessary]

odd <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 1)
even <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 0)



plt.r1 <- plt[,odd] 
plt.r2 <- plt[,even]
stopifnot ( gsub("^(...).*$","\\1",colnames(plt.r1)) == gsub("^(...).*$","\\1",colnames(plt.r2)) )
rm(plt)


plt <- log2(plt.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
            plt.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) )



plt <- rbind(plt,
  dplyr::full_join(
    data.frame(sid = colnames(plt.r1)) %>%
      dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
      dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
      dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna) %>%
      dplyr::mutate(sid = NULL) ,
    data.frame(sid = colnames(plt.r2)) %>%
      dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
      dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))%>%
      dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) %>%
      dplyr::mutate(sid = NULL) , by = c('pid'='pid'))  %>%
    dplyr::mutate(`tumor-% DNA` = log2(tumour.percentage.dna.R1 / tumour.percentage.dna.R2)) %>%
    dplyr::mutate(tumour.percentage.dna.R1 = NULL, tumour.percentage.dna.R2 = NULL) %>%
    tibble::column_to_rownames('pid') %>%
    t()
  )

stopifnot(labels$hugo_symbol == rownames(plt))


labels <- labels %>% 
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels, 6, 1)




#ggsave("output/figures/paper_dge_corrplot_logFc_gene_per_pair.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_logFc_gene_per_pair.pdf",width = 1200 * 2.4, height = 1200 * 2.4 ,units="px" )







### corrplot 'm + C4 + Pericytes [log2(r1/r2)] ----



Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(Pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(C4 = hugo_symbol %in% c(C4A,C4B)) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | Pericyte | C4 | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2")) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>% # 
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>% 
  dplyr::mutate(Pericyte = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>% 
  dplyr::mutate(C4 = ifelse(hugo_symbol %in% c(C4A,C4B), NA , F)) %>% 
  dplyr::select(c("gid","hugo_symbol",
                  "Pericyte",
                  "endothelial",
                  "tumor-% DNA",
                  "chr7 gained (tumor)",
                  "C4",
                  "microglia/TAM",
                  "TIL / T-cell",
                  'MES subtype',
                  'Extracellular Matrix',
                  'CL subtype',
                  "astrocyte",
                  "oligodendrocyte",
                  'PN subtype',
                  "neuron",
  ))





plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    
    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



# Add NMF 1 - 3
# Add ECM signature
# Add EPIC macrophage score [not necessary]

odd <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 1)
even <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 0)



plt.r1 <- plt[,odd] 
plt.r2 <- plt[,even]
stopifnot ( gsub("^(...).*$","\\1",colnames(plt.r1)) == gsub("^(...).*$","\\1",colnames(plt.r2)) )
rm(plt)


plt <- log2(plt.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
              plt.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) )



plt <- rbind(plt,
             dplyr::full_join(
               data.frame(sid = colnames(plt.r1)) %>%
                 dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
                 dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
                 dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna) %>%
                 dplyr::mutate(sid = NULL) ,
               data.frame(sid = colnames(plt.r2)) %>%
                 dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
                 dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))%>%
                 dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) %>%
                 dplyr::mutate(sid = NULL) , by = c('pid'='pid'))  %>%
               dplyr::mutate(`tumor-% DNA` = log2(tumour.percentage.dna.R1 / tumour.percentage.dna.R2)) %>%
               dplyr::mutate(tumour.percentage.dna.R1 = NULL, tumour.percentage.dna.R2 = NULL) %>%
               tibble::column_to_rownames('pid') %>%
               t()
)

stopifnot(labels$hugo_symbol == rownames(plt))


labels <- labels %>% 
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels)




ggsave("output/figures/paper_dge_corrplot_C4_logFc_per_pair.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_C4_logFc_per_pair.pdf",width = 1200 * 2, height = 900 * 2 ,units="px" )



### corrplot 'm + C4 + Pericytes [vst expr] ----
Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(Pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(C4 = hugo_symbol %in% c(C4A,C4B)) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | Pericyte | C4 | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2")) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>% # 
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>% 
  dplyr::mutate(Pericyte = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>% 
  dplyr::mutate(C4 = ifelse(hugo_symbol %in% c(C4A,C4B), NA , F)) %>% 
  dplyr::select(c("gid","hugo_symbol",
                  "endothelial",
                  'CL subtype',
                  "astrocyte",
                  "tumor-% DNA",
                  "chr7 gained (tumor)",
                  "oligodendrocyte",
                  "neuron",
                  'PN subtype',
                  "Pericyte",
                  'MES subtype',
                  'Extracellular Matrix',
                  "TIL / T-cell",
                  "microglia/TAM",
                  "C4"
  ))



# labels <- labels %>% dplyr::filter(!is.na(endothelial))



plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    
    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')) %>% dplyr::rename(`tumor-% DNA`='tumour.percentage.dna'), by=c('sid'='sid')) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



labels <- labels %>%
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels)




ggsave("output/figures/paper_dge_corrplot_C4_vst_per_resection.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_C4_vst_per_resection.pdf",width = 1200 * 2, height = 900 * 2 ,units="px" )



### corrplot 'm + C5 + Pericytes [log2(r1/r2)] ----


Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(Pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(C5 = hugo_symbol %in% C5) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | Pericyte | C5 | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2")) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>% # 
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>% 
  dplyr::mutate(Pericyte = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>% 
  dplyr::mutate(C5 = ifelse(hugo_symbol %in% C5, NA , F)) %>% 
  dplyr::select(c("gid","hugo_symbol",
                  "oligodendrocyte",
                  'PN subtype',
                  "neuron",
                  'CL subtype',
                  "astrocyte",
                  "tumor-% DNA",
                  "chr7 gained (tumor)",
                  "TIL / T-cell",
                  "microglia/TAM",
                  'MES subtype',
                  'Extracellular Matrix',
                  "C5",
                  "Pericyte",
                  "endothelial"
  ))





plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    
    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



# Add NMF 1 - 3
# Add ECM signature
# Add EPIC macrophage score [not necessary]

odd <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 1)
even <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 0)



plt.r1 <- plt[,odd] 
plt.r2 <- plt[,even]
stopifnot ( gsub("^(...).*$","\\1",colnames(plt.r1)) == gsub("^(...).*$","\\1",colnames(plt.r2)) )
rm(plt)


plt <- log2(plt.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
              plt.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) )



plt <- rbind(plt,
             dplyr::full_join(
               data.frame(sid = colnames(plt.r1)) %>%
                 dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
                 dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
                 dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna) %>%
                 dplyr::mutate(sid = NULL) ,
               data.frame(sid = colnames(plt.r2)) %>%
                 dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
                 dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))%>%
                 dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) %>%
                 dplyr::mutate(sid = NULL) , by = c('pid'='pid'))  %>%
               dplyr::mutate(`tumor-% DNA` = log2(tumour.percentage.dna.R1 / tumour.percentage.dna.R2)) %>%
               dplyr::mutate(tumour.percentage.dna.R1 = NULL, tumour.percentage.dna.R2 = NULL) %>%
               tibble::column_to_rownames('pid') %>%
               t()
)

stopifnot(labels$hugo_symbol == rownames(plt))


labels <- labels %>% 
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels)




ggsave("output/figures/paper_dge_corrplot_C5_logFc_per_pair.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_C5_logFc_per_pair.pdf",width = 1200 * 2, height = 900 * 2 ,units="px" )



### corrplot 'm + C5 + Pericytes [vst expr] ----



Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(Pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(C5 = hugo_symbol %in% C5) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | Pericyte | C5 | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2")) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>% # 
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>% 
  dplyr::mutate(Pericyte = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>% 
  dplyr::mutate(C5 = ifelse(hugo_symbol %in% C5, NA , F)) %>% 
  dplyr::select(c("gid","hugo_symbol",
                  "chr7 gained (tumor)",
                  "tumor-% DNA",
                  "astrocyte",
                  'CL subtype',
                  "endothelial",
                  
                  "oligodendrocyte",
                  "neuron",
                  'PN subtype',
                  "microglia/TAM",
                  "TIL / T-cell",
                  
                  'Extracellular Matrix',
                  'MES subtype',
                  "C5",
                  "Pericyte"
                  
  ))


# labels <- labels %>% dplyr::filter(!is.na(endothelial))



plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    
    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')) %>% dplyr::rename(`tumor-% DNA`='tumour.percentage.dna'), by=c('sid'='sid')) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



labels <- labels %>%
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels)




ggsave("output/figures/paper_dge_corrplot_C5_vst_per_resection.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_C5_vst_per_resection.pdf",width = 1200 * 2, height = 900 * 2 ,units="px" )



### corrplot 'm + C6 + Pericytes [log2(r1/r2)] ----


C6 <- c('CRABP2', 'CILP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
        'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
        "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
        "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
        "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")
Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(Pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(C6 = hugo_symbol %in% C6) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | Pericyte | C6 | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2")) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>% # 
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>% 
  dplyr::mutate(Pericyte = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>% 
  dplyr::mutate(C6 = ifelse(hugo_symbol %in% C6, NA , F)) %>% 
  dplyr::select(c("gid","hugo_symbol",
                  "tumor-% DNA",
                  "chr7 gained (tumor)",
                  "endothelial",
                  
                  "oligodendrocyte",
                  'PN subtype',
                  "neuron",
                  "astrocyte",
                  'CL subtype',
                  "C6",
                  'Extracellular Matrix',
                  'MES subtype',
                  "Pericyte",
                  "TIL / T-cell",
                  "microglia/TAM"
  ))





plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    
    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



# Add NMF 1 - 3
# Add ECM signature
# Add EPIC macrophage score [not necessary]

odd <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 1)
even <- 1:ncol(plt) %>% purrr::keep(~ . %% 2 == 0)



plt.r1 <- plt[,odd] 
plt.r2 <- plt[,even]
stopifnot ( gsub("^(...).*$","\\1",colnames(plt.r1)) == gsub("^(...).*$","\\1",colnames(plt.r2)) )
rm(plt)


plt <- log2(plt.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
              plt.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) )



plt <- rbind(plt,
             dplyr::full_join(
               data.frame(sid = colnames(plt.r1)) %>%
                 dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
                 dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
                 dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna) %>%
                 dplyr::mutate(sid = NULL) ,
               data.frame(sid = colnames(plt.r2)) %>%
                 dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
                 dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))%>%
                 dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) %>%
                 dplyr::mutate(sid = NULL) , by = c('pid'='pid'))  %>%
               dplyr::mutate(`tumor-% DNA` = log2(tumour.percentage.dna.R1 / tumour.percentage.dna.R2)) %>%
               dplyr::mutate(tumour.percentage.dna.R1 = NULL, tumour.percentage.dna.R2 = NULL) %>%
               tibble::column_to_rownames('pid') %>%
               t()
)

stopifnot(labels$hugo_symbol == rownames(plt))


labels <- labels %>% 
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels)




ggsave("output/figures/paper_dge_corrplot_C6_logFc_per_pair.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_C6_logFc_per_pair.pdf",width = 1200 * 2, height = 900 * 2 ,units="px" )



### corrplot 'm + C6 + Pericytes [vst expr] ----


C6 <- c('CRABP2', 'CILP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
        'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
        "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
        "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
        "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")
Pericyte <- c("RGS5","PDGFRB","CD248","HEYL","CFH")


labels <- results.out %>% 
  dplyr::mutate(Pericyte = hugo_symbol %in% Pericyte) %>% 
  dplyr::mutate(C6 = hugo_symbol %in% C6) %>% 
  dplyr::mutate(primary.marker.genes = ifelse( primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) %>%
  dplyr::filter((!is.na(primary.marker.genes) | Pericyte | C6 | hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2")) & !is.na(gid)) %>%
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) %>% # 
  dplyr::select(c('gid', 'hugo_symbol' , 'primary.marker.genes')) %>%
  rbind(data.frame('gid'="NMF:V[1] (MES)", hugo_symbol = "NMF:V[1] (MES)", primary.marker.genes = 'MES subtype')) %>%
  rbind(data.frame('gid'="NMF:V[2] (CL)", hugo_symbol = "NMF:V[2] (CL)", primary.marker.genes = 'CL subtype')) %>%
  rbind(data.frame('gid'="NMF:V[3] (PN)", hugo_symbol = "NMF:V[3] (PN)", primary.marker.genes = 'PN subtype')) %>%
  rbind(data.frame('gid'='tumor-% DNA', hugo_symbol = 'tumor-% DNA', primary.marker.genes = 'tumor-% DNA')) %>%
  reshape2::dcast (gid + hugo_symbol ~ primary.marker.genes, fill = F,fun.aggregate = as.logical,value.var= 'primary.marker.genes') %>%
  dplyr::mutate('Extracellular Matrix'=F) %>% 
  dplyr::mutate(`Extracellular Matrix` = ifelse(hugo_symbol %in% c("MFAP2","COL1A1","COL5A1","ADAMTS2"), NA , F)) %>% 
  dplyr::mutate(Pericyte = ifelse(hugo_symbol %in% Pericyte, NA , F)) %>% 
  dplyr::mutate(C6 = ifelse(hugo_symbol %in% C6, NA , F)) %>% 
  dplyr::select(c("gid","hugo_symbol",
                  "tumor-% DNA",
                  "chr7 gained (tumor)",
                  "astrocyte",
                  'CL subtype',
                  "endothelial",
                  
                  "oligodendrocyte",
                  "neuron",
                  'PN subtype',
                  "C6",
                  'Extracellular Matrix',
                  'MES subtype',
                  "Pericyte",
                  "TIL / T-cell",
                  "microglia/TAM"
  ))


# labels <- labels %>% dplyr::filter(!is.na(endothelial))



plt <- labels %>%
  dplyr::filter(gid != 'tumor-% DNA') %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::select(c('gid','hugo_symbol',colnames(gsam.gene.expression.all.paired))) %>% # <3
  dplyr::filter(gid %in% labels$gid) %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::filter(hugo_symbol %in% c('NMF:V[1] (MES)', 'NMF:V[2] (CL)', 'NMF:V[3] (PN)') == F)  %>% # lines below add the NMF matrices
  tibble::column_to_rownames('hugo_symbol') %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(
    
    readRDS("tmp/combi.gbm_nmf_150.new.Rds") %>%
      purrr::pluck('123456') %>%
      purrr::pluck('H') %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(c('NMF:V[1] (MES)','NMF:V[2] (CL)','NMF:V[3] (PN)')) %>% 
      tibble::rownames_to_column('sid') %>%
      dplyr::mutate(sid = gsub('^GSAM-','',sid))
    
    , by=c('sid' = 'sid')
  ) %>%
  dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')) %>% dplyr::rename(`tumor-% DNA`='tumour.percentage.dna'), by=c('sid'='sid')) %>%
  tibble::column_to_rownames('sid') %>%
  t() %>%
  as.data.frame()



labels <- labels %>%
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hugo_symbol')


cor_cor_plot(plt, labels)




ggsave("output/figures/paper_dge_corrplot_C6_vst_per_resection.png",width = 1200 * 2, height = 900 * 2 ,units="px" )
ggsave("output/figures/paper_dge_corrplot_C6_vst_per_resection.pdf",width = 1200 * 2, height = 900 * 2 ,units="px" )





## old random stuff? ----

a.cl = results.out %>%
  dplyr::filter(!is.na(TCGA.subtype.marker) | !is.na(GliTS.reduxsubtype.marker)) %>%
  dplyr::filter(TCGA.subtype.marker == "TCGA-CL" | GliTS.reduxsubtype.marker == "GliTS-CL") %>%
  dplyr::select(c(hugo_symbol, gid,TCGA.subtype.marker,  GliTS.reduxsubtype.marker)) 


a.mes = results.out %>%
  dplyr::filter(!is.na(TCGA.subtype.marker) | !is.na(GliTS.reduxsubtype.marker)) %>%
  dplyr::filter(TCGA.subtype.marker == "TCGA-MES" | GliTS.reduxsubtype.marker == "GliTS-MES") %>%
  dplyr::select(c(hugo_symbol, gid,TCGA.subtype.marker,  GliTS.reduxsubtype.marker)) 


a.pn = results.out %>%
  dplyr::filter(!is.na(TCGA.subtype.marker) | !is.na(GliTS.reduxsubtype.marker) | hugo_symbol %in% c('NCAM1', 'NRCAN', 'OLIG1', 'SOX4', 'SOX6') ) %>%
  dplyr::filter(TCGA.subtype.marker == "TCGA-PN" | GliTS.reduxsubtype.marker == "GliTS-PN" | hugo_symbol %in% c('NCAM1', 'NRCAN', 'OLIG1', 'SOX4', 'SOX6')) %>%
  dplyr::select(c(hugo_symbol, gid,TCGA.subtype.marker,  GliTS.reduxsubtype.marker)) 


a <- rbind(a.cl, a.mes, a.pn)
a <- a.mes

b = gsam.gene.expression.all.vst %>% 
  as.data.frame %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% a$gid) %>%
  dplyr::left_join(a, by=c('gid'='gid')) %>%
  dplyr::mutate(hugo_symbol =  ifelse(is.na(TCGA.subtype.marker) , paste0(hugo_symbol, "     "), paste0(hugo_symbol, " [T]")) ) %>%
  dplyr::mutate(hugo_symbol =  ifelse(is.na(GliTS.reduxsubtype.marker) , paste0(hugo_symbol, "     "), paste0(hugo_symbol, "[G]")) ) %>%
  tibble::column_to_rownames('hugo_symbol')  %>%
  dplyr::mutate(gid = NULL, TCGA.subtype.marker= NULL, GliTS.reduxsubtype.marker= NULL)


corrplot::corrplot(cor(t(b), method="pearson") %>% `colnames<-`(NULL), method = "circle",tl.cex=1, order='hclust')



# beste voor MES: COL1A1, LAMB1, RUNX1, S100A11
# beste voor CL:  ELOVL2, VAV3, SOCS2, MLC1
# beste voor PN: GABRB3, MAST1 GNG4, ERBB3, SOX10, 



plt <- data.frame(
  hugo_symbol = c(

  "CREB5",	"TRIM24",	"ETV1", "COA1", # tumor/chr7 gain
  "CD45",

  "CACHD1","AHCYL1","GPR37L1","BMPR1B", # astroctyes
  #"RBFOX3", "GABRB2", "SLC17A7","SST", # neuron
  "RBFOX3", "GABRB2","GABRA1","GABRG2","GABBR2",

  #"SSTR1","SSTR2","SSTR3","SSTR5", # Antibodies
  #"GABRA1","GABRA2","GABRB1","GABRB2",
  #"TNNT1", "TNNT2", "TNNT3",

  "PLP1", "OPALIN", "TMEM144","CLCA4", # oligodendrocyte
  #"PDGFA", "PDGFRA", "OLIG1", "OLIG2", "OLIG3",

  "TIE1","PEAR1","RGS5","NOSTRIN",  # endothelial
  "CD163",  "CD14", "C1QA","THEMIS2" # TAM/MG

  ,"CD4","CD2", "CD3D", "CD3E","CD8A" # t-cell? TIL

  ###, "MME", "ERG", "FCER2", "EPCAM", "EREG" << !!
  ,"EGFR"
  ,"EREG","AREG", "BTC","HBEGF","NGF","TGFA","EGF","EPGN",

  #"ARHGAP28","RHGEF26","BVES-AS1","BVES","CACNA2D1","CDH4","CNGA3","COL11A1","ELOVL2","ETV4","EVA1A","FGFR3","GNAI1","LFNG","LHFPL6","POPDC3","PPARGC1A","RFX4","RNF180","ROBO2","SLC24A3","SOCS2-AS1","SOCS2","SOX9","SULF1","TACR1","TAP1","VAV3"

  "TOP2A", "CDK1", "DTL", "CCNB1", "XRCC2", "CCNE2", "DSN1", "TIMELESS", # Cell Cycle genes Patel/Bernstein
  "VEGFA", "ADM", "TREM1", "ENO1", "LDHA", "NRN1", "UBC", "GBE1", "MIF" # Hypoxia genes Patel Bernstein
  )) %>%
  dplyr::left_join(results.out %>% 
      dplyr::select(c('gid','hugo_symbol','McKenzie_celltype_top_human_specificity','show.marker.chr7',TCGA.subtype.marker,  GliTS.reduxsubtype.marker)) %>%
      dplyr::rename(type=McKenzie_celltype_top_human_specificity) ,
    by=c('hugo_symbol'='hugo_symbol'))  %>%
  dplyr::mutate(type = ifelse(hugo_symbol %in% c('SLC17A7','SST'), 'neuron', type) ) %>%
  dplyr::mutate(type = ifelse(is.na(type), 'NA', type) ) %>%
  dplyr::mutate(type = ifelse(show.marker.chr7 & type != 'astrocyte' , 'chr7/gain' , type)) %>%
  dplyr::mutate(show.marker.chr7 = NULL) %>%
  dplyr::filter(!is.na(gid))  %>%
  dplyr::mutate(TCGA.subtype.marker = NULL) %>%
  dplyr::mutate(GliTS.reduxsubtype.marker = NULL)



tmp <- gsam.gene.expression.all.vst %>% 
  as.data.frame() %>%
  dplyr::select(colnames(gsam.gene.expression.all.paired)) %>% # <3
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% plt$gid) %>%
  tibble::column_to_rownames('gid')

odd <- 1:ncol(tmp) %>% purrr::keep(~ . %% 2 == 1)
even <- 1:ncol(tmp) %>% purrr::keep(~ . %% 2 == 0)

tmp.r1 <- tmp[,odd] 
tmp.r2 <- tmp[,even]
stopifnot ( gsub("^(...).*$","\\1",colnames(tmp.r1)) == gsub("^(...).*$","\\1",colnames(tmp.r2)) )
rm(tmp)


tmp <- plt %>% dplyr::left_join(
  log2(tmp.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
       tmp.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) ) %>%
    `rownames<-`(gsub("^ENSG.+\\|(.+)\\|chr.+$","\\1",rownames(.))) %>%
    tibble::rownames_to_column('hugo_symbol')
  , by = c('hugo_symbol'='hugo_symbol')) %>%
  dplyr::mutate(type = case_when(
    type == "chr7/gain" ~ "chr7", 
    type == "astrocyte" ~ "astr",
    type == "neuron" ~ "neur",
    type == "endothelial" ~ "endo",
    type == "microglia/TAM" ~ "TAM",
    type == "oligodendrocyte" ~ "olig",
    T ~ "?" )) %>%
  dplyr::mutate(hugo_symbol = paste0(hugo_symbol , " [" , type, "]") ) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(   gid    = NULL,  type=NULL)



tmp.2 <- dplyr::full_join(
  data.frame(sid = colnames(tmp.r1)) %>%
    dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
    dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
    dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna) %>%
    dplyr::mutate(sid = NULL) ,
  data.frame(sid = colnames(tmp.r2)) %>%
    dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
    dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))%>%
  dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) %>%
  dplyr::mutate(sid = NULL) , by = c('pid'='pid'))  %>%
  dplyr::mutate(`tumor-% DNA` = log2(tumour.percentage.dna.R1 / tumour.percentage.dna.R2)) %>%
  dplyr::mutate(tumour.percentage.dna.R1 = NULL, tumour.percentage.dna.R2 = NULL) %>%
  tibble::column_to_rownames('pid') %>%
  t()


stopifnot(colnames(tmp) == colnames(tmp.2))
tmp <- rbind(tmp.2, tmp)

# #"ward.D", "single", "complete", "average", "mcquitty",  "median", "centroid", "ward.D2"
# # volgorde voor plotten:
# #h = hclust(as.dist(1 - abs(cor(t(tmp)))),method = "ward.D")
# h = hclust(as.dist(1 - cor(t(tmp))))
# plot(h)
# o = h$labels[h$order]
# 
# tmp <- t(tmp) %>% as.data.frame %>% select(o) %>% t()


tmp.3 <- tmp %>% 
          t() %>%
        cor(method="pearson") %>%
      `colnames<-`(NULL)


#png(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient.png", width = 1200 * 0.8, height = 900 *0.8 )
corrplot::corrplot(tmp.3, method = "circle",tl.cex=1, order='hclust')
#dev.off()






## 2.16 Corrected LFC + individual tophits ----

# stap 1 - maak de selectie (en visualiseer deze)


plt <- results.out %>%
  dplyr::filter( !is.na(log2FoldChange.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(log2FoldChange.glass.tpc.res)  ) %>%
  dplyr::filter( !is.na(padj.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(padj.glass.tpc.res)  ) %>%
  
  dplyr::left_join(readRDS("tmp/clusters.Rds") %>%
                     tibble::rownames_to_column('hugo_symbol'), by=c('hugo_symbol'='hugo_symbol')) %>%
  
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2.5, 2.5 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2.5, -2.5 , log2FoldChange.gsam.tpc.res)) %>%
  
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(log2FoldChange.glass.tpc.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 2.5, 2.5 , log2FoldChange.glass.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res < -2.5, -2.5 , log2FoldChange.glass.tpc.res)) %>%
  
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") ) %>%
  
  dplyr::mutate(significant = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 & #padj.glass.tpc.res < 0.05
                  abs(log2FoldChange.glass.tpc.res) > 0.5 & 
                  direction.gsam.tpc.res == direction.glass.tpc.res #lfcSE.gsam.tpc.res < 0.3 &
                ) %>%
  dplyr::mutate(show.label = significant == T) 
  #%>%
  #dplyr::mutate(show.label = hugo_symbol %in% c("SLC35E3", "DACT2", "TLX2", "SEPT12", 
  #                                              "TLR10", "SCIN", "HAMP", "CXCR2",
  #                                              "GLI1", "PRB2", "MAFA", "HAPLN1")) %>%
  #dplyr::mutate(show.label = !is.na(C5) & C5 == T)


# a = plt %>%
#   dplyr::filter(significant == T) %>% 
#   dplyr::filter(log2FoldChange.gsam.tpc.res < 0 )



p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)") +
  youri_gg_theme + xlim(-3, 3)

p2 <- ggplot(plt, aes(x=log2FoldChange.glass.tpc.res ,
                      y=statistic.glass.cor.tpc  ,
                      col=show.label,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt,  show.label == F),cex=0.35) +
  geom_point(data=subset(plt, show.label == T),col="red",cex=0.65) +
  
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,                   nudge_x = 3.1, direction = "y", hjust = "left" ) +
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,                  nudge_x = -3.1, direction = "y", hjust = "right"  ) +
  
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 GLASS",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)

p1 / p2



ggsave("output/figures/figure-DGE_all_corrected.png",  width = 5.7 * 1.4 , height = 4 * 1.4 * 2)





### corrplot 'm (Real Figure 5?) ----




labels <- results.out %>%
  dplyr::filter( !is.na(log2FoldChange.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(log2FoldChange.glass.tpc.res)  ) %>%
  dplyr::filter( !is.na(padj.gsam.tpc.res)   ) %>%
  dplyr::filter( !is.na(padj.glass.tpc.res)  ) %>%
  
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
  dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") ) %>%
  
  dplyr::mutate(direction_up = ifelse(log2FoldChange.gsam.tpc.res > 0 , T, F) ) %>% 
  dplyr::mutate(direction_down = ifelse(log2FoldChange.gsam.tpc.res < 0 , T, F) ) %>% 

  dplyr::mutate(TCGA.CL = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-CL" , T, F) ) %>% 
  dplyr::mutate(TCGA.MES = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-MES" , T, F) ) %>% 
  dplyr::mutate(TCGA.PN = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-PN" , T, F) ) %>% 
  
  dplyr::mutate(Patel.Cell.cycle = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Cell cycle" , T, F) ) %>% 
  dplyr::mutate(Patel.Complete.Immune.response = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Complete/Immune response" , T, F) ) %>% 
  dplyr::mutate(Patel.Hypoxia = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Hypoxia" , T, F) ) %>% 
  
  dplyr::rename(Neftel.AC = neftel.meta.module.AC) %>% 
  dplyr::rename(Neftel.NPC1 = neftel.meta.module.NPC1) %>% 
  dplyr::rename(Neftel.NPC2 = neftel.meta.module.NPC2) %>% 
  dplyr::rename(Neftel.OPC = neftel.meta.module.OPC) %>% 
  dplyr::rename(Neftel.MES1 = neftel.meta.module.MES1) %>% 
  dplyr::rename(Neftel.MES2 = neftel.meta.module.MES2) %>% 
  
  dplyr::rename(`Extracellular Matrix` = EM.struct.constituent) %>%
  
  dplyr::mutate(significant = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 &
                  abs(log2FoldChange.glass.tpc.res) > 0.5 & 
                  direction.gsam.tpc.res == direction.glass.tpc.res) %>%
  
  dplyr::filter(significant == T) %>%
  dplyr::select(c('gid', 'hugo_symbol' , 'McKenzie_celltype_top_human_specificity', 'direction_up', 'direction_down',
                  TCGA.CL, TCGA.MES, TCGA.PN ,
                  Patel.Cell.cycle, Patel.Complete.Immune.response, Patel.Hypoxia, 
                  Neftel.AC, Neftel.OPC, Neftel.MES1, Neftel.MES2, Neftel.NPC1, Neftel.NPC2,
                  `Extracellular Matrix`
                  
                  ,significant,direction.gsam.tpc.res ,log2FoldChange.gsam.tpc.res,  direction.glass.tpc.res ,log2FoldChange.glass.tpc.res
                  ))
  #%>% dplyr::slice_head(n=50)
  # dplyr::filter(hugo_symbol %in%c("CRABP2",'CILP2','DPT','FGF7')) 





# add checks to check dcast
n_astrocyte <- labels %>%
  dplyr::filter(!is.na(`McKenzie_celltype_top_human_specificity`)) %>%
  dplyr::filter(`McKenzie_celltype_top_human_specificity` == 'astrocyte') %>% 
  nrow()
n_endothelial <- labels %>%
  dplyr::filter(!is.na(`McKenzie_celltype_top_human_specificity`)) %>%
  dplyr::filter(`McKenzie_celltype_top_human_specificity` == 'endothelial') %>% 
  nrow()
n_neuron <- labels %>%
  dplyr::filter(!is.na(`McKenzie_celltype_top_human_specificity`)) %>%
  dplyr::filter(`McKenzie_celltype_top_human_specificity` == 'neuron') %>% 
  nrow()
n_oligodendrocyte <- labels %>%
  dplyr::filter(!is.na(`McKenzie_celltype_top_human_specificity`)) %>%
  dplyr::filter(`McKenzie_celltype_top_human_specificity` == 'oligodendrocyte') %>% 
  nrow()
  
#labels %>%
  ##reshape2::dcast(hugo_symbol + direction_up + direction_down + gid ~ `McKenzie_celltype_top_human_specificity`, fill = NA,fun.aggregate = as.logical) %>%
  #dplyr::filter(grepl("GAB", hugo_symbol))


labels <- labels %>%
  reshape2::dcast(hugo_symbol +
                    direction_up + direction_down + 
                    TCGA.CL + TCGA.MES + TCGA.PN +
                    Neftel.AC + Neftel.NPC1 + Neftel.NPC2 + Neftel.OPC + Neftel.MES1 + Neftel.MES2 + 
                    Patel.Cell.cycle + Patel.Complete.Immune.response + Patel.Hypoxia +
                    `Extracellular Matrix` +
                    gid ~ `McKenzie_celltype_top_human_specificity`, fill = NA,fun.aggregate = as.logical) %>%
  dplyr::mutate('NA'=NULL) %>% 
  
  dplyr::mutate(direction_up = ifelse(direction_up , F, NA) ) %>% # F functions as T, NA as F, for dcast2
  dplyr::mutate(direction_down = ifelse(direction_down , F, NA) ) %>% # F functions as T, NA as F, for dcast2
  
  dplyr::mutate(TCGA.CL = ifelse(TCGA.CL , F, NA) ) %>%
  dplyr::mutate(TCGA.MES = ifelse(TCGA.MES , F, NA) ) %>%
  dplyr::mutate(TCGA.PN = ifelse(TCGA.PN , F, NA) ) %>%
  
  dplyr::mutate(Neftel.AC = ifelse(Neftel.AC , F, NA) ) %>%
  dplyr::mutate(Neftel.NPC1 = ifelse(Neftel.NPC1 , F, NA) ) %>%
  dplyr::mutate(Neftel.NPC2 = ifelse(Neftel.NPC2 , F, NA) ) %>%
  dplyr::mutate(Neftel.OPC = ifelse(Neftel.OPC , F, NA) ) %>%
  dplyr::mutate(Neftel.MES1 = ifelse(Neftel.MES1 , F, NA) ) %>%
  dplyr::mutate(Neftel.MES2 = ifelse(Neftel.MES2 , F, NA) ) %>%
  
  dplyr::mutate(`Extracellular Matrix` = ifelse(`Extracellular Matrix` , F, NA) ) %>%

  dplyr::mutate(Patel.Cell.cycle = ifelse(Patel.Cell.cycle , F, NA) ) %>%
  dplyr::mutate(Patel.Complete.Immune.response = ifelse(Patel.Complete.Immune.response , F, NA) ) %>%
  dplyr::mutate(Patel.Hypoxia = ifelse(Patel.Hypoxia , F, NA) ) %>%

  dplyr::mutate(TCGA.CL = NULL , TCGA.MES = NULL , TCGA.PN = NULL ,
                Patel.Cell.cycle = NULL , Patel.Complete.Immune.response = NULL , Patel.Hypoxia = NULL,
                C1 = NULL, C2 = NULL, C3 = NULL, C4 = NULL, C5 = NULL, C6 = NULL)



plt <- labels %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL)



labels <- labels %>% 
  `rownames<-`(NULL) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL) %>%
  dplyr::mutate_all(function(arg) { return (!is.na(arg)) }) %>%
  tibble::rownames_to_column('hugo_symbol') %>%
  dplyr::left_join(readRDS("tmp/clusters.Rds") %>%
    tibble::rownames_to_column('hugo_symbol'), by=c('hugo_symbol'='hugo_symbol')) %>%
  tibble::column_to_rownames('hugo_symbol')
  



stopifnot(labels %>% dplyr::filter(astrocyte == T) %>% nrow == n_astrocyte)
stopifnot(labels %>% dplyr::filter(endothelial == T) %>% nrow == n_endothelial)
stopifnot(labels %>% dplyr::filter(neuron == T) %>% nrow == n_neuron)
stopifnot(labels %>% dplyr::filter(oligodendrocyte == T) %>% nrow == n_oligodendrocyte)


# 
# cor_cor_plot(plt, labels, method="average")
# ggsave("output/figures/cor_ward.average.png", height=30, width=30) # then average, then ward.D1?
# 
# cor_cor_plot(plt, labels, method="complete")
# ggsave("output/figures/cor_ward.complete.png", height=30, width=30) # then average, then ward.D1?

# cor_cor_plot <- function(normalised_correlated_data, labels, font_scale , legend_scale , method="ward.D2")
h <- cor_cor_plot(plt, labels, 3,6, method="ward.D2")
#ggsave("output/figures/cor_ward.D2x.pdf", height=30, width=30) # then average, then ward.D1?
ggsave("output/figures/cor_ward.D2x.svg", height=30, width=30) # then average, then ward.D1?


export <- data.frame(hugo_symbol = rev(h$labels[h$order]),cluster=NA) 

c <- NA
for(i in 1:nrow(export)) {
  slice <- export[i,]
  if(slice$hugo_symbol == "PTPRN") {
    c <- "C1"
  } else if(slice$hugo_symbol == "TLX2") {
    c <- NA
  } else if(slice$hugo_symbol == "RASGRF1") {
    c <- "C2"
  } else if(slice$hugo_symbol == "VWF") {
    c <- "C3"
  } else if(slice$hugo_symbol == "NPY2R") {
    c <- NA
  } else if(slice$hugo_symbol == "SOD3") {
    c <- "C4"
  } else if(slice$hugo_symbol == "PRF1") {
    c <- "C5"
  } else if(slice$hugo_symbol == "CRABP2") {
    c <- "C6"
  }
  
  export[i,]$cluster <- c
}

stopifnot(nrow(export) == 422)


export <- export %>%  dplyr::left_join(
  results.out %>%
    dplyr::filter( !is.na(log2FoldChange.gsam.tpc.res)   ) %>%
    dplyr::filter( !is.na(log2FoldChange.glass.tpc.res)  ) %>%
    dplyr::filter( !is.na(padj.gsam.tpc.res)   ) %>%
    dplyr::filter( !is.na(padj.glass.tpc.res)  ) %>%
    
    dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) %>%
    dplyr::mutate(direction.glass.tpc.res = ifelse(log2FoldChange.glass.tpc.res > 0 , "up", "down") ) %>%
    
    dplyr::mutate(direction_up = ifelse(log2FoldChange.gsam.tpc.res > 0 , T, F) ) %>% 
    dplyr::mutate(direction_down = ifelse(log2FoldChange.gsam.tpc.res < 0 , T, F) ) %>% 
    
    dplyr::mutate(TCGA.CL = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-CL" , T, F) ) %>% 
    dplyr::mutate(TCGA.MES = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-MES" , T, F) ) %>% 
    dplyr::mutate(TCGA.PN = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-PN" , T, F) ) %>% 
    
    dplyr::mutate(Patel.Cell.cycle = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Cell cycle" , T, F) ) %>% 
    dplyr::mutate(Patel.Complete.Immune.response = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Complete/Immune response" , T, F) ) %>% 
    dplyr::mutate(Patel.Hypoxia = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Hypoxia" , T, F) ) %>% 
    
    dplyr::rename(Neftel.AC = neftel.meta.module.AC) %>% 
    dplyr::rename(Neftel.NPC1 = neftel.meta.module.NPC1) %>% 
    dplyr::rename(Neftel.NPC2 = neftel.meta.module.NPC2) %>% 
    dplyr::rename(Neftel.OPC = neftel.meta.module.OPC) %>% 
    dplyr::rename(Neftel.MES1 = neftel.meta.module.MES1) %>% 
    dplyr::rename(Neftel.MES2 = neftel.meta.module.MES2) %>% 
    
    dplyr::rename(`Extracellular Matrix` = EM.struct.constituent) %>%
    
    dplyr::mutate(significant = 
                    padj.gsam.tpc.res < 0.01 &
                    abs(log2FoldChange.gsam.tpc.res) > 0.5 &
                    abs(log2FoldChange.glass.tpc.res) > 0.5 & 
                    direction.gsam.tpc.res == direction.glass.tpc.res) %>%
    
    dplyr::filter(significant == T) %>%
    dplyr::select(c('gid', 'hugo_symbol' , 'McKenzie_celltype_top_human_specificity', 'direction_up', 'direction_down',
                    TCGA.CL, TCGA.MES, TCGA.PN ,
                    Patel.Cell.cycle, Patel.Complete.Immune.response, Patel.Hypoxia, 
                    Neftel.AC, Neftel.OPC, Neftel.MES1, Neftel.MES2, Neftel.NPC1, Neftel.NPC2,
                    `Extracellular Matrix`
                    
                    ,significant,
                    
                    
                    padj.gsam.tpc.res,
                    direction.gsam.tpc.res ,log2FoldChange.gsam.tpc.res, 
                    
                    padj.glass.tpc.res,
                    direction.glass.tpc.res ,log2FoldChange.glass.tpc.res
    ))
  , by=c('hugo_symbol'='hugo_symbol')) %>% 
  dplyr::mutate(significant=NULL)
stopifnot(nrow(export) == 422)


xlsx::write.xlsx(export, "output/tables/corcorplot.xlsx", sheetName = "Sheet1")



# cor_cor_plot(plt, labels, method="ward")
# ggsave("output/figures/cor_ward.D1x.png", height=30, width=30) # then average, then ward.D1?







## 2.17a Corrected LFC + individual tophits ----



plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2, 2 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2, -2 , log2FoldChange.gsam.tpc.res)) %>%
  
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  
  #dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res <= 0.01 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  dplyr::mutate(show.label = 
                  !is.na(padj.gsam.res) &  padj.gsam.res < 0.01 &
                  !is.na(padj.glass.res) & padj.glass.res < 0.01)
                plt %>% filter(show.label) %>% pull(hugo_symbol)


p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)") +
  youri_gg_theme + xlim(-3, 3)

p2 <- ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc  ,
                      col=show.label,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt,  show.label == F),cex=0.35) +
  geom_point(data=subset(plt, show.label == T),col="red",cex=0.65) +
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,                   nudge_x = 3.1, direction = "y", hjust = "left" ) +
  #geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 ,                  nudge_x = -3.1, direction = "y", hjust = "right"  ) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 GLASS",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-3, 3)

p1 + p2


#plot(plt$padj.glass.res , plt$padj.gsam.res, pch=19,cex=0.1)
plot(plt$stat.glass.res, plt$stat.gsam.res, col=plt$show.label + 1, pch = 19, cex = 0.4)



## 2.17 Corrplot DE + marker genes (Figure 5?) ----



plt.genes <- data.frame(
  hugo_symbol = c(
    # "CREB5",	"TRIM24",	"ETV1", "COA1", # tumor [chr7]
    # "CACHD1","AHCYL1","GPR37L1","BMPR1B", # astroctyes
    # "RBFOX3", "GABRB2", "SLC17A7", "SST", # neuron; significant anyway
    # "PLP1", "OPALIN", "TMEM144","CLCA4", # oligodendrocyte [    #"PDGFA", "PDGFRA", "OLIG1", "OLIG2", "OLIG3",]
    # "TIE1","PEAR1","RGS5","NOSTRIN",  # endothelial
    # "CD163",  "CD14", "C1QA","THEMIS2", # TAM/MG
    # 
    # #"SSTR1","SSTR2","SSTR3","SSTR5", # Antibodies
    # #"GABRA1","GABRA2","GABRB1","GABRB2",
    # #"TNNT1", "TNNT2", "TNNT3",
    # 
    # ###, "MME", "ERG", "FCER2", "EPCAM", "EREG" << !!
    # #,"EGFR"
    # #,"EREG","AREG", "BTC","HBEGF","NGF","TGFA","EGF","EPGN"
    
    #"NDRG1","VEGFA",
    
    (results.out %>% dplyr::filter(!is.na(padj.gsam.res) &  padj.gsam.res < 0.01 & !is.na(padj.glass.res) & padj.glass.res < 0.01 ) %>% pull(hugo_symbol)))) %>%
  dplyr::filter(!duplicated(hugo_symbol)) %>%
  dplyr::left_join(results.out %>% 
                     dplyr::select(c('gid','hugo_symbol','McKenzie_celltype_top_human_specificity','show.marker.chr7','ensembl_id')) %>%
                     dplyr::rename(type=McKenzie_celltype_top_human_specificity) ,
                   by=c('hugo_symbol'='hugo_symbol'))  %>%
  dplyr::mutate(type = ifelse(hugo_symbol %in% c('SLC17A7','SST'), 'neuron', type) ) %>%
  dplyr::mutate(type = ifelse(is.na(type), 'NA', type) ) %>%
  dplyr::mutate(type = ifelse(show.marker.chr7 & type != 'astrocyte' , 'chr7/gain' , type)) %>%
  dplyr::mutate(show.marker.chr7 = NULL) %>%
  dplyr::filter(!is.na(gid))



### G-SAM ----


tmp <- gsam.gene.expression.all.vst %>% 
  as.data.frame() %>%
  dplyr::select(colnames(gsam.gene.expression.all.paired)) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% plt.genes$gid) %>%
  tibble::column_to_rownames('gid')
tmp.r1 <- tmp[,1:ncol(tmp) %>% purrr::keep(~ . %% 2 == 1)] # odd
tmp.r2 <- tmp[,1:ncol(tmp) %>% purrr::keep(~ . %% 2 == 0)] # even
stopifnot ( gsub("^(...).*$","\\1",colnames(tmp.r1)) == gsub("^(...).*$","\\1",colnames(tmp.r2)) )
rm(tmp)



tmp <- plt.genes %>% dplyr::left_join(
  log2(tmp.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
       tmp.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) ) %>%
    `rownames<-`(gsub("^ENSG.+\\|(.+)\\|chr.+$","\\1",rownames(.))) %>%
    tibble::rownames_to_column('hugo_symbol')
  , by = c('hugo_symbol'='hugo_symbol')) %>%
  dplyr::mutate(type = case_when(
    type == "chr7/gain" ~ "chr7", 
    type == "astrocyte" ~ "astr",
    type == "neuron" ~ "neur",
    type == "endothelial" ~ "endo",
    type == "microglia/TAM" ~ "TAM",
    type == "oligodendrocyte" ~ "olig",
    T ~ "?" )) %>%
  dplyr::mutate(hugo_symbol = paste0(hugo_symbol , " [" , type, "]") ) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL, type=NULL, ensembl_id=NULL)
tmp.2 <- dplyr::full_join(
  data.frame(sid = colnames(tmp.r1)) %>%
    dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
    dplyr::mutate(pid = gsub("^(...).*$","\\1",sid)) %>%
    dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna) %>%
    dplyr::mutate(sid = NULL) ,
  data.frame(sid = colnames(tmp.r2)) %>%
    dplyr::left_join(gsam.metadata.all.paired %>% dplyr::select(c('sid','tumour.percentage.dna')), by=c('sid'='sid')) %>%
    dplyr::mutate(pid = gsub("^(...).*$","\\1",sid))%>%
    dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) %>%
    dplyr::mutate(sid = NULL) , by = c('pid'='pid'))  %>%
  dplyr::mutate(`tumor-% DNA` = log2(tumour.percentage.dna.R1 / tumour.percentage.dna.R2)) %>%
  dplyr::mutate(tumour.percentage.dna.R1 = NULL, tumour.percentage.dna.R2 = NULL) %>%
  tibble::column_to_rownames('pid') %>%
  t()
stopifnot(colnames(tmp) == colnames(tmp.2))
lfc.gsam <- rbind(tmp.2, tmp) # genes plus tumor percentage?!
rm(tmp.2, tmp, tmp.r1, tmp.r2)



#"ward.D", "single", "complete", "average", "mcquitty",  "median", "centroid", "ward.D2"
# volgorde voor plotten:
#h = hclust(as.dist(1 - abs(cor(t(tmp)))),method = "ward.D")
#h = hclust(as.dist(1 - cor(t(tmp),method="pearson", use = "pairwise.complete.obs")))
#h = hclust(dist(t(tmp)))
#h = hclust(as.dist(1 - cor(t(tmp),method="pearson", use = "pairwise.complete.obs")),method="average")
#h = hclust(dist(t(scale(t(tmp)))), method="average")
#h = hclust(dist(t(scale(t(tmp)))) * 2, method="cen")
#h = hclust(as.dist(1 - cor(t(tmp),method="pearson", use = "pairwise.complete.obs")))# ,method="average"
#plot(h, hang = -1)
#h2 <- cutree(h, k = 4)
#colors = c("red", "blue", "green", "black")
#plot(h, hang = -1)
#d = as.dendrogram(h) %>%  set("branches_k_color", k=8) 
#plot(d)
#p1 <- as.ggdend(d)


h <- hclust(as.dist(1 - abs(cor(t(lfc.gsam),method="pearson", use = "pairwise.complete.obs")))) # nicest
order <- h$labels[h$order]
rm(h)
plt.gsam <- lfc.gsam %>% as.matrix %>% t() %>% as.data.frame %>% dplyr::select(order) %>%  cor()
stopifnot(colnames(plt.gsam) == rownames(plt.gsam))
plt <- plt.gsam
labels <- rownames(plt)
labels <- gsub("^[^ ]+ ","",labels)
labels <- gsub("[?]","",labels,fixed=T)
labels <- gsub("[","",labels,fixed=T)
labels <- gsub("]","",labels,fixed=T)
rownames(plt) <- labels
corrplot::corrplot(plt.gsam, method = "circle",tl.cex=0.75) # , order="hclust", addrect=4)


# pheatmap(tmp,scale="row",
#          clustering_distance_rows = "correlation",
#          clustering_distance_cols = "euclidean",
#          cluster_cols = T,
#         cutree_cols = 3,
#   cutree_rows = 4,
# fontsize_row = 2.1 )
# 



png(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient_and_DE_genes.png", width = 1200 * 2.8, height = 900 * 2.8 )
corrplot::corrplot(cor(t(tmp), method="pearson"), method = "circle",tl.cex=2.25)
dev.off()


#pheatmap::pheatmap(t(tmp), scale="column",clustering_distance_rows="correlation")

# a <- data.frame(a= c(1,2,5,6,7,8,1,2,5,6,7,8,1,2,5,6,7,8), b= c(4,4,6,8,8,9,4,4,6,8,8,9,4,4,6,8,8,9))
# a <- cbind(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)
# colnames(a) <- make.names(colnames(a),unique = T)
# ggcorrplot(cor(a)) # incredibly slow


### GLASS ----


tmp.ids <- data.frame(R1 = colnames(glass.gene.expression.all.vst)) %>%
  dplyr::filter(grepl("-TP-",R1)) %>%
  dplyr::mutate(unified = gsub("-TP-.+$","",R1,fixed=F)) %>%
  dplyr::full_join(
    data.frame(R2 = colnames(glass.gene.expression.all.vst)) %>%
      dplyr::filter(grepl("-TP-",R2) == F) %>%
      dplyr::mutate(unified = gsub("-R[1-4]-.+$","",R2,fixed=F))
    , by=c('unified' = 'unified')) %>%
  dplyr::filter(! is.na(R1) & !is.na(R2))

tmp <- glass.gene.expression.all.vst %>%
  as.data.frame() %>%
  dplyr::select(sort(colnames(.))) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% plt.genes$ensembl_id) %>%
  tibble::column_to_rownames('gid')
tmp.r1 <- tmp %>% dplyr::select(tmp.ids$R1) # R1
tmp.r2 <- tmp %>% dplyr::select(tmp.ids$R2) # R2,R3,R4
stopifnot (gsub("^(.{12}).+$","\\1",colnames(tmp.r1)) == gsub("^(.{12}).+$","\\1",colnames(tmp.r2)) )
rm(tmp.ids, tmp)


# data, pre correlation
lfc.glass <- plt.genes %>% dplyr::left_join(
  log2(tmp.r1 %>% `colnames<-`(gsub("^(.{12}).*$","\\1",colnames(.))) /
         tmp.r2 %>% `colnames<-`(gsub("^(.{12}).*$","\\1",colnames(.))) ) %>%
    tibble::rownames_to_column('ensembl_id')
  , by = c('ensembl_id'='ensembl_id')) %>%
  dplyr::mutate(type = case_when(
    type == "chr7/gain" ~ "chr7", 
    type == "astrocyte" ~ "astr",
    type == "neuron" ~ "neur",
    type == "endothelial" ~ "endo",
    type == "microglia/TAM" ~ "TAM",
    type == "oligodendrocyte" ~ "olig",
    T ~ "?" )) %>%
  dplyr::mutate(hugo_symbol = paste0(hugo_symbol , " [" , type, "]") ) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(   gid    = NULL,  type=NULL, ensembl_id = NULL)
rm(tmp.r1, tmp.r2)


# h = hclust(as.dist(1 - abs(cor(t(tmp),method="pearson", use = "pairwise.complete.obs")))) # nicest
# #h = hclust(as.dist(1 - abs(cor(t(tmp))))) # nicest
# h = hclust(as.dist(1 - cor(t(tmp)))) # nicest
# h = hclust(dist(t(scale(t(tmp),center=F)))) # nicest
# d = t(scale(t(tmp), center = F))
# pc <- prcomp(d)
# plot(pc$x[,1],pc$x[,2],)
# text(pc$x[,1],pc$x[,2],rownames(pc$x),cex=0.6,pos=3)
# plot(pc$x[,3],pc$x[,2],)
# plot(pc$x[,4],pc$x[,2],)
# plot(pc$x[,5],pc$x[,2],)
# h = hclust(Dist(pc$x[,c(1,3,4,5)],method="euclidean"),method=m[6]) #m6 was nice with eucl/corr on full data
# a = pc$x %>% as.data.frame %>% dplyr::mutate(PC1.abs = abs(PC1)) %>% dplyr::arrange('PC1.abs') %>% rownames_to_column('gid') %>% pull('gid')
#h = hclust(Dist(t(scale(t(tmp), center = F)), method="abscorrelation"),method=m[ 1]) #m6 was nice
#plot(h, hang = -1)
# h2 <- cutree(h, k = 4)
# colors = c("red", "blue", "green", "black")
# plot(h, hang = -1)
# 
# d = as.dendrogram(h) %>%  set("branches_k_color", k=8) 
# #plot(d)
# #p1 <- as.ggdend(d)
#plt <- cor(t(tmp %>% as.matrix %>% t() %>% as.data.frame %>% dplyr::select(h$labels[h$order]) %>% t()))
#plt <- cor(t(t(scale(t(tmp),center=F)) %>% as.matrix %>% t() %>% as.data.frame %>% dplyr::select(a) %>% t()))


h = hclust(Dist(t(scale(t(lfc.glass), center = F)), method="euclidean"),method=m[6]) #m6 was nice

plt.order <- data.frame(id = h$labels[h$order]) %>%
  dplyr::mutate(i = nrow(.):1) %>%
  dplyr::arrange(i) %>%
  dplyr::mutate(cluster = as.factor(case_when(
               i <= 4 ~  'glcm1', # glioblastoma longtitudinal correlation metafeature
    i >= 4   & i <= 30 ~ 'glcm2',
    i >= 30  & i <= 53 ~ 'glcm3',
    i >= 53  & i <= 89 ~ 'glcm4',
    i >= 89  & i <= 223 ~ 'glcm5',
    i >= 223 & i <= 238 ~ 'glcm6',
    i >= 238 & i <= 257 ~ 'glcm7',
    i >= 257 & i <= 263 ~ 'glcm8',
    i >= 263            ~ 'glcm9'
    ))) %>%
  dplyr::mutate(cluster = factor(cluster, levels = 
       c("glcm8","glcm5","glcm4", "glcm6","glcm2","glcm7","glcm3","glcm1","glcm9") )) %>%
  dplyr::arrange(cluster, i) %>%
  dplyr::mutate(k = 1:nrow(.))

#plt.order %>% pull(id)


plt.glass <- cor(t(t(scale(t(lfc.glass),center=F)) %>% as.matrix %>% t() %>% as.data.frame %>% dplyr::select(h$labels[h$order]) %>% t()))
plt <- plt.glass
head(plt[1:6,1:6])
#plt <- cor(t(t(scale(t(tmp),center=F)) %>% as.matrix %>% t() %>% as.data.frame %>% dplyr::select(plt.order %>% dplyr::arrange(cluster, i) %>% pull(id)) %>% t()))


labels <- colnames(plt)
labels <- gsub("^[^ ]+ ","",labels)
labels <- gsub("[?]","",labels,fixed=T)
labels <- gsub("[","",labels,fixed=T)
labels <- gsub("]","",labels,fixed=T)
colnames(plt) <- labels

labels <- rownames(plt)
labels <- gsub(" .+$","",labels)
labels <- gsub("TNNT1","***        TNNT1",labels)
labels <- gsub("TNNT2","***        TNNT2",labels)
labels <- gsub("^GABR","***        GABR",labels)
rownames(plt) <- labels


lines <- plt.order %>% dplyr::group_by(cluster) %>% summarise(max = max(k)) %>% dplyr::mutate(max = max) %>% pull(max) %>% c(0)



png(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient_and_DE_genes_GLASS.png", width = 1200 * 2.8, height = 900 * 2.8 )
corrplot::corrplot(plt, method = "circle",tl.cex=0.75) # , order="hclust") #, addrect=4
for(line in lines) {
  lines(c(nrow(plt) - line,nrow(plt) - line) + 0.5 ,c(0,nrow(plt)) + 0.5, lwd=3.5, col="black" )
  lines(c(0,nrow(plt)) + 0.5 , c( line ,line) + 0.5, lwd=3.5, col="black" )
}
dev.off()




svg(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient_and_DE_genes_GLASS.svg", width = 120, height = 90 )
corrplot::corrplot(plt, method = "circle",tl.cex=2.25)

for(line in lines) {
  lines(c(nrow(plt) - line,nrow(plt) - line),c(0,nrow(plt) + 0.5), lwd=1, col="black" )
  #lines(c(0.5,nrow(plt) + 0.5) , c( line + 0.5,line + 0.5), lwd=1, col="black" )
}
dev.off()


# nos2: https://pure.uva.nl/ws/files/2782854/178651_04.pdf # page25



plt <- plt.gsam %>%
  as.data.frame() %>%
  dplyr::select(h$labels[h$order]) %>%
  t() %>% as.data.frame() %>% dplyr::select(h$labels[h$order]) %>% as.matrix
plt[1:6,1:6]
stopifnot(rownames(plt) == colnames(plt))
head(plt[1:6,1:6])



### combined ----

lfc.combi <- rbind(
    lfc.glass %>% t() , 
    lfc.gsam %>% t() %>% as.data.frame() %>% dplyr::mutate(`tumor-% DNA` = NULL)
  )

# lines = c(0, 30, 49 , 66, 68, 86, 113, 172, 233, nrow(plt.combi))
# h <- lfc.combi %>%
#   t() %>%
#   #scale() %>%
#   scale(center = F) %>%
#   Dist(method = "correlation") %>%
#   hclust(method=m[1])
# order <- h$labels[h$order] %>% purrr::discard(. %in% c("tumor-% DNA"))
# rm(h)


h <- lfc.combi %>%
  t() %>%
  #scale() %>%
  #scale(center = F) %>%
  amap::Dist(method = "pearson") %>%
  hclust(method=m[1])
order <- h$labels[h$order] %>% purrr::discard(. %in% c("tumor-% DNA"))
rm(h)

plt.combi <- lfc.combi %>% t() %>%
  scale(center=F) %>%
  as.matrix %>% 
  t() %>% 
  as.data.frame %>%
  dplyr::select(order) %>%
  cor()


labels <- colnames(plt.combi)
labels <- gsub("^[^ ]+ ","",labels)
labels <- gsub("[?]","",labels,fixed=T)
labels <- gsub("[","",labels,fixed=T)
labels <- gsub("]","",labels,fixed=T)
colnames(plt.combi) <- labels

labels <- rownames(plt.combi)
labels <- gsub(" .+$","",labels)
labels <- gsub("TNNT1","***        TNNT1",labels)
labels <- gsub("TNNT2","***        TNNT2",labels)
labels <- gsub("NOS2","***        NOS2",labels)
labels <- gsub("VAV3","***        VAV3",labels)
labels <- gsub("^GABR","***        GABR",labels)
rownames(plt.combi) <- labels


lines = c(0, 16, 45, 72, 88, 108, 126, 128,  161, 173, 239, nrow(plt.combi))
png(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient_and_DE_genes_combi.png", width = 1200 * 2.8, height = 900 * 2.8 )
corrplot::corrplot(plt.combi, method = "circle",tl.cex=0.75) # , order="hclust") #, addrect=4
for(line in lines) {
  lines(c(nrow(plt.combi) - line,nrow(plt.combi) - line) + 0.5 ,c(0,nrow(plt.combi)) + 0.5, lwd=3.5, col="black" )
  lines(c(0,nrow(plt.combi)) + 0.5 , c( line ,line) + 0.5, lwd=3.5, col="black" )
}
dev.off()




## 2.18 per patient clustering ----

#mat_breaks <- seq(min(mat), max(mat), length.out = 20)

# quantile_breaks <- function(xs, n = 10) {
#   breaks <- quantile(xs$values, probs = seq(0, 1, length.out = n))
#   breaks[!duplicated(breaks)]
# }
# mat_breaks <- quantile_breaks(mat, n = 11)


# the expected value of a logFc is 0, and thus scaling needs to be adapted to that
scale0 <- function(m) {
    return ( apply(m, 2,  function(v) {
      svar <- sum(v^2) / (length(v) - 1)
      return(v / sqrt(svar))
    }     ) )
}

# m <- matrix( c(-4,-4,-3,-3,-2,-2,2,2,3,3,4,4) ,nrow=2)  %>%
#   t()
# m %>% scale() %>% t() 
# m %>% scale0() %>% t() 


# calc per gene var
h <- lfc.combi %>%
  as.matrix %>%
  t() %>%
  scale0 %>%
  t() %>%
  amap::Dist(method='correlation') %>%
  hclust(method = m[2])
order <- h$labels[h$order]

gltc <- data.frame(pid = h$labels[h$order],
                   gltc = paste0('GLTC',cutree(h,3)[h$order])) %>%
          dplyr::mutate(dataset = ifelse(grepl("GLSS|TCGA",pid),"GLASS","GSAM"))

plt <-  lfc.combi %>% t() %>% as.data.frame %>% dplyr::select(gltc$pid) %>% as.matrix() %>% cor()


png(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient_and_DE_genes_pat.png", width = 1200 * 2.8, height = 900 * 2.8 )
corrplot::corrplot(plt)

lines  <- gltc %>%
  dplyr::mutate(i = 1:nrow(.)) %>%
  group_by(gltc) %>%
  summarise(max = max(i)) %>% dplyr::mutate(max = max) %>% pull(max) %>% c(0) %>% sort()
for(line in lines) {
  lines(c(nrow(lfc.combi) - line,nrow(lfc.combi) - line) + 0.5 ,c(0, nrow(lfc.combi)) + 0.5, lwd=2.5, col="black" )
  lines(c(0, nrow(lfc.combi)) + 0.5 , c(line ,line) + 0.5, lwd=2.5, col="black" )
}

dev.off()




# 
# 
# plt <- cor(t(lfc.combi))
# colnames(plt) <- data.frame(cur.pid = colnames(plt)) %>%
#   dplyr::left_join(gltc, by=c('cur.pid'='pid')) %>%
#   dplyr::pull(gltc)
# 
# png(file = "output/figures/paper_dge_corrplot_logFc_gene_per_patient_and_DE_genes_pat.png", width = 1200 * 2.8, height = 900 * 2.8 )
# corrplot(plt, order='hclust')
# for(line in lines) {
#   lines(c(nrow(lfc.combi) - line,nrow(lfc.combi) - line) + 0.5 ,c(0, nrow(lfc.combi)) + 0.5, lwd=2.5, col="black" )
#   lines(c(0, nrow(lfc.combi)) + 0.5 , c(line ,line) + 0.5, lwd=2.5, col="black" )
# }
# dev.off()
# 


### SVVL!? ----

library(survival)
library(survminer)


plt <- gltc %>%
  dplyr::filter(dataset == "GSAM") %>%
  dplyr::left_join(gsam.patient.metadata, by=c('pid'='studyID')) %>%
  dplyr::left_join(gsam.rna.metadata %>%
                     dplyr::filter(resection == 'r1' & !is.na(resection_pair)) %>%
                     dplyr::select(c('pid','NMF.123456.PCA.SVM.class')) %>%
                     dplyr::rename('subtype.R1'='NMF.123456.PCA.SVM.class')
                   , by=c('pid'='pid')) %>%
  dplyr::left_join(gsam.rna.metadata %>%
                     dplyr::filter(resection == 'r2' & !is.na(resection_pair)) %>%
                     dplyr::select(c('pid','NMF.123456.PCA.SVM.class')) %>%
                     dplyr::rename('subtype.R2'='NMF.123456.PCA.SVM.class')
                   , by=c('pid'='pid')) %>%
  dplyr::mutate(subtype.stable = subtype.R1 == subtype.R2)


stopifnot(duplicated(plt$pid) == F)


ggplot(plt, aes(x = gltc, fill = tumorLocation)) +
  geom_bar(position = "stack")
  

# dit is wel next level
ggplot(plt, aes(x = gltc, fill = treatedWithTMZ)) +
  geom_bar(position = "stack")


ggplot(plt, aes(x = gltc, y = sizeLargestLesionSecondSurgery)) +
  geom_point()


ggplot(plt, aes(x = gltc, y = daysToLastCycleTMZ)) +
  geom_point()


ggplot(plt, aes(x = gltc, y = nMutationsRecurrent)) +
  geom_point()


ggplot(plt, aes(x = gltc, fill = subtype.R1)) +
  geom_bar(position = "stack")
ggplot(plt, aes(x = gltc, fill = subtype.R2)) + # :)
  geom_bar(position = "stack")
ggplot(plt, aes(x = gltc, fill = subtype.stable)) +
  geom_bar(position = "stack")






surv.model <- Surv(time = plt$progressionFreeDays)
surv.model <- Surv(time = plt$daysToSecondSurgery)
surv.model <- Surv(time = plt$survivalFromSecondSurgeryDays , event = plt$survival.event)
surv.model <- Surv(time = plt$survivalDays , event = plt$survival.event)

surv.fit <- survfit(surv.model ~ gltc, data = plt)
ggsurvplot(surv.fit , data = plt, pval = TRUE)

surv.fit <- survfit(surv.model ~ subtype.R2, data = plt)
ggsurvplot(surv.fit , data = plt, pval = TRUE)




## 2.18b clustering on VST values? ----

plt.genes <- data.frame(
  hugo_symbol = c((results.out %>% dplyr::filter(!is.na(padj.gsam.res) &  padj.gsam.res < 0.01 & !is.na(padj.glass.res) & padj.glass.res < 0.01 ) %>% pull(hugo_symbol)))) %>%
  dplyr::filter(!duplicated(hugo_symbol)) %>%
  dplyr::left_join(results.out %>% 
                     dplyr::select(c('gid','hugo_symbol','McKenzie_celltype_top_human_specificity','show.marker.chr7','ensembl_id')) %>%
                     dplyr::rename(type=McKenzie_celltype_top_human_specificity) ,
                   by=c('hugo_symbol'='hugo_symbol'))  %>%
  dplyr::mutate(type = ifelse(hugo_symbol %in% c('SLC17A7','SST'), 'neuron', type) ) %>%
  dplyr::mutate(type = ifelse(is.na(type), 'NA', type) ) %>%
  dplyr::mutate(type = ifelse(show.marker.chr7 & type != 'astrocyte' , 'chr7/gain' , type)) %>%
  dplyr::mutate(show.marker.chr7 = NULL) %>%
  dplyr::filter(!is.na(gid))

tmp <- gsam.gene.expression.all.vst %>%
  as.data.frame %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% plt.genes$gid) %>%
  dplyr::left_join(plt.genes, by=c('gid'='gid')) %>%
  dplyr::mutate(type = case_when(
    type == "chr7/gain" ~ "chr7", 
    type == "astrocyte" ~ "astr",
    type == "neuron" ~ "neur",
    type == "endothelial" ~ "endo",
    type == "microglia/TAM" ~ "TAM",
    type == "oligodendrocyte" ~ "olig",
    T ~ "?" )) %>%
  dplyr::mutate(hugo_symbol = paste0(hugo_symbol , " [" , type, "]") ) %>%
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL, type = NULL, ensembl_id = NULL)

# why are these different?!
#h <- hclust(as.dist(1 - abs(cor(t(tmp),method="pearson", use = "pairwise.complete.obs")))) # nicest

h <- tmp %>%
  t() %>%
  scale(center = T) %>%
  t() %>%
  Dist(method = "pearson") %>%
  hclust(method=m[1])

order <- h$labels[h$order]
rm(h)
plt <- tmp %>% as.matrix %>% t() %>% as.data.frame %>% dplyr::select(order) %>%  cor()

labels <- colnames(plt)
labels <- gsub("^[^ ]+ ","",labels)
labels <- gsub("[?]","",labels,fixed=T)
labels <- gsub("[","",labels,fixed=T)
labels <- gsub("]","",labels,fixed=T)
colnames(plt) <- labels

labels <- rownames(plt)
labels <- gsub(" .+$","",labels)
labels <- gsub("TNNT1","***        TNNT1",labels)
labels <- gsub("TNNT2","***        TNNT2",labels)
labels <- gsub("NOS2","***        NOS2",labels)
labels <- gsub("VAV3","***        VAV3",labels)
labels <- gsub("^GABR","***        GABR",labels)
rownames(plt) <- labels

png(file = "output/figures/paper_dge_corrplot_vst_and_DE_genes_combi.png", width = 1200 * 2.8, height = 900 * 2.8 )
lines = c(0,3, 22, 30, 51, 75, 96,98, 139, 155, 222,223, nrow(plt.combi))
corrplot::corrplot(plt.combi, method = "circle",tl.cex=0.75) # , order="hclust") #, addrect=4
for(line in lines) {
  lines(c(nrow(plt.combi) - line,nrow(plt.combi) - line) + 0.5 ,c(0,nrow(plt.combi)) + 0.5, lwd=2.5, col="black" )
  lines(c(0,nrow(plt.combi)) + 0.5 , c( line ,line) + 0.5, lwd=2.5, col="black" )
}
dev.off()


## 2.19 pathology antibody genes ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.glass.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2, 2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2, -2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 3, 3 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -3, -3 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(show.label = hugo_symbol %in% IHC_diagnostiek_IVD_antibodies$hugo_symbol & !is.na(McKenzie_celltype_top_human_specificity))


plt %>% filter(!is.na(McKenzie_celltype_top_human_specificity) & show.label)



ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_label_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, fill="white",label.size = 0) + 
  geom_label_repel(data=subset(plt, show.label == T  & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, fill="white",label.size = 0) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Neuron marker genes") +
  youri_gg_theme + xlim(-2.5, 2.5)



## g:profiler upregulated ----

# significantly UO regulated genes are related to synaptic signalling and neurons / neuronal development / cell junctions
results.out %>%
  dplyr::filter(log2FoldChange.gsam.tpc.res  >= 0.5 ) %>%
  dplyr::filter(padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3) %>%
  #dplyr::top_n(500, -padj.tpc.res) %>%
  dplyr::pull(hugo_symbol) 

#%>%   length()


## g:profiler downregulated top1000 ----

# removing LFC cut-off makes signal stronger

# significantly DOWN regulated genes are related to blood vessels / angiogenesis
results.out %>%
  dplyr::filter(log2FoldChange.gsam.tpc.res <= -0.5) %>%
  dplyr::filter(padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3) %>%
  #dplyr::top_n(1000, -padj.gsam.tpc.res) %>%
  dplyr::arrange(padj.gsam.tpc.res) %>%
  dplyr::pull(hugo_symbol)



## g:profiler all top1000 ----

# DE regulated genes are related to
gsam.gene.res.combined %>%
  dplyr::filter(padj.tpc.res < 0.01 & p.value >= 0.05) %>%
  dplyr::top_n(1000, -padj.tpc.res) %>%
  #dplyr::arrange(padj.tpc.res) %>%
  dplyr::pull(hugo_symbol.tpc.res)



## plot pharmaco genes ----


p1 <- ggplot(gsam.gene.res.combined, aes(x=log2FoldChange.res ,
                                         y=statistic,
                                         #shape=significant.tpc.res,
                                         col = pharma.relation,
                                         label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == F ), cex=0.05, pch=19) +
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #scale_shape_manual(values = c('TRUE'=1, 'FALSE' = 19) ) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.res == T & abs(log2FoldChange.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme

p2 <- ggplot(gsam.gene.res.combined, aes(x=log2FoldChange.tpc.res ,
                                         y= statistic,
                                         #shape=significant.tpc.res,
                                         col = pharma.relation,
                                         label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == F ), cex=0.05, pch=19) +
  geom_point(data=subset(gsam.gene.res.combined, pharma.relation == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  #scale_shape_manual(values = c('TRUE'=1, 'FALSE' = 19) ) +
  #geom_text_repel(data = subset(gsam.gene.res.combined, significant.tpc.res == T & abs(log2FoldChange.tpc.res) > 1.5), size = 2.4 )  +
  geom_text_repel(data = subset(gsam.gene.res.combined, significant.tpc.res == T & pharma.relation), size = 5,
                  #box.padding = 1.5,
                  nudge_x = 3.1, direction = "y", hjust = "left" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


## plot synaptic signalling GO:0099537 ----

tss <- read.delim('data/gProfiler_trans-synaptic_signaling_GO.0099537.tsv',header=T, sep=",")
tss <- read.delim('data/gProfiler_Morphine_addiction_KEGG.05032.tsv',header=T, sep=",")


#summary(as.factor(gsub("^.+(chr[^:]+):.+$","\\1",plt %>% dplyr::filter(significant.tpc.res) %>% pull(gid.res) )))



plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = ensembl_id.res %in% tss$converted_alias ) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label == F & significant.res == T ), cex=0.15, pch=19, alpha=1) +
  geom_point(data=subset(plt, show.label == F & significant.res == F ), cex=0.15, pch=19, alpha=0.4) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label == F & significant.tpc.res == T ), cex=0.15, pch=19, alpha=1) +
  geom_point(data=subset(plt, show.label == F & significant.tpc.res == F ), cex=0.15, pch=19, alpha=0.4) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot glioma stem cells (GSC) markers ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("PROM1|CD44|FUT4|CD70|S100A4|ALDH1A3|NANOG$|^SOX2$|NES", hugo_symbol.tpc.res)) %>%
  #grepl("MGMT|(TP53)$|(EGFR)$|ABCB(1|11|4|5)$|CCN2|(HTRA3)$|GLI1|MTOR$", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot PI3K/Akt pathway  ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("AKT[0-9]$|TSC[1-9]$|^PTEN$|CDKN2[AB]$|^PIP[0-9]$|OAF|^MTOR|IRS1$|NFKB1", hugo_symbol.tpc.res)) %>%
  # grepl("MPG$|HMGA2|^POLB|PARP1$", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


#
gid <- 'FLT1|KDR'
rownames(gsam.gene.expression.all)[grepl(gid,rownames(gsam.gene.expression.all))]
plt$gid.res[grepl(gid,plt$gid.res )]



## plot BER mechanism  ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("MPG$|HMGA2|^POLB|PARP1$|HDAC", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


## plot MMR mechanism  ----
# https://cdrjournal.com/article/view/3779

plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  grepl("MLH1|MLH3|MSH2|MSH3|MSH6|PMS1|^PMS2$|^POLD[1-4]$|PCNA|^RPA$|HMGB1$|RFC|LIG1", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


## plot MGMT ----

# https://www.sciencedirect.com/science/article/pii/S2352304216300162?via%3Dihub


plt <-gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
  #                grepl("PROM1|CD44|FUT4|CD70|S100A4|ALDH1A3|NANOG$|^SOX2$|NES", hugo_symbol.tpc.res)) %>%
  grepl("MGMT|(TP53)$|(EGFR)$|ABCB(1|11|4|5)$|CCN2|(HTRA3)$|GLI1|MTOR$", hugo_symbol.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  

p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2


gid <- 'SOX2'
rownames(gsam.gene.expression.all)[grepl(gid,rownames(gsam.gene.expression.all))]
plt$gid.res[grepl(gid,plt$gid.res )]



## plot transporter genes


## hh signalling ? -----

## subtype genes ----
### subtype genes [CL] ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Classical") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Classical") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2

### subtype genes [PN] ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Proneural") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Proneural") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2




### subtype genes [MES] ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot ABCB1 ----


plt <- data.frame(ABCB1 = as.numeric(gsam.gene.expression.all.vst %>%
                                       as.data.frame() %>%
                                       dplyr::filter(grepl("ABCB1\\|", rownames(.)))) ,
                  res = gsub("^...(.).*$","R\\1",colnames( gsam.gene.expression.all.vst)) )

ggplot(plt, aes(x = res, y = ABCB1)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y = "ABCB1 relative expression")

rm(plt)

## plot MGMT ----


plt <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  dplyr::filter(grepl("MGMT\\|", rownames(.))) %>% 
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pid') %>%
  dplyr::rename(MGMT.vst = 'ENSG00000170430.10_4|MGMT|chr10:131265454-131569247(+)') %>%
  dplyr::mutate(res = gsub("^...(.).*$","R\\1", pid)) %>%
  dplyr::filter(pid %in% colnames(gsam.gene.expression.all))


#                  res = gsub("^...(.).*$","R\\1",colnames( gsam.rnaseq.expression.vst))) %>%
#dplyr::select(colnames(gsam.gene.expression.all))

ggplot(plt, aes(x = res, y = MGMT.vst)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y = "MGMT relative expression")


wilcox.test(
  plt %>% dplyr::filter(res == "R1") %>% dplyr::pull(MGMT.vst),
  plt %>% dplyr::filter(res == "R2") %>% dplyr::pull(MGMT.vst)
)




## plot MGMT + TP53 + p21 + Chk1 + Chk2 + ATM ----
## "Increased levels of MGMT or loss of the mismatch repair capacity confer resistance to temozolomide (3)" - doi:10.1158/1535-7163.MCT-05-0428


## plot GLI1 ----

plt <- gsam.gene.res.combined %>%
  dplyr::mutate(show.label = 
                  ensembl_id.tpc.res  %in%  (wang.glioma.intrinsic.genes %>%
                                               dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                               dplyr::pull(ENSG.short)) | 
                  hugo_symbol.tpc.res %in% (wang.glioma.intrinsic.genes %>%
                                              dplyr::filter(Subtyping_Signature_Gene. == "Mesenchymal") %>%
                                              dplyr::pull(Gene_Symbol))) %>%
  dplyr::mutate(log2FoldChange.res = ifelse(log2FoldChange.res > 3, 3, log2FoldChange.res) ) %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res))  


p1 <- ggplot(plt, aes(x=log2FoldChange.res , y=statistic, col = show.label , label=hugo_symbol.tpc.res)) + 
  geom_point(data=subset(plt, show.label == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res > 0 & significant.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  #geom_text_repel(data = subset(plt, show.label & log2FoldChange.res < 0 & significant.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2)") +
  youri_gg_theme

p2 <- ggplot(plt, aes(x=log2FoldChange.tpc.res , y= statistic, col = show.label  , label=hugo_symbol.tpc.res )) + 
  geom_point(data=subset(plt, show.label  == F ), cex=0.05, pch=19) +
  geom_point(data=subset(plt, show.label  == T ), cex=0.5, pch=19) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res > 0 & significant.tpc.res ), size = 5, nudge_x = 3.1, direction = "y", hjust = "left" )  +
  geom_text_repel(data = subset(plt, show.label & log2FoldChange.tpc.res < 0 & significant.tpc.res ), size = 5, nudge_x = -3.1, direction = "y", hjust = "right" )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)", y="Correlation t-statistic with tumour percentage", col="Difference significant (R1 ~ R2; tumour percentage corrected)") +
  youri_gg_theme

p1 + p2



## plot Akt pathway [TMZ resist?] ----


## plot [BMPR1B, CTGF, CYP4F2, EDNRB, ELL, EPHA3, FOS, GJA1, GPM6A, HTRA3, IGFBP2, IGFBP7, LEF1, RAI, RGL, SEC3L1, SRP72, SSB1, ZNF436] doi:10.1158/1535-7163.MCT-05-0428 ----






# - - - - -

# 
# EnhancedVolcano(res,
#                 lab = res$hugo_symbol,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = 0.01)





#save.image(file = "dge_supervised_expression_analysis.RData")
#load(file = "dge_supervised_expression_analysis.RData")

# TP53 is interesting

ggid <- 'SAA3'
gid <- res %>% dplyr::filter(hugo_symbol == ggid) %>% dplyr::pull('gid')
plt <- data.frame(sid = colnames(vst),  expr = as.numeric(as.data.frame(vst) %>% dplyr::filter(rownames(.) == gid))) %>%
  dplyr::left_join(metadata, by = c ('sid' = 'sid') )
ggplot(plt, aes(x = tumour.percentage.dna, y= expr,col=resection)) + #)) + # , 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  labs(y = paste0("Expression: ", ggid))




EnhancedVolcano(glass.gene.res.res,
                lab = glass.gene.res.res$gene_symbol ,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01)


## dplyr inner join (!)

combined.gene.res.res <- dplyr::inner_join(
  gsam.gene.res.res , glass.gene.res.res , by=c('ensembl_id.res'='ensembl_id.glass.res') ) %>%
  dplyr::inner_join( gsam.gene.res.tpc.res, by=c('ensembl_id.res'='ensembl_id.tpc.res') ) %>%
  dplyr::left_join(gsam.gene.expression.all.cor, by=c('ensembl_id.res' = 'ensembl_id')) %>%
  dplyr::mutate(chr = gsub("^.+(chr[^:]+):.+$","\\1",gid))



head(gsam.gene.res.res)
head(glass.gene.res.res)


dim(gsam.gene.res.res)

dim(combined.gene.res.res)


## plot GSAM x GLASS uncorrected ----

plt <- combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) # plotting limit


p1 <- ggplot(plt, aes(x=log2FoldChange.res,
                                        y= log2FoldChange.glass.res
                                        #,col=significant.tpc.res
                      )) +
  #geom_smooth (method = "lm",lty=2,lwd=0.7,se = F) +
  geom_vline(xintercept = 0, col="red", lwd=0.2) +
  geom_hline(yintercept = 0, col="red", lwd=0.2) +
  geom_point(pch=19,cex=0.3) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "log2FC R1 vs. R2 G-SAM (unpaired)",
       y= "log2FC R1 vs. R2 GLSS (unpaired)",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-2,3) +
  ylim(-2,3)

p2 <- ggplot(plt,  aes(x=log2FoldChange.tpc.res,
                       y= log2FoldChange.glass.res,
                       col = significant.tpc.res,
                       label=hugo_symbol.tpc.res)) +
  geom_point(pch=19,cex=0.3) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(combined.gene.res.res, significant.tpc.res == T & (log2FoldChange.tpc.res > 1.5 | log2FoldChange.tpc.res < 0.85)), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme


mx <- plt %>% dplyr::filter(chr == "chr19") %>% dplyr::pull(log2FoldChange.tpc.res) %>% median()
my <- plt %>% dplyr::filter(chr == "chr19") %>% dplyr::pull(log2FoldChange.glass.res) %>% median()

p3 <- ggplot(plt,  aes(x= log2FoldChange.tpc.res,
                       y= log2FoldChange.glass.res,
                       col = chr == "chr19",
                       label=hugo_symbol.tpc.res)) +
  geom_point(data = dplyr::filter(plt , chr != "chr19") , pch=19,cex=0.3) +
  geom_point(data = dplyr::filter(plt , chr == "chr19") , pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "FC GLASS",
       y="FC G-SAM",
       col="Difference significant (R1 ~ R2)"
  ) +
  geom_hline(yintercept=my) +
  geom_vline(xintercept=mx) +
  youri_gg_theme


p1 + p2 + p3





## de COL genen van het MES type gaan aanzienlijk omlaag


## plot GLASS x GSAM signi ----

plt <-  a = combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) %>% # plotting limit
  dplyr::mutate(significant = significant.tpc.res & padj.glass.res < 0.05 ) %>%
  dplyr::mutate(show.label = significant ) %>%
  dplyr::filter(significant == T) %>%
  dplyr::pull(hugo_symbol.tpc.res) 
  
wang.glioma.intrinsic.genes%>% dplyr::filter (Gene_Symbol %in% a )


ggplot(plt,  aes(x= log2FoldChange.tpc.res,
                 y= log2FoldChange.glass.res,
                 col = show.label,
                 label=hugo_symbol.tpc.res)) +
  geom_point(data = dplyr::filter(plt , show.label == F) , pch=19,cex=0.3) +
  geom_point(data = dplyr::filter(plt , show.label == T) , pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  labs(x = "FC GLASS",
       y="FC G-SAM",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme



## plot chr [chr19]


combined.gene.res.res %>%
  dplyr::filter(significant.glass.res) %>%
  dplyr::pull(chr) %>%
  as.factor %>%
  summary()

plt <- combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) # plotting limit


# clearly chr19 is more to the right
chr.candidate = "chr19"
ggplot(combined.gene.res.res, aes(x= log2FoldChange.glass.res,
                                  y= statistic,
                                  label=hugo_symbol.tpc.res,
                                  #col=significant.glass.res,
                                  #col = ensembl_id.res %in% (wang.glioma.intrinsic.genes %>% dplyr::filter( Subtyping_Signature_Gene. == "Mesenchymal") %>% dplyr::pull(ENSG.short))
                                  col = chr == chr.candidate
)) +
  geom_point(data = subset(combined.gene.res.res, chr != chr.candidate), pch=19,cex=0.3) +
  geom_point(data = subset(combined.gene.res.res, chr == chr.candidate), pch=19,cex=0.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(combined.gene.res.res, significant.glass.res == T & abs(log2FoldChange.glass.res) > 1.5), size = 2.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme



## plot GLI1 ----



plt <- combined.gene.res.res %>%
  dplyr::mutate(log2FoldChange.tpc.res = ifelse(log2FoldChange.tpc.res > 3, 3, log2FoldChange.tpc.res)) %>% # plotting limit
  dplyr::mutate(show.label = grepl("GLI1", hugo_symbol.tpc.res))


ggplot(plt, aes(x= log2FoldChange.glass.res,
                                  y= statistic,
                                  label=hugo_symbol.tpc.res,
                                  col = show.label )) +
  geom_point(data = subset(plt, show.label == F), pch=19,cex=0.3) +
  geom_point(data = subset(plt, show.label == T), pch=19,cex=1.5) +
  scale_color_manual(values = c('TRUE'='red', 'FALSE' = 'black') ) +
  geom_text_repel(data = subset(plt, show.label), size = 3.4 )  +
  labs(x = "log2FC R1 vs. R2 (unpaired)",
       y="Correlation t-statistic with tumour percentage",
       col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme




# c("ENSG00000168137.16_5|SETD5|chr3:9439299-9520924(+)", "ENSG00000167785.9_3|ZNF558|chr19:8916846-8942990(-)", "ENSG00000171649.12_4|ZIK1|chr19:58089824-58105258(+)", "ENSG00000198521.7|ZNF43|chr19:21987752-22034927(-)", "ENSG00000204514.10_5|ZNF814|chr19:58360099-58400442(-)", "ENSG00000175414.7_4|ARL10|chr5:175792490-175828866(+)", "ENSG00000066117.14_4|SMARCD1|chr12:50478755-50494495(+)", "ENSG00000106344.8_5|RBM28|chr7:127937738-127983962(-)", "ENSG00000188321.13_4|ZNF559|chr19:9434448-9461838(+)", "ENSG00000078177.14_4|N4BP2|chr4:40058470-40159872(+)", "ENSG00000197299.12_4|BLM|chr15:91260577-91359396(+)", "ENSG00000008311.15_3|AASS|chr7:121713603-121784309(-)", "ENSG00000197128.11_3|ZNF772|chr19:57978031-57988938(-)", "ENSG00000138443.16_5|ABI2|chr2:204192962-204312451(+)", "ENSG00000105708.9_3|ZNF14|chr19:19821281-19843921(-)", "ENSG00000197928.10_3|ZNF677|chr19:53738634-53758151(-)", "ENSG00000198799.12_4|LRIG2|chr1:113615820-113674882(+)", "ENSG00000131115.16_4|ZNF227|chr19:44711700-44741421(+)", "ENSG00000103037.11_3|SETD6|chr16:58549383-58555085(+)", "ENSG00000204519.11_4|ZNF551|chr19:58193337-58228669(+)", "ENSG00000198551.10_5|ZNF627|chr19:11670189-11729976(+)", "ENSG00000189079.17_5|ARID2|chr12:46123489-46301823(+)", "ENSG00000007392.16_4|LUC7L|chr16:238968-279462(-)", "ENSG00000135164.18_5|DMTF1|chr7:86781677-86825653(+)", "ENSG00000105866.15_5|SP4|chr7:21467661-21554440(+)", "ENSG00000164828.18_7|SUN1|chr7:856252-936072(+)", "ENSG00000105486.14_6|LIG1|chr19:48618702-48673860(-)", "ENSG00000167637.17_4|ZNF283|chr19:44331473-44356169(+)", "ENSG00000254004.7_3|ZNF260|chr19:37001589-37019173(-)", "ENSG00000177839.6_5|PCDHB9|chr5:140566701-140571114(+)", "ENSG00000162086.15_4|ZNF75A|chr16:3355406-3368852(+)", "ENSG00000134744.14_6|TUT4|chr1:52873954-53019159(-)", "ENSG00000120784.16_5|ZFP30|chr19:38104650-38183238(-)", "ENSG00000147274.14_4|RBMX|chrX:135930163-135962923(-)", "ENSG00000129351.17_3|ILF3|chr19:10764937-10803093(+)", "ENSG00000151612.17_6|ZNF827|chr4:146678779-146859975(-)", "ENSG00000166704.11_3|ZNF606|chr19:58488421-58514717(-)", "ENSG00000106443.16_3|PHF14|chr7:11012963-11209257(+)", "ENSG00000004139.14_4|SARM1|chr17:26691378-26731067(+)", "ENSG00000122779.18_4|TRIM24|chr7:138145004-138274741(+)", "ENSG00000160961.12_3|ZNF333|chr19:14800613-14844558(+)", "ENSG00000167380.16_4|ZNF226|chr19:44669226-44682534(+)", "ENSG00000122970.16_4|IFT81|chr12:110562140-110656598(+)", "ENSG00000071575.11_2|TRIB2|chr2:12857015-12882860(+)", "ENSG00000263001.6_5|GTF2I|chr7:74064571-74175022(+)", "ENSG00000236609.4_3|ZNF853|chr7:6655241-6663921(+)", "ENSG00000197647.11_5|ZNF433|chr19:12125547-12146556(-)", "ENSG00000104885.18_3|DOT1L|chr19:2163932-2232577(+)", "ENSG00000113387.12_4|SUB1|chr5:32531739-32604185(+)", "ENSG00000167384.10_2|ZNF180|chr19:44978645-45004576(-)")
# c("ENSG00000111885.7_4|MAN1A1|chr6:119498370-119670926(-)", "ENSG00000141506.13_4|PIK3R5|chr17:8782228-8869029(-)", "ENSG00000166272.18_5|WBP1L|chr10:104503705-104594273(+)", "ENSG00000198624.13_3|CCDC69|chr5:150560613-150603654(-)", "ENSG00000197746.14_4|PSAP|chr10:73576055-73611082(-)", "ENSG00000151726.14_4|ACSL1|chr4:185676749-185747972(-)", "ENSG00000107968.10_3|MAP3K8|chr10:30722950-30750762(+)", "ENSG00000266412.5_3|NCOA4|chr10:51565108-51590734(+)", "ENSG00000131370.16_5|SH3BP5|chr3:15295860-15382875(-)", "ENSG00000100365.16_5|NCF4|chr22:37257030-37274059(+)", "ENSG00000082397.17_7|EPB41L3|chr18:5392380-5630699(-)", "ENSG00000204161.14_6|TMEM273|chr10:50362770-50396630(-)", "ENSG00000108639.7_5|SYNGR2|chr17:76164639-76169608(+)", "ENSG00000248905.8_6|FMN1|chr15:33057746-33486934(-)", "ENSG00000084070.12_5|SMAP2|chr1:40810522-40888998(+)", "ENSG00000110079.18_5|MS4A4A|chr11:59953175-60085417(+)", "ENSG00000155252.13_3|PI4K2A|chr10:99400443-99436191(+)", "ENSG00000134996.11_2|OSTF1|chr9:77703459-77762181(+)", "ENSG00000130775.16_3|THEMIS2|chr1:28199054-28213196(+)", "ENSG00000122359.18_5|ANXA11|chr10:81910645-81965328(-)", "ENSG00000110324.10_3|IL10RA|chr11:117857063-117873752(+)", "ENSG00000165457.14_4|FOLR2|chr11:71927645-71932994(+)", "ENSG00000155926.14_4|SLA|chr8:134048973-134115156(-)", "ENSG00000134516.17_6|DOCK2|chr5:169064272-169510386(+)", "ENSG00000173372.17_4|C1QA|chr1:22963121-22966171(+)", "ENSG00000163131.11_4|CTSS|chr1:150702664-150738268(-)", "ENSG00000153071.15_5|DAB2|chr5:39371777-39462402(-)", "ENSG00000101336.14_6|HCK|chr20:30639991-30689659(+)", "ENSG00000175155.9_4|YPEL2|chr17:57409016-57479090(+)", "ENSG00000128805.14_4|ARHGAP22|chr10:49654077-49864310(-)", "ENSG00000148180.19_5|GSN|chr9:123970072-124095121(+)", "ENSG00000142185.16_6|TRPM2|chr21:45770046-45862964(+)", "ENSG00000111912.20_7|NCOA7|chr6:126102307-126253180(+)", "ENSG00000155659.15_5|VSIG4|chrX:65241580-65259967(-)", "ENSG00000137462.8_4|TLR2|chr4:154605222-154626854(+)", "ENSG00000180353.11_3|HCLS1|chr3:121350246-121379774(-)", "ENSG00000147459.18_4|DOCK5|chr8:25042204-25275598(+)", "ENSG00000136250.11_5|AOAH|chr7:36552557-36764154(-)", "ENSG00000235568.6_2|NFAM1|chr22:42776416-42828401(-)", "ENSG00000198879.11_3|SFMBT2|chr10:7200586-7453448(-)", "ENSG00000197142.10_2|ACSL5|chr10:114133776-114188138(+)", "ENSG00000143119.14_4|CD53|chr1:111413810-111442544(+)", "ENSG00000138964.17_5|PARVG|chr22:44568836-44615413(+)", "ENSG00000167613.16_5|LAIR1|chr19:54862991-54882165(-)", "ENSG00000183484.12_5|GPR132|chr14:105515732-105531782(-)", "ENSG00000130830.15_5|MPP1|chrX:154006959-154049282(-)", "ENSG00000101160.14_3|CTSZ|chr20:57570240-57582309(-)", "ENSG00000205744.10_4|DENND1C|chr19:6467218-6482568(-)", "ENSG00000100368.14_6|CSF2RB|chr22:37309670-37336481(+)", "ENSG00000012779.11_3|ALOX5|chr10:45869624-45941567(+)")

results.out %>%
  dplyr::arrange(-statistic.gsam.cor.tpc) %>%
  dplyr::select(c(gid, ensembl_id, hugo_symbol, statistic.gsam.cor.tpc)) %>%
  dplyr::slice_head(n=50) %>%
  dplyr::pull(ensembl_id)

results.out %>%
  dplyr::arrange(statistic.gsam.cor.tpc) %>%
  dplyr::select(c(gid, ensembl_id, hugo_symbol, statistic.gsam.cor.tpc)) %>%
  dplyr::slice_head(n=50) %>%
  dplyr::pull(ensembl_id)




df %>%
  tibble::rownames_to_column('loc') %>% 
  dplyr::mutate(chr = as.character(gsub("^(chr[^:]):.+$","\\1",loc) ))  %>%
  dplyr::mutate(start = as.numeric(   gsub("^chr[^:]:([0-9]+).+$","\\1",loc)   )) %>%
  dplyr::arrange(chr, start)


## N NPC1 NPC2 N :: Figure S6 A-B ----
# a. eerst losse PCA bepaling, dan correlatie daar tussen?


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


labels <- results.out %>% 
  dplyr::mutate(C1 = hugo_symbol %in% neuron.genes) %>%
  dplyr::filter(!is.na(gid) ) %>% 
  dplyr::filter(C1 | neftel.meta.module.NPC1 | neftel.meta.module.NPC2  ) %>% 
  dplyr::select(gid, ensembl_id, hugo_symbol, C1, neftel.meta.module.NPC1 , neftel.meta.module.NPC2, McKenzie_celltype_top_human_specificity) %>% 
  dplyr::mutate(col = "") %>% 
  dplyr::mutate(col = ifelse(C1, "C1", col)) %>% 
  dplyr::mutate(col = ifelse(neftel.meta.module.NPC1, paste0(col,",NPC1"), col)) %>% 
  dplyr::mutate(col = ifelse(neftel.meta.module.NPC2, paste0(col,",NPC2"), col)) %>% 
  dplyr::mutate(col = as.factor(gsub("^[,]+","",col))) %>% 
  dplyr::mutate(label = ifelse(neftel.meta.module.NPC1 | neftel.meta.module.NPC2 | McKenzie_celltype_top_human_specificity == "neuron",hugo_symbol, ""))

#factoextra::fviz_eig(res.pca)
# 
# factoextra::fviz_pca_biplot(res.pca, repel = TRUE,
#                 col.var = labels$col, # Variables color
#                 col.ind = "#696969"  # Individuals color 
# )




plt <- labels %>% 
  dplyr::select('gid', 'hugo_symbol') %>% 
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::mutate(gid = NULL) %>% 
  tibble::column_to_rownames('hugo_symbol')
res.pca <- prcomp(t(plt), scale = TRUE)
factoextra::fviz_pca_var(res.pca,
                         col.var = labels$col,
                         repel = TRUE   ) +
  youri_gg_theme +
  scale_color_manual(values=c('C1'='#ff5f68',
                              'C1,NPC2'='purple',
                              'NPC1'='#6ba6e5',
                              'NPC2'='#6be5d9',
                              'NPC1,NPC2'='#6bcee5'
                              ))


ggsave("output/figures/figure_S6_C1_NPC1_NPC2_biplot.pdf", width=7, height=6)





plt <- labels %>% 
  dplyr::select('gid') %>% 
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  tibble::column_to_rownames('gid')
res.pca <- prcomp(t(plt), scale = TRUE)
stopifnot(rownames(res.pca$rotation) == labels$gid)

tmp <- res.pca$rotation %>% 
  as.data.frame %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::left_join(labels, by=c('gid'='gid'))

aa <- 0.7
ggplot(tmp, aes(x=PC1, y = PC2, fill=col, label=label)) +
  geom_point(size=2.75, pch=21, col="black") +
  geom_text_repel(size=2.25) + 
  labs(x = "PC1 loadings", y = "PC2 loadings", title="PCA on G-SAM VST expression of C1 and NPC1 & NPC2 (Neftel) genes") + 
  youri_gg_theme +
  scale_fill_manual(values=c('C1'=alpha('#ff5f68',aa),
                              'C1,NPC2'=alpha('purple',aa),
                              'NPC1'=alpha('#6ba6e5',aa),
                              'NPC2'=alpha('#6be5d9',aa),
                              'NPC1,NPC2'=alpha('#6bcee5',aa)))


ggsave("output/figures/figure_S6_C1_NPC1_NPC2_PC-loadings_scatterplot.pdf", width=7, height=7.5)



plt <- plt %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::left_join(labels %>% dplyr::select('gid','hugo_symbol'), by=c('gid'='gid')) %>% 
  dplyr::mutate(gid=NULL) %>% 
  tibble::column_to_rownames('hugo_symbol')


labels2 <- labels %>% 
  dplyr::select(c('hugo_symbol','C1','McKenzie_celltype_top_human_specificity','neftel.meta.module.NPC1','neftel.meta.module.NPC2')) %>% 
  
  dplyr::mutate(`Neuron Marker` = ifelse(is.na(McKenzie_celltype_top_human_specificity), "NA", McKenzie_celltype_top_human_specificity)) %>% 
  dplyr::mutate(`Neuron Marker` = `Neuron Marker` == "neuron") %>% 
  dplyr::mutate(`Neuron Marker` = ifelse(`Neuron Marker` == T, NA , F)) %>% 
  dplyr::mutate(McKenzie_celltype_top_human_specificity = NULL) %>% 
  
  dplyr::rename(`Neftel.NPC1` = neftel.meta.module.NPC1) %>% 
  dplyr::mutate(Neftel.NPC1 = ifelse(Neftel.NPC1 == T, NA , F)) %>% 

  dplyr::rename(`Neftel.NPC2` = neftel.meta.module.NPC2) %>% 
  dplyr::mutate(Neftel.NPC2 = ifelse(Neftel.NPC2 == T, NA , F)) %>% 
  
  dplyr::mutate(C1 = ifelse(C1 == T, NA , F)) %>% 
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(`Neuron Marker` = NULL)


cor_cor_plot(plt, labels2)
ggsave("output/figures/figure_S6_C1_NPC1_NPC2_corr-corr-plot.pdf", width=7*3, height=7.5*3)
ggsave("output/figures/figure_S6_C1_NPC1_NPC2_corr-corr-plot.svg", width=7*3, height=7.5*3)




## OD OPC :: Figure S6 C-D ----
# a. eerst losse PCA bepaling, dan correlatie daar tussen?


oligodendrocyte.genes <- c("RASGRF1","FBXO2","PPP1R16B","TPPP","SEC14L5","TMEM151A","LGI3","TMCC2","HHATL","RAB11FIP4","PDIA2",
                           "HCN2","KCNJ9","DNAJC6","TUBB4A","ADCY5","CSDC2","AC118754.1","PLIN4","HSPA2","PI16","PTGDS","CDK18",
                           "FA2H","AATK","NKX6-2","MAG","PLCH2","FAM131C","PPP1R14A","CHADL","TMEM88B","ABCA2","PLPP2","CERCAM",
                           "BOK","ACP7","CYS1","ANXA3","TNFSF9","AVPI1","MYH11","ADH1B","TUBA4A","CORO6","IL12RB2","TESPA1",
                           "MPP7","RSPO3","KCNJ12","OPN4","MKX","FRAS1","CPNE7")

labels <- results.out %>% 
  dplyr::mutate(C2 = hugo_symbol %in% oligodendrocyte.genes) %>%
  dplyr::filter(!is.na(gid) ) %>% 
  dplyr::filter(C2 | neftel.meta.module.OPC  ) %>% 
  dplyr::select(gid, ensembl_id, hugo_symbol, C2, neftel.meta.module.OPC, McKenzie_celltype_top_human_specificity) %>% 
  dplyr::mutate(col = "") %>% 
  dplyr::mutate(col = ifelse(C2, "C2", col)) %>% 
  dplyr::mutate(col = ifelse(neftel.meta.module.OPC, paste0(col,",OPC"), col)) %>% 
  dplyr::mutate(col = as.factor(gsub("^[,]+","",col))) %>% 
  dplyr::mutate(label = ifelse(neftel.meta.module.OPC | McKenzie_celltype_top_human_specificity == "oligodendrocyte",
                               hugo_symbol, ""))





plt <- labels %>% 
  dplyr::select('gid', 'hugo_symbol') %>% 
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  dplyr::mutate(gid = NULL) %>% 
  tibble::column_to_rownames('hugo_symbol')
res.pca <- prcomp(t(plt), scale = TRUE)
factoextra::fviz_pca_var(res.pca,
                         col.var = labels$col,
                         repel = TRUE   ) +
  youri_gg_theme + 
  scale_color_manual(values=c('C2'='#ff5f68',
                              'OPC'='#6ba6e5'
  ))



ggsave("output/figures/figure_S6_C2_OPC_biplot.pdf", width=7, height=6)





plt <- labels %>% 
  dplyr::select('gid') %>% 
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  tibble::column_to_rownames('gid')
res.pca <- prcomp(t(plt), scale = TRUE)
stopifnot(rownames(res.pca$rotation) == labels$gid)

tmp <- res.pca$rotation %>% 
  as.data.frame %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::left_join(labels, by=c('gid'='gid'))

ggplot(tmp, aes(x=PC1, y = PC2, fill=col, label=hugo_symbol)) +
  geom_point(size=2.75, pch=21, col="black") +
  geom_text_repel(size=2.25) + 
  labs(x = "PC1 loadings", y = "PC2 loadings", title="PCA on G-SAM VST expression of C2 and OPC (Neftel) genes") + 
  youri_gg_theme +
  scale_fill_manual(values=c('C2'=alpha('#ff5f68',aa),
                             'OPC'=alpha('#6ba6e5',aa)
                    ))

ggsave("output/figures/figure_S6_C2_OPC_PC-loadings_scatterplot.pdf", width=7, height=7.5)




plt <- plt %>% 
  tibble::rownames_to_column('gid') %>% 
  dplyr::left_join(labels %>% dplyr::select('gid','hugo_symbol'), by=c('gid'='gid')) %>% 
  dplyr::mutate(gid=NULL) %>% 
  tibble::column_to_rownames('hugo_symbol')


labels2 <- labels %>% 
  dplyr::select(c('hugo_symbol','C2','McKenzie_celltype_top_human_specificity','neftel.meta.module.OPC')) %>% 

  dplyr::mutate(`Oligodendrocyte Marker` = ifelse(is.na(McKenzie_celltype_top_human_specificity), "NA", McKenzie_celltype_top_human_specificity)) %>% 
  dplyr::mutate(`Oligodendrocyte Marker` = `Oligodendrocyte Marker` == "oligodendrocyte") %>% 
  dplyr::mutate(`Oligodendrocyte Marker` = ifelse(`Oligodendrocyte Marker` == T, NA , F)) %>% 
  dplyr::mutate(McKenzie_celltype_top_human_specificity = NULL) %>% 

  dplyr::rename(`Neftel.OPC` = neftel.meta.module.OPC) %>% 
  dplyr::mutate(Neftel.OPC = ifelse(Neftel.OPC == T, NA , F)) %>% 
  
  dplyr::mutate(C2 = ifelse(C2 == T, NA , F)) %>% 
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(`Oligodendrocyte Marker` = NULL)


cor_cor_plot(plt, labels2)
ggsave("output/figures/figure_S6_C2_OPC_corr-corr-plot.pdf", width=7*1.5, height=7.5*1.5)


