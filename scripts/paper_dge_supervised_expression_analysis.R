#!/usr/bin/env R

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



# load data ----

#source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
#ensembl_genes <- get_ensembl_hsapiens_gene_ids()

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


source("scripts/R/ligands.R")
source("scripts/R/subtype_genes.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")

source('scripts/R/wang_glioma_intrinsic_genes.R')

source("scripts/R/glass_expression_matrix.R") # glass & tcga validation set


# prepare data ----

## GSAM ----


# gsam.metadata.all
# gsam.metadata.r1
# gsam.metadata.r2

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
wilcox.test(gsam.metadata.all %>%
              dplyr::filter(resection == "r1") %>%
              dplyr::pull(tumour.percentage.dna) ,
            gsam.metadata.all %>%
              dplyr::filter(resection == "r2") %>%
              dplyr::pull(tumour.percentage.dna)
            , alternative = "two.sided")



# gsam.gene.expression.all
# gsam.bfg.expression.all

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



glass.metadata.all  <- glass.gbm.rnaseq.metadata %>%
  dplyr::mutate(condition = ifelse(resection == "TP","Primary","NotPrimary")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Primary","NotPrimary") )) %>%
  dplyr::filter(idh_status == "IDHwt") %>% # exclude the three IDH mutants according to Synapse WGS/WES VCF files
  dplyr::filter(pid %in%  (glass.gbm.rnaseq.metadata %>%
                             dplyr::filter(grade %in% c('II','III')) %>%
                             dplyr::pull(pid)) == F )



glass.gene.expression.all <- glass.gbm.rnaseq.expression %>%
  dplyr::select(glass.metadata.all$sid)


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)




## per-gene results table ----


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
  dplyr::slice_head(n=100) %>%
  dplyr::filter(gene %in% unique(.[["gene"]][duplicated(.[["gene"]])]) == F) %>% # only those that match 1 cell type
  dplyr::ungroup() %>%
  dplyr::add_row(gene = "CD163",Celltype = "microglia/TAM") %>%
  dplyr::add_row(gene = "RBFOX3",Celltype = "neuron") %>%
  dplyr::add_row(gene = "FLT1",Celltype = "endothelial") %>%
  dplyr::mutate(show.marker = ifelse(gene %in% c(
    
    "TBX3","ERG", "FLT1","RERGL","VCAM1" , # endothelial
    
    "BCAS1","BLNK","ERBB3","FGFR3","GFAP","GABRA1","GABRB2","OPALIN","PLP1","RBFOX3","RUNX1","SOX8","SOX9","TGFA","VIP") | grepl("^CD[0-9]", gene) , T , F) ) %>%
  dplyr::rename(McKenzie_celltype_top_human_specificity = Celltype) %>%
  dplyr::mutate(grand_mean = NULL)
#dim(tmp)

results.out$McKenzie_celltype_top_human_specificity <- NULL
results.out$show.marker <- NULL
results.out <- results.out %>%
  dplyr::left_join(tmp, by=c('hugo_symbol'='gene'))
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



# corr TPC [GSAM] ----


stopifnot(colnames(gsam.gene.expression.all.vst) == gsam.metadata.all$sid)
gsam.gene.expression.all.cor.estimate <- data.frame(apply(gsam.gene.expression.all.vst, 1, function (x)  cor.test(x, gsam.metadata.all$tumour.percentage.dna) %>% purrr::pluck('estimate') ), stringsAsFactors = F) %>%
  `colnames<-`("estimate") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor.statistic <- data.frame(apply(gsam.gene.expression.all.vst, 1, function (x)  cor.test(x, gsam.metadata.all$tumour.percentage.dna) %>% purrr::pluck('statistic') ), stringsAsFactors = F) %>%
  `colnames<-`("statistic") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor.p.value <- data.frame(apply(gsam.gene.expression.all.vst, 1, function (x)  cor.test(x, gsam.metadata.all$tumour.percentage.dna) %>% purrr::pluck('p.value') ), stringsAsFactors = F) %>%
  `colnames<-`("p.value") %>%
  tibble::rownames_to_column('gid')

gsam.gene.expression.all.cor <- gsam.gene.expression.all.cor.estimate %>%
  dplyr::left_join(gsam.gene.expression.all.cor.statistic , by=c('gid' = 'gid') ) %>%
  dplyr::left_join(gsam.gene.expression.all.cor.p.value , by=c('gid' = 'gid') ) %>%
  `colnames<-`(paste0(colnames(.), ".gsam.cor.tpc")) %>%
  dplyr::rename(gid = gid.gsam.cor.tpc)

rm(gsam.gene.expression.all.cor.estimate , gsam.gene.expression.all.cor.statistic , gsam.gene.expression.all.cor.p.value)



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






# DE unpaired all [G-SAM] ----



gsam.gene.res.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                              dplyr::filter(rowSums(.) > ncol(.) * 3)
                                            , gsam.metadata.all, ~resection ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>%
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  `colnames<-`(paste0(colnames(.),".gsam.res")) %>%
  dplyr::rename(gid = gid.gsam.res)

results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.res , by = c('gid' = 'gid'))



# TODO remove .paired twice:
gsam.gene.res.tpc.res <- DESeqDataSetFromMatrix(gsam.gene.expression.all %>%
                                                  dplyr::filter(rowSums(.) > ncol(.) * 3),
                                                gsam.metadata.all, ~tpc + resection ) %>% # + resection; corrected for tpc
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5 ) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('gid') %>%
  #dplyr::mutate(ensembl_id = gsub('\\..+\\|.+$','',gid) ) %>%
  #dplyr::mutate(hugo_symbol = gsub("^.+\\|([^\\]+)\\|.+$","\\1",gid) ) %>%
  `colnames<-`(paste0(colnames(.),".gsam.tpc.res")) %>%
  dplyr::rename(gid = gid.gsam.tpc.res)

results.out <- results.out %>%
  dplyr::left_join(gsam.gene.res.tpc.res , by = c('gid' = 'gid'))



# DE unpaired all [GLASS] ----

# no per-sample tumor percentage available in Synapse portal


stopifnot(colnames(glass.gene.expression.all) == glass.metadata.all$sid)

glass.gene.res.res <- DESeqDataSetFromMatrix(glass.gene.expression.all %>%
                                               dplyr::filter(rowSums(.) > ncol(.) * 3)
                                             , glass.metadata.all, ~condition ) %>% # + resection
  DESeq(parallel = T) %>%
  results() %>% 
  as.data.frame() %>%
  dplyr::mutate(significant = padj <= 0.01 & abs(log2FoldChange) > 0.5) %>%
  dplyr::arrange(padj) %>%
  tibble::rownames_to_column('ensembl_id') %>% 
  `colnames<-`(paste0(colnames(.),".glass.res")) %>%
  dplyr::rename(ensembl_id = ensembl_id.glass.res)


results.out <- results.out %>%
  dplyr::left_join(glass.gene.res.res, by = c('ensembl_id' = 'ensembl_id'))


# hypeR enrichment ---- 

# ensembl to entrez?
# data(examplePathways)
# data(exampleRanks)
# # https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html


## https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

genesets <- list()

y <- results.out %>% 
  dplyr::filter(!is.na(stat.gsam.tpc.res)) %>%
  dplyr::filter(!is.na(hugo_symbol)) %>%
  dplyr::filter(lfcSE.gsam.tpc.res < 0.3) %>%
  #dplyr::filter(baseMean >= 10) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(hugo_symbol) %>% 
  dplyr::summarize(stat.gsam.tpc.res=mean(stat.gsam.tpc.res))

y.ordered <- (y %>% dplyr::select(hugo_symbol, stat.gsam.tpc.res) %>% dplyr::arrange(stat.gsam.tpc.res)) %>% tibble::deframe()
y.ordered.abs <- (y %>% dplyr::select(hugo_symbol, stat.gsam.tpc.res) %>% dplyr::arrange(abs(stat.gsam.tpc.res))) %>% tibble::deframe()


for(a in names(genesets)) {
  for(b in names(genesets[[a]])) {
    if(grepl("blood vessel",b,ignore.case = T) |   grepl("adasd vascul",b,ignore.case = T) ) {
      print(a)
      print(b)
      print("")
    }
  }
}


grepl("asd","AsD",ignore.case = T)


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

## 2.1 :: TPC violin ----


m1 <- gsam.metadata.all %>%
  dplyr::filter(!is.na(tumour.percentage.dna) & resection == "r1") %>%
  dplyr::pull(tumour.percentage.dna) %>%
  median()

m2 <- gsam.metadata.all %>%
  dplyr::filter(!is.na(tumour.percentage.dna) & resection == "r2") %>%
  dplyr::pull(tumour.percentage.dna) %>%
  median()


ggplot(gsam.metadata.all %>% dplyr::mutate(resection = ifelse(resection == 'r1', "Resection 1", "Resection 2"))
       , aes(x = resection, y = tumour.percentage.dna, col=resection)) +
  geom_violin() +
  geom_segment(aes(x=0.55, y=m1, xend=1.45, yend=m1),lty=1,lwd=0.15, col="gray40") +
  geom_segment(aes(x=1.55, y=m2, xend=2.45, yend=m2),lty=1,lwd=0.15, col="gray40") +
  geom_jitter( position=position_jitter(0.2), size=0.9) +
  ylim(0, 100) +
  labs(x = NULL, col = NULL, y = "WES estimate tumour cell percentage" ) +
  job_gg_theme
rm(m1,m2)


ggsave('output/figures/tpc_estimate_wes.png', width = 5.7, height = 4.3)



## 2.3 :: FC res x corr TPC G-SAM  [diag; chr7 ; chr10 ] ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::mutate(is.limited.gsam.res = as.character(log2FoldChange.gsam.res > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2, 2 , log2FoldChange.gsam.res)) %>%
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
  xlim(-2,2) +
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
  xlim(-2,2) +
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
  xlim(-2,2) +
  ylim(-16.5,16.5)

p1 / p2 / p3


ggsave('output/figures/paper_dge_log2FoldChange.res_x_statistic.cor.tpc.png', width = 5.7 * 1.4 , height = 4 * 1.4 * 3)




## 2.5 :: G-SAM.res & GLASS.res ----


### a :: GLASS res x G-SAM TPC corr ----


plt <- results.out %>%
  dplyr::filter(!is.na(log2FoldChange.glass.res) & !is.na(statistic.gsam.cor.tpc)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(log2FoldChange.glass.res) > 2.5)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res > 2.5, 2.5 , log2FoldChange.glass.res)) %>%
  dplyr::mutate(log2FoldChange.glass.res = ifelse(log2FoldChange.glass.res < -2.5, -2.5 , log2FoldChange.glass.res)) 

p1 <- ggplot(plt, aes(x = log2FoldChange.glass.res ,
                      y =  statistic.gsam.cor.tpc,
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



### b :: GLASS ----


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




### A :: Corrected ----


p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res ,
                y=statistic.gsam.cor.tpc ) ) + 
  geom_point(cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.05 &  is.limited.gsam.tpc.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  
  scale_color_manual(values = c('remaining'='black','chr7'='red') ) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       #,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-2.5,2.5)



### B :: Un-corrected ----


p2 <- ggplot(plt, aes(x=log2FoldChange.gsam.res ,
                y=statistic.gsam.cor.tpc ) ) + 
  geom_point(cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.res > 0.05 &  is.limited.gsam.res == "FALSE"),method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="red" , size=0.8) +
  scale_color_manual(values = c('remaining'='black','chr7'='red') ) +
  labs(x = "log2FC R1 vs. R2",
       y="Correlation t-statistic with tumour percentage"
       #,col="Difference significant (R1 ~ R2)"
  ) +
  youri_gg_theme + 
  xlim(-2.5,2.5)


p1 + p2


ggsave('output/figures/paper_dge_gsam_corrected_uncorrected.png', width = 5.7 * 1.1 * 2 , height = 4 * 1.6 )



## 2.9 Corrected LFC + GAB(R)A labels ----


plt <- results.out %>%
  #dplyr::filter(is.na(padj.glass.res) | (!is.na(padj.glass.res) & padj.glass.res  < 0.01) ) %>%
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  #dplyr::mutate(show.label = significant & abs(log2FoldChange.gsam.tpc.res) > 1.5) %>%
  dplyr::mutate(show.label.neuronal.label = hugo_symbol %in% c("RBFOX3", "SATB2", "SLC17A7", "RORB", "GAD1", "GAD2", "SST", "LHX6","ADARB2","VIP")) %>%
  dplyr::mutate(show.label.neuronal = hugo_symbol %in% c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3","GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST","VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5","NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK","ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A","RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87","ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2","DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2","CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4","ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2","CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1","GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2","RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")) %>%
  dplyr::mutate(show.label.gaba = hugo_symbol %in% c("GABBR2","GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6","GABRB1","GABRB2","GABRB3","GABRD","GABRE","GABRG1","GABRG2","GABRG3","GABRP","GABRQ","GABRR1","GABRR2","GABRR3")) %>%
  dplyr::mutate(show.label = as.factor ( case_when(
    show.label.gaba == T ~ 'GABA',
    show.label.neuronal == T & show.label.gaba == F ~ 'neuronal',
    significant == T ~ 'significant',
    T ~ 'not significant')))

p1 <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res ,
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



ggsave("/tmp/gabra.png",height=10 * 1.3,width=4.5 * 1.3)




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
  dplyr::filter(!is.na(log2FoldChange.gsam.res) & !is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  dplyr::mutate(show.label.gains = hugo_symbol %in% c(
    # "MGMT",
    # "HER2","RAS","PTEN","MTOR","TP53","GFAP","TAP1","TGFB1","RB1","NF1","ATRX",
    # "CDK2","PIK3CA", "BRAF","PIK3R1","TERT","TERC",
    # "MCM","MCM2", "MCM4",
    # "SOX2", "SOX9", "FGFR3","CCNE1", "AKT3","LSAMP","CCND1","CCND2", "PIK3C2B","CDK6","CDK4","MDM2","MDM4",
    "AKT1","AKT3","EGFR", "PDGFRA","MET", "PIK3C2B", "MDM2","MDM4", "CDK4","CDK6","SOX2","FGFR3","MYCN","MYC","CCND1","CCND2","BMI1" # gains
  )) %>%
  dplyr::mutate(show.label.other = hugo_symbol %in% c( 
    "TP53", "PIK3CA","PIK3R1","NF1","SPTA1","GABRA6","ABCC6","CXCL12",
    #"SRC","PTK2",
    #CXCL12=https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5193023/ = strongly associated w/ non tumour cells
    # 'CDK4','MDM2',"EGFR","PDGFRA","OCT4","MGMT"# "MET","HER2","RAS","PTEN","MTOR","TP53","GFAP","TAP1","TGFB1","RB1","NF1","ATRX",
    # "CDK2","PIK3CA", "BRAF","PIK3R1","TERT","TERC",
    # "MCM","MCM2", "MCM4",
    "LTBP4", "TGFB1","PREX1","MSH6", "MSH2", "MLH1","VEGFA", "PTPN11","STAT3","DCC", "MGMT"
    # "MSH6" =- Recurrent Glioblastoma: From Molecular Landscape to New Treatment Perspectives
    # DCC = progression marker? - https://www.intechopen.com/books/neurooncology-newer-developments/genetic-alterations-of-glioblastoma + DOI:10.1371/journal.pone.0025408
  )) %>%
  dplyr::mutate(show.label.losses = hugo_symbol %in% c(
    "CDKN2A","CDKN2B","PTEN","RB1","NF1","QKI","CDKN2C","TP53","MTAP","ELAVL2","PARK2"
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3194060/
  ))


p1a <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
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
p1b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.other == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.other == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.other == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.other == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Other genes of interest in GBM") +
  youri_gg_theme + xlim(-3, 3)


p2a <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
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
p2b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.gains == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.gains == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.gains == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.gains == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Commonly gained genes in GBM") +
  youri_gg_theme + xlim(-3, 3)



p3a <- ggplot(plt, aes(x=log2FoldChange.gsam.tpc.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
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
p3b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=significant, label = hugo_symbol )) + 
  geom_point(data=subset(plt, significant != T),cex=0.35) +
  geom_point(data=subset(plt, significant == T),cex=0.45) +
  geom_point(data=subset(plt, show.label.losses == T & significant == F ),col="black",cex=0.65) +
  geom_point(data=subset(plt, show.label.losses == T & significant == T ),col="red",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label.losses == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label.losses == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Commonly lost genes in GBM") +
  youri_gg_theme + xlim(-3, 3)


(p2a + p1a + p3a) / (p2b + p1b + p3b) 



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
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Neuron marker genes") +
  youri_gg_theme + xlim(-2, 2)
p1b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Neuron marker genes") +
  youri_gg_theme + xlim(-3, 3)



plt <- plt %>% dplyr::mutate(show.label = show.label.endothelial)
p2a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Endothelial marker genes") +
  youri_gg_theme + xlim(-2, 2)
p2b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.glass.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.glass.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Endothelial marker genes") +
  youri_gg_theme + xlim(-3, 3)



plt <- plt %>% dplyr::mutate(show.label = show.label.endothelial)
p3a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="TAM/Microglia marker genes") +
  youri_gg_theme + xlim(-2, 2)
p3b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="TAM/Microglia marker genes") +
  youri_gg_theme + xlim(-3, 3)



plt <- plt %>% dplyr::mutate(show.label = show.label.oligodendrocytes)
p4a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Oligodendrocyte marker genes") +
  youri_gg_theme + xlim(-2, 2)
p4b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Oligodendrocyte marker genes") +
  youri_gg_theme + xlim(-3, 3)


plt <- plt %>% dplyr::mutate(show.label = show.label.astrocytes)
p5a <- ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.gsam.tpc.res > 0.1 &  is.limited.gsam.res == "FALSE"),method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5) , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & show.marker == T & log2FoldChange.gsam.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 [G-SAM]",
       y="Corr t-stat tumour-%"
       ,col="Astrocyte marker genes") +
  youri_gg_theme + xlim(-2, 2)
p5b <- ggplot(plt, aes(x=log2FoldChange.glass.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_smooth(data = subset(plt, padj.glass.res > 0.1 &  is.limited.gsam.res == "FALSE"), method="lm", se = FALSE,  formula=y ~ x, orientation="y", col=rgb(1,0,0,0.5)  , size=0.4) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2/3/4 [GLASS]",
       y="Corr t-stat tumour-%"
       ,col="Astrocyte marker genes") +
  youri_gg_theme + xlim(-3, 3)



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
(p5a + p5b)


ggsave("output/figures/paper_dge_cell-type_genes.png",width=2*6, height=6*6)




## 2.14 Corrected LFC + individual tophits ----



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
  
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2, 2 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2, -2 , log2FoldChange.gsam.tpc.res)) %>%
  
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) %>%

  dplyr::mutate(significant = padj.gsam.tpc.res < 0.01 & lfcSE.gsam.tpc.res < 0.3 & abs(log2FoldChange.gsam.tpc.res) > 0.5 ) %>%
  dplyr::mutate(show.label = hugo_symbol %in% c(
    #'GLIS1',"CPEB1",'GLI1',
    
    #  plt  %>% arrange(pvalue.glass.res  * pvalue.gsam.tpc.res) %>% filter (significant) %>% head(n=25) %>% pull(hugo_symbol)
    
    #"MME","CORO6","CCBE1","SCUBE3","MMP11","CYSLTR2","RBP4","HS3ST2","TNNT2","SYN1","LARP6","CERCAM","PPP4R4","AK5","SLIT3","CABP1","RYR2","CAMKK1","SNCG","ATP8A2","REPS2","RSPO3","PTGFR","GPR83","PHYHIP"
    
    # SYN1 = synaptic vesicle & neurotransmitter (serotonin) [related to DSCAM DCC and EGFR?][Multiple epidermal growth factor-like domains protein 5]
    # MME = in pathway met VIP (sterke neuron marker)
    # SLIT3 = axon guidance    
    # HS3ST2, brain specific, cancer and invasion related
    # CERCAM = cerebral endothelial cell adhesion molecule << while endothelial markers go down?
    
    #"CD3D", "CD3E", "CD3G", "CD163"
    #"CD274", "PDCD1LG2", "HLA-C", "B2M" # PDL1, PDL2, MHC, B2M
    
    #"PLP1", "CLDN11", "CNTN2", "CTNNA3", ""
    #"CCL4", "CCL33", "CX3CD1", "CD53", "RUNX1","IL1B","IL6R","PTPRC","BLNK"
    "CCL4L1"
  ))


ggplot(plt, aes(x=log2FoldChange.gsam.res , y=statistic.gsam.cor.tpc , col=show.label, label = hugo_symbol )) + 
  geom_point(data=subset(plt, show.label == F),cex=0.35) +
  geom_point(data=subset(plt, show.label == T ),col="black",cex=0.65) +
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35) + 
  geom_text_repel(data=subset(plt, show.label == T & log2FoldChange.gsam.tpc.res < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35) + 
  scale_color_manual(values = c('TRUE'=rgb(0,0,0,0.35),'FALSE'='gray60')) +
  labs(x = "log2FC R1 vs. R2 (Tumor cell percentage corrected)",
       y="Correlation t-statistic with tumour percentage"
       ,col="Difference significant (R1 ~ R2)") +
  youri_gg_theme + xlim(-3, 3)



ggplot(plt, aes(x=log2FoldChange.glass.res ,
                      y=statistic.gsam.cor.tpc  ,
                      col=show.label,
                      label = hugo_symbol )) + 
  geom_point(data=subset(plt,  show.label == F),cex=0.35) +
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


plt <- data.frame(MGMT = as.numeric(gsam.gene.expression.all.vst %>%
                                      as.data.frame() %>%
                                      dplyr::filter(grepl("TP53\\|", rownames(.)))) ,
                  res = gsub("^...(.).*$","R\\1",colnames( gsam.gene.expression.all.vst)))

ggplot(plt, aes(x = res, y = MGMT)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y = "MGMT relative expression")

rm(plt)



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

