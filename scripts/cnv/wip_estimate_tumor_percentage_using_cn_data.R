#!/usr/bin/env R

"
Goal of this script is to determine the tumour percentage

This is technically possible with CNVKit + bam files + loads of other analysis,
but may most likely take weeks


https://cnvkit.readthedocs.io/en/stable/heterogeneity.html

Requires VCF Files, which require running mutect, which require the bam files, of which ~25% are missing
"


setwd("~/projects/gsam")

# load libs & funcs ----


library(tidyverse)
library(ggplot2)
library(limma)


# load data----


chrs_hg19 <- {{}}
chrs_hg19['chr1'] <- 249250621
chrs_hg19['chr2'] <- 243199373
chrs_hg19['chr3'] <- 198022430
chrs_hg19['chr4'] <- 191154276
chrs_hg19['chr5'] <- 180915260
chrs_hg19['chr6'] <- 171115067
chrs_hg19['chr7'] <- 159138663
chrs_hg19['chr8'] <- 146364022
chrs_hg19['chr9'] <- 141213431
chrs_hg19['chr10'] <- 135534747
chrs_hg19['chr11'] <- 135006516
chrs_hg19['chr12'] <- 133851895
chrs_hg19['chr13'] <- 115169878
chrs_hg19['chr14'] <- 107349540
chrs_hg19['chr15'] <- 102531392
chrs_hg19['chr16'] <- 90354753
chrs_hg19['chr17'] <- 81195210
chrs_hg19['chr18'] <- 78077248
chrs_hg19['chr19'] <- 59128983
chrs_hg19['chr20'] <- 63025520
chrs_hg19['chr21'] <- 48129895
chrs_hg19['chr22'] <- 51304566
chrs_hg19['chrX'] <- 155270560
chrs_hg19['chrY'] <- 59373566
n_hg19 = sum(chrs_hg19)


chrs_hg19_c <- {{}}
chrs_hg19_c['chr1'] <- 0
chrs_hg19_c['chr2'] <- chrs_hg19['chr1'] 
chrs_hg19_c['chr3'] <- chrs_hg19['chr1'] + chrs_hg19['chr2']
chrs_hg19_c['chr4'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3']
chrs_hg19_c['chr5'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] 
chrs_hg19_c['chr6'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] 
chrs_hg19_c['chr7'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] 
chrs_hg19_c['chr8'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] 
chrs_hg19_c['chr9'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] 
chrs_hg19_c['chr10'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9']
chrs_hg19_c['chr11'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10']
chrs_hg19_c['chr12'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11']
chrs_hg19_c['chr13'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12']
chrs_hg19_c['chr14'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13']
chrs_hg19_c['chr15'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14']
chrs_hg19_c['chr16'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15']
chrs_hg19_c['chr17'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16']
chrs_hg19_c['chr18'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17']
chrs_hg19_c['chr19'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17'] + chrs_hg19['chr18']
chrs_hg19_c['chr20'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17'] + chrs_hg19['chr18'] + chrs_hg19['chr19']
chrs_hg19_c['chr21'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17'] + chrs_hg19['chr18'] + chrs_hg19['chr19'] + chrs_hg19['chr20']
chrs_hg19_c['chr22'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17'] + chrs_hg19['chr18'] + chrs_hg19['chr19'] + chrs_hg19['chr20'] + chrs_hg19['chr21']
chrs_hg19_c['chrX'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17'] + chrs_hg19['chr18'] + chrs_hg19['chr19'] + chrs_hg19['chr20'] + chrs_hg19['chr21'] + chrs_hg19['chr22']
chrs_hg19_c['chrY'] <- chrs_hg19['chr1'] + chrs_hg19['chr2'] + chrs_hg19['chr3'] + chrs_hg19['chr4'] + chrs_hg19['chr5'] + chrs_hg19['chr6'] + chrs_hg19['chr7'] + chrs_hg19['chr8'] + chrs_hg19['chr9'] + chrs_hg19['chr10'] + chrs_hg19['chr11'] + chrs_hg19['chr12'] + chrs_hg19['chr13'] + chrs_hg19['chr14'] + chrs_hg19['chr15'] + chrs_hg19['chr16'] + chrs_hg19['chr17'] + chrs_hg19['chr18'] + chrs_hg19['chr19'] + chrs_hg19['chr20'] + chrs_hg19['chr21'] + chrs_hg19['chr22'] + chrs_hg19['chrX']




# load data ----

#source("scripts/R/gsam_metadata.R") - auto loaded in script below
source("scripts/R/cnv_matrix.R")
source("scripts/R/youri_gg_theme.R")
source("scripts/R/job_gg_theme.R")


cnv_segments <- read.delim('data/gsam/output/tables/cnv_copynumber-ratio.cns_all.txt', stringsAsFactors = F) %>%
  dplyr::mutate(gene = NULL) %>%
  dplyr::mutate(segment.length = end - start)



# some test code ----
# 5,6 = pair without further gains and losses


i = 1
j = 2

off <- median(cnv_matrix [, i] - cnv_matrix[,j])
s1 <- cnv_matrix [, i ]
s2 <- (cnv_matrix[,j] + off)

plot(c(-2.2, 2.2),c(-2.2, 2.2), type="n",xlab=colnames(cnv_matrix)[i],ylab=colnames(cnv_matrix)[j])
points(s1, s2, pch=19,cex=0.2,col=rgb(0,0,0,0.05))



# per sample example ----

## IMPORTANT - FROM CNVKit:
## The log2 ratio values of CNAs in a tumor sample correspond to
## integer copy numbers in tumor cells, and in aggregate these
## log2 values will cluster around values that indicate subclone
## populations, each with a given ploidy and clonality. For
## example, a single-copy loss in a 50% pure tumor sample will
## have 3/4 the coverage of a neutral site (2/2 normal copies,
## 1/2 tumor copies), for a log2 value of log2(.75) = -0.415.
## This calculation can also be generalized to other copy number
## states.


data.frame(n = 2, lfc = c(-0.415, 0.415,    -0.559452, 0.835, log2(1.25))) %>%
  dplyr::mutate(fc = 2^(lfc) ) %>%
  dplyr::mutate(pct = case_when(
    fc < 1 ~     (1 - fc) * n * 100.0 ,
    fc > 1 & fc < 1.5 ~ (fc * 2) - 2 * 100.0
    ))



pct = c(1.0, 0.75,0.5,0.25,0.1,0)
print(pct)
fc = ((1-pct)*2 + pct*3) / n
print(fc)
pct.2 = (fc * 2) - 2
pct.2 = ((fc * 2) - 2)/2
print(pct.2)



#pct = 0.835
#(pct*2 + pct*3) / (2 + 2) = 
#(pct*2 + pct*3) = fc  / (2 + 2)

# naar 0.835

1.04375 = (pct*2 + pct*3) / (2 + 2)
1.04375 * (2 + 2) = pct*2 + pct*3
1.04375 * 4 / 5
0.835 * 4 / 5




# AAF1 kan beter ~= 0.74%
# AAS2

##  estiamte them


tpc.estimate = data.frame()
for(k in 1:ncol(cnv_matrix)) {
  print(k)
  #k = 7
  #k = 115 # AAF
  
  
  #k = which(colnames(cnv_matrix) == "KAE1")
  
  if(colnames(cnv_matrix)[k]  %in% cnv_segments$patient.id) {
    
    plt <- cnv_matrix %>%
      tibble::rownames_to_column('segment.loc') %>%
      dplyr::mutate(cur.sample = colnames(cnv_matrix)[k] ) %>% # define sid
      dplyr::filter(row_number() %% 4 == 1) %>%
      dplyr::mutate(segment.chr = gsub(':.+$', '', segment.loc)) %>%
      dplyr::mutate(start = as.numeric(gsub("^.+:(.+)\\-.+$","\\1", segment.loc))) %>%
      dplyr::mutate(end = as.numeric(gsub("^.+:.+\\-(.+)$","\\1", segment.loc))) %>%
      dplyr::mutate(pos = chrs_hg19_c[segment.chr] + ((start + end) / 2)) %>%
      dplyr::mutate(pos.end = NA) %>%
      dplyr::rename(log2 = unique(.$cur.sample)) %>%
      dplyr::mutate(type = 'point') %>%
      dplyr::filter(segment.chr %in% c('chrX', 'chrY') == F)


    tmp <- cnv_segments %>%
      dplyr::mutate(cur.sample = unique(plt$cur.sample)) %>%
      dplyr::filter(patient.id == cur.sample) %>%
      dplyr::mutate(chromosome = paste0("chr", chromosome)) %>%
      dplyr::rename(segment.chr = chromosome) %>%
      dplyr::mutate(pos = start + chrs_hg19_c[segment.chr] ) %>%
      dplyr::mutate(pos.end = end + chrs_hg19_c[segment.chr] ) %>%
      dplyr::mutate(type = 'segment') %>%
      dplyr::filter(segment.chr %in% c('chrX', 'chrY') == F) %>%
      dplyr::filter(segment.length >= 5000000) %>%
      dplyr::filter(log2 < 1.1)


    plt <- rbind(
      plt %>% dplyr::select(c('segment.chr', 'pos', 'log2', 'type', 'pos.end')),
      tmp %>% dplyr::select(c('segment.chr', 'pos', 'log2', 'type', 'pos.end')))


    out <- data.frame()
    for(frac in  1:100 / 100) {
      fc.p.4 <- ((1 - frac) * 2 + frac * 4) / 2
      fc.p.3 <- ((1 - frac) * 2 + frac * 3) / 2
      fc.n.1 <-((1 - frac) * 2 + frac * 1) / 2
      
      lfc.p.4 <- log2(fc.p.4)
      lfc.p.3 <- log2(fc.p.3)
      lfc.n.1 <- log2(fc.n.1)
      
      
      dists <- c()
      dist <- 0
      for(i in 1:nrow(tmp)) {
        # these sample are typically clonally different at particular chr's
        if(colnames(cnv_matrix)[k] == "AAF1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr10', 'chr5', 'chr8','chr18','chr1') == F) }
        else if(colnames(cnv_matrix)[k] == "EAD2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1', 'chr2', 'chr14','chr18','chr22') == F) }
        else if(colnames(cnv_matrix)[k] == "EAZ1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2', 'chr3', 'chr4','chr5','chr6','chr7','chr8','chr10','chr14') == T) }
        else if(colnames(cnv_matrix)[k] == "EBW1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr2', 'chr3', 'chr4','chr7','chr10','chr11','chr13','chr14') == T) }
        else if(colnames(cnv_matrix)[k] == "AAS2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr10', 'chr13','chr14') == T) }
        else if(colnames(cnv_matrix)[k] == "AAT1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr4', 'chr15','chr16') == F) }
        else if(colnames(cnv_matrix)[k] %in% c("AAU1","AAU2")) { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1', 'chr2','chr3','chr9','chr10','chr11','chr12','chr18','chr19') == T) }
        else if(colnames(cnv_matrix)[k] == "ABA1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr3', 'chr4','chr5', 'chr7','chr10','chr16') == T) }
        else if(colnames(cnv_matrix)[k] == "ACA1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr3') == F) }
        else if(colnames(cnv_matrix)[k] == "AFA1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr2','chr3') == T) }
        else if(colnames(cnv_matrix)[k] == "AZB2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr4','chr5','chr6','chr7','chr10','chr13') == T) }
        else if(colnames(cnv_matrix)[k] == "AZH1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr6','chr8','chr16','chr17','chr18') == F) }
        else if(colnames(cnv_matrix)[k] == "BAC1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr13','chr15','chr19') == F) }
        else if(colnames(cnv_matrix)[k] == "BAW1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr2','chr4','chr10','chr13','chr14','chr15') == T) }
        else if(colnames(cnv_matrix)[k] == "CBT1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "CDA1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr3','chr4','chr9','chr17','chr18','chr19') == F) }
        else if(colnames(cnv_matrix)[k] == "DAB1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr7','chr10','chr13') == T) }
        else if(colnames(cnv_matrix)[k] == "EAO2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr7','chr10','chr13') == T) }
        else if(colnames(cnv_matrix)[k] %in% c("ECA1","ECA2")) { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr2','chr7','chr10','chr22') == T) }
        else if(colnames(cnv_matrix)[k] == "FAF1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr3','chr7','chr10','chr13','chr22') == T) }
        else if(colnames(cnv_matrix)[k] == "FAF2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr3','chr7','chr10','chr13','chr22') == T) }
        else if(colnames(cnv_matrix)[k] == "FAJ2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr7','chr10','chr13','chr15') == T) }
        else if(colnames(cnv_matrix)[k] == "FAP1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr14','chr15') == T) }
        else if(colnames(cnv_matrix)[k] == "FAP2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2',"chr7",'chr10','chr14','chr15') == T) }
        else if(colnames(cnv_matrix)[k] == "GAI1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2',"chr7",'chr9','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "GAI2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2',"chr7",'chr9','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "HAB2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2',"chr3",'chr7','chr10','chr15') == T) }
        else if(colnames(cnv_matrix)[k] == "HAF1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2',"chr3",'chr7','chr10','chr15') == T) }
        else if(colnames(cnv_matrix)[k] == "HAK1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr5',"chr10",'chr13','chr17') == T) }
        else if(colnames(cnv_matrix)[k] == "JAB1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr3',"chr4",'chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "JAE2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr7',"chr12",'chr13') == T) }
        else if(colnames(cnv_matrix)[k] == "JAI1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr7','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "JAL2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr7','chr8','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "JAN1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr7','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "JAN2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr7','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "KAB2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr7','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "KAC2") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr2','chr13','chr22') == T) }
        else if(colnames(cnv_matrix)[k] == "KAD1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr7','chr10') == T) }
        else if(colnames(cnv_matrix)[k] == "KAE1") { tmp <- tmp %>% dplyr::filter(segment.chr %in% c('chr1','chr7','chr10') == T) }
        
        
        
        e <- tmp[i,]
        #if(e$log2 < lfc.p.4 * 1.2) { # extreme outlier; skip
          d <- c(
            (((e$log2 - 0) * e$weight)^2),
            (((e$log2 - lfc.p.4) * e$weight)^2) * 1.1, # penalize, do never prefer
            (((e$log2 - lfc.p.3) * e$weight)^2),
            (((e$log2 - lfc.n.1) * e$weight)^2)
          )
          d <- min(d)
          dists <- c(dists, d)
          
          dist <- dist + d
        #}
      }
      
      out <- rbind(out, data.frame(pct = frac * 100,
                                   lfc.3p = lfc.p.3,
                                   lfc.4p = lfc.p.4,
                                   lfc.n = lfc.n.1,
                                   dist = dist)
                   )
    }
    
    
    
    lfc <- out %>%
      dplyr::arrange(dist) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(sample = unique(tmp$cur.sample))
    
    #print(lfc)
    tpc.estimate <- rbind(tpc.estimate, lfc)
    #print(tpc.estimate)
    #print("")
    
    
    ggplot(plt, aes(x = pos, y = log2, col=segment.chr) ) + 
      geom_hline(yintercept = lfc$lfc.3p , col="red") +
      geom_hline(yintercept = lfc$lfc.4p , col="red", lwd=0.5, alpha=0.5) +
      geom_hline(yintercept = lfc$lfc.n , col="red") +
      geom_hline(yintercept = 0, col="gray", lwd=0.5, lty=2) +
      geom_point(pch=19, size=0.2, alpha=0.5, data = subset(plt, type == "point") ) +
      theme_void() + 
      job_gg_theme +
      theme(legend.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank() ) +
      ylim(-2, 2) +
      geom_segment(aes( xend = pos.end, yend = log2, colour = "segment"), data = subset(plt, type == "segment") , col="white", lwd=2 ) + 
      geom_segment(aes( xend = pos.end, yend = log2, colour = "segment"), data = subset(plt, type == "segment") , col="black" , lwd=1.35) +
      labs(subtitle = paste0(unique(tmp$cur.sample), ": log2Fc: ", round(lfc$lfc.3p,3), " ~ ", round(lfc$lfc.n,3), "  est. tumor: ", lfc$pct, "%    dist:", round(lfc$dist,1)))
    
    
    
    fn <- paste0("output/figures/cnv/",unique(tmp$cur.sample),".png")
    #print(fn)
    ggsave(fn,width=10,height=6)
    
  }
}

write.table(tpc.estimate , "output/tables/cnv/tumor.percentage.estimate.txt")




#data.frame(lfc = c(2, 1,0.75,0.5,0.25,0.07,0)) %>%
#  dplyr::mutate(cnv = 2^lfc) %>%
#  dplyr::mutate(pct = (cnv - 1) * 100)


data.frame(chr.n = 2, lfc = c(2, 1, 0.75, 0.5849, 0.5,0.25, 0.1, 0.07,0)) %>%
  dplyr::mutate(fc = 2^lfc) %>%
  dplyr::mutate(cnv = chr.n * fc) %>% # if the foldChange is 3, we expect 6 chromosomes (n=2)
  dplyr::mutate(pct = (cnv - chr.n) * 100.0)

#chr.n = 2
#lfc = 0.5849

pct = ((chr.n * (2^lfc)) - chr.n) * 100

n = 2
pct = c(100, 10, 5, 0)
pct = 60
lfc = log2(((pct / 100) + n) / n)
# lfc[0.5849625] -> 100%
# lfc[0.07038933] -> 10%
# lfc[0.03562391] -> 5%




values <- c()

plot(c(-1,2), c(0, max(abs(tmp$segment.length))), type="n")
for(i in 1:nrow(tmp)) {
  e <- tmp[i,]
  
  lines(c(abs(e$log2), abs(e$log2)), c(0, e$segment.length))
  
  
  values <- c(values,
    plt %>%
      dplyr::filter(type == 'point' &
                      segment.chr == e$segment.chr
                      & pos >= e$pos
                      & pos <= e$pos.end
      ) %>% 
    dplyr::pull(log2))

}




# there are insufficient snp calls for snp based

# (re-)estimate using RNA-seq

tpc.estimate <- read.table("output/tables/cnv/tumor.percentage.estimate.txt")

source('scripts/R/gsam_rna-seq_expression.R')
source('data/wang/msig.library.12.R') # no license w/ code provided, can't include it in source tree


expression <- gsam.rnaseq.expression.vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% (gsam.rnaseq.expression %>%
                            dplyr::filter(rowSums(.) >= ncol(.) * 10) %>%
                            tibble::rownames_to_column('gid') %>%
                            dplyr::pull(gid))) %>%
  tibble::column_to_rownames('gid') %>%
  as.matrix()



expression.pca <- expression %>%
  prcomp() %>%
  purrr::pluck('rotation') %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(pid = gsub("^(....).*$","\\1",sid) ) %>%
  dplyr::left_join(tpc.estimate , by=c('pid'='sample') ) %>%
  dplyr::mutate(dr_ratio = pct / PC3) %>%
  dplyr::mutate(outlier = abs(dr_ratio) > 2000) %>%
  dplyr::mutate(resection = gsub("^...(.).*$","R\\1",pid))



ggplot(expression.pca, aes(x = pct, y= PC3, col = resection , label = sid , group=1 ) )  + 
  geom_smooth(method='lm', col="gray60",se = FALSE, lty=2,lwd=0.65) +
  geom_point(cex = 0.8) +
  xlim(0,100) +
  youri_gg_theme +
  labs(x = "WES estimated tumour percentage", y = "RNA-seq PC-3")

ggsave("output/figures/cnv/tumour_cell_percentage_RNA_PC-3.png")



fit <- lm(PC3 ~ pct, data = expression.pca)


expression.pca <- expression.pca %>%
  dplyr::mutate(PC3.scaled = (PC3 - fit$coefficients[1]) /  fit$coefficients[2] ) %>%
  dplyr::mutate(outlier = abs(PC3.scaled - pct) > 28) %>%
  tibble::rownames_to_column('sid')


lm(PC3.scaled ~ pct, data = expression.pca)



ggplot(expression.pca, aes(x = pct, y= PC3.scaled, col = outlier , label = sid  ) ) + 
  geom_point(cex = 0.6) +
  geom_text_repel(data = subset(expression.pca, outlier == T) )
  


# does not seem to reveal clear improvements?


