#!/usr/bin/env R

# read GSE103224 ----

if(file.exists("tmp/GSE103224.scRNA.counts.Rds")) {
  GSE103224 <- readRDS("tmp/GSE103224.scRNA.counts.Rds")
} else {
  a <- read.delim("data/scRNA/GSE103224/GSM2758471_PJ016.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758471.cell.",(1:(ncol(.)-2))+2) ))
  b <- read.delim("data/scRNA/GSE103224/GSM2758472_PJ017.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758472.cell.",(1:(ncol(.)-2))+2) ))
  c <- read.delim("data/scRNA/GSE103224/GSM2758473_PJ018.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758473.cell.",(1:(ncol(.)-2))+2) ))
  d <- read.delim("data/scRNA/GSE103224/GSM2758474_PJ025.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758474.cell.",(1:(ncol(.)-2))+2) ))
  e <- read.delim("data/scRNA/GSE103224/GSM2758475_PJ030.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758475.cell.",(1:(ncol(.)-2))+2) ))
  f <- read.delim("data/scRNA/GSE103224/GSM2758476_PJ032.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758476.cell.",(1:(ncol(.)-2))+2) ))
  g <- read.delim("data/scRNA/GSE103224/GSM2758477_PJ035.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758477.cell.",(1:(ncol(.)-2))+2) ))
  
  # brain meta?
  #h <- read.delim("data/scRNA/GSE103224/GSM2940098_PJ048.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  #  `colnames<-`(c("ENS", "HGNC", paste0("GSM2940098.cell.",(1:(ncol(a)-2))+2) ))
  
  
  
  GSE103224 <- a %>%
  dplyr::left_join(b %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(c %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(d %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(e %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(f %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(g %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS'))
  
  #saveRDS(GSE103224, file="tmp/GSE103224.scRNA.counts.Rds")
  
  rm(a,b,c,d,e,f,g)
}



# read GSE117891, GSE131928 ----

# read GSE117891, GSE131928 ----

# further processing ----