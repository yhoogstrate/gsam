#!/usr/bin/env R 

# load  ----
## libs ----


library(tidyverse)


## gencode 22 ----

if(!exists('gencode.22')) {
  rm.gencode <- T
  
  if(file.exists('tmp/gencode.22.Rds')) {
    gencode.22 <- readRDS('tmp/gencode.22.Rds')
  }
  else {
    gencode.22 <- read.table("data/tcga/read-counts/gencode.v22.annotation.gtf",sep="\t",comment.char="#",stringsAsFactors = F) %>%
      dplyr::filter(V3 == "gene") %>%
      dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
      dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
      dplyr::mutate(V9 = NULL) %>%
      dplyr::mutate(ENSG.short = gsub('\\..*?$','',ENSG) )
    
    saveRDS(gencode.22, file='tmp/gencode.22.Rds')
  }
} else {
  rm.gencode <- F
}


## Wang paper genes ----

wang.glioma.intrinsic.genes <- read.table('data/wang/table_s1.txt',
                                          stringsAsFactors = F,
                                          sep="\t",
                                          header=T)


isct <- readRDS("tmp/isct.Rds")
if(!exists('isct')) {
  
  tcga.gtf <- read.table('data/wang/Homo_sapiens.GRCh37.64.gtf', sep="\t", header=F, stringsAsFactors = F) %>% # head() %>%
    dplyr::mutate(ENSG.short = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
    dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
    dplyr::select(c('ENSG.short', 'GENE')) %>% # This is non-unique matching - multiple ENSEMBL entries per GENE
    dplyr::distinct(ENSG.short, GENE) 
  
  isct <- data.frame('ENSG' = c(), 'GENE'= c(), 'ENSG.short'= c())
  for(gene in wang.glioma.intrinsic.genes$Gene_Symbol) {
    a <- gencode.22 %>% 
      dplyr::filter(GENE == gene) %>%
      dplyr::select(c('ENSG', 'GENE', 'ENSG.short'))
    
    b <- tcga.gtf %>%
      dplyr::filter(GENE == gene)
    
    if(nrow(a) == 1) { # direct match in gencode
      m <- a
    
    } else if(nrow(a) > 1) { # multiple matches in gencode
      c <- a %>%
        dplyr::filter(ENSG.short %in% b$ENSG.short)
      
      if(nrow(c) == 0) {
        # pick the one from a with the lexigraphically smallest ENSG
        # seems to happen only with: gene <- 'C17orf100'
        
        m <- a %>%
          dplyr::arrange(ENSG.short) %>%
          dplyr::top_n(1)
      } else if(nrow(c) == 1) {
        # pick the intersect
        
        m <- c
      } else {
        # pick the one from c with the lexigraphically smallest ENSG
        #print("C")
        print("shoud not happen")
        stopifnot(FALSE)
      }
    
    } else if(nrow(a) == 0) { # direct gene name does not match gencode.22 - let's see if we can translate it through tcga.gtf
  
      if(nrow(b) == 1) {
        m <- gencode.22 %>%
          dplyr::filter(ENSG.short %in% b$ENSG.short) %>%
          dplyr::select(c('ENSG', 'GENE', 'ENSG.short'))
      } else if(nrow(b) >= 2) {
        m <- gencode.22 %>%
          dplyr::filter(ENSG.short %in% b$ENSG.short) %>%
          dplyr::select(c('ENSG', 'GENE', 'ENSG.short')) %>%
          dplyr::arrange(ENSG.short) %>%
          dplyr::top_n(1)
      } else {
        #print("contr")
        #print(gene)
        m <- data.frame(ENSG = NA , GENE=gene, ENSG.short = NA)
      }
    
    }
    
    if(! is.null(m)) {
      if(nrow(m) == 0) {
        m <- data.frame(ENSG = NA , GENE=gene, ENSG.short = NA)
      } else if(m$GENE != gene) {
        m <- data.frame(ENSG = NA , GENE=gene, ENSG.short = NA)
      }
      isct <- rbind(isct, m)
    }
    
    rm(m)
  }
  
  stopifnot( nrow(isct) == nrow(wang.glioma.intrinsic.genes))
  stopifnot(sum(duplicated(isct %>% dplyr::filter(!is.na(ENSG.short)) %>% dplyr::pull(ENSG.short))) == 0)

  saveRDS(isct, "tmp/isct.Rds")
  
  rm(a,b,c, tcga.gtf)
}




wang.glioma.intrinsic.genes <- wang.glioma.intrinsic.genes %>% 
  dplyr::left_join(isct, by=c('Gene_Symbol' = 'GENE'))


if(rm.gencode) {
  rm(gencode.22)
}


rm(isct, rm.gencode)



