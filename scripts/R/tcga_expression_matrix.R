#!/usr/bin/env R 

# load  ----
## libs ----

library(tidyverse)
library(NMF)

## cfg ----

source('scripts/R/youri_gg_theme.R')

## identifiers ----


tcga.identifiers.gbm <- read.delim('data/tcga/tables/rna-seq/gdc_sample_sheet.2019-01-28.tsv',stringsAsFactors = F) %>%
  dplyr::mutate( tcga.htseq.id = paste0('tcga.',gsub("-",".",gsub(".htseq.counts.gz","",File.Name,fixed=T))) ) %>%
  dplyr::mutate( Data.Category = NULL ) %>%
  dplyr::mutate(Data.Type = NULL) %>%
  dplyr::mutate(TCGA.pat.ID = gsub("^(....)-[^-]+-[0]*([^-]+).+$","\\1.\\2",Sample.ID)) %>%
  dplyr::filter(tcga.htseq.id != "tcga.18fbf794.a85b.4b0d.9a5c.c1c0e0ede3d8") %>% # replicate of TCGA-06-0156
  dplyr::filter(tcga.htseq.id != "tcga.f67cf68d.56b3.451d.ba21.566664213691") %>% # replicate of TCGA-06-0211
  dplyr::filter((Case.ID %in% c("TCGA-06-0211", "TCGA-06-0210", "TCGA-14-1034", "TCGA-06-0125", "TCGA-06-0190", "TCGA-19-4065" ) &
                   Sample.Type == "Recurrent Tumor") == F ) # From these, exclude the secondary - as both prim and secondary are present: "TCGA-06-0211" "TCGA-06-0210" "TCGA-14-1034" "TCGA-06-0125" "TCGA-06-0190" "TCGA-19-4065"
  #dplyr::filter(Sample.Type == "Primary Tumor") %>% Recurrent tumours were also used?


## Wang / GlioVis paper data ----


TCGA_GBM <- readRDS('data/gliovis/data/datasets/TCGA_GBM.Rds')
TCGA_GBM$sample_rnaseq <- TCGA_GBM$rseq %>%
  dplyr::mutate(Sample = gsub('\\.', '-', Sample, fixed=F )) %>%
  dplyr::select(c('EGFR', 'Sample', 'Histology', 'Recurrence')) %>%
  dplyr::mutate(Sample.RNA = case_when( # kan niet met NA values dealen
    is.na(EGFR) ~ "Not RNA-sequenced" ,
    Histology != "GBM" ~ "No GBM",
    Sample %in% tcga.identifiers.gbm$Case.ID == F ~ "no identifier",
    TRUE ~ "Match")) 
TCGA_GBM$sample_rnaseq %>% dplyr::filter(Sample.RNA %in% c('Not RNA-sequenced', 'Match') == F )



TCGA_GBM <- readRDS('data/gliovis/data/datasets/TCGA_GBM.Rds')
TCGA_GBM$sample_rnaseq <- TCGA_GBM$rseq %>%
  dplyr::mutate(Sample = gsub('\\.', '-', Sample, fixed=F )) %>%
  dplyr::select(c('EGFR', 'Sample', 'Histology', 'Recurrence')) %>%
  dplyr::mutate(Sample.RNA = case_when( # kan niet met NA values dealen
    is.na(EGFR) ~ "NA" ,
    Histology != "GBM" ~ "NA",
    Sample %in% tcga.identifiers.gbm$Case.ID == F ~ "NA",
    TRUE ~ Sample)) %>%
  dplyr::mutate(Sample.RNA = ifelse(Sample.RNA == "NA" , NA , Sample.RNA)) %>%
  dplyr::select(Sample, Sample.RNA)

# Not included:
# 3   TCGA-06-0675 Non-tumor        No GBM
# 4   TCGA-06-0678 Non-tumor        No GBM
# 5   TCGA-06-0680 Non-tumor        No GBM
# 6   TCGA-06-0681 Non-tumor        No GBM
# 7   TCGA-06-5415       GBM no identifier  not in TCGA?


# Included - but are actually recurrent samples (!) :
# 1   TCGA-06-0152       GBM no identifier  recurrent?
# 2   TCGA-06-0171       GBM no identifier  recurrent?
# 8   TCGA-14-0736       GBM no identifier  recurrent?
# 9   TCGA-14-1402       GBM no identifier  recurrent?
# 10  TCGA-19-0957       GBM no identifier  recurrent?
# 11  TCGA-19-1389       GBM no identifier  recurrent?


stopifnot( (TCGA_GBM$sample_rnaseq %>%
  dplyr::filter(!is.na(Sample.RNA)) %>%
  dplyr::pull(Sample.RNA)) %in% tcga.identifiers.gbm$Case.ID)


## actual TCGA counts and match w/ Wang ----


tcga.expression.gbm <-read.delim('data/tcga/tables/rna-seq/merged-expression_TCGA-GBM.htseq.counts.txt') %>%
  `colnames<-`(gsub("^X","tcga.",colnames(.))) %>%
  `colnames<-`(gsub("^([a-f])","tcga.\\1",colnames(.))) %>%
  `colnames<-`(gsub(".htseq.counts.gz","",colnames(.),fixed=T)) %>%
  dplyr::filter(gene.id %in% c('__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned', '__too_low_aQual') == F) %>%
  tibble::column_to_rownames('gene.id') %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('sample.id') %>%
  dplyr::left_join(tcga.identifiers.gbm %>% dplyr::select(c('tcga.htseq.id', "Case.ID" )) , by=c('sample.id'='tcga.htseq.id')) %>%
  dplyr::filter(!is.na(Case.ID)) %>%
#a = tcga.expression.gbm  %>%
  #dplyr::select('sample.id' , 'Case.ID', 'ENSG00000281884.1','ENSG00000281887.1') %>%
  dplyr::mutate(sample.id = NULL) %>%
  tibble::column_to_rownames('Case.ID') %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(gsub(".","-",colnames(.),fixed=T)) %>%
  dplyr::select( TCGA_GBM$sample_rnaseq$Sample.RNA %>% purrr::discard(is.na)  ) # only take those also used in gliovis
  

stopifnot(TCGA_GBM$sample_rnaseq$Sample.RNA %>% purrr::discard(is.na) == colnames(tcga.expression.gbm))



## gencode 22 ----


gencode.22 <- read.table("data/tcga/read-counts/gencode.v22.annotation.gtf",sep="\t",comment.char="#",stringsAsFactors = F) %>%
  dplyr::filter(V3 == "gene") %>%
  dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(V9 = NULL) %>%
  dplyr::mutate(ENSG.short = gsub('\\..*?$','',ENSG) )


stopifnot(rownames(tcga.expression.gbm) %in% gencode.22$ENSG)


#sum(duplicated(gencode.22$GENE))
#sum(duplicated(gencode.22$ENSG))



## Wang paper genes ----

wang.glioma.intrinsic.genes <- read.table('data/wang/table_s1.txt', stringsAsFactors = F, sep="\t", header=T)

tcga.gtf <- read.table('data/wang/Homo_sapiens.GRCh37.64.gtf', sep="\t", header=F, stringsAsFactors = F) %>% # head() %>%
  dplyr::mutate(ENSG.short = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
  dplyr::select(c('ENSG.short', 'GENE')) %>% # This is non-unique matching - multiple ENSEMBL entries per GENE
  dplyr::distinct(ENSG.short, GENE) 


isct <- readRDS("tmp/isct.Rds")
if(!exists('isct')) {
   
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
  
  rm(a,b,c)
}


wang.glioma.intrinsic.genes <- wang.glioma.intrinsic.genes %>% 
  dplyr::left_join(isct, by=c('Gene_Symbol' = 'GENE'))


#a = wang.glioma.intrinsic.genes %>% dplyr::filter(!is.na(ENSG)) %>% dplyr::pull('ENSG')
#b = gencode.22 %>% dplyr::filter(!is.na(ENSG)) %>% dplyr::pull('ENSG')



# NMF all genes ----

samples <- TCGA_GBM$sample_rnaseq %>%
  dplyr::filter(!is.na(Sample.RNA)) %>%
  dplyr::pull(Sample.RNA)


# expression of all genes with 3 reads on average per sample OR from the Wang list
expression.matrix.tcga <- tcga.expression.gbm %>%
  dplyr::select(samples) %>%
  dplyr::filter(
    rowSums(.) > 0 &
    (
      rownames(.) %in% (wang.glioma.intrinsic.genes %>% dplyr::filter(!is.na(ENSG)) %>% dplyr::pull('ENSG')) | 
      rowSums(.) >= 3 * ncol(.)
    )
  )



metadata.tcga <- data.frame(sample = colnames(expression.matrix.tcga)) %>%
  dplyr::left_join(
    TCGA_GBM$rseq %>%
      dplyr::select(c('Sample', 'Histology', 'Grade', 'Recurrence', 'Subtype', 'CIMP_status', 'survival', 'status')) %>%
      dplyr::mutate(Sample = gsub('\\.','-', Sample) )
    ,
    by = c('sample' = 'Sample')
    )

# these become NA because of zero counts, while being in the list of Wang
# "ENSG00000100146.15" "ENSG00000100219.15" "ENSG00000105371.8"  "ENSG00000139574.8"  "ENSG00000164708.5"  "ENSG00000213689.8"  "ENSG00000256453.1"  "ENSG00000275221.1" 
expression.matrix.tcga.vst <- expression.matrix.tcga %>%
  DESeq2::DESeqDataSetFromMatrix(S4Vectors::DataFrame(tmp=as.factor(1:ncol(expression.matrix.tcga))), ~tmp) %>%
  DESeq2::varianceStabilizingTransformation(blind=T) %>%
  SummarizedExperiment::assay() %>%
  data.frame(stringsAsFactors = F) %>%
  dplyr::mutate(mad = rapply(data.frame(t(.)) , stats::mad) )

#    0.5786233          0.6078969          0.3170584


#
wang.glioma.intrinsic.genes %>%
  dplyr::filter(U133A_Gene. == "YES") %>%
  dplyr::filter(is.na(ENSG)) %>%
  dim()



expression.matrix.tcga.vst.BFG.U133A <- expression.matrix.tcga.vst %>%
  dplyr::filter(rownames(.) %in% 
                  (wang.glioma.intrinsic.genes %>%
                     dplyr::filter(U133A_Gene. == "YES") %>%
                     dplyr::filter(!is.na(ENSG)) %>%
                     dplyr::pull(ENSG)))



expression.matrix.tcga.vst.150 <- expression.matrix.tcga.vst %>%
  dplyr::filter(rownames(.) %in% 
                  (wang.glioma.intrinsic.genes %>%
                     dplyr::filter(Subtyping_Signature_Gene. != "") %>%
                     dplyr::filter(!is.na(ENSG)) %>%
                     dplyr::pull(ENSG)))



# ---- PCA: TCGA [BFG/U133A] -----

plt <- prcomp(expression.matrix.tcga.vst.BFG.U133A %>% t(), scale=F)$x %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(sid = gsub("\\.","-", sid) ) %>%
  dplyr::left_join(metadata.tcga , by = c('sid' = 'sample'))


ggplot2::ggplot(plt, aes(x = PC1, y=PC3, col = Subtype)) + 
  geom_point() +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype")

# ---- PCA: TCGA [150] -----

plt <- prcomp(expression.matrix.tcga.vst.150 %>% t(), scale=F)$x %>%
  data.frame(stringsAsFactors = F) %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::mutate(sid = gsub("\\.","-", sid) ) %>%
  dplyr::left_join(metadata.tcga , by = c('sid' = 'sample'))


ggplot2::ggplot(plt, aes(x = PC1, y=PC2, col = Subtype)) + 
  geom_point() +
  youri_gg_theme +
  labs(colour="GlioVis GBM subtype")



# ---- NMF: TCGA [BFG/U133A] -----


nmf.input <- expression.matrix.tcga.vst.BFG.U133A %>%
  t() %>%
  #scale() %>%
  - min(data.frame(.)) %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(  gsub(".","-",colnames(.) , fixed=T ))



set.seed(12345) # same seed as Verhaak / Wang
#nmf.fit.2 <- nmf(df, 2 , seed = 12345)
nmf.fit.3 <- nmf(nmf.input, 3 , seed = 12345)
fit <- nmf.fit.3



p <- data.frame(sid = colnames(coef(fit)),
                nmf.1 = coef(fit)[1,],
                nmf.2 = coef(fit)[2,]
                #,nmf.3 = fit.3@fit@H[3,]
) %>%
  dplyr::left_join(metadata.tcga , by=c('sid' = 'sample'))

ggplot(p, aes(x = nmf.1 , y = nmf.2 , col= Subtype) ) + 
  geom_point()




#nmf.fit.4 <- nmf(nmf.input, 4 , seed = 12345)
#nmf.fit.5 <- nmf(nmf.input, 5 , seed = 12345)
#nmf.fit.6 <- nmf(nmf.input, 6 , seed = 12345)





ar <- metadata.tcga %>%
  dplyr::select('sample', 'Subtype') %>%
  tibble::column_to_rownames('sample')

coefmap(minfit(fit), info = TRUE, annCol = ar, hclustfun='complete', distfun='correlation')


basismap(fit, info = TRUE)
consensusmap(fit)

#basiscor(fit)
#profcor(fit)




# ---- NMF: TCGA [150s] -----


nmf.input <- expression.matrix.tcga.vst.150 %>%
  t() %>%
  scale() %>%
  - min(data.frame(.)) %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(  gsub(".","-",colnames(.) , fixed=T ))



set.seed(12345) # same seed as Verhaak / Wang
#nmf.fit.2 <- nmf(df, 2 , seed = 12345)
nmf.fit.3 <- nmf(nmf.input, 3 , seed = 12345)
fit <- nmf.fit.3



p <- data.frame(sid = colnames(coef(fit)),
                nmf.1 = coef(fit)[1,],
                nmf.2 = coef(fit)[2,]
                #,nmf.3 = fit.3@fit@H[3,]
) %>%
  dplyr::left_join(metadata.tcga , by=c('sid' = 'sample'))

ggplot(p, aes(x = nmf.1 , y = nmf.2 , col= Subtype) ) + 
  geom_point()




#nmf.fit.4 <- nmf(nmf.input, 4 , seed = 12345)
#nmf.fit.5 <- nmf(nmf.input, 5 , seed = 12345)
#nmf.fit.6 <- nmf(nmf.input, 6 , seed = 12345)





ar <- metadata.tcga %>%
  dplyr::select('sample', 'Subtype') %>%
  tibble::column_to_rownames('sample')

coefmap(minfit(fit), info = TRUE, annCol = ar, hclustfun='complete', distfun='correlation')


basismap(fit, info = TRUE)
consensusmap(fit)

#basiscor(fit)
#profcor(fit)


# ---- NMF: TCGA Arrays ----

# 364 / 369 TCGA GBMs? numbers are inconsistent w/ the paper

samples <- TCGA_GBM$pData %>%
  dplyr::filter(IDH1_status == "Wild-type") %>% 
  dplyr::filter(CIMP_status == "NON G-CIMP") %>%
  dplyr::mutate(Sample = gsub('\\.+', '-' , Sample) ) %>%
  dplyr::pull(Sample) %>%
  unique()


subtype <- TCGA_GBM$expr %>%
  dplyr::mutate(Sample = gsub('\\.+', '-' , Sample) ) %>%
  dplyr::filter(Sample %in% samples) %>%
  dplyr::select(c('Subtype', 'Sample'))


expr.arr <- TCGA_GBM$expr %>%
  dplyr::mutate(Sample = gsub('\\.+', '-' , Sample) ) %>%
  tibble::column_to_rownames('Sample') %>%
  dplyr::select(-c('Histology', 'Grade', 'Recurrence', 'Subtype', 'CIMP_status', 'survival', 'status')) %>%
  dplyr::filter(rownames(.) %in% samples) %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`( gsub("\\.", "-" , colnames(.) ) ) %>%
  `rownames<-`( gsub("\\.", "-" , rownames(.) ) )



genes <- wang.glioma.intrinsic.genes %>%
  dplyr::filter( U133A_Gene. == "YES") %>%
  dplyr::pull(Gene_Symbol)


expr.arr.sub <- expr.arr %>%
  dplyr::filter(rownames(.) %in% genes )


# dze doen t niet meer!?
#genes[genes %in% rownames(expr.arr) == F]


## perform actual NMF w/ lib ----


set.seed(12345) # same seed as Verhaak / Wang

nmf.input.scaled <- expr.arr.sub %>%
  t() %>%
  scale() %>%
  - min(data.frame(.)) %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`(  gsub(".","-",colnames(.) , fixed=T ))


nmf.fit <- nmf( expr.arr.sub, 3 , seed = 12345)
nmf.fit.scaled <- nmf( nmf.input.scaled , 3 , seed = 12345)



p <- data.frame(sid = colnames(coef(nmf.fit)),
                nmf.1 = coef(nmf.fit)[1,],
                nmf.2 = coef(nmf.fit)[2,]
) %>%
  dplyr::left_join(subtype , by=c('sid' = 'Sample'))

ggplot(p, aes(x = nmf.1 , y = nmf.2 , col= Subtype) ) + 
  geom_point()




coefmap(minfit(fit), info = TRUE, annCol = ar, hclustfun='complete', distfun='correlation')


basismap(fit, info = TRUE)
consensusmap(fit)



## perform verhaak/wang CNMF ----

source('scripts/msig.library.12.R')

# NMF / NMF.div difference w/ normalisation?

a.3 = NMF( as.matrix(expr.arr.sub)  , 3, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10)
#a.4 = NMF( as.matrix(expr.arr.sub)  , 4, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10)
b.3 = NMF.div( as.matrix(expr.arr.sub)  , 3, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10)
#b.4 = NMF.div( as.matrix(expr.arr.sub)  , 4, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10)

# near identical convergence
plot(a.3$H[1,] , b.3$H[1,])
plot(a.3$H[2,] , b.3$H[2,])
plot(a.3$H[3,] , b.3$H[3,])


# make connectivity matrix? - assign kan veel efficienter
assign <- matrix(0, nrow = 3, ncol = ncol(expr.arr.sub) ) # assigned to which Metagene
for (i in 1:3) {
  for (j in 1:ncol(expr.arr.sub)  ) { # Find membership
    class <- order(a.3$H[,j], decreasing=T)
    assign[i, j] <- class[1]
  }
}


# assignment as in their code:
plot( as.factor(assign[1,]) ,  p$Subtype)



# ---- test normalised ----

a = t(scale(t(expr.arr.sub))) 
a = a - min(a) + 0.000001

a.3.nm = NMF( as.matrix(a)  , 3, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10)


# near identical convergence
plot(a.3.nm$H[1,] , a.3$H[1,])
plot(a.3.nm$H[2,] , a.3$H[2,])
plot(a.3.nm$H[3,] , a.3$H[3,])


# make connectivity matrix? - assign kan veel efficienter
assign <- matrix(0, nrow = 3, ncol = ncol(expr.arr.sub) ) # assigned to which Metagene
for (i in 1:3) {
  for (j in 1:ncol(expr.arr.sub)  ) { # Find membership
    class <- order(a.3.nm$H[,j], decreasing=T)
    assign[i, j] <- class[1]
  }
}


# assignment as in their code:
plot( as.factor(assign[1,]) ,  p$Subtype)




# ---- plot ----

plt.a <- a.3$H %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  `colnames<-`( paste0('NMF.', 1:ncol(.) ) ) %>%
  dplyr::mutate(sid = p$sid) %>%
  dplyr::mutate(subtype = p$Subtype)




ggplot(plt.a, aes(x = NMF.1, y = NMF.2, col=subtype)) +
  geom_point()




plot_ly(x=plt.a$NMF.1, y = plt.a$NMF.2 , z=plt.a$NMF.3 , type="scatter3d", color=plt.a$subtype, size=5,
        colors = "Dark2")



