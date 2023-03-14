#!/usr/bin/env R

#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
#Load gene annotation data
gencode.31 <- read.delim("data/gsam/ref/star-hg19/gencode.v31lift37.annotation.gtf", comment.char="#",stringsAsFactors = F,header=F) %>%
  dplyr::filter(V3 == "gene") %>%
  dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(V9 = NULL)

# Load actual expression matrix
# New refers to the run containing the addtional ~25? samples re-done
expression_matrix_full_new <- read.delim("data/gsam/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#") %>%
  `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1", colnames(.) ,fixed=F)) %>%
  dplyr::left_join(gencode.31, by=c('Geneid' = 'ENSG') ) %>%
  dplyr::mutate(rn = paste0 ( Geneid,"|", GENE , "|", unlist(lapply(str_split(Chr,";") , unique)), ':', unlist(lapply(str_split(Start,";") , min)), '-', unlist(lapply(str_split(End,";") , max)) , '(', unlist(lapply(str_split(Strand,";") , unique)), ')' ) ) %>%
  dplyr::mutate(Chr=NULL, Start = NULL, End = NULL, Strand = NULL, Length=NULL, Geneid=NULL, GENE=NULL) %>%
  dplyr::select(-paste0("V",1:8)) %>%
  tibble::column_to_rownames("rn") %>%
  `colnames<-`(gsub ('.', '-', colnames(.) , fixed=T)) # bla.new -> bla-new




# RPKM = count_assigned_to_gene / (gene_length/1000 * totalReads /1e6)
# choose totalReads carefully as rRNA type outliers can mess things tramendously up



# unstranded counts, useful to determine possible DNA contamination ratios
# only force parsing this file when a flag has been set, and remove flag accordingly
if(exists("PARSE_EXPRESSION_MATRIX_UNSTRANDED")) {
  expression_matrix_full_new_unstranded <- read.delim("data/gsam/output/tables/gsam_featureCounts_readcounts.unstranded_new.txt",stringsAsFactors = F,comment="#") %>%
    `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1", colnames(.) ,fixed=F)) %>%
    dplyr::left_join(gencode.31, by=c('Geneid' = 'ENSG') ) %>%
    dplyr::mutate(rn = paste0 ( Geneid,"|", GENE , "|", unlist(lapply(str_split(Chr,";") , unique)), ':', unlist(lapply(str_split(Start,";") , min)), '-', unlist(lapply(str_split(End,";") , max)) , '(', unlist(lapply(str_split(Strand,";") , unique)), ')' ) ) %>%
    dplyr::mutate(Chr=NULL, Start = NULL, End = NULL, Strand = NULL, Length=NULL, Geneid=NULL, GENE=NULL) %>% dplyr::select(-paste0("V",1:8)) %>%
    tibble::column_to_rownames("rn")

  rm(PARSE_EXPRESSION_MATRIX_UNSTRANDED)
}


