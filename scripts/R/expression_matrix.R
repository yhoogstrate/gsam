#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
#Load data
gencode.31 <- read.delim("data/ref/star-hg19/gencode.v31lift37.annotation.gtf", comment.char="#",stringsAsFactors = F,header=F) %>%
  dplyr::filter(V3 == "gene") %>%
  dplyr::mutate(ENSG = gsub("^.+(ENSG[^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(GENE = gsub("^.+gene_name ([^;]+);.+$","\\1",V9)) %>%
  dplyr::mutate(V9 = NULL)


# file is currently not at server and @ wrong permissions @ ccbc
#expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#")
expression_matrix_full_new <- read.delim("data/output/tables/gsam_featureCounts_readcounts_new.txt",stringsAsFactors = F,comment="#") %>%
  `colnames<-`(gsub("^.+RNA.alignments\\.(.+)\\.Aligned.sortedByCoord.+$","\\1", colnames(.) ,fixed=F)) %>%
  dplyr::left_join(gencode.31, by=c('Geneid' = 'ENSG') ) %>%
  dplyr::mutate(rn = paste0 ( Geneid,"|", GENE , "|", unlist(lapply(str_split(Chr,";") , unique)), ':', unlist(lapply(str_split(Start,";") , min)), '-', unlist(lapply(str_split(End,";") , max)) , '(', unlist(lapply(str_split(Strand,";") , unique)), ')' ) ) %>%
  dplyr::mutate(Chr=NULL, Start = NULL, End = NULL, Strand = NULL, Length=NULL, Geneid=NULL) %>% dplyr::select(-paste0("V",1:8)) %>%
  tibble::column_to_rownames("rn")