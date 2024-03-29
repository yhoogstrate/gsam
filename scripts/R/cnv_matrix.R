#!/usr/bin/env R

# load ----

cnv_matrix <- read.table("data/gsam/output/tables/cnv_copynumber-ratio.cnr_log2_all.txt",stringsAsFactors = F,header=T) %>%
  dplyr::rename(chr = chromosome) %>%
  dplyr::mutate(chr = paste0('chr', chr) ) %>%
  dplyr::mutate(loc = paste0(chr,":",start, "-", end)) %>%
  dplyr::arrange(loc)



#cnv_matrix[ncol(cnv_matrix)-4:ncol(cnv_matrix)]



## make a separate table for the gene information ----
cnv_matrix_genes <- cnv_matrix %>%
  dplyr::select(c('chr','start','end','gene', 'loc')) %>%
  tibble::column_to_rownames('loc')
  
cnv_matrix <- cnv_matrix  %>%
  dplyr::select( - c('chr','start','end','gene')) %>%
tibble::column_to_rownames('loc')



#cnv_matrix[ncol(cnv_matrix)-4:ncol(cnv_matrix)]



stopifnot(rownames(cnv_matrix) == rownames(cnv_matrix_genes))


## 


# remove rows with NA values ----

naPerRow <- apply(cnv_matrix, 1, function(tmp) sum(is.na(tmp)))
tmp <- nrow(cnv_matrix)
sel <- naPerRow <= 4
cnv_matrix <- cnv_matrix[sel,]
cnv_matrix_genes <- cnv_matrix_genes[sel,]
print(paste0("Removal of ",tmp - nrow(cnv_matrix) , " CNV regions only present in a subset of samples, because of different batches?" ))



# remove columns with NA values ----

naPerCol <- apply(cnv_matrix, 2, function(x) sum(is.na(x)))
tmp <- ncol(cnv_matrix)
cnv_matrix <- cnv_matrix[,naPerCol == 0]
print(paste0("Removal of ",tmp - ncol(cnv_matrix) , " CNV regions due to missing/corrupt and excluded files" ))
stopifnot(sum(is.na(head(cnv_matrix))) == 0) # no NA values may exist


# add check on colnames ----

colnames(cnv_matrix) <- gsub("(.b1|.b2)$","",colnames(cnv_matrix))


# change PD ids (sequencing ids) to study ids ----
if(!exists('gsam.cnv.metadata')) {
  source("scripts/R/gsam_metadata.R")
}


#
tmp <- match(colnames(cnv_matrix), gsam.cnv.metadata$PD_ID) # order
stopifnot(sum(is.na(tmp)) == 0)
colnames(cnv_matrix) <- gsam.cnv.metadata[tmp,]$donor_ID


# reorder to match metadata order
#cnv_matrix <- cnv_matrix[,match(gsam.cnv.metadata$donor_ID , colnames(cnv_matrix))]
#rm(tmp, naPerRow, naPerCol, sel)



# 〰 © Dr. Youri Hoogstrate 〰 ----


