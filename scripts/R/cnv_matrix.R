#!/usr/bin/env R


cnv_matrix <- read.table("data/output/tables/cnv_copynumber-ratio.cnr_log2_all.txt",stringsAsFactors = F,header=T)

# make a separate table for the gene information
cnv_matrix_genes <- cnv_matrix[,1:4]
cnv_matrix_genes$chr <- paste0('chr',cnv_matrix_genes$chromosome)
cnv_matrix_genes$loc <- paste0(cnv_matrix_genes$chr,":",cnv_matrix_genes$start, "-", cnv_matrix_genes$end)
cnv_matrix <- cnv_matrix[,-(1:4)] # exclude gene details
rownames(cnv_matrix) <- cnv_matrix_genes$loc

# reorder by sorted colnames
cnv_matrix <- cnv_matrix[,order(colnames(cnv_matrix))]



# remove rows with NA values
naPerRow <- apply(cnv_matrix, 1, function(tmp) sum(is.na(tmp)))
tmp <- nrow(cnv_matrix)
sel <- naPerRow <= 4
cnv_matrix <- cnv_matrix[sel,]
cnv_matrix_genes <- cnv_matrix_genes[sel,]
print(paste0("Removal of ",tmp - nrow(cnv_matrix) , " CNV regions only present in a subset of samples, because of different batches?" ))



# remove columns with NA values
naPerCol <- apply(cnv_matrix, 2, function(x) sum(is.na(x)))
tmp <- ncol(cnv_matrix)
cnv_matrix <- cnv_matrix[,naPerCol == 0]
print(paste0("Removal of ",tmp - ncol(cnv_matrix) , " CNV regions due to missing/corrupt and excluded files" ))



rm(tmp, naPerRow, naPerCol, sel)

