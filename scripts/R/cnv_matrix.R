#!/usr/bin/env R


cnv_matrix <- read.table("data/output/tables/cnv_copynumber-ratio.cnr_log2_all.txt",stringsAsFactors = F,header=T)
naPerRow <- apply(cnv_matrix, 1, function(tmp) sum(is.na(tmp)))
tmp <- nrow(cnv_matrix)
# there are four samples of which the CNV data is missing at this moment, and are always NA
cnv_matrix <- cnv_matrix[naPerRow <= 4,]
print(paste0("Removal of ",tmp - nrow(cnv_matrix) , " CNV regions only present in a subset of samples, because of different batches?" ))

naPerCol <- apply(cnv_matrix, 2, function(x) sum(is.na(x)))
tmp <- ncol(cnv_matrix)
cnv_matrix <- cnv_matrix[,naPerCol == 0]
print(paste0("Removal of ",tmp - ncol(cnv_matrix) , " CNV regions due to missing/corrupt and excluded files" ))

rm(tmp)

