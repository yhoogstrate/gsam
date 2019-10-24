#!/usr/bin/env R


setwd("~/projects/gsam")

# ---- load libs ----

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(limma)


# ---- load functions ----


# ---- load data----

source("scripts/R/gsam_metadata.R")
source("scripts/R/cnv_matrix.R")

source("scripts/R/chrom_sizes.R")

# there is a difference in size of results files of CNVKit
# they either have 27493 lines (all in ExtraSequencing_CN_data_29Mar2018/)
# or have 27287 lines (195/374) = 50%



# ---- PCA ----

"
We would like to exclude chrX and chrY regions unless we would like to sample swaps (with we also need to do)
"

cnv_matrix_autosomes <- cnv_matrix[!cnv_matrix_genes$chr %in% c("chrX", "chrY"),]
cnv_matrix_genes_autosomes <- cnv_matrix_genes[!cnv_matrix_genes$chr %in% c("chrX", "chrY"),]
print(paste0("Removal of ",nrow(cnv_matrix) - nrow(cnv_matrix_autosomes) , " CNV regions at chrX & chrY" ))





#dd <- dist(d)


ntop <- 5000
variances <- rowVars(as.matrix(d))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- d[select,]

dim(high_variance_genes)


dd <- dist(high_variance_genes)
dd2 <- dist(t(high_variance_genes))

fit <- cmdscale(dd,eig=TRUE, k=2)
fit2 <- cmdscale(dd2,eig=TRUE, k=2)


x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",pch=19,cex=0.5)







#chrs_hg19_s[genes$chr]


# ---- t-SNE ----


# ---- hierarchical clustering ---- 

e <- d[,1:50]

# normalise? median should be 0?
# colMedians(as.matrix(d))

ntop <- 2500
variances <- rowVars(as.matrix(e))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e[select,]


dist_mat <- dist(t(e), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)



# random correlation

plot(c(-3, 7), c(-3, 7), type="n")
points(d[,3] , d[,4], pch=19,cex=0.15)

lines(c(-2,10),c(-2,10),col="red")


# ---- pca ----


e <- d[,]#1:250]

# normalise? median should be 0?
# colMedians(as.matrix(d))

ntop <- 2500
variances <- rowVars(as.matrix(e))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e[select,]

# scree plot = pca variance plot
# screeplot(pc)

pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19, main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"))#,col=as.numeric(cond)+1)

#((pc$sdev / sum(pc$sdev) * 100)


# ---- Pheatmap ----

