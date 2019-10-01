#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)
library(ggplot2)

# ---- load data ----
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)
source("scripts/R/ensembl_to_geneid.R")

source("scripts/R/dna_idh.R")

# ---- unsupervised DESeq2 PCA [ all genes, all samples ] ----
e <- expression_matrix

colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id

# EGFR
#egfr_high <- colnames(e.vst) [ e.vst[match("ENSG00000146648" ,gsub("\\..+$","",rownames(e.vst))),] > 12.5 ]
#cond <- as.factor(gsub("FALSE","EGFR.low",gsub("TRUE","EGFR.high",as.character(colnames(e) %in% egfr_high),fixed=T),fixed=T))

# cdk4
#cdk4_high <- colnames(e.vst) [ e.vst[match("ENSG00000135446" ,gsub("\\..+$","",rownames(e.vst))),] > 8.64 ]
#cond <- as.factor(colnames(e) %in% cdk4_high)


# exon read count
#cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$featurecounts.Assigned < 2000000))
#cond <- as.factor(gsub("TRUE", "low.Q", gsub("FALSE", "high.Q", cond, fixed=T), fixed=T))

# resection
cond <- factor(paste0("resection",gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$resection))

# gender
cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$donor_sex))
cond[is.na(cond)] <- "NA"
cond <- as.factor(cond)




# ---- PCA 1 & 2  x  resection ----

# MAKE SURE ALL IDs MATCH OR QUIT - otherwise metadata was swapped
stopifnot(
  sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id)
  ==
    ncol(e) )

# resection
cond <- factor(paste0("resection",gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$resection))
colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))
p <- plotPCA(vst(dds), intgroup=c("cond"),ntop=500)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black", nudge_x = 3, nudge_y = 2, size=2.75 )

rm(dds, e.vst, p, cond)


ggsave("output/figures/unspervised_expression_analysis_vst_pca_x_resection.png",width=7,height=7,scale=1.5)


# ---- PCA 3 & 4  x  gender ----
stopifnot(
  sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id)
  ==
    ncol(e) )

cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$donor_sex))
cond[is.na(cond)] <- "NA"
cond <- as.factor(cond)


png("output/figures/unspervised_expression_analysis_vst_pca_x_gender.png",width=480*2,height=480*2,res=72*2)

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))

ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- prcomp(t(high_variance_genes))
pc1 <- 3
pc2 <- 4
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)

legend("bottomleft", unique(as.character(cond)),col=unique(as.numeric(cond) + 1),pch=19)

# draw labels of outliers
select <- select <- pc$x[,pc1] < 4 & cond == "Female"
text( pc$x[select,pc1] , pc$x[select,pc1] , rownames(pc$x)[select], cex=0.7, pos=4)

rm(cond, dds, e.vst, ntop, variances, select, high_variance_genes, pc, pc1, pc2)


dev.off()


# ---- PCA 1 & 2  x  vIII ----
png("output/figures/unspervised_expression_analysis_vst_pca_x_vIII.png",width=480*2,height=480*2,res=72*2)


cond <- paste0("vIII.",gsub("FALSE","neg",gsub("TRUE","pos",as.character(gsam.metadata$v3.rnaseq.v3.percentage > 1.0),fixed=T),fixed=T))
cond[is.na(cond) | cond == "vIII.NA"] <- "NA"
cond <- as.factor(cond)

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))


stopifnot(
  sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id)
  ==
    ncol(e) )


ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)

legend("bottomleft", unique(as.character(cond)),col=unique(as.numeric(cond) + 1),pch=19)


dev.off()

# ---- PCA 1 & 2  x  IDH(1+2) ----
png("output/figures/unspervised_expression_analysis_vst_pca_x_IDH_1_2.png",width=480*2,height=480*2,res=72*2)

cond <- as.factor(gsub("FALSE","DNA.wt", gsub("TRUE","DNA.IDH",as.character(colnames(e) %in% is_idh))))

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))

stopifnot(
  sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id)
  ==
    ncol(e) )


ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)


# draw labels of outliers
select <- colnames(e) %in% is_idh
text( pc$x[select,pc1] , pc$x[select,pc2] , rownames(pc$x)[select], cex=0.7, pos=4)


legend("bottomleft", unique(as.character(cond)),col=unique(as.numeric(cond) + 1),pch=19)


dev.off()



# ---- do DE ----

outliers <- rownames(pc$x)[ pc$x[,1] > 14 ] 
cond <- as.factor(colnames(e) %in% outliers)


dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)

# Keep only on average 3 reads per sample
print(dim(dds))
keep <- (rowSums(counts(dds)) >= ncol(e) * 3)
dds <- dds[keep,]
print(dim(dds))

dds <- DESeq(dds)
res <- results(dds)
dim(res)
res <- res[!is.na(res$padj),] # remove those that have failed, probably 0,0,0,0, high, high, high, high?
res <- res[order(res$padj),]
dim(res)

res <- res[order(res$padj),]
res$gid <- ens_ids[match(gsub("\\..+","",rownames(res)), gsub("\\..+","",ens_ids$ens.id)),]$gene.id
sig <- res[res$padj < 0.01,]

o <- match(rownames(sig) , gene_matrix$Geneid)
sig$chr <- gene_matrix[o,]$Chr
sig$start <- gene_matrix[o,]$Start

write.table(sig, "/tmp/sig-right-area.txt")


chr19 <- sig[sig$chr == "chr10",]
plot(as.integer(gsub("M","",chr19$start)), chr19$log2FoldChange, pch=19,cex=0.6)




# ---- de group having some IDH1 mutants ----

outliers2 = rownames(pc$x)[ pc$x[,1] < 0 & pc$x[,2] > 12.5] 
cond <- as.factor(colnames(e) %in% outliers2)


dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)

# Keep only on average 3 reads per sample
print(dim(dds))
keep <- (rowSums(counts(dds)) >= ncol(e) * 3)
dds <- dds[keep,]
print(dim(dds))

dds <- DESeq(dds)
res <- results(dds)
res$gid <- ens_ids[match(gsub("\\..+","",rownames(res)), gsub("\\..+","",ens_ids$ens.id)),]$gene.id
dim(res)
res <- res[!is.na(res$padj),] # remove those that have failed, probably 0,0,0,0, high, high, high, high?
res <- res[order(res$padj),]
dim(res)

res <- res[order(res$padj),]
sig <- res[res$padj < 0.01,]

o <- match(rownames(sig) , gene_matrix$Geneid)
sig$chr <- gene_matrix[o,]$Chr
sig$start <- gene_matrix[o,]$Start


chr19 <- sig[sig$chr == "chr9",]
plot(as.integer(gsub("M","",chr19$start)), chr19$log2FoldChange, pch=19,cex=0.6)



#o <- match(rownames(sig) , gene_matrix$Geneid)
#sig$chr <- gene_matrix[o,]$Chr
#sig$start <- paste0(round(as.numeric(gene_matrix[o,]$Start) / 1000000), 'M')

write.table(sig, "/tmp/sig-upper-corner-cluster.txt")





