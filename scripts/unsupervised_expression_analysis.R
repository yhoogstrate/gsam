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
source("scripts/R/chrom_sizes.R")

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
#cond <- factor(paste0("resection",gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$resection))

# gender
#cond <- as.character((gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$donor_sex))
#cond[is.na(cond)] <- "NA"
#cond <- as.factor(cond)




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
rm(pc, pc1, pc2)



dev.off()



# ---- MDS ----
ntop <- 5000
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

#pc <- prcomp(t(high_variance_genes))
sampleDists <- dist(t(high_variance_genes))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists)


# ---- do DE - rechtse wolk----

ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- prcomp(t(high_variance_genes))
#sampleDists <- dist(t(high_variance_genes))

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
res <- res[order(res$padj),]
o <- match(rownames(res) , gene_matrix$Geneid)
res$chr <- gene_matrix[o,]$Chr
res$start <- as.numeric(gene_matrix[o,]$Start)

sig <- res[res$padj < 0.01,]


cdkn2ab_r <- res[res$chr == "chr9" & res$start > 21600000 & res$start < 22300000,]




write.table(sig, "/tmp/sig-right-area.txt")



png("output/figures/unsupervised_expression_analysis__de_right_pca_cloud.png", width=480*3.5 , height=480*2 , res=72*2)

#plot(c(chrs_hg19_s["chr6"] + 10000000, chrs_hg19_s["chr6"] + 40000000 ), c(min(sig$log2FoldChange), max(sig$log2FoldChange)), pch=19, cex=0.2 , main="DE Rechter wolk PCA",type="n")
plot(c(chrs_hg19_s["chr1"] + 10000000, chrs_hg19_e["chr22"] ), c(min(sig$log2FoldChange), max(sig$log2FoldChange)), pch=19, cex=0.2 , main="DE Rechter wolk PCA",type="n",xlab="Chromosomal location",xaxt="n",ylab="LogFC DE Genes (padj<0.01)")
points(chrs_hg19_s[sig$chr] + as.numeric(sig$start) , sig$log2FoldChange, pch=19, cex=0.2)
abline(v=chrs_hg19_s)

#chr19 <- sig[sig$chr == "chr10",]
#plot(as.integer(gsub("M","",chr19$start)), chr19$log2FoldChange, pch=19,cex=0.6)

for(chr in names(chrs_hg19)) {
  text( (chrs_hg19_s[chr] + chrs_hg19_e[chr]) / 2 , -6 , chr, srt=90,cex=0.7)
}

#abline(v=chrs_hg19_s["chr6"] + 30000000, col="red", lwd=2)
#abline(v=chrs_hg19_s["chr6"] + 34000000, col="red", lwd=2)

dev.off()






# ---- de group having some IDH1 mutants ----

# create VST transformed data
cond <- as.factor(round(runif(ncol(e))) + 1)
dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))
rm(cond, dds)

# do PCA on VST transformed data
ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]
pc <- prcomp(t(high_variance_genes))

# create subgroups based on that PCA
outliers2 = rownames(pc$x)[ pc$x[,1] < 0 & pc$x[,2] > 12.5] 
cond <- as.factor(colnames(e) %in% outliers2)

# do DE analyses based on the subgroups
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
o <- match(rownames(res) , gene_matrix$Geneid)
res$chr <- gene_matrix[o,]$Chr
res$start <- as.numeric(gene_matrix[o,]$Start)

sig <- res[res$padj < 0.01,]


cdkn2ab_tl <- res[res$chr == "chr9" & res$start > 21600000 & res$start < 22300000,]



png("output/figures/unsupervised_expression_analysis__de_topleft_idh_corner.png", width=480*3.5 , height=480*2 , res=72*2)

plot(c(chrs_hg19_s["chr1"] + 10000000, chrs_hg19_e["chr22"] ), c(min(sig$log2FoldChange), max(sig$log2FoldChange)), pch=19, cex=0.2 , main="DE top-left corner w/IDH samples in PCA",type="n",xlab="Chromosomal location",xaxt="n",ylab="LogFC DE Genes (padj<0.01)")

points(chrs_hg19_s[sig$chr] + as.numeric(sig$start) , sig$log2FoldChange, pch=19, cex=0.2)
abline(v=chrs_hg19_s)

#chr19 <- sig[sig$chr == "chr10",]
#plot(as.integer(gsub("M","",chr19$start)), chr19$log2FoldChange, pch=19,cex=0.6)

for(chr in names(chrs_hg19)) {
  text( (chrs_hg19_s[chr] + chrs_hg19_e[chr]) / 2 , -5.25 , chr, srt=90,cex=0.7)
}

#abline(v=chrs_hg19_s["chr6"] + 30000000, col="red", lwd=2)
#abline(v=chrs_hg19_s["chr6"] + 34000000, col="red", lwd=2)

dev.off()




write.table(sig, "/tmp/sig-upper-corner-cluster.txt")


# fgsea
getLeadingSymbols <- function(x){
  z <- lapply(x$leadingEdge, function(y){ paste(unique(na.omit(clusterProfiler::bitr(y, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = org.Hs.eg.db::org.Hs.eg.db)$SYMBOL)), collapse = ', ')})
  return(unlist(z))
}



ens <- gsub("\\..+$","",rownames(res))
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c("hgnc_symbol","entrezgene_id","ensembl_gene_id","chromosome_name","start_position","end_position"),
  mart = mart)
entrez <- genes[ match(ens, genes$ensembl_gene_id) ,]$entrezgene_id
stopifnot(length(entrez) == length(ens))

res$entrez <- entrez
gseaDat <- res[!is.na(res$entrez),]
gseaDat <- gseaDat[!duplicated(gseaDat$entrez),] # conversion leads to two duplicated ENTREZ ids, remove them

pws <- reactomePathways(as.character(gseaDat$entrez))
pws <- fgsea::gmtPathways('misc/h.all.v6.2.entrez.gmt')
pws <- fgsea::gmtPathways('misc/msigdb.v6.2.entrez.gmt')
#names(pws) <- gsub("HALLMARK_","",names(pws),fixed=T)


gseaDat <- gseaDat [order(gseaDat$stat) , ]
ranks <- gseaDat$stat
#gseaDat <- gseaDat [order(gseaDat$log2FoldChange) , ]
#ranks <- gseaDat$log2FoldChange

names(ranks) <- gseaDat$entrez

fgseaRes <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results

fgseaRes <- fgseaRes[order(fgseaRes$pval),]


fgseaResa <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results
fgseaResb <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results
 
plot(fgseaResa$pval , fgseaResb$pval, pch=19,cex=0.1)



#plotEnrichment(pws, ranks, gseaParam = 1, ticksSize = 0.2)
plotEnrichment(pws$`Neuronal System`, ranks, gseaParam = 1, ticksSize = 0.2)

plotEnrichment(pws$`SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6`, ranks, gseaParam = 1, ticksSize = 0.2)

#pdf("gsea.msigdb.hallmark.pdf")
plotGseaTable(pws, ranks, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 0.8, 1.2, 1.2), render = TRUE)
#dev.off()


