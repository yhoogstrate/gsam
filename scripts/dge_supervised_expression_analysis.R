# Set wd

# ---- initialization setup ----

wd <- "/home/youri/projects/gsam"
wd <- "/home/yhoogstrate/projects/gsam"

setwd(wd)

rm(wd)


# ---- load: libs ----

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(limma)


# ---- load: functions ----
"
get_ensembl_hsapiens_gene_ids

Loads a function that allows to translate ENSEMBL IDs into HUGO symbols (handy for sharing tables) and Entrez IDs (handy for GSEa like analysis)
"
ensembl_genes <- get_ensembl_hsapiens_gene_ids()


# ---- load: data ----
# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)

source("scripts/R/dna_idh.R")


source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
ensembl_genes <- get_ensembl_hsapiens_gene_ids()

source("scripts/R/chrom_sizes.R")
source("scripts/R/plotDiffAnalysis.R")

# ---- DGE: rechtse wolk ----
# De wolk aan de rechterkant
# DGE / DESeq2 part
"
@todo factor out gender
"

e <- expression_matrix

stopifnot(sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id) == ncol(e))
cond <- as.factor(gsub("[ \\+\\.\\/\\-]{1,}",".",gsam.metadata$primary710,fixed=F)) # chr7gain + chr10 loss
cond2 <- as.factor(gsub("FALSE","DNA.wt", gsub("TRUE","DNA.IDH",as.character(colnames(e) %in% is_idh)))) # idh

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds, blind=T))


ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- prcomp(t(high_variance_genes))
outliers <- rownames(pc$x)[ pc$x[,1] > 14 ] 
cond <- as.factor(colnames(e) %in% outliers)


dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))

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
res <- res[order(res$padj),]
res$gid <- ens_ids[match(gsub("\\..+","",rownames(res)), gsub("\\..+","",ens_ids$ens.id)),]$gene.id
res <- res[order(res$padj),]
# Add things like chr/pos , HGNC gene symbol & Entrez gene id
o <- match(rownames(res) , gene_matrix$Geneid)
res$chr <- gene_matrix[o,]$Chr
res$start <- as.numeric(gene_matrix[o,]$Start)
rownames(res) <- gsub("\\..+$","",rownames(res)) 
#
write.table(res, "output/tables/de_right_cloud.txt")



# subtract significant ones (not for GSEA, which requires all of them)
sig <- res[res$padj < 0.01,]



# take a look at some specific genes/regions
#cdkn2ab_r <- res[res$chr == "chr9" & res$start > 21600000 & res$start < 22300000,]
#mgmt




##write.table(sig, "/tmp/sig-right-area.txt")



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



## GSEA / fgsea part

```{r}
res$entrez <- ensembl_genes[ match(gsub("\\..+$","",rownames(res)), ensembl_genes$ensembl_gene_id) ,]$entrezgene_id
stopifnot(length(res$entrez) == nrow(res))

dim(res)
#[1] 26729    10

gseaDat <- res[!is.na(res$entrez),]
dim(gseaDat)
#[1] 16366    10

gseaDat <- gseaDat[!duplicated(gseaDat$entrez),] # conversion leads to two duplicated ENTREZ ids, remove them
dim(gseaDat)
#[1] 16346    10


pws <- reactomePathways(as.character(gseaDat$entrez))
#pws <- fgsea::gmtPathways('misc/h.all.v6.2.entrez.gmt')
#pws <- fgsea::gmtPathways('misc/msigdb.v6.2.entrez.gmt')


#gseaDat <- gseaDat [order(gseaDat$log2FoldChange) , ]
#ranks <- gseaDat$log2FoldChange
gseaDat <- gseaDat [order(gseaDat$stat) , ]
ranks <- gseaDat$stat
names(ranks) <- gseaDat$entrez

# do GSEA and reorder table
fgseaRes <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results
fgseaRes <- fgseaRes[order(fgseaRes$pval),]


# High scoring gene set:
plotEnrichment(pws$`REACTOME_ADAPTIVE_IMMUNE_SYSTEM`, ranks, gseaParam = 1, ticksSize = 0.2)

# Gene set scoring high in other comparison
plotEnrichment(pws$`SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6`, ranks, gseaParam = 1, ticksSize = 0.2)

#pdf("gsea.msigdb.hallmark.pdf")
plotGseaTable(pws, ranks, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 0.8, 1.2, 1.2), render = TRUE)
#dev.off()


# DGE: top left corner containing IDH mutants
## DGE / DESeq2
e <- expression_matrix

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
rownames(res) <- gsub("\\..+$","",rownames(res))

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


## GSEA
genes <- get_ensembl_hsapiens_gene_ids()
res$entrez <- genes[ match(ens, genes$ensembl_gene_id) ,]$entrezgene_id
stopifnot(length(res$entrez) == length(ens))

dim(res)
gseaDat <- res[!is.na(res$entrez),]
dim(gseaDat)
gseaDat <- gseaDat[!duplicated(gseaDat$entrez),] # conversion leads to two duplicated ENTREZ ids, remove them
dim(gseaDat)

pws <- reactomePathways(as.character(gseaDat$entrez))
#pws <- fgsea::gmtPathways('misc/h.all.v6.2.entrez.gmt') # 50 gene sets
pws <- fgsea::gmtPathways('misc/msigdb.v6.2.entrez.gmt')


gseaDat <- gseaDat [order(gseaDat$stat) , ]
ranks <- gseaDat$stat
#gseaDat <- gseaDat [order(gseaDat$log2FoldChange) , ]
#ranks <- gseaDat$log2FoldChange

names(ranks) <- gseaDat$entrez

fgseaRes <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results

fgseaRes <- fgseaRes[order(fgseaRes$pval),]


#fgseaResa <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results
#fgseaResb <- fgsea(pws, ranks, minSize=15, maxSize = 500, nperm=10000)# higher number of permutations for more stable results

plot(fgseaResa$pval , fgseaResb$pval, pch=19,cex=0.1)



#plotEnrichment(pws, ranks, gseaParam = 1, ticksSize = 0.2)
plotEnrichment(pws$`Neuronal System`, ranks, gseaParam = 1, ticksSize = 0.2)

plotEnrichment(pws$`SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6`, ranks, gseaParam = 1, ticksSize = 0.2)

#pdf("gsea.msigdb.hallmark.pdf")
plotGseaTable(pws, ranks, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 0.8, 1.2, 1.2), render = TRUE)
#dev.off()





# Survival analysis on groups

# First make subgroups, append to table, and add to metadata and plot to confirm if it fits
e <- expression_matrix

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
outliers_topleft_idh <- rownames(pc$x)[ pc$x[,1] < 0 & pc$x[,2] > 12.5] 
outliers_right <- rownames(pc$x)[ pc$x[,1] > 14 ] 
cond <- rep("main", ncol(e))
cond[ colnames(e) %in% outliers_topleft_idh ] <- "idh"
cond[ colnames(e) %in% outliers_right ] <- "right"
cond <- as.factor(cond)

plot(pc$x[,1], pc$x[,2], col = as.numeric(cond) + 1,pch=19,cex=0.5)


survival.data <- data.frame(sid = colnames(e),pid=gsub("[0-9]$","",colnames(e)), group=cond)
survival.data$studyID <- tmp[match(survival.data$pid , tmp$studyID),]$studyID
survival.data$survivalDays <- tmp[match(survival.data$pid , tmp$studyID),]$survivalDays
survival.data$status <- tmp[match(survival.data$pid , tmp$studyID),]$status
survival.data <- survival.data[!is.na(survival.data$survivalDays),] # exclude 2 NA's


s2 <- survfit(Surv(survival.data$survivalDays) ~ survival.data$group, data=survival.data)
#survplot(s2)
ggsurvplot(s2, data=survival.data, pval=T)

# geen verschil

```


