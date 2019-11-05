# Set wd

# ---- initialization setup ----

setwd("~/projects/gsam")


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
source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
ensembl_genes <- get_ensembl_hsapiens_gene_ids()


# ---- load: data ----
# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)

source("scripts/R/dna_idh.R")

source("scripts/R/chrom_sizes.R")




# ---- PCA: pc1 & pc2  x  IDH(1+2) ----


e <- expression_matrix
stopifnot(sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id) == ncol(e))

# resection as condition
#cond <- factor(paste0("resection",gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$resection))
cond <- as.factor(paste0("idh.",dna_idh[match(colnames(e), dna_idh$donor_ID),]$IDH.mut != "-"))
colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds,blind=T))

# remove effect of gender
#e.vst <- removeBatchEffect(e.vst, cond)




ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]



pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)


# now with full data
#e2 <- expression_matrix_full
e2 <- expression_matrix_full[,match(colnames(e), colnames(expression_matrix_full))]
cond <- as.factor(paste0("idh.",dna_idh[match(colnames(e2), dna_idh$donor_ID),]$IDH.mut != "-"))
dds <- DESeqDataSetFromMatrix(e2, DataFrame(cond), ~cond)
e.vst2 <- assay(vst(dds,blind=T))
#high_variance_genes2 <- e.vst2[match(rownames(high_variance_genes), rownames(e.vst2)),]
high_variance_genes2 <- e.vst2[match(rownames(high_variance_genes), rownames(e.vst2)),]

pc <- prcomp(t(high_variance_genes2))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)


# re-apply
e3 <- expression_matrix_full
cond <- as.factor(paste0("idh.",dna_idh[match(colnames(e3), dna_idh$donor_ID),]$IDH.mut != "-"))
dds <- DESeqDataSetFromMatrix(e3, DataFrame(cond), ~cond)
e.vst3 <- assay(vst(dds,blind=T))
e.vst3 <- e.vst3[match(rownames(high_variance_genes), rownames(e.vst3)),]
pcn <- scale(t(e.vst3), pc$center, pc$scale) %*% pc$rotation

plot(pcn[,1], pcn[,2], cex=0.7,pch=19,col=as.numeric(cond)+1)







# ---- PCA: pc1 & pc2  x  loss10   ----
png("output/figures/unspervised_expression_analysis_vst_pca_x_loss10.png",width=480*2,height=480*2,res=72*2)

e <- expression_matrix

stopifnot(sum(colnames(e) == gsam.metadata[match(colnames(e) , gsam.metadata$sample.id),]$sample.id) == ncol(e))
#cond <- as.factor(gsub("[ \\+\\.\\/\\-]{1,}",".",gsam.metadata$primary710,fixed=F)) # chr7gain + chr10 loss
#cond2 <- as.factor(gsub("FALSE","DNA.wt", gsub("TRUE","DNA.IDH",as.character(colnames(e) %in% is_idh)))) # idh
cond <- as.factor(paste0("chr10.",gsam.metadata$primary10))

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds))

ntop <- 500
variances <- rowVars(e.vst)
high_variance_genes <- e.vst[order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))],]

pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2

plot( c(min(pc$x[,pc1]), max(pc$x[,pc1])), c(min(pc$x[,pc2]) - 15,max(pc$x[,pc2])), type="n" , main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),xlab=paste0("PCA pc-",pc1),ylab=paste0("PCA pc-",pc2))
points(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,col=as.numeric(cond)+2)


legend("bottomright", unique(as.character(cond)),
       col=unique(as.numeric(cond) + 2),pch=19)







rm(pc, pc1, pc2, pch, cond , cond2, ntop, e, e.vst, dds, high_variance_genes, variances)


dev.off()



