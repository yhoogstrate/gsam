# Set wd

# ---- initialization setup ----

setwd("~/projects/gsam")


# ---- load: libs ----

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(fgsea)
library(limma)
library(EnhancedVolcano)


# ---- load: functions ----

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

source("scripts/R/job_gg_theme.R")




# ---- |--------------| ----

# read new table

expression_matrix_full <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt",stringsAsFactors = F,comment="#")
colnames(expression_matrix_full) <- gsub("^[^\\.]+\\.([^\\.]+)\\..+$","\\1",colnames(expression_matrix_full),fixed=F)

#gene_matrix <- expression_matrix[,1:6]

rownames(expression_matrix_full) <- expression_matrix_full$Geneid
expression_matrix_full$Geneid <- NULL
expression_matrix_full$Chr <- NULL
expression_matrix_full$Start <- NULL
expression_matrix_full$End <- NULL
expression_matrix_full$Strand <- NULL
expression_matrix_full$Length <- NULL

# old and new values correlate very high
#em <- expression_matrix[,20:30]
#tmp <- expression_matrix_full[,match(colnames(em) , colnames(expression_matrix_full))]
#colnames(tmp) <- paste0(colnames(tmp), ".f")
#tmp <- cbind(em,tmp)
#c <- cor(tmp,method="spearman")
#corrplot(c)


# resection as condition
e <- expression_matrix_full
cond <- factor(paste0("r",gsub("^.+([0-9])$","\\1",colnames(e))))

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds,blind=T))


ntop <- 500
variances <- rowVars(e.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

idh = as.factor(gsam.rna.metadata$IDH.mut != "-")

pc <- prcomp(t(high_variance_genes))
pc1 <- 3
pc2 <- 4
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=idh)


legend("bottomleft", unique(as.character(cond)),col=unique(as.numeric(cond) + 1),pch=19)




# https://stats.stackexchange.com/questions/2592/how-to-project-a-new-vector-onto-pca-space

#e <- expression_matrix_full



e2 <- expression_matrix_full #[,match(colnames(e), colnames(expression_matrix_full))]
cond <- as.character( dna_idh[match(colnames(e2), dna_idh$donor_ID),]$IDH.mut != "-" )
cond[is.na(cond)] <- "NA"
cond <- as.factor(cond)

#pcn <- scale(e2, pc$center, pc$scale) %*% pc$rotation

select <- match(rownames(high_variance_genes), rownames(e.vst2))
high_variance_genes2 <- e.vst2[select,]

pc <- prcomp(t(high_variance_genes2))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)


pc1 <- 1
pc2 <- 2
plot(pcn[,c(pc1,pc2)],  cex=0.7,pch=19,
     main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"),col=as.numeric(cond)+1)


# ---- * DE: res1 <=> res2 ----


# load VST data
e <- expression_matrix_full
m <- gsam.rna.metadata
e <- e[match(gsub('-','.',m$sid,fixed=T), colnames(e))] # match with metadata; those that are of suff quali
stopifnot( sum(gsub(".","-",colnames(e),fixed=T) != m$sid) == 0 ) # throw error if id's do not match with metadata

# remove low QC samples [remove blacklist.pca == T]
e <- e[,m$blacklist.pca == F]
e <- e[rowSums(e) / ncol(e) > 4,]
m <- m[m$blacklist.pca == F,]
stopifnot( sum(gsub(".","-",colnames(e),fixed=T) != m$sid) == 0 ) # throw error if id's do not match with metadata

# simple test, not correct for gender nor patient
cond <- as.factor(gsub("r","resection", m$resection,fixed=T))
dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
dds <- DESeq(dds)

res <- results(dds)
ens <- ensembl_genes[match(gsub("\\..+$","",rownames(res)), ensembl_genes$ensembl_gene_id),]$hgnc_symbol
ens[ens == ""] <- NA
ens[is.na(ens)] <- gsub("\\..+$","",rownames(res)[is.na(ens)])
res$gene <- ens
rm(ens)
res <- res[order(res$padj),]


EnhancedVolcano(res,
                lab = res$gene,
                x = 'log2FoldChange',
                y = 'pvalue')

#save.image(file = "dge_supervised_expression_analysis.RData")
#load(file = "dge_supervised_expression_analysis.RData")



