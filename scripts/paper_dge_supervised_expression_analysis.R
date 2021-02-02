# Set wd

# load ----
## libs ----

library(tidyverse)

#library(DESeq2)
#library(ggplot2)
#library(ggrepel)
#library(pheatmap)
#library(fgsea)
#library(limma)
#library(EnhancedVolcano)

## functions & templates ----

#source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
#ensembl_genes <- get_ensembl_hsapiens_gene_ids()

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")


## data ----
# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")
source("scripts/R/subtype_genes.R")

source("scripts/R/gsam_metadata.R")
source("scripts/R/gsam_rna-seq_expression.R")
source("scripts/R/dna_idh.R") # IDH mutants must be excluded

source("scripts/R/glass_metadata.R")

source('scripts/R/wang_glioma_intrinsic_genes.R')


#source("scripts/R/chrom_sizes.R")


# data transformation(s) ----

# - exclude IDH mutants
# - exclude low tumour percentage


# some code? ----


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


# DE ----
## TP ~ R1 (~15k genes) ----


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



