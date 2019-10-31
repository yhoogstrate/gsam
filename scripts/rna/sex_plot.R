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


source("scripts/R/job_gg_theme.R")

# ---- unsupervised DESeq2 PCA / MDS / clustering / t-SNE ----

## PCA: pc1 & pc2  x  resection
e <- expression_matrix
# make sure order metadata aligns with expression data (!)
stopifnot(sum(colnames(e) == gsam.rna.metadata[match(colnames(e) , gsam.rna.metadata$sample.id),]$sample.id) == ncol(e))

# resection as condition
cond <- factor(paste0("resection",gsam.rna.metadata[match(colnames(e) , gsam.rna.metadata$sample.id),]$resection))
colnames(e) == gsam.rna.metadata[match(colnames(e) , gsam.rna.metadata$sample.id),]$sample.id

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds,blind=T))

# remove effect of gender
e.vst <- removeBatchEffect(e.vst, cond)

plot(e.vst[rownames(e.vst) == "ENSG00000229807.12_5",],
     col=as.numeric(as.factor(gsam.rna.metadata$gender)) , ylab="ENSG00000229807 / XIST VST nromalised expression" ,
     xlab="Sample", pch=19,cex=0.7)
abline(h=8.5,col="gray",lty=2)
# True gender = > 8.5

plot(e.vst2[rownames(e.vst2) == "ENSG00000229807.12_5",],
     col=as.numeric(as.factor(gsam.rna.metadata$gender)) , ylab="ENSG00000229807 / XIST VST nromalised expression" ,
     xlab="Sample", pch=19,cex=0.7)
#abline(h=8.5,col="gray",lty=2)
# True gender = > 8.5


variances <- rowVars(e.vst)

#plot(sort(variances,decreasing=T))
plot(c(1,1000),c(0,max(variances)),type="n")
points(sort(variances,decreasing=T))

ntop <- sum(variances>2)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]
high_variance_genes <- rownames(high_variance_genes)
high_variance_genes <- gene_matrix[match(high_variance_genes, gene_matrix$Geneid),c(1,2)]

xist <- high_variance_genes[high_variance_genes$Chr == "chrX",]$Geneid[1:1]
genes_y <- high_variance_genes[high_variance_genes$Chr  == "chrY",]$Geneid[1:2]

xist <- e.vst[rownames(e.vst) %in% xist,]
genes_y <- e.vst[rownames(e.vst) %in% genes_y,]
genes_y <- colSums(genes_y)


tmp <- data.frame(
  x = xist,
  y = genes_y,
  gender = as.factor(gsam.patient.metadata[match(gsub("[12]$","",colnames(e.vst)), gsam.patient.metadata$studyID),]$gender) , 
  sid = colnames(e.vst)
)


mislabels <- subset(tmp, 
                    (x < 8.5) & (gender == "Female")
                    |
                    (x >= 8.5) & (gender == "Male"))


gg <- ggplot(tmp, aes(x=x, y=y, label=sid)) + 
  geom_point(aes(col=gender)) + 
  geom_encircle(aes(x=x, y=y), 
                data=subset(tmp, x < 8.5),
                color="grey50", 
                size=0.75, 
                expand=0.08) +
  geom_encircle(aes(x=x, y=y), 
                data=subset(tmp, x >= 8.5),
                color="grey50", 
                size=0.75, 
                expand=0.08) +
  geom_text_repel(
    nudge_y       = 23 - mislabels$y,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x",
    data = mislabels
  ) +
  job_gg_theme +
  ylim(NA, 23) + 
  labs(x = "XIST (vst transformed expression)",y="chrY (sum of normalised expression)") +
  ggtitle("Sex plot RNA-seq data [1st 96 samples]")

plot(gg)
ggsave("output/figures/rna/sex-plot_dna_rna.pdf")





