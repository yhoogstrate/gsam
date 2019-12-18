# Set wd

# ---- initialization setup ----

setwd("~/projects/gsam")


# ---- load: libs ----

library(DESeq2)
library(ggplot2)
library(ggalt)
library(limma)
library(ggrepel)


# ---- load: functions ----

# ---- load: data ----
# Load the data (expression values; metadata) into data frames


source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)

source("scripts/R/dna_idh.R")

source("scripts/R/job_gg_theme.R")



# ---- sex plot [all samples] ----



e <- expression_matrix_full
# make sure order metadata aligns with expression data (!)
#stopifnot(sum(colnames(e) == gsam.rna.metadata[match(colnames(e) , gsam.rna.metadata$sample.id),]$sample.id) == ncol(e))

# resection as condition
cond <- as.factor(gsub("^.+([0-9])$","resection\\1",colnames(e)))
#colnames(e) == gsam.rna.metadata[match(colnames(e) , gsam.rna.metadata$sample.id),]$sample.id

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds,blind=T))

variances <- rowVars(e.vst)
names(variances) <- rownames(e.vst)
variances <- sort(variances,decreasing=T)

ntop <- sum(variances>2)
select <- variances[variances > 2]
select <- names(select)
select <- match(select, rownames(e.vst))

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
                    (y >= 12.5) & (gender == "Female")
                    |
                      (y < 12.5) & (gender == "Male"))


gg <- ggplot(tmp, aes(x=x, y=y, label=sid)) + 
  geom_point(aes(col=gender)) + 
  geom_point(aes(col=gender),data=mislabels) + 
  geom_encircle(aes(x=x, y=y), 
                data=subset(tmp, y < 12.5),
                color="grey50", 
                size=0.75, 
                expand=0.08) +
  geom_encircle(aes(x=x, y=y), 
                data=subset(tmp, y >= 12.5),
                color="grey50", 
                size=0.75, 
                expand=0.08) +
  geom_text_repel(
    nudge_y       = 22 - mislabels$y,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x",
    data = mislabels
  ) +
  job_gg_theme +
  ylim(NA, 23) + 
  labs(x = "XIST (vst transformed expression)",y="chrY (sum of normalised expression)") +
  ggtitle(paste0("Sex plot RNA-seq data [all ",ncol(e)," samples]"))

plot(gg)







ggsave("output/figures/rna/sex-plot_dna_rna.pdf")
ggsave("output/figures/rna/sex-plot_dna_rna.png")






