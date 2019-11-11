# Set wd

# ---- initialization setup ----

setwd("~/projects/gsam")


# ---- load: libs ----

library(DESeq2)
library(ggplot2)
library(ggalt)
library(pheatmap)
library(fgsea)
library(limma)


# ---- load: functions ----
"
get_ensembl_hsapiens_gene_ids

Loads a function that allows to translate ENSEMBL IDs into HUGO symbols (handy for sharing tables) and Entrez IDs (handy for GSEa like analysis)
"
#source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
#ensembl_genes <- get_ensembl_hsapiens_gene_ids()


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
pcn <- scale(t(e.vst3), pc$center, pc$scale) %*% pc$rotation # rescale new data onto pca of first 96 samples

cond.idh <- as.factor(paste0("idh.",dna_idh[match(colnames(e3), dna_idh$donor_ID),]$IDH.mut != "-"))
cond.chr10 <- as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn))   , gsam.patient.metadata$studyID   ),]$primary10)
cond <- as.factor(paste0(cond.idh,".",cond.chr10))
plot(pcn[,1], pcn[,2], cex=0.7,pch=19,col=as.numeric(cond)+1)



tmp <- data.frame(
  sid = colnames(e3),
  pc1 = pcn[,1],  pc2 = pcn[,2],  pc3 = pcn[,3],  pc4 = pcn[,4],  pc5 = pcn[,5],  pc6 = pcn[,6],
  resection = as.factor(gsub("^.+([1-2])$","\\1",colnames(e3))),
  idh = as.factor(gsub("TRUE","mut",gsub("FALSE","wt",paste0("",dna_idh[match(colnames(e3), dna_idh$donor_ID),]$IDH.mut != "-")))),
  chr10 = as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn)),
                                                 gsam.patient.metadata$studyID   ),]$primary10),
  chr7 = as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn)),
                                                gsam.patient.metadata$studyID   ),]$primary7),
  atrx = as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn)),
                                                gsam.patient.metadata$studyID   ),]$ATRX),
  gender = as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn))   , gsam.patient.metadata$studyID   ),]$gender) ,
  
  mgmt = as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn)),
                                                gsam.patient.metadata$studyID   ),]$initialMGMT),
  
  tlc  = as.factor(gsam.patient.metadata[match( gsub("[12]$","",rownames(pcn))                              , gsam.patient.metadata$studyID   ),]$tumorLocation))
tmp$chr7.wt.chr10.wt <- tmp$chr10 == "Normal" & tmp$chr7 == "Normal"


# pca x co-amplificatie chr7&chr10
ggplot(tmp, aes(x=pc1, y=pc2, col=chr7.wt.chr10.wt)) +
  geom_point(size=2) +
  labs(
      title = paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)"),
              fill = "Chr7-wt & Chr10-wt") +
  job_gg_theme 
ggsave("output/figures/rna/pca_rna_x_chr7_chr10_ampl.png")


ggplot(tmp, aes(x=pc1, y=pc2, col=idh)) +
  geom_point(size=2) +
  ggtitle(paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)")) +
  labs(fill = "IDH-mutation",
       subtitle="*Only high confidence IDH mutations included") +
  job_gg_theme 
ggsave("output/figures/rna/pca_rna_x_idh.png")


ggplot(tmp, aes(x=pc1, y=pc2, col=atrx)) +
  geom_point(size=2) +
  geom_point(size=2, data=subset(tmp, atrx != "Wildtype")) +
  ggtitle(paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)")) +
  labs(fill = "ATRX") +
  job_gg_theme 
ggsave("output/figures/rna/pca_rna_x_atrx.png")



ggplot(tmp, aes(x=pc1, y=pc2, col=atrx)) +
  geom_point(size=2, data=subset(tmp, atrx == "Wildtype")) +
  geom_point(size=2, data=subset(tmp, atrx != "Wildtype")) +
  ggtitle(paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)")) +
  labs(fill = "ATRX") +
  job_gg_theme 
ggsave("output/figures/rna/pca_rna_x_atrx.png")


#geom_point(size=2, data=subset(tmp, atrx != "Wildtype")) +

ggplot(tmp, aes(x=pc1, y=pc2, col=mgmt)) +
  geom_point(size=2, data=subset(tmp, is.na(mgmt))) +
  geom_point(size=2, data=subset(tmp, mgmt != "NA")) +
  ggtitle(paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)")) +
  labs(fill = "MGMT") +
  job_gg_theme 
ggsave("output/figures/rna/pca_rna_x_mgmt.png")



ggplot(tmp, aes(x=pc1, y=pc2, col=mgmt)) +
  geom_point(size=2, data=subset(tmp, is.na(mgmt)),col="gray90") +
  geom_point(size=2, data=subset(tmp, mgmt != "NA")) +
  ggtitle(paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)")) +
  labs(fill = "MGMT") +
  geom_encircle(
                size=0.75,
                expand=0.04,
                aes(col=mgmt),
                data = subset(tmp, !is.na(tmp$mgmt))
                ) +
  scale_x_continuous(expand = c(0.08, 0)) +
  scale_y_continuous(expand = c(0.124, 0)) +
  job_gg_theme 

ggsave("output/figures/rna/pca_rna_x_mgmt.png")





ggplot(tmp, aes(x=pc1, y=pc2, col=resection)) +
  geom_point(size=2) +
  ggtitle(paste0("PCA GSAM RNA normalised expression (",nrow(tmp)," samples)")) +
  labs(fill = "MGMT") +
  geom_encircle(
    size=0.75,
    expand=0.04,
    aes(col=mgmt),
    data = subset(tmp, !is.na(tmp$mgmt))
  ) +
  scale_x_continuous(expand = c(0.08, 0)) +
  scale_y_continuous(expand = c(0.124, 0)) +
  job_gg_theme 



# location and lowgrade is correlated
#plot(tmp$tlc, tmp$chr7)
#plot(tmp$tlc, tmp$idh)




df <- data.frame(
  c10=as.factor(gsam.patient.metadata$primary10),
  c710=as.factor(gsam.patient.metadata$primary710
                )
  )
df <- df[!is.na(df$c10) & !is.na(df$c710),]




