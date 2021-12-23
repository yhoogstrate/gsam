#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(infercnv)
library(AnnotationHub)
library(ensembldb)


# cluster genes ----

C3 <- c('VWF', 'TIE1', 'HIGD1B', 'MMRN1', 'CYSLTR2', 'MMP25','FLT4', 'BCL6B', 'GRAP', 'LAMC3', 'DPEP1', 'PXDNL', 'ANGPT2',
        'PALD1', 'ADGRD1', 'GBP6', 'SLC52A3', 'CLDN5', 'VWA2', 'ABCB1', 'THSD7B', 'SPINK8', 'FOXQ1', 'ZIC3', 'NODAL')

C4A <- c('SOD3', "FSTL3", "FAM180A", "OSGIN1", "NDRG1", "AC010327.1","TRIM29", "HSPB7", "TNNT1", "CCN5", "MICAL2", "GLIS1", "SLIT3",
        "CYP26B1", "NPR3", "FGF5", "CCBE1", "GPR68", "SH3RF2")
C4B <- c("WNT11", "SCUBE3", "KRT17", "GPR78","CPZ","GLI1", "PRB2","MAFA","HAPLN1")

C5 <- c("PRF1", "ARHGAP9", "FCMR","LXN","KCNE3", "NR5A2","FPR2", "CCL13", "MMP7", "CALCR", "LRG1", "SAA2", "PI3", "LIF", "HSPA6")

C6 <- c('CRABP2', 'CLIP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
        'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
        "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
        "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
        "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")




# generate annotation table ----

# obtain gene annotation file for gene location index
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("homo sapiens", "EnsDb"),  ignore.case = TRUE)
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]
annotation.genes <- genes(edb, return.type = "data.frame")
#annotations.sub <- annotation.genes[, c("gene_id","seq_name","gene_seq_start","gene_seq_end")]
annotations.sub <- annotation.genes[, c("symbol","seq_name","gene_seq_start","gene_seq_end")] %>%
  dplyr::mutate(seq_name = paste0('chr',seq_name)) %>%
  dplyr::filter(symbol != "") %>%
  dplyr::filter(!duplicated(symbol)) %>%
  dplyr::filter(seq_name %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23','chrX','chrY'))

write.table(annotations.sub, file = "data/scRNA/GSE103224_Yuan/annotation_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote=F)


# GSE103224 :: Yuan J. et. al ----

if(file.exists("tmp/GSE103224.scRNA.counts.Rds")) {
  GSE103224 <- readRDS("tmp/GSE103224.scRNA.counts.Rds")
} else {
  # diagnosis: Glioblastoma, WHO grade IV - idh1 status: R132H
  #a <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758471_PJ016.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  #`colnames<-`(c("ENS", "HGNC", paste0("GSM2758471.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-ampli
  b <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758472_PJ017.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758472_PJ017.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-non ampli
  c <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758473_PJ018.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758473_PJ018.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-ampli
  d <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758474_PJ025.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758474_PJ025.cell.",(1:(ncol(.)-2))+2) ))
  
  # Anaplastic Astrocytoma, WHO grade III, IDH-wt, EGFR-non ampli
  #e <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758475_PJ030.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  #`colnames<-`(c("ENS", "HGNC", paste0("GSM2758475.cell.",(1:(ncol(.)-2))+2) ))
  
  # GBM RECURRENT - IDH-wt, EGFR ampli
  f <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758476_PJ032.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758476_PJ032.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, recurrent, dh1 status: wt, egfr status: amplified in initial resection
  g <- read.delim("data/scRNA/GSE103224_Yuan/GSM2758477_PJ035.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758477_PJ035.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, IDH-wt, EGFR-non ampli
  h <- read.delim("data/scRNA/GSE103224_Yuan/GSM2940098_PJ048.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
    `colnames<-`(c("ENS", "HGNC", paste0("GSM2940098_PJ048.cell.",(1:(ncol(.)-2))+2) ))
  
  

  
  
  GSE103224 <- b %>%
  dplyr::left_join(c %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(d %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(f %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS')) %>%
  dplyr::left_join(g %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS'))%>%
  dplyr::left_join(h %>% dplyr::mutate(HGNC=NULL), by=c('ENS'='ENS'))
  
  saveRDS(GSE103224, file="tmp/GSE103224.scRNA.counts.Rds")
  
  rm(b,c,d,f,g,h)
}


## convert to seurat/10x ----


GSE103224.genes <- GSE103224 %>%
  dplyr::select(ENS, HGNC)

# make it a numerical-only data.frame
GSE103224 <- GSE103224 %>%
  tibble::column_to_rownames('ENS') %>% 
  dplyr::mutate(HGNC = NULL)

# remove empty features
GSE103224 <- GSE103224 %>%
  dplyr::filter(rowSums(.) > 0)
gc()

# sync gene data frame
GSE103224.genes  <- GSE103224.genes %>% 
  dplyr::filter(ENS %in% rownames(GSE103224))

# double check if both tables are still in-sync
stopifnot(rownames(GSE103224) == GSE103224.genes$ENS)


# store identifiers of the cells
GSE103224.cell.ids <- colnames(GSE103224)



# # convert to dgTMarix
# GSE103224 <- GSE103224 %>%
#   as.matrix() %>%
#   as("dgTMatrix")
# gc()



# export each
for(sample in unique(gsub(".cell.+$","",GSE103224.cell.ids ))) {
  print(paste0("Exporting ",sample))
  
  GSE103224.subset <- GSE103224 %>%
    as.data.frame() %>%
    dplyr::select(contains(sample))
  print(dim(GSE103224.subset))

  GSE103224.cell.ids.subset <- colnames(GSE103224.subset)
  
  write10xCounts(path        = paste0("/home/youri/projects/gsam/data/scRNA/GSE103224_Yuan/", sample),
                 x           = GSE103224.subset %>% as.matrix() %>% as("dgTMatrix"),
                 barcodes    = GSE103224.cell.ids.subset , 
                 gene.id     = GSE103224.genes$ENS,
                 gene.symbol = GSE103224.genes$HGNC)
  
  rm(GSE103224.subset, GSE103224.cell.ids.subset, sample)
  gc()
}


rm(GSE103224, GSE103224.cell.ids, GSE103224.genes)
gc()


## B :: GSM2758472_PJ017 :: T,MG,TC ---- 
# Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-ampli

rm(object_1)
gc()


sid <- 'GSM2758472_PJ017'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224_Yuan/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & 
                     nFeature_RNA <4000 & nCount_RNA >200 &
                     nCount_RNA <10000 & percent.mito < 0.025)
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot2

### scaling of data
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1)

### cluster the cells
d <- 25

object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1)
head(Idents(object_1), 50)


object_1 <- RunTSNE(object_1, dims = 1:d)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^0$","Tumor/AC 0",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^8$","Tumor/AC 8",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^2$","Tumor/AC 2",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^7$","Tumor/AC 7",levels(object_1$seurat_clusters))

levels(object_1$seurat_clusters) <- gsub("^1$","TAM/MG 1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^3$","TAM/MG 3",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^4$","TAM/MG 4",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^5$","TAM/MG 5",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^6$","TAM/MG 6",levels(object_1$seurat_clusters))


levels(object_1$seurat_clusters) <- gsub("^9$","T-cells",levels(object_1$seurat_clusters))

# DimPlot(object_1, reduction = "tsne", label = TRUE, pt.size = .8, group.by = "seurat_clusters")
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters") 


object_1@reductions

dge <- FindMarkers(object_1, ident.1 = '9')
head(dge,15)

rm(dge)

### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 1)

object_1$youri_clusters <- as.character(object_1$seurat_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(1,3,4,5,6),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0,2,8,7),"Tumor", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(9),"T-Cell", object_1$youri_clusters)

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=0, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] < 0 &
  object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=3, linetype="dashed", color = "red") +
  geom_vline(xintercept=10, linetype="dashed", color = "red") +
  geom_hline(yintercept=-9, linetype="dashed", color = "red") +
  geom_hline(yintercept=-2, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
    object_1@reductions$umap@cell.embeddings[,1] >= 3 &
    object_1@reductions$umap@cell.embeddings[,1] <= 10 &
    object_1@reductions$umap@cell.embeddings[,2] >= -9 &
    object_1@reductions$umap@cell.embeddings[,2] <= -2 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

object_1.PJ017 <- object_1



#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "STMN2") # Tumor


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB", "CDK4", "GFAP")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M



# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")


#### 2A. TAM/mg/monocytes (+) ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))

#### 2B. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "CD2")


#### 3. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
# FeaturePlot(object = object_1, features = "GABRB2")
# FeaturePlot(object = object_1, features = "RBFOX1")
# FeaturePlot(object = object_1, features = "DDN")
# FeaturePlot(object = object_1, features = "TNNT2")
# FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")


#### 4. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")


#### 5A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 5B. Pericytes (3 cells ?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


FeaturePlot(object = object_1, features = "FOXJ1") # ependymal cells? - https://www.cell.com/cell/pdf/S0092-8674(18)30395-7.pdf


#### 6. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "TRBC2")

# Deze werken nergens op T-cellen (?):
# FeaturePlot(object = object_1, features = "CACHD1")
# FeaturePlot(object = object_1, features = "BMPR1B")
# FeaturePlot(object = object_1, features = "GPR37L1")

#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(object = object_1, features = C4A[1:4])
FeaturePlot(object = object_1, features = C4A[5:8])
FeaturePlot(object = object_1, features = C4A[9:12])
FeaturePlot(object = object_1, features = C4A[13:16])
FeaturePlot(object = object_1, features = C4A[17:19])

FeaturePlot(object = object_1, features = C4B[1:4])
FeaturePlot(object = object_1, features = C4B[5:9])



#### C5 (down) ----

f <- C5
DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5[1:4])
FeaturePlot(object = object_1, features = C5[5:8])
FeaturePlot(object = object_1, features = C5[9:12])
FeaturePlot(object = object_1, features = C5[13:16])


#### C6 (up) ----


DotPlot(object = object_1, features = C6, group.by = "seurat_clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


VlnPlot(object = object_1, features = C6, stack = T, sort = T)
VlnPlot(object = object_1, features = C6, stack = T, sort = T)

FeaturePlot(object = object_1, features = C6)


FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CLIP2" )
# FeaturePlot(object = object_1, features =  "DPT" )
# FeaturePlot(object = object_1, features =  "FGF7" )
# FeaturePlot(object = object_1, features =  "COL10A1" )
FeaturePlot(object = object_1, features =  "FBN1" )
# FeaturePlot(object = object_1, features =  "GLT8D2" )
FeaturePlot(object = object_1, features =  "IRX3" )
FeaturePlot(object = object_1, features =  "MFAP5" )
FeaturePlot(object = object_1, features =  "MFAP4" ) # tumor + MG
FeaturePlot(object = object_1, features =  "COL8A2" )
# FeaturePlot(object = object_1, features =  "FNDC1" )
FeaturePlot(object = object_1, features =  "MMP11" )
FeaturePlot(object = object_1, features =  "MFAP2" )
FeaturePlot(object = object_1, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL5A1" )
FeaturePlot(object = object_1, features =  "ADAMTS2" )
# FeaturePlot(object = object_1, features =  "TPSB2" )
FeaturePlot(object = object_1, features =  "KRT8" )
# FeaturePlot(object = object_1, features =  "OMD" )
# FeaturePlot(object = object_1, features =  "OGN" )
# FeaturePlot(object = object_1, features =  "MME" )
FeaturePlot(object = object_1, features =  "MLPH" )
# FeaturePlot(object = object_1, features =  "MRC1L1" )
# FeaturePlot(object = object_1, features =  "PTGFR" )
# FeaturePlot(object = object_1, features =  "TWIST2" )
# FeaturePlot(object = object_1, features =  "C5orf46" )
FeaturePlot(object = object_1, features =  "TNNT3" )
FeaturePlot(object = object_1, features =  "ASS1" )
FeaturePlot(object = object_1, features =  "PERP" )
# FeaturePlot(object = object_1, features =  "KLHDC7B" )
FeaturePlot(object = object_1, features =  "CCL8" )



## C :: GSM2758473_PJ018 :: T,MG,OL,EN,PE ----
# Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-non ampli

rm(object_1, sid)
gc()

sid <- 'GSM2758473_PJ018'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224_Yuan/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & 
                     nFeature_RNA <10000 &
                     nCount_RNA >200 &
                     nCount_RNA < 20000 & 
                     percent.mito < 0.02) # aan de hand van de violin plot
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#CombinePlots(plots = list(plot1, plot2))     

### scaling of data
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1)

### cluster the cells

object_1 <- FindNeighbors(object_1, dims = 1:25)
object_1 <- FindClusters(object_1, resolution = 1) # 0.5 - 1 ?
#head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:25)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(7|3|0|6|2|1|4|5)$","Tumor.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^8$","Oligodendrocytes",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^9$","TAM/Microglia",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^10$","Endothelial + Pericytes",levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 1)

object_1$youri_clusters <- as.character(object_1$seurat_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(10) & object_1@reductions$umap@cell.embeddings[,1] < -2.2 ,"Pericytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(10) & object_1@reductions$umap@cell.embeddings[,1] >= -2.2 ,"Endothelial", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(8),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(9),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(00,1,2,3,4,5,6,7),"Tumor", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-2.2, linetype="dashed", color = "red")

object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -4 &
    object_1@reductions$umap@cell.embeddings[,1] <= -2.2 &
    object_1@reductions$umap@cell.embeddings[,2] >= 13 &
    object_1@reductions$umap@cell.embeddings[,2] <= 16 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-4, linetype="dashed", color = "red") +
  geom_vline(xintercept=-2.2, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 13, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 16, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -2.2 &
    object_1@reductions$umap@cell.embeddings[,1] <= -0 &
    object_1@reductions$umap@cell.embeddings[,2] >= 13 &
    object_1@reductions$umap@cell.embeddings[,2] <= 16 &
    object_1$youri_clusters != "Endothelial", "ambiguous", object_1$youri_clusters)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-4, linetype="dashed", color = "red") +
  geom_vline(xintercept=-2.2, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 13, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 16, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 17 &
    object_1@reductions$umap@cell.embeddings[,1] <= 21.5 &
    object_1@reductions$umap@cell.embeddings[,2] >= -5 &
    object_1@reductions$umap@cell.embeddings[,2] <= 0 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=17, linetype="dashed", color = "red") +
  geom_vline(xintercept=21.5, linetype="dashed", color = "red") + 
  geom_hline(yintercept = -5, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -9 &
    object_1@reductions$umap@cell.embeddings[,1] <= -4.5 &
    object_1@reductions$umap@cell.embeddings[,2] >= 3 &
    object_1@reductions$umap@cell.embeddings[,2] <= 5.2 &
    object_1$youri_clusters != "Oligodendrocytes", "ambiguous", object_1$youri_clusters)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-9, linetype="dashed", color = "red") +
  geom_vline(xintercept=-4.5, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 3, linetype="dashed", color = "red") + 
  geom_hline(yintercept = 5.2, linetype="dashed", color = "red")


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

object_1.PJ018 <- object_1



#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "STMN2") # Tumor


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB", "CDK4", "GFAP")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M



# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")


#### 2A. TAM/mg/monocytes (+) ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 2B. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "CD2")

#### 3. Neurons (-) ----

FeaturePlot(object = object_1, features = "RBFOX3")
# FeaturePlot(object = object_1, features = "GABRG2")

# FeaturePlot(object = object_1, features = "RBFOX1")
# FeaturePlot(object = object_1, features = "DDN")
# FeaturePlot(object = object_1, features = "TNNT2")
# FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")


#### 4. Oligodendrocytes (+) ----

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


#### 5A. Endothelial (+) ----

f <- c("ABCB1","CD34","FLT4","TIE1","ITGA1","RGS5","PDGFRB","CD248")
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)



FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 5B. Pericytes (+ heel dicht met endothelial?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


# FeaturePlot(object = object_1, features = "FOXJ1") # ependymal cells? - https://www.cell.com/cell/pdf/S0092-8674(18)30395-7.pdf



#### 6. Til/T-cell (?) ----

# FeaturePlot(object = object_1, features = "CD3D")
# FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "TRBC2")

#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


f <- C5
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


f <- C6
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6)

FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CLIP2" )
FeaturePlot(object = object_1, features =  "DPT" )
FeaturePlot(object = object_1, features =  "FGF7" )
# FeaturePlot(object = object_1, features =  "COL10A1" )
FeaturePlot(object = object_1, features =  "FBN1" )
FeaturePlot(object = object_1, features =  "GLT8D2" )
# FeaturePlot(object = object_1, features =  "IRX3" )
FeaturePlot(object = object_1, features =  "MFAP5" )
FeaturePlot(object = object_1, features =  "MFAP4" ) # tumor + MG
FeaturePlot(object = object_1, features =  "COL8A2" )
# FeaturePlot(object = object_1, features =  "FNDC1" )
FeaturePlot(object = object_1, features =  "MMP11" )
FeaturePlot(object = object_1, features =  "MFAP2" )
FeaturePlot(object = object_1, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL5A1" )
FeaturePlot(object = object_1, features =  "ADAMTS2" )
FeaturePlot(object = object_1, features =  "TPSB2" )
FeaturePlot(object = object_1, features =  "KRT8" )
# FeaturePlot(object = object_1, features =  "OMD" )
# FeaturePlot(object = object_1, features =  "OGN" )
# FeaturePlot(object = object_1, features =  "MME" )
FeaturePlot(object = object_1, features =  "MLPH" )
# FeaturePlot(object = object_1, features =  "MRC1L1" )
# FeaturePlot(object = object_1, features =  "PTGFR" )
# FeaturePlot(object = object_1, features =  "TWIST2" )
# FeaturePlot(object = object_1, features =  "C5orf46" )
FeaturePlot(object = object_1, features =  "TNNT3" )
FeaturePlot(object = object_1, features =  "ASS1" )
FeaturePlot(object = object_1, features =  "PERP" )
# FeaturePlot(object = object_1, features =  "KLHDC7B" )
FeaturePlot(object = object_1, features =  "CCL8" )



## D :: GSM2758474_PJ025 :: T,MG,TC ----
# Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-ampli

rm(object_1, sid)

sid <- 'GSM2758474_PJ025'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224_Yuan/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & 
                     nFeature_RNA <10000 &
                     nCount_RNA >200 &
                     nCount_RNA < 20000 & 
                     percent.mito < 0.1) # aan de hand van de violin plot
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))     

### scaling of data 
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

### cluster the cells

object_1 <- FindNeighbors(object_1, dims = 1:25)
object_1 <- FindClusters(object_1, resolution = 1) # 0.5 - 1 ?
#head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:25)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^10$","Pericytes",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^14","TAM/mg/monocytes",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^12$","Endothelial",levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")



#tmp2 = FindMarkers(object = object_1, ident.1 = 10)

### Prepare for integration ----


object_1 <- FindClusters(object_1, resolution = 1)

object_1$youri_clusters <- as.character(object_1$seurat_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0:9,11,13),"Tumor", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(14),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(12),"Endothelial", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(10),"Pericytes", object_1$youri_clusters)


# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=3.5, linetype="dashed", color = "red") +
  geom_vline(xintercept=7, linetype="dashed", color = "red") +
  geom_hline(yintercept=12, linetype="dashed", color = "red") +
  geom_hline(yintercept=17, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 3.5 &
    object_1@reductions$umap@cell.embeddings[,1] <= 7 &
    object_1@reductions$umap@cell.embeddings[,2] >= 12 &
    object_1@reductions$umap@cell.embeddings[,2] <= 17 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)

# Tam
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=12, linetype="dashed", color = "red") +
  geom_vline(xintercept=16, linetype="dashed", color = "red") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_hline(yintercept=5, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 12 &
    object_1@reductions$umap@cell.embeddings[,1] <= 16 &
    object_1@reductions$umap@cell.embeddings[,2] >= 0 &
    object_1@reductions$umap@cell.embeddings[,2] <= 5 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)

# endo
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=13, linetype="dashed", color = "red") +
  geom_vline(xintercept=17, linetype="dashed", color = "red") +
  geom_hline(yintercept=-8, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 13 &
    object_1@reductions$umap@cell.embeddings[,1] <= 17 &
    object_1@reductions$umap@cell.embeddings[,2] >= -8 &
    object_1@reductions$umap@cell.embeddings[,2] <= -3 &
    object_1$youri_clusters != "Endothelial", "ambiguous", object_1$youri_clusters)


# Tumor
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-5, linetype="dashed", color = "red") +
  geom_vline(xintercept=4, linetype="dashed", color = "red") +
  geom_hline(yintercept=-8, linetype="dashed", color = "red") +
  geom_hline(yintercept=10, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 0 &
    object_1@reductions$umap@cell.embeddings[,1] <= 3 &
    object_1@reductions$umap@cell.embeddings[,2] >= 13 &
    object_1@reductions$umap@cell.embeddings[,2] <= 17 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

sum(object_1$youri_clusters == "ambiguous")



object_1.PJ025 <- object_1




#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "STMN2") # Tumor


# succes met vinden van een marker

FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM", "AURKB", "CDK4", "GFAP")) # Tumor


FeaturePlot(object = object_1, features = c("TMPO")) # G2/M



# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")


#### 2A. TAM/mg/monocytes (+) ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))

#### 2B. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "CD2")


#### 3. Neurons (-) ----

FeaturePlot(object = object_1, features = "RBFOX3")
# FeaturePlot(object = object_1, features = "RBFOX1")
# FeaturePlot(object = object_1, features = "DDN")
# FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")


#### 4. Oligodendrocytes (-) ----

# FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


#### 5A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 5B. Pericytes (+ very strong) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### 6. Til/T-cell (-) ----


FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "TRBC2")

# FeaturePlot(object = object_1, features = "CACHD1")
# FeaturePlot(object = object_1, features = "BMPR1B")
# FeaturePlot(object = object_1, features = "GPR37L1")

#### 7. EGFR ligands etc. ----

FeaturePlot(object = object_1, features = c("SOCS2"))
FeaturePlot(object = object_1, features = c("EGFR", "HBEGF"))
FeaturePlot(object = object_1, features = c("EGFR", "SOCS2"))
FeaturePlot(object = object_1, features = c("EGFR", "MEOX2"))
FeaturePlot(object = object_1, features = c("EGFR", "EGF"))
FeaturePlot(object = object_1, features = c("EGFR", "AREG", "EREG"))
FeaturePlot(object = object_1, features = c("EGFR", "AREG", "EREG"))
FeaturePlot(object = object_1, features = c("EGFR", "TGFA"))
# FeaturePlot(object = object_1, features = c("EGFR", "BTC"))
# FeaturePlot(object = object_1, features = c("EGFR", "EPGN"))
# FeaturePlot(object = object_1, features = c("EGFR", "NGF"))
FeaturePlot(object = object_1, features = c("EGFR", "FGFR3"))



#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


f <- C5
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----

f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CLIP2" )
# FeaturePlot(object = object_1, features =  "DPT" )
FeaturePlot(object = object_1, features =  "FGF7" )
# FeaturePlot(object = object_1, features =  "COL10A1" )
FeaturePlot(object = object_1, features =  "FBN1" )
FeaturePlot(object = object_1, features =  "GLT8D2" )
FeaturePlot(object = object_1, features =  "IRX3" )
# FeaturePlot(object = object_1, features =  "MFAP5" )
FeaturePlot(object = object_1, features =  "MFAP4" ) # tumor + MG
FeaturePlot(object = object_1, features =  "COL8A2" )
# FeaturePlot(object = object_1, features =  "FNDC1" )
FeaturePlot(object = object_1, features =  "MMP11" )
FeaturePlot(object = object_1, features =  "MFAP2" ) # endo?
FeaturePlot(object = object_1, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL5A1" )
# FeaturePlot(object = object_1, features =  "ADAMTS2" )
# FeaturePlot(object = object_1, features =  "TPSB2" )
FeaturePlot(object = object_1, features =  "KRT8" )
# FeaturePlot(object = object_1, features =  "OMD" )
# FeaturePlot(object = object_1, features =  "OGN" )
FeaturePlot(object = object_1, features =  "MME" )
# FeaturePlot(object = object_1, features =  "MLPH" )
# FeaturePlot(object = object_1, features =  "MRC1L1" )
# FeaturePlot(object = object_1, features =  "PTGFR" )
FeaturePlot(object = object_1, features =  "TWIST2" )
# FeaturePlot(object = object_1, features =  "C5orf46" )
# FeaturePlot(object = object_1, features =  "TNNT3" )
FeaturePlot(object = object_1, features =  "ASS1" )
FeaturePlot(object = object_1, features =  "PERP" )
FeaturePlot(object = object_1, features =  "KLHDC7B" )
# FeaturePlot(object = object_1, features =  "CCL8" )






#### C2 ~ neuronal ----

FeaturePlot(object = object_1, features = "PLPP2")


#### OPC netfel ----
# lijkt wel in tumor te zitten(!)
FeaturePlot(object = object_1, features = "PLP1") # wel/ook in tumor [ is dit niet een Neftel OPC gene?]




#### Astrocytes ----
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")
FeaturePlot(object = object_1, features = "GFAP")




## F :: GSM2758476_PJ032 :: T,MG ----
# age: 63
# gender: Female
# location: left temporal
# diagnosis: Glioblastoma, recurrent
# idh1 status: wt
# egfr status: amplified in recurrent tumor

rm(object_1, sid)
gc()

sid <- "GSM2758476_PJ032"
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224_Yuan/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1,
                   subset = nFeature_RNA > 700 &
                     nFeature_RNA < 10000 &
                     nCount_RNA > 200 &
                     nCount_RNA < 9000 &
                     percent.mito <0.015)
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))     

### scaling of data
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1)

### cluster the cells

object_1 <- FindNeighbors(object_1, dims = 1:25)
object_1 <- FindClusters(object_1, resolution = 1)
head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:25)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(0|2|4)$","Tumor.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(1|3|5)","TAM/MG.\\1",levels(object_1$seurat_clusters))



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")



### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 1)

object_1$youri_clusters <- as.character(object_1$seurat_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(1,3,5),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0,2,4),"Tumor", object_1$youri_clusters)


object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 1 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=1, linetype="dashed", color = "red")

object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] < 1 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=1, linetype="dashed", color = "red")


object_1.PJ032 <- object_1



#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "STMN2") # Tumor


# succes met vinden van een marker

FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB", "CDK4", "GFAP")) # Tumor

# 
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M



# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")


#### 2A. TAM/mg/monocytes (+) ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))

#### 2B. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")

# FeaturePlot(object = object_1, features = "CD2")


# FeaturePlot(object = object_1, features = "CACHD1")
# FeaturePlot(object = object_1, features = "BMPR1B")
# FeaturePlot(object = object_1, features = "GPR37L1")




#### 3. Neurons (-) ----

# FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
# FeaturePlot(object = object_1, features = "DDN")
# FeaturePlot(object = object_1, features = "TNNT2")
# FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")


#### 4. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


#### 5A. Endothelial (-) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 5B. Pericytes (?small cluster?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


f <- C5
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


f <- C6
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CLIP2" )
FeaturePlot(object = object_1, features =  "DPT" )
FeaturePlot(object = object_1, features =  "FGF7" )
# FeaturePlot(object = object_1, features =  "COL10A1" )
FeaturePlot(object = object_1, features =  "FBN1" )
FeaturePlot(object = object_1, features =  "GLT8D2" )
FeaturePlot(object = object_1, features =  "IRX3" )
# FeaturePlot(object = object_1, features =  "MFAP5" )
FeaturePlot(object = object_1, features =  "MFAP4" ) # tumor + MG
FeaturePlot(object = object_1, features =  "COL8A2" )
# FeaturePlot(object = object_1, features =  "FNDC1" )
FeaturePlot(object = object_1, features =  "MMP11" )
FeaturePlot(object = object_1, features =  "MFAP2" )
FeaturePlot(object = object_1, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL5A1" )
FeaturePlot(object = object_1, features =  "ADAMTS2" )
FeaturePlot(object = object_1, features =  "TPSB2" )
FeaturePlot(object = object_1, features =  "KRT8" )
# FeaturePlot(object = object_1, features =  "OMD" )
# FeaturePlot(object = object_1, features =  "OGN" )
# FeaturePlot(object = object_1, features =  "MME" )
# FeaturePlot(object = object_1, features =  "MLPH" )
# FeaturePlot(object = object_1, features =  "MRC1L1" )
# FeaturePlot(object = object_1, features =  "PTGFR" )
FeaturePlot(object = object_1, features =  "TWIST2" )
# FeaturePlot(object = object_1, features =  "C5orf46" )
# FeaturePlot(object = object_1, features =  "TNNT3" )
FeaturePlot(object = object_1, features =  "ASS1" )
FeaturePlot(object = object_1, features =  "PERP" )
# FeaturePlot(object = object_1, features =  "KLHDC7B" )
FeaturePlot(object = object_1, features =  "CCL8" )






## G :: GSM2758477_PJ035 :: T,MG,EN,PE ----
# Glioblastoma, IDH-wt, recurrent, egfr status: amplified in initial resection


rm(object_1, sid)
gc()


sid <- 'GSM2758477_PJ035'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224_Yuan/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & 
                     nFeature_RNA < 2400 & 
                     nCount_RNA > 750 & 
                     nCount_RNA < 3750 &
                     percent.mito < 0.015)
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))     


### scaling of data
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1)

### cluster the cells


object_1 <- FindNeighbors(object_1, dims = 1:15)
object_1 <- FindClusters(object_1, resolution = .5)
head(Idents(object_1), 50)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:15)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(0|1|4|2|5)$","Tumor.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^3$","TAM/Microglia",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^6$","Endothelial",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^7$","Pericytes",levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")



### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 0.5)

object_1$youri_clusters <- as.character(object_1$seurat_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0,1,2,4,5),"Tumor", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(3),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(6),"Endothelial", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(7),"Pericytes", object_1$youri_clusters)


# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-11, linetype="dashed", color = "red") +
  geom_vline(xintercept=-8, linetype="dashed", color = "red") +
  geom_hline(yintercept=-15.5, linetype="dashed", color = "red") +
  geom_hline(yintercept=-13, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -11 &
    object_1@reductions$umap@cell.embeddings[,1] <= -8 &
    object_1@reductions$umap@cell.embeddings[,2] >= -15.5 &
    object_1@reductions$umap@cell.embeddings[,2] <= -13 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)

# Tam
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-15, linetype="dashed", color = "red") +
  geom_vline(xintercept=-11, linetype="dashed", color = "red") +
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  geom_hline(yintercept=8, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -15 &
    object_1@reductions$umap@cell.embeddings[,1] <= -11 &
    object_1@reductions$umap@cell.embeddings[,2] >= 2 &
    object_1@reductions$umap@cell.embeddings[,2] <= 8 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)

# endo
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  geom_vline(xintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=13, linetype="dashed", color = "red") +
  geom_hline(yintercept=17, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 0 &
  object_1@reductions$umap@cell.embeddings[,1] <= 3 &
  object_1@reductions$umap@cell.embeddings[,2] >= 13 &
  object_1@reductions$umap@cell.embeddings[,2] <= 17 &
  object_1$youri_clusters != "Endothelial", "ambiguous", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

object_1.PJ035 <- object_1





#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "STMN2") # Tumor


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB", "CDK4", "GFAP")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M



# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")


#### 2A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 2B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "TRBC2")

# FeaturePlot(object = object_1, features = "CACHD1")
# FeaturePlot(object = object_1, features = "BMPR1B")
# FeaturePlot(object = object_1, features = "GPR37L1")


#### 3. Neurons (-) ----

FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")


#### 4. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


#### 5A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 5B. Pericytes (++) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### 6. other cluster(s) ----

FeaturePlot(object = object_1, features = "COL1A1")
FeaturePlot(object = object_1, features = "VEGFA")
FeaturePlot(object = object_1, features = "NDRG1")


#### C3 (down) :: endothelial ----

DotPlot(object = object_1, features = c(C3), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C3)


#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


f <- C5
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


f <- C6
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CLIP2" )
FeaturePlot(object = object_1, features =  "DPT" )
FeaturePlot(object = object_1, features =  "FGF7" )
# FeaturePlot(object = object_1, features =  "COL10A1" )
FeaturePlot(object = object_1, features =  "FBN1" )
FeaturePlot(object = object_1, features =  "GLT8D2" )
FeaturePlot(object = object_1, features =  "IRX3" )
# FeaturePlot(object = object_1, features =  "MFAP5" )
FeaturePlot(object = object_1, features =  "MFAP4" ) # tumor + MG
FeaturePlot(object = object_1, features =  "COL8A2" )
# FeaturePlot(object = object_1, features =  "FNDC1" )
FeaturePlot(object = object_1, features =  "MMP11" )
FeaturePlot(object = object_1, features =  "MFAP2" )
FeaturePlot(object = object_1, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL5A1" )
FeaturePlot(object = object_1, features =  "ADAMTS2" )
FeaturePlot(object = object_1, features =  "TPSB2" )
FeaturePlot(object = object_1, features =  "KRT8" )
# FeaturePlot(object = object_1, features =  "OMD" )
# FeaturePlot(object = object_1, features =  "OGN" )
# FeaturePlot(object = object_1, features =  "MME" )
# FeaturePlot(object = object_1, features =  "MLPH" )
# FeaturePlot(object = object_1, features =  "MRC1L1" )
# FeaturePlot(object = object_1, features =  "PTGFR" )
FeaturePlot(object = object_1, features =  "TWIST2" )
# FeaturePlot(object = object_1, features =  "C5orf46" )
# FeaturePlot(object = object_1, features =  "TNNT3" )
FeaturePlot(object = object_1, features =  "ASS1" )
FeaturePlot(object = object_1, features =  "PERP" )
# FeaturePlot(object = object_1, features =  "KLHDC7B" )
FeaturePlot(object = object_1, features =  "CCL8" )


### Find DE genes voor 5 ----

dge <- FindMarkers(object_1, ident.1 = '5')
head(dge,15)

C4


### InferCNV plot ----


# write.table(object_1@assays$RNA@counts %>% as.data.frame(stringsAsFactors=F),
#   file = paste0("data/scRNA/GSE103224_Yuan/",sid,'/matrix.only-included-cells.txt'),
#   row.names = T, col.names = T, sep = "\t",quote = FALSE )


infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=paste0("data/scRNA/GSE103224_Yuan/",sid,'/matrix.only-included-cells.txt'),
                                     annotations_file=paste0("data/scRNA/GSE103224_Yuan/",sid,'/seurat_cluster_annotations.txt'),
                                     gene_order_file='data/scRNA/GSE103224_Yuan/annotation_genes.txt',
                                     #annotations_file=system.file("extdata", "data/scRNA/GSE103224_Yuan/GSM2758477_PJ035/seurat_cluster_annotations.txt", package = "infercnv"),
                                     #gene_order_file=system.file("extdata", "data/scRNA/GSE103224_Yuan/annotation_genes.txt", package = "infercnv"),
                                     ref_group_names=c("SeuratCluster7","SeuratCluster8","SeuratCluster4")) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='infercnv', 
                             cluster_by_groups=TRUE,
                             denoise=F,
                             HMM=F
)


## H :: GSM2940098_PJ048 :: T,OD,EN,PE ----
# Glioblastoma, WHO grade IV, IDH-wt, EGFR-non ampli


rm(object_1, sid)
gc()

sid <- 'GSM2940098_PJ048'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224_Yuan/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
# object_1 <- subset(x = object_1,
#                    subset = nFeature_RNA >= 500 &
#                      nFeature_RNA < 3500 &
#                      nCount_RNA >= 700 &
#                      nCount_RNA < 6000 &
#                      percent.mito <0.015)

object_1 <- subset(x = object_1, subset = nFeature_RNA > 500 & nFeature_RNA <35000 & nCount_RNA >200 & nCount_RNA <6000 & percent.mito <0.1)

object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 


top10 <- head(VariableFeatures(object_1), 10)


# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2

### scaling of data
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ---

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1)

### cluster the cells

object_1 <- FindNeighbors(object_1, dims = 1:25)
object_1 <- FindClusters(object_1, resolution = 0.4, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:25)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")

levels(object_1$seurat_clusters) <- gsub("^0$","Tumor/CL~EGFR",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^4$","Tumor/Neuronal",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^1$","Tumor/OPC",levels(object_1$seurat_clusters))

levels(object_1$seurat_clusters) <- gsub("^7$","Endothelial",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^8$","Pericytes",levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


### Prepare for integration ----


object_1 <- FindClusters(object_1, resolution = 0.4)


object_1$youri_clusters <- as.character(object_1$seurat_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_hline(yintercept=-9, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0:6) & object_1@reductions$umap@cell.embeddings[,2] > -9,"Tumor", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0:6) & object_1@reductions$umap@cell.embeddings[,2] <= -9,"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(7),"Endothelial", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(8),"Pericytes", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")



# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-12, linetype="dashed", color = "red") +
  geom_vline(xintercept=-8, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -12 &
    object_1@reductions$umap@cell.embeddings[,1] <= -8 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)




# Endo
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-21, linetype="dashed", color = "red") +
  geom_vline(xintercept=-17, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -21 &
    object_1@reductions$umap@cell.embeddings[,1] <= -17 &
    object_1$youri_clusters != "Endothelial", "ambiguous", object_1$youri_clusters)



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

sum(object_1$youri_clusters == "ambiguous")



object_1.PJ048 <- object_1





#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = c("GFAP","SOX4", "ANXA2","AKAP12","TNR"    , "TMPO" )) # Tumor/AC

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "FIBIN") # Tumor/OPC
FeaturePlot(object = object_1, features = "GPR17") # Tumor/OPC
FeaturePlot(object = object_1, features = "ANXA2") # Tumor/MES
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "ANXA1") # Tumor/MES
FeaturePlot(object = object_1, features = "NEU4") # Tumor/NPC
FeaturePlot(object = object_1, features = "SOX4") # Tumor/NPC
FeaturePlot(object = object_1, features = "SOX11") # Tumor/NPC
FeaturePlot(object = object_1, features = "ETV1") # Tumor/NPC
FeaturePlot(object = object_1, features = "STMN2") # Tumor

# AC
FeaturePlot(object = object_1, features = c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1"))
FeaturePlot(object = object_1, features = c("GFAP","FABP7","BCAN","PON2","METTL7B","SPARC","GATM","RAMP1")) # Tumor
FeaturePlot(object = object_1, features = c("PMP2","AQP4","DBI","EDNRB","PTPRZ1","CLU","PMP22","ATP1A2")) # Tumor
FeaturePlot(object = object_1, features = c("S100A16","HEY1","PCDHGC3","TTYH1","NDRG2","PRCP","ATP1B2","AGT","PLTP","GPM6B"))
FeaturePlot(object = object_1, features = c("F3","RAB31","PPAP2B","ANXA5","TSPAN7")) # Tumor

# MES1
FeaturePlot(object = object_1, features = c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2")) # Tumor

# MES2
FeaturePlot(object = object_1, features = c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3","NAMPT")) # Tumor

# OPC
FeaturePlot(object = object_1, features = c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4")) # Tumor

# NPC2
FeaturePlot(object = object_1, features = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4")) # Tumor



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M



# GFAP en ANXA1 astro markers?

FeaturePlot(object = object_1, features = "ANXA1") # Tumor / Neuronal? weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "ANXA2") # Tumor

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTPZ1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "OLIG1")


#### 2A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))

#### 2B. Til/T-cell (-) ----

# FeaturePlot(object = object_1, features = "CD3D")
# FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "TRBC2")


#### 3. Neurons (-) ----

FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 4. Oligodendrocytes (+) ----

FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")


#### 5A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 5B. Pericytes (+) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")

#### C3 - C6 ----


DotPlot(object = object_1, features = c(C3,C4A,C4B,C5,C6), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### C3 (down) :: endothelial ----

DotPlot(object = object_1, features = c(C3), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C3)

#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


f <- C5
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


f <- C6
DotPlot(object_1, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6)

FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CLIP2" )
#FeaturePlot(object = object_1, features =  "DPT" )
FeaturePlot(object = object_1, features =  "FGF7" )
#FeaturePlot(object = object_1, features =  "COL10A1" )
FeaturePlot(object = object_1, features =  "FBN1" )
FeaturePlot(object = object_1, features =  "GLT8D2" )
FeaturePlot(object = object_1, features =  "IRX3" )
FeaturePlot(object = object_1, features =  "MFAP5" )
FeaturePlot(object = object_1, features =  "MFAP4" )
# FeaturePlot(object = object_1, features =  "COL8A2" )
# FeaturePlot(object = object_1, features =  "FNDC1" )
FeaturePlot(object = object_1, features =  "MMP11" )
FeaturePlot(object = object_1, features =  "MFAP2" )
FeaturePlot(object = object_1, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_1, features =  "COL5A1" )
# FeaturePlot(object = object_1, features =  "ADAMTS2" )
# FeaturePlot(object = object_1, features =  "TPSB2" )
# FeaturePlot(object = object_1, features =  "KRT8" )
# FeaturePlot(object = object_1, features =  "OMD" )
FeaturePlot(object = object_1, features =  "OGN" )
FeaturePlot(object = object_1, features =  "MME" )
# FeaturePlot(object = object_1, features =  "MLPH" )
# FeaturePlot(object = object_1, features =  "MRC1L1" )
# FeaturePlot(object = object_1, features =  "PTGFR" )
# FeaturePlot(object = object_1, features =  "TWIST2" )
# FeaturePlot(object = object_1, features =  "C5orf46" )
# FeaturePlot(object = object_1, features =  "TNNT3" )
FeaturePlot(object = object_1, features =  "ASS1" )
FeaturePlot(object = object_1, features =  "PERP" )
# FeaturePlot(object = object_1, features =  "KLHDC7B" )
#FeaturePlot(object = object_1, features =  "CCL8" )




## Combined ----
# All sampels can be included
# some have few cell types but plenty TAM/MG

#dim(object_1.PJ017)
#dim(object_1.PJ025)
object_c <- merge(object_1.PJ017, object_1.PJ018)
object_c <- merge(object_c, object_1.PJ025)
object_c <- merge(object_c, object_1.PJ032)
object_c <- merge(object_c, object_1.PJ035)
object_c <- merge(object_c, object_1.PJ048)
object_c$seurat_clusters <- NULL
object_c <- NormalizeData(object = object_c, normalization.method = "LogNormalize", scale.factor = 1e4)
object_c <- FindVariableFeatures(object = object_c, selection.method = "vst", nfeatures = 2000)
object_c[["state"]] <- as.factor(gsub(".cell.+$","",colnames(object_c)) )
# top10 <- head(VariableFeatures(object_c), 10)
plot1 <- VariableFeaturePlot(object_c)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2

object_c <- ScaleData(object_c, features = rownames(object_c))
object_c <- RunPCA(object_c, features = VariableFeatures(object = object_c))
VizDimLoadings(object_c, dims = 1:2, reduction = "pca")
DimPlot(object_c, reduction = "pca")
ElbowPlot(object_c, ndims=45)

d <- 25


object_c <- FindNeighbors(object_c, dims = 1:d)
object_c <- FindClusters(object_c, resolution = 1.2, algorithm=1)
object_c <- RunUMAP(object_c, dims = 1:d)

object_c@meta.data$pt = sapply(strsplit(rownames(object_c@meta.data), "[.]"), "[", 1)

levels(object_c$seurat_clusters) <- gsub("^(18)$","Pericytes.\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(20)$","Endothelial.\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(10|16|24)$","TAM/MG.\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(23)$","Oligodendrocytes.\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(26)$","T-Cells.\\1",levels(object_c$seurat_clusters))


levels(object_c$seurat_clusters) <- gsub("^(15)$","Tumor [PJ017].\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(2|7|17|25)$","Tumor [PJ018].\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(0|1|11|12|9|22)$","Tumor [PJ025].\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(14)$","Tumor [PJ032].\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(4|8|19|21)$","Tumor [PJ035].\\1",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(3|5|6|13)$","Tumor [PJ048].\\1",levels(object_c$seurat_clusters))


DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")




DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "state")



#### C4 (up) ----


f <- c(C4A,C4B)
DotPlot(object_c, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


DotPlot(object_c, features=f, group.by = "youri_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(object = object_c, features = C4A)
FeaturePlot(object = object_c, features = C4B)


#### C5 (down) ----


f <- C5
DotPlot(object_c, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



DotPlot(object_c, features=f, group.by = "youri_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


f <- C6

f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
DotPlot(object_c, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


VlnPlot(object = object_c, features = c(C6), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_c, features = c(C6), group.by = "seurat_clusters",stack=T)

VlnPlot(object = object_c, features = c(C6), group.by = "youri_clusters",stack=T)
RidgePlot(object = object_c, features = c(C6), group.by = "youri_clusters",stack=T)


FeaturePlot(object = object_c, features = C6)

FeaturePlot(object = object_c, features =  "CRABP2" )
FeaturePlot(object = object_c, features =  "CLIP2" )
FeaturePlot(object = object_c, features =  "DPT" )
FeaturePlot(object = object_c, features =  "FGF7" )
#FeaturePlot(object = object_c, features =  "COL10A1" )
FeaturePlot(object = object_c, features =  "FBN1" ) # J
FeaturePlot(object = object_c, features =  "GLT8D2" )
FeaturePlot(object = object_c, features =  "IRX3" )
FeaturePlot(object = object_c, features =  "MFAP5" )
FeaturePlot(object = object_c, features =  "MFAP4" ) # tumor
# FeaturePlot(object = object_c, features =  "COL8A2" )
# FeaturePlot(object = object_c, features =  "FNDC1" )
FeaturePlot(object = object_c, features =  "MMP11" ) # ja
FeaturePlot(object = object_c, features =  "MFAP2" ) # endo
FeaturePlot(object = object_c, features =  "COL1A2" ) # high in pericytes
FeaturePlot(object = object_c, features =  "COL1A1" ) # high in pericytes
FeaturePlot(object = object_c, features =  "COL5A1" )
FeaturePlot(object = object_c, features =  "ADAMTS2" )
FeaturePlot(object = object_c, features =  "TPSB2" )
FeaturePlot(object = object_c, features =  "KRT8" )
# FeaturePlot(object = object_c, features =  "OMD" )
FeaturePlot(object = object_c, features =  "OGN" )
FeaturePlot(object = object_c, features =  "MME" )


