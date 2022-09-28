#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
#library(infercnv)
#library(AnnotationHub)
#library(ensembldb)


# cluster genes ----

C3 <- c('VWF', 'TIE1', 'HIGD1B', 'MMRN1', 'CYSLTR2', 'MMP25','FLT4', 'BCL6B', 'GRAP', 'LAMC3', 'DPEP1', 'PXDNL', 'ANGPT2',
        'PALD1', 'ADGRD1', 'GBP6', 'SLC52A3', 'CLDN5', 'VWA2', 'ABCB1', 'THSD7B', 'SPINK8', 'FOXQ1', 'ZIC3', 'NODAL')

C4A <- c('SOD3', "FSTL3", "FAM180A", "OSGIN1", "NDRG1", "AC010327.1","TRIM29", "HSPB7", "TNNT1", "CCN5", "MICAL2", "GLIS1", "SLIT3",
        "CYP26B1", "NPR3", "FGF5", "CCBE1", "GPR68", "SH3RF2")
C4B <- c("WNT11", "SCUBE3", "KRT17", "GPR78","CPZ","GLI1", "PRB2","MAFA","HAPLN1")

C5 <- c("PRF1", "ARHGAP9", "FCMR","LXN","KCNE3", "NR5A2","FPR2", "CCL13", "MMP7", "CALCR", "LRG1", "SAA2", "PI3", "LIF", "HSPA6")

C6 <- c('CRABP2', 'CILP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
        'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
        "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
        "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
        "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")






# EGAS00001004422 :: Couturier  ----
# https://www.nature.com/articles/s41467-020-17186-5
# https://www.frontiersin.org/articles/10.3389/fonc.2021.683007/full



## BT322 [95-100% tumor] -----

rm(object_1)
gc()

sid <- "BT322.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1900,col="red") +
  geom_hline(yintercept = 8000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 55000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1900 &
                     nFeature_RNA < 8000 &
                     nCount_RNA > 500 &
                     nCount_RNA < 55000 &
                     percent.mito < 0.175)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 55)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("Oligodendrocyte"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


tmp.14 <- FindMarkers(object_1, ident.1 = 14)


#### 5. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = c("TMEM125","TMEM144"))


## BT324-GSC :: T,MG,OD ----

rm(object_1, sid)
gc()

sid <- "BT324-GSC.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1250,col="red") +
  geom_hline(yintercept = 6750,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 35000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1250 &
                     nFeature_RNA < 6750 &
                     nCount_RNA > 500 &
                     nCount_RNA < 35000 &
                     percent.mito < 0.15)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(4|9)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(15)$",paste0("Mono/Leukocyte?"),levels(object_1$seurat_clusters)) # zelfde cluster als in 363-GSC

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.13.14 <- FindMarkers(object_1, ident.1 = c(13,14)) 

#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


#### 2. Astrocyte (???) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")


#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor
FeaturePlot(object = object_1, features = "HBG2") # Tumor


#### 3D. ? Mono/Leukocyte ?? ----


# These are cluster-13 DE genes, of which some at genecards seem related to leukocytes?
FeaturePlot(object = object_1, features = c("LAMP3","IRF4","NCCRP1","CRIP1","SYNPO2","CCR7","EHF","CCL22","VTN","LSP1","CDX2"))


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")


#### 5. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = "TMEM144")


#### 6A. Endothelial (-) ----

FeaturePlot(object = object_1, features = "CD34")


#### 6B. Pericytes (?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "CD248")


#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)



FeaturePlot(object = object_1, features = C4A[1:4])
FeaturePlot(object = object_1, features = C4A[5:8])
FeaturePlot(object = object_1, features = C4A[9:12])
FeaturePlot(object = object_1, features = C4A[13:16])
FeaturePlot(object = object_1, features = C4A[17:19])

FeaturePlot(object = object_1, features = C4B[1:4])
FeaturePlot(object = object_1, features = C4B[5:9])



#### C5 (down) ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


FeaturePlot(object = object_1, features = C5[1:4])
FeaturePlot(object = object_1, features = C5[5:8])
FeaturePlot(object = object_1, features = C5[9:12])
FeaturePlot(object = object_1, features = C5[13:16])


#### C6 (up) ----


DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)




VlnPlot(object = object_1, features = C6, stack = T, sort = T)
VlnPlot(object = object_1, features = C6, stack = T, sort = T)

FeaturePlot(object = object_1, features = C6)


## BT326-GSC [95-100% tumor] ----

rm(object_1, sid)
gc()

sid <- "BT326-GSC.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 2350,col="red") +
  geom_hline(yintercept = 8250,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 80000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 2350 &
                     nFeature_RNA < 8250 &
                     nCount_RNA > 500 &
                     nCount_RNA < 80000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 50)

object_1 <- FindNeighbors(object_1, dims = 1:35)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:35)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


#### 1. Tumor (+) ----


FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


#### 2. Astrocyte (???) ----

FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor


#### 5. Oligodendrocytes (+) ----

FeaturePlot(object = object_1, features = "TMEM144")


## BT333-GSC [95-100% tumor?] ----


rm(object_1)
gc()

sid <- 'BT333-GSC.filtered_gene_matrices'
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1950,col="red") +
  geom_hline(yintercept = 7000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 40000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1950 &
                     nFeature_RNA < 7000 &
                     nCount_RNA > 1000 &
                     nCount_RNA < 40000 &
                     percent.mito < 0.15)

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


# scaling of data
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

# cluster the cells


object_1 <- FindNeighbors(object_1, dims = 1:35)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:35)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|6|7|8|9|10|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))

levels(object_1$seurat_clusters) <- gsub("^(11)$",paste0("TAM/Microglia"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2","ETNPPL",   "AURKB")) # Tumor


#### 3A. TAM/mg/monocytes (-)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))

#### 5. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = "TMEM144")


## BT333 Total [95-100% tumor?] ----


rm(object_1)
gc()

sid <- 'BT333.filtered_gene_matrices'
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
#object_1 <- Read10X(data.dir = "data/scRNA/EGAS00001004422_Couturier/filtered/BT333-GSC.filtered_gene_matrices/")
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 6000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 35000,col="red") # + scale_y_log10()



# object_1 <- subset(x = object_1, subset = 
#                      nFeature_RNA > 1000 &
#                      nFeature_RNA < 6000 & 
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 35000 &
#                      percent.mito < 0.15)

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


# scaling of data
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

# cluster the cells


object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


#levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|9|12|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


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



#### 2. Astrocyte (???) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor



#### 3A. TAM/mg/monocytes (-)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")


#### 5. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = "TMEM144")


#### 6A. Endothelial (-) ----

FeaturePlot(object = object_1, features = "CD34")


#### 6B. Pericytes (?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "CD248")

#### C4 (up) ----

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)

DotPlot(object = object_1, features = c(C4A, C4B))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)


#### C5 (down) ----


FeaturePlot(object = object_1, features = C5)


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


#### C6 (up) ----

FeaturePlot(object = object_1, features = C6)


DotPlot(object = object_1, features = c(C6), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)



## BT338 [1+2/2] :: T+,MG-,TC-,OD,PE++ ----


rm(object_1,sid)
gc()

sid <- "BT338_1of2.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

# RenameCells(object_1, new.names = paste0("Z_",colnames(object_1)))
# RenameCells(object_1,add.cell.id = "Zz")
# head(x = colnames(x = object_1))

sid <- "BT338_2of2.filtered_gene_matrices"
object_1.tmp <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1.tmp <- CreateSeuratObject(counts = object_1.tmp, min.cells = 3, min.features = 200, project="Couturier")

object_1.m <- merge(object_1, y=object_1.tmp, add.cell.ids = c("1of2","2of2"), project="Couturier")

rm(object_1, object_1.tmp)
object_1 <- object_1.m
rm(object_1.m)
gc()



mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1250,col="red") +
  geom_hline(yintercept = 7000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 50000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1250 &
                     nFeature_RNA < 7000 &
                     nCount_RNA > 500 &
                     nCount_RNA < 50000 &
                     percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- as.factor(paste(gsub("_.+$","",sid),gsub("_.+$","",colnames(object_1)),sep="-"))
object_1$dataset <- as.character(object_1$state)


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:35)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----


object_1 <- RunUMAP(object_1, dims = 1:35)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(1,4,11),"PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(8),"OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(13),"TAM", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(6)  & object_1@reductions$umap@cell.embeddings[,2] < -14,"TC", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0,2,3,5,6,7,9,10,12)  & object_1@reductions$umap@cell.embeddings[,2] >= -14,"T", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","2. T","3. T","5. T","6. T","7. T","9. T","10. T","12. T",
  "8. OD",
  "1. PE", "11. PE","4. PE",
  "13. TAM",
  "6. TC"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.png"),width=10,height=8)




median(object_1[,object_1$seurat_clusters == "1. PE"]$nCount_RNA) # -> strongest C6 cluster
median(object_1[,object_1$seurat_clusters == "1. PE"]$nFeature_RNA)
mean(object_1[,object_1$seurat_clusters == "1. PE"]$nCount_RNA)
mean(object_1[,object_1$seurat_clusters == "1. PE"]$nFeature_RNA)

median(object_1[,object_1$seurat_clusters == "4. PE"]$nCount_RNA)
median(object_1[,object_1$seurat_clusters == "4. PE"]$nFeature_RNA)
mean(object_1[,object_1$seurat_clusters == "4. PE"]$nCount_RNA)
mean(object_1[,object_1$seurat_clusters == "4. PE"]$nFeature_RNA)

median(object_1[,object_1$seurat_clusters == "11. PE"]$nCount_RNA)
median(object_1[,object_1$seurat_clusters == "11. PE"]$nFeature_RNA)
mean(object_1[,object_1$seurat_clusters == "11. PE"]$nCount_RNA)
mean(object_1[,object_1$seurat_clusters == "11. PE"]$nFeature_RNA)


### Prepare for integration ----


object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_hline(yintercept=-14, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(1,4,11),"Pericytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(8),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(13),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(6)  & object_1@reductions$umap@cell.embeddings[,2] < -14,"T-Cells", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0,2,3,5,6,7,9,10,12)  & object_1@reductions$umap@cell.embeddings[,2] >= -14,"Tumor", object_1$youri_clusters)

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")


# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=10, linetype="dashed", color = "red") +
  geom_vline(xintercept=15, linetype="dashed", color = "red") +
  geom_hline(yintercept=-4, linetype="dashed", color = "red") +
  geom_hline(yintercept=5, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 10 &
    object_1@reductions$umap@cell.embeddings[,1] <= 15 &
    object_1@reductions$umap@cell.embeddings[,2] >= -4 &
    object_1@reductions$umap@cell.embeddings[,2] <= 5 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)


# T-cells
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-2.05, linetype="dashed", color = "red") +
  geom_vline(xintercept=-1, linetype="dashed", color = "red") +
  geom_hline(yintercept=-20, linetype="dashed", color = "red") +
  geom_hline(yintercept=-14, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -2.05 &
    object_1@reductions$umap@cell.embeddings[,1] <= -1 &
    object_1@reductions$umap@cell.embeddings[,2] >= -20 &
    object_1@reductions$umap@cell.embeddings[,2] <= -14 &
    object_1$youri_clusters != "T-Cells", "ambiguous", object_1$youri_clusters)


# TAM/MG
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_vline(xintercept=-2.05, linetype="dashed", color = "red") +
  geom_hline(yintercept=-20, linetype="dashed", color = "red") +
  geom_hline(yintercept=-14, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -3 &
    object_1@reductions$umap@cell.embeddings[,1] <= -2.05 &
    object_1@reductions$umap@cell.embeddings[,2] >= -20 &
    object_1@reductions$umap@cell.embeddings[,2] <= -14 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)

# OD
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-11, linetype="dashed", color = "red") +
  geom_vline(xintercept=-7, linetype="dashed", color = "red") +
  geom_hline(yintercept=13, linetype="dashed", color = "red") +
  geom_hline(yintercept=18, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -11 &
    object_1@reductions$umap@cell.embeddings[,1] <= -7 &
    object_1@reductions$umap@cell.embeddings[,2] >= 13 &
    object_1@reductions$umap@cell.embeddings[,2] <= 18 &
    object_1$youri_clusters != "Oligodendrocyte", "ambiguous", object_1$youri_clusters)


# Tumor
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-10, linetype="dashed", color = "red") +
  geom_vline(xintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-10, linetype="dashed", color = "red") +
  geom_hline(yintercept=9, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -10 &
    object_1@reductions$umap@cell.embeddings[,1] <= 3 &
    object_1@reductions$umap@cell.embeddings[,2] >= -10 &
    object_1@reductions$umap@cell.embeddings[,2] <= 9 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)




DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")
sum(object_1$youri_clusters == "ambiguous")
object_1.BT338 <- object_1




#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = "OLIG2") # Tumor/OPC+NPC1

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


#### 2. Astrocyte (-) ----




FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor


#### 3D. ? Leukocyte ?? ----


# These are cluster-13 DE genes, of which some at genecards seem related to leukocytes?
FeaturePlot(object = object_1, features = c("LAMP3","IRF4","NCCRP1","CRIP1","SYNPO2","CCR7","EHF","CCL22","VTN","LSP1","CDX2"))


#### 4. Neurons (-) ----

DotPlot(object = object_1, features = c("EGFR", "GFAP","MOG", "PLP1", "TMEM144", 
                                        "RBFOX1", "RBFOX2", "RBFOX3", "CD2",
                                        "CD3D", "P2RY12", "CD163", "ABCB1", "RGS5"
))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


FeaturePlot(object = object_1, features = c("SOX4", "RBFOX3"))


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 5. Oligodendrocytes + OPC (+) ----

FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")

DotPlot(object = object_1, features = c("MOG","PLP1","TMEM144"),group.by = "seurat_clusters")



##### Figure S7b ----


tmp.c3 <- results.out |>
  dplyr::filter(!is.na(.data$C3.2022)) |> 
  dplyr::filter(.data$C3.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.opc <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.OPC)) |> 
  dplyr::filter(.data$neftel.meta.modules.OPC == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.c3.opc <- intersect(tmp.c3, tmp.opc)
tmp.c3 <- setdiff(tmp.c3, tmp.c3.opc)
tmp.opc <- setdiff(tmp.opc, tmp.c3.opc)


sid_print <- sid |> 
  stringr::str_replace(".filtered_gene_matrices","") |> 
  stringr::str_replace("_2of2"," (1 & 2 of 2)")


DotPlot(object = object_1, features =list('C3'=tmp.c3, 'OPC'=tmp.opc, 'C3+OPC'=tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC] in: ",sid_print, " (Couturier dataset)"))



ggsave(paste0("output/figures/2022_figure_S7b.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)




#DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
#  labs(x = paste0("Features [C2/OPC] in: ",sid))

#ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
#ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)



#### 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "ESM1")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = c("RGS5","PDGFRB","CD248","PEAR1", "HEYL" , "CFH"))

FeaturePlot(object = object_1, features = c("RGS5"))
FeaturePlot(object = object_1, features = c("PDGFRB"))
FeaturePlot(object = object_1, features = c("CD248"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))




#### C3 (down) :: endothelial ----


endo <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_specificity') %>%
  dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
  dplyr::filter(Celltype == 'end') %>% 
  dplyr::arrange(desc(grand_mean)) %>%
  dplyr::filter(gene %in% all.genes ) %>%
  dplyr::slice_head(n=25) %>%
  dplyr::mutate(grand_mean = NULL) %>% 
  dplyr::pull(gene)


C3.only <- setdiff(C3, endo)
C3.and.endo <- intersect(endo, C3)
endo.only <- setdiff(endo, C3)



DotPlot(object = object_1, features = list('C3'=C3.only, 'C3+endo'= C3.and.endo, 'endo'=endo.only,'pericyte'=c('PDGFRB','CD248','RGS5')), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C3 & top25 McKenzy endothelial cell markers] in: ",sid))


ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.pdf"),width=7.5, height=3,scale=2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.png"),width=7.5, height=3,scale=2)



RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C3)


#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)



FeaturePlot(object = object_1, features = C4A[1:4])
FeaturePlot(object = object_1, features = C4A[5:8])
FeaturePlot(object = object_1, features = C4A[9:12])
FeaturePlot(object = object_1, features = C4A[13:16])
FeaturePlot(object = object_1, features = C4A[17:19])

FeaturePlot(object = object_1, features = C4B[1:4])
FeaturePlot(object = object_1, features = C4B[5:9])



#### C5 (down) ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)



FeaturePlot(object = object_1, features = C5[1:4])
FeaturePlot(object = object_1, features = C5[5:8])
FeaturePlot(object = object_1, features = C5[9:12])
FeaturePlot(object = object_1, features = C5[13:16])


#### C6 (up) ----


DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)


VlnPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))


VlnPlot(object = object_1, features = C6, stack = T, sort = F, group.by = "seurat_clusters")
VlnPlot(object = object_1, features = 'PC_1', group.by = "seurat_clusters") # PE
VlnPlot(object = object_1, features = 'PC_2', group.by = "seurat_clusters")
VlnPlot(object = object_1, features = 'PC_3', group.by = "seurat_clusters")
VlnPlot(object = object_1, features = 'PC_4', group.by = "seurat_clusters")
VlnPlot(object = object_1, features = 'PC_5', group.by = "seurat_clusters")
VlnPlot(object = object_1, features = 'PC_6', group.by = "seurat_clusters")


plot(object_1@reductions$pca@feature.loadings[,1], 
     1:length(object_1@reductions$pca@feature.loadings[,1]),
     col=rownames(object_1@reductions$pca@feature.loadings) %in% C6  + 1 , pch=19)

wilcox.test(object_1@reductions$pca@feature.loadings[rownames(object_1@reductions$pca@feature.loadings) %in% C6 ,1],
            object_1@reductions$pca@feature.loadings[rownames(object_1@reductions$pca@feature.loadings) %in% C6 == F ,1])

wilcox.test(object_1@reductions$pca@feature.loadings[rownames(object_1@reductions$pca@feature.loadings) %in% C5 ,1],
            object_1@reductions$pca@feature.loadings[rownames(object_1@reductions$pca@feature.loadings) %in% C5 == F ,1])

wilcox.test(object_1@reductions$pca@feature.loadings[rownames(object_1@reductions$pca@feature.loadings) %in% c(C4A,C4B) ,1],
            object_1@reductions$pca@feature.loadings[rownames(object_1@reductions$pca@feature.loadings) %in% c(C4A,C4B) == F ,1])


plot(object_1@reductions$pca@feature.loadings[,3], 
     1:length(object_1@reductions$pca@feature.loadings[,1]),
     col=rownames(object_1@reductions$pca@feature.loadings) %in% C6  + 1 , pch=19)



FeaturePlot(object = object_1, features = C6)



#### CC-2022 (up) ----


DotPlot(object = object_1, features =list('C1'=
      results.out |> dplyr::filter(C1.2022) |> dplyr::pull(hugo_symbol)
                                             , 'Peri'=c("RGS5", "PDGFRB", "CD248")), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C1] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.png"),width=7.5, height=4,scale=1.2)


#### Wang Sub-types ----


object_t <- subset(object_1, class == "T")
object_t <- RunPCA(object_t, features = subtype.classical$symbol, reduction.name = 'pca.subtype.cl', reduction.key = 'PCcl')
object_t <- RunPCA(object_t, features = subtype.mesenchymal$symbol, reduction.name = 'pca.subtype.mes', reduction.key = 'PCmes')
object_t <- RunPCA(object_t, features = subtype.proneural$symbol, reduction.name = 'pca.subtype.pn', reduction.key = 'PCpn')

if(sum(object_t@reductions$pca.subtype.mes@feature.loadings[,1] < 0) > sum(object_t@reductions$pca.subtype.mes@feature.loadings[,1] > 0)) {
  object_t@reductions$pca.subtype.mes@cell.embeddings[,1] <- object_t@reductions$pca.subtype.mes@cell.embeddings[,1] * -1
}
if(sum(object_t@reductions$pca.subtype.cl@feature.loadings[,1] < 0) > sum(object_t@reductions$pca.subtype.cl@feature.loadings[,1] > 0)) {
  object_t@reductions$pca.subtype.cl@cell.embeddings[,1] <- object_t@reductions$pca.subtype.cl@cell.embeddings[,1] * -1
}
if(sum(object_t@reductions$pca.subtype.pn@feature.loadings[,1] < 0) > sum(object_t@reductions$pca.subtype.pn@feature.loadings[,1] > 0)) {
  object_t@reductions$pca.subtype.pn@cell.embeddings[,1] <- object_t@reductions$pca.subtype.pn@cell.embeddings[,1] * -1
}



DotPlot(object = object_t, features = list('MES (Wang)'=subtype.mesenchymal$symbol,
                                           'CL (Wang)'=subtype.classical$symbol,
                                           'PN (Wang)'=subtype.proneural$symbol
                                           
), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


object_sty <- merge(subset(object_1, class != "T"), object_t)
object_sty@reductions <- object_t@reductions
object_sty@reductions$umap <- object_1@reductions$umap
rm(object_t)


FeaturePlot(object = object_sty, features = c("PCcl_1","PCmes_1","PCpn_1"), max.cutoff = 5)



#rm(object_sty)



#### Neftel Sub-types ----

source('scripts/R/neftel_meta_modules.R')

object_t <- subset(object_1, class == "T")
object_t <- RunPCA(object_t, features = neftel.meta.modules.MES2.tt2, reduction.name = 'pca.subtype.mes2', reduction.key = 'PCmes2')
object_t <- RunPCA(object_t, features = neftel.meta.modules.MES1.tt2, reduction.name = 'pca.subtype.mes1', reduction.key = 'PCmes1')
object_t <- RunPCA(object_t, features = neftel.meta.modules.AC.tt2, reduction.name = 'pca.subtype.ac', reduction.key = 'PCac')
object_t <- RunPCA(object_t, features = neftel.meta.modules.OPC.tt2, reduction.name = 'pca.subtype.opc', reduction.key = 'PCopc')
object_t <- RunPCA(object_t, features = neftel.meta.modules.NPC1.tt2, reduction.name = 'pca.subtype.npc1', reduction.key = 'PCnpc1')
object_t <- RunPCA(object_t, features = neftel.meta.modules.NPC2.tt2, reduction.name = 'pca.subtype.npc2', reduction.key = 'PCnpc2')

if(sum(object_t@reductions$pca.subtype.mes2@feature.loadings[,1] < 0) > sum(object_t@reductions$pca.subtype.mes2@feature.loadings[,1] > 0)) {
  object_t@reductions$pca.subtype.mes2@cell.embeddings[,1] <- object_t@reductions$pca.subtype.mes2@cell.embeddings[,1] * -1
}

# DotPlot(object = object_t, features = list('MES2 (Neftel)'=neftel.meta.modules.MES2.tt2,
#                                            'CL (Wang)'=subtype.classical$symbol,
#                                            'PN (Wang)'=subtype.proneural$symbol
#                                            
# ), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
# 



object_sty <- merge(subset(object_1, class != "T"), object_t)
object_sty@reductions <- object_t@reductions
object_sty@reductions$umap <- object_1@reductions$umap
rm(object_t)


FeaturePlot(object = object_sty, features = c("PCmes2_1","PCmes1_1",
                                              "PCac_1","PCopc_1",
                                              "PCnpc1_1","PCnpc2_2"))



rm(object_sty)


## BT346 [poor separation or high tumor?] ----


sid <- "BT346.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 5750,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 25000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1000 &
                     nFeature_RNA < 5750 &
                     nCount_RNA > 500 &
                     nCount_RNA < 25000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



#### 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### C6 (up) :: PRESENT! <Levi> ----


DotPlot(object = object_1, features = C6) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(object = object_1, features = C6, stack = T, sort = T)
VlnPlot(object = object_1, features = C6, stack = T, sort = T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])



## BT363-GSC :: T,MG,OD,unknown type ----


rm(sid, object_1)
gc()

sid <- "BT363-GSC.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4500,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(0|1|2|4|5|6|7|8|10|11|12)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(3)$",paste0("TAM/microglia next to T-cell"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("?Leukocyte?.\\1"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




#tmp.13 <- FindMarkers(object_1, ident.1 = 13)
#head(tmp.13, 20) # LAMP3,IRF4,NCCRP1,CRIP1,SYNPO2,CCR7,EHF,CCL22,VTN,LSP1,CDX2,IDO1,TRAF3IP3,RP11-290F5.1,ANXA3,RP1-28O10.1,RP11-66B24.4,AP003774.1,NR4A3,CCL17
#FeaturePlot(object = object_1, features = "LAMP3")


#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


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



#### 2. Astrocyte (???) ----




FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")




#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor


#### 3D. ? Leukocyte ?? ----


# These are cluster-13 DE genes, of which some at genecards seem related to leukocytes?
FeaturePlot(object = object_1, features = c("LAMP3","IRF4","NCCRP1","CRIP1","SYNPO2","CCR7","EHF","CCL22","VTN","LSP1","CDX2"))



# FeaturePlot(object = object_1, features = c("CD19", "CD20", "CD34", "CD38")) # b-cell markers?



#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")


#### 5. Oligodendrocytes (-) ----

FeaturePlot(object = object_1, features = c("TMEM125","TMEM144"))


#### 6A. Endothelial (-) ----

FeaturePlot(object = object_1, features = "CD34")


#### 6B. Pericytes (-) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "CD248")


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



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6


DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)

FeaturePlot(object = object_1, features = C6[1:8])
FeaturePlot(object = object_1, features = C6[9:16])
FeaturePlot(object = object_1, features = C6[17:24])
FeaturePlot(object = object_1, features = C6[25:33])




## BT363 [1+2/2] :: T,MG,OD,PE,T(mitotic) ----


rm(sid, object_1)
gc()

sid <- "BT363_1of2.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

sid <- "BT363_2of2.filtered_gene_matrices"
object_1.tmp <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1.tmp <- CreateSeuratObject(counts = object_1.tmp, min.cells = 3, min.features = 200, project="Couturier")

object_1.m <- merge(object_1, y=object_1.tmp, add.cell.ids = c("1of2","2of2"), project="Couturier")

rm(object_1, object_1.tmp)
object_1 <- object_1.m
rm(object_1.m)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 700,col="red") +
  geom_hline(yintercept = 5250,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset = 
                     nFeature_RNA > 700 &
                     nFeature_RNA < 5250 & 
                     nCount_RNA > 1000 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.15)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- as.factor(paste( gsub("_.+$","",sid) , gsub("_.+$","",colnames(object_1)), sep='-'))
object_1[["state"]]
object_1$dataset <- as.character(object_1$state)


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2

# scaling of data
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


# cluster the cells
object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|8|10|11|12|14)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(7|13)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6|9)$","Tam/MG.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^15$","Endothelial",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^16$","Pericytes",levels(object_1$seurat_clusters))
object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(6,9),"TAM", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(7,13),"OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0:5,8,10:12,14) ,"T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(15,16) & 
                           object_1@reductions$umap@cell.embeddings[,2] > 0.2
                         ,"PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(16) & 
                           object_1@reductions$umap@cell.embeddings[,2] < 0.2 # -0.75
                         ,"EN|PE?", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(15) & 
                           object_1@reductions$umap@cell.embeddings[,2] < 0.2 # -0.75
                         ,"EN", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)



object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","1. T","2. T","3. T","4. T","5. T","8. T","10. T","11. T","12. T","14. T",
  "7. OD","13. OD",
  "15. EN",
  "16. EN|PE?",
  "16. PE",
  "6. TAM","9. TAM"
))



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters") +
  labs(subtitle=sid)



ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.png"),width=10,height=8)


median(object_1[,object_1$seurat_clusters == "16. PE"]$nCount_RNA)
median(object_1[,object_1$seurat_clusters == "16. PE"]$nFeature_RNA)
mean(object_1[,object_1$seurat_clusters == "16. PE"]$nCount_RNA)
mean(object_1[,object_1$seurat_clusters == "16. PE"]$nFeature_RNA)

median(object_1[,object_1$seurat_clusters == "16. EN|PE?"]$nCount_RNA)
median(object_1[,object_1$seurat_clusters == "16. EN|PE?"]$nFeature_RNA)
mean(object_1[,object_1$seurat_clusters == "16. EN|PE?"]$nCount_RNA)
mean(object_1[,object_1$seurat_clusters == "16. EN|PE?"]$nFeature_RNA)



### Prepare for integration ----


object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-0.75, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(6,9),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(7,13),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0:5,8,10:12,14) ,"Tumor", object_1$youri_clusters)
#object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(15,16) ,"Endothelial", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(15,16) & 
                                    object_1@reductions$umap@cell.embeddings[,2] > 0.2
                                  ,"Pericytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(15,16) & 
                                    object_1@reductions$umap@cell.embeddings[,2] < 0.2 # -0.75
                                  ,"Endothelial", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")


# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=12, linetype="dashed", color = "red") +
  geom_vline(xintercept=16, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_hline(yintercept=1, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 12 &
    object_1@reductions$umap@cell.embeddings[,1] <= 16 &
    object_1@reductions$umap@cell.embeddings[,2] >= 0.2 &
    object_1@reductions$umap@cell.embeddings[,2] <= 1 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)

# # Peri|Endo
# DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
#   geom_vline(xintercept=12, linetype="dashed", color = "red") +
#   geom_vline(xintercept=16, linetype="dashed", color = "red") +
#   geom_hline(yintercept=-0.75, linetype="dashed", color = "red") +
#   geom_hline(yintercept=0.2, linetype="dashed", color = "red")
# object_1$youri_clusters <- ifelse(
#   object_1@reductions$umap@cell.embeddings[,1] >= 12 &
#     object_1@reductions$umap@cell.embeddings[,1] <= 16 &
#     object_1@reductions$umap@cell.embeddings[,2] >= -0.75 &
#     object_1@reductions$umap@cell.embeddings[,2] <= 0.2 &
#     object_1$youri_clusters != "Pericytes|Endothelial", "ambiguous", object_1$youri_clusters)


# Endo
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=12, linetype="dashed", color = "red") +
  geom_vline(xintercept=16, linetype="dashed", color = "red") +
  geom_hline(yintercept=-2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-0.75, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 12 &
    object_1@reductions$umap@cell.embeddings[,1] <= 16 &
    object_1@reductions$umap@cell.embeddings[,2] >= -2 &
    object_1@reductions$umap@cell.embeddings[,2] <= 0.2 & #-0.75 &
    object_1$youri_clusters != "Endothelial", "ambiguous", object_1$youri_clusters)



# TAM/MG
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-2, linetype="dashed", color = "red") +
  geom_vline(xintercept=4, linetype="dashed", color = "red") +
  geom_hline(yintercept=-18, linetype="dashed", color = "red") +
  geom_hline(yintercept=-8, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -2 &
    object_1@reductions$umap@cell.embeddings[,1] <= 4 &
    object_1@reductions$umap@cell.embeddings[,2] >= -18 &
    object_1@reductions$umap@cell.embeddings[,2] <= -8 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)


# OD
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-15, linetype="dashed", color = "red") +
  geom_vline(xintercept=-8, linetype="dashed", color = "red") +
  geom_hline(yintercept=-9, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -15 &
    object_1@reductions$umap@cell.embeddings[,1] <= -8 &
    object_1@reductions$umap@cell.embeddings[,2] >= -9 &
    object_1@reductions$umap@cell.embeddings[,2] <= -3 &
    object_1$youri_clusters != "Oligodendrocytes", "ambiguous", object_1$youri_clusters)



# Tumor
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-8, linetype="dashed", color = "red") +
  geom_vline(xintercept=10, linetype="dashed", color = "red") +
  geom_hline(yintercept=-5, linetype="dashed", color = "red") +
  geom_hline(yintercept=10, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -8 &
    object_1@reductions$umap@cell.embeddings[,1] <= 10 &
    object_1@reductions$umap@cell.embeddings[,2] >= -5 &
    object_1@reductions$umap@cell.embeddings[,2] <= 10 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)

sum(object_1$youri_clusters == "ambiguous")



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")



object_1.BT363 <- object_1



#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


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



#### 2. Astrocyte (-) ----




FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")





#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor



#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")
FeaturePlot(object = object_1, features = "ELAVL4")# ? https://www.biorxiv.org/content/10.1101/775197v1



#### 5. Oligodendrocytes + OPC (+) ----

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")




##### Figure S7c ----


tmp.c3 <- results.out |>
  dplyr::filter(!is.na(.data$C3.2022)) |> 
  dplyr::filter(.data$C3.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.opc <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.OPC)) |> 
  dplyr::filter(.data$neftel.meta.modules.OPC == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.c3.opc <- intersect(tmp.c3, tmp.opc)
tmp.c3 <- setdiff(tmp.c3, tmp.c3.opc)
tmp.opc <- setdiff(tmp.opc, tmp.c3.opc)


sid_print <- sid |> 
  stringr::str_replace(".filtered_gene_matrices","") |> 
  stringr::str_replace("_2of2"," (1 & 2 of 2)")


DotPlot(object = object_1, features =list('C3'=tmp.c3, 'OPC'=tmp.opc, 'C3+OPC'=tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC] in: ",sid_print, " (Couturier dataset)"))



ggsave(paste0("output/figures/2022_figure_S7c.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)




#### 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248",pt.size=0.04)


#### C3 endo (down) ----



endo <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_specificity') %>%
  dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
  dplyr::filter(Celltype == 'end') %>% 
  dplyr::arrange(desc(grand_mean)) %>%
  dplyr::filter(gene %in% all.genes ) %>%
  dplyr::slice_head(n=25) %>%
  dplyr::mutate(grand_mean = NULL) %>% 
  dplyr::pull(gene)


C3.only <- setdiff(C3, endo)
C3.and.endo <- intersect(endo, C3)
endo.only <- setdiff(endo, C3)


DotPlot(object = object_1, features = list('C3'=C3.only, 'C3+endo'= C3.and.endo, 'endo'=endo.only,'pericyte'=c('PDGFRB','CD248','RGS5')), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C3 & top25 McKenzy endothelial cell markers] in: ",sid))


ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.pdf"),width=7.5, height=3,scale=2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.png"),width=7.5, height=3,scale=2)




#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)




FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)

DotPlot(object = object_1, features = c(C4A, C4B))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)




FeaturePlot(object = object_1, features = C5)


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)



#### C6 (up) :: some cor w/ pericytes ----


DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))

ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)



RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)



FeaturePlot(object = object_1, features =  "CRABP2" )
FeaturePlot(object = object_1, features =  "CILP2" )
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



#### CC-2022 (up) ----



DotPlot(object = object_1, features =list('C1'=
                                            results.out |> dplyr::filter(C1.2022) |> dplyr::pull(hugo_symbol)
                                          , 'Peri'=c("RGS5", "PDGFRB", "CD248")), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C1] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.png"),width=7.5, height=4,scale=1.2)



## BT364 [1+2/2] :: T,MG,OD ----


rm(sid, object_1)
gc()

sid <- "BT364_1of2.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

sid <- "BT364_2of2.filtered_gene_matrices"
object_1.tmp <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1.tmp <- CreateSeuratObject(counts = object_1.tmp, min.cells = 3, min.features = 200, project="Couturier")

object_1.m <- merge(object_1, y=object_1.tmp, add.cell.ids = c("1of2","2of2"), project="Couturier")

rm(object_1, object_1.tmp)
object_1 <- object_1.m
rm(object_1.m)



mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 6000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 35000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset = 
                     nFeature_RNA > 1000 &
                     nFeature_RNA < 6000 & 
                     nCount_RNA > 1000 &
                     nCount_RNA < 35000 &
                     percent.mito < 0.15)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- as.factor(paste( gsub("_.+$","",sid) , gsub("_.+$","",colnames(object_1)), sep='-'))
object_1[["state"]]
object_1$dataset <- as.character(object_1$state)


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



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
ElbowPlot(object_1, ndims = 45)


### cluster the cells

object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|9|12|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^15$","Tumor\nimmune cell recruiting?",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(10|14)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6|11|8)$","TAM/MG\\1",levels(object_1$seurat_clusters))

object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(10,14) &
                           object_1@reductions$umap@cell.embeddings[,2] > 11.5, 
                         "OD", object_1$class)
object_1$class <- ifelse( object_1@reductions$umap@cell.embeddings[,1] >= -3 &
                          object_1@reductions$umap@cell.embeddings[,1] <= 0 &
                          object_1@reductions$umap@cell.embeddings[,2] >= 10 &
                          object_1@reductions$umap@cell.embeddings[,2] <= 11.5   , 
                         "PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(13) &
                           object_1@reductions$umap@cell.embeddings[,1] < -10.25,
                         "T ?", object_1$class) # NPC2 cluster
object_1$class <- ifelse(object_1$seurat_clusters %in% c(6,8,11),"TAM", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0:5,7,9,12,13,15) &
                           object_1@reductions$umap@cell.embeddings[,1] >= -10.25 ,"T", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","1. T","2. T","3. T","4. T","5. T","7. T","9. T","12. T","13. T","15. T",  "13. T ?",
  "10. OD","14. OD",
  "14. PE",
  "6. TAM","8. TAM","11. TAM"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") +
  labs(subtitle=sid)



ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.png"),width=10,height=8)



### Prepare for integration ----


object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_hline(yintercept=11.5, linetype="dashed", color = "red") +
  geom_vline(xintercept=-10.25, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(10,14) &
                                    object_1@reductions$umap@cell.embeddings[,2] > 11.5, 
                                  "Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(14) &
                                    object_1@reductions$umap@cell.embeddings[,2] <= 11.5, 
                                  "Pericytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(13) &
                                    object_1@reductions$umap@cell.embeddings[,1] < -10.25,
                                  "Tumor", object_1$youri_clusters) # NPC2 cluster
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(6,8,11),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0:5,7,9,12,13,15) &
                                    object_1@reductions$umap@cell.embeddings[,1] >= -10.25 ,"Tumor", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")


# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  geom_hline(yintercept=10, linetype="dashed", color = "red") +
  geom_hline(yintercept=11.5, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -3 &
    object_1@reductions$umap@cell.embeddings[,1] <= 0 &
    object_1@reductions$umap@cell.embeddings[,2] >= 10 &
    object_1@reductions$umap@cell.embeddings[,2] <= 11.5 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)



# TAM/MG
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=9, linetype="dashed", color = "red") +
  geom_vline(xintercept=17, linetype="dashed", color = "red") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_hline(yintercept=7, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 9 &
    object_1@reductions$umap@cell.embeddings[,1] <= 17 &
    object_1@reductions$umap@cell.embeddings[,2] >= 0 &
    object_1@reductions$umap@cell.embeddings[,2] <= 7 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)


# OD
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-6, linetype="dashed", color = "red") +
  geom_vline(xintercept=1, linetype="dashed", color = "red") +
  geom_hline(yintercept=11.5, linetype="dashed", color = "red") +
  geom_hline(yintercept=17, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -6 &
    object_1@reductions$umap@cell.embeddings[,1] <= 1 &
    object_1@reductions$umap@cell.embeddings[,2] >= 11.5 &
    object_1@reductions$umap@cell.embeddings[,2] <= 17 &
    object_1$youri_clusters != "Oligodendrocytes", "ambiguous", object_1$youri_clusters)



# Tumor
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=-10.25, linetype="dashed", color = "red") +
  geom_vline(xintercept=6.5, linetype="dashed", color = "red") +
  geom_hline(yintercept=-9.5, linetype="dashed", color = "red") +
  geom_hline(yintercept=5, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= -10.25 &
    object_1@reductions$umap@cell.embeddings[,1] <= 6.5 &
    object_1@reductions$umap@cell.embeddings[,2] >= -9.5 &
    object_1@reductions$umap@cell.embeddings[,2] <= 5 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)




DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")
sum(object_1$youri_clusters == "ambiguous")

object_1.BT364 <- object_1




#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))


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



#### 2. Astrocyte (+) ----

# FeaturePlot(object = object_1, features = "STMN2") # Tumor
# FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor



#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 5. Oligodendrocytes + OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)




##### Figure S7d ----


tmp.c3 <- results.out |>
  dplyr::filter(!is.na(.data$C3.2022)) |> 
  dplyr::filter(.data$C3.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.opc <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.OPC)) |> 
  dplyr::filter(.data$neftel.meta.modules.OPC == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.c3.opc <- intersect(tmp.c3, tmp.opc)
tmp.c3 <- setdiff(tmp.c3, tmp.c3.opc)
tmp.opc <- setdiff(tmp.opc, tmp.c3.opc)


sid_print <- sid |> 
  stringr::str_replace(".filtered_gene_matrices","") |> 
  stringr::str_replace("_2of2"," (1 & 2 of 2)")


DotPlot(object = object_1, features =list('C3'=tmp.c3, 'OPC'=tmp.opc, 'C3+OPC'=tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC] in: ",sid_print, " (Couturier dataset)"))



ggsave(paste0("output/figures/2022_figure_S7d.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)





#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")




#### C3 endo (down) ----


endo <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_specificity') %>%
  dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
  dplyr::filter(Celltype == 'end') %>% 
  dplyr::arrange(desc(grand_mean)) %>%
  dplyr::filter(gene %in% all.genes ) %>%
  dplyr::slice_head(n=25) %>%
  dplyr::mutate(grand_mean = NULL) %>% 
  dplyr::pull(gene)


C3.only <- setdiff(C3, endo)
C3.and.endo <- intersect(endo, C3)
endo.only <- setdiff(endo, C3)


DotPlot(object = object_1, features = list('C3'=C3.only, 'C3+endo'= C3.and.endo, 'endo'=endo.only,'pericyte'=c('PDGFRB','CD248','RGS5')), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C3 & top25 McKenzy endothelial cell markers] in: ",sid))


ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.pdf"),width=7.5, height=3,scale=2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.png"),width=7.5, height=3,scale=2)




#### C4 (up) ----

DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)



FeaturePlot(object = object_1, features = C4A[1:4])
FeaturePlot(object = object_1, features = C4A[5:8])
FeaturePlot(object = object_1, features = C4A[9:12])
FeaturePlot(object = object_1, features = C4A[13:16])
FeaturePlot(object = object_1, features = C4A[17:19])

FeaturePlot(object = object_1, features = C4B[1:4])
FeaturePlot(object = object_1, features = C4B[5:9])



#### C5 (down) ----

f <- C5
DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


FeaturePlot(object = object_1, features = C5[1:4])
FeaturePlot(object = object_1, features = C5[5:8])
FeaturePlot(object = object_1, features = C5[9:12])
FeaturePlot(object = object_1, features = C5[13:16])

#### C6 (up) ----

# DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248","HEYL","CFH") ), group.by = "seurat_clusters") +
DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))

ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)



RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)

FeaturePlot(object = object_1, features = C6[1:8])
FeaturePlot(object = object_1, features = C6[9:16])
FeaturePlot(object = object_1, features = C6[17:24])
FeaturePlot(object = object_1, features = C6[25:33])


#### CC-2022 (up) ----


DotPlot(object = object_1, features =list('C1'=
                                            
                                            results.out |> dplyr::filter(C1.2022) |> dplyr::pull(hugo_symbol)
                                            , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C1] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.png"),width=7.5, height=4,scale=1.2)




## BT368-GSC [95-100% tumor] ----


rm(object_1)
gc()

sid <- "BT368-GSC.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



## BT368 Total [95-100% tumor?] <geen scheiding, per pericytes, vraag levi voor hulp?> ----

rm(object_1)
gc()


sid <- "BT368.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 6000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 35000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1000 &
                     nFeature_RNA < 6000 &
                     nCount_RNA > 1000 &
                     nCount_RNA < 35000 &
                     percent.mito < 0.15)

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


### scaling of data ----
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

### cluster the cells ----


object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


#levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|9|12|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


#### 1. Tumor (+) ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))


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


#### 2. Astrocyte (+) ----


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")

#### 4. Neurons (?) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 5. Oligodendrocytes (?) ----

FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "PLP1")



#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) ----

DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)



#### C6 (up) ----

#f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T) 
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6)




## BT389 [95-100% tumor, low res] ----

rm(object_1)
gc()

sid <- "BT389.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 

# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 35000,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 1000 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 35000 &
#                      percent.mito < 0.15)

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


### scaling of data ----
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

### cluster the cells ----


object_1 <- FindNeighbors(object_1, dims = 1:40)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


#levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|9|12|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


#levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|9|12|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


#### 1. Tumor (+) ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))


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


#### 2. Astrocyte (?) ----


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")

#### 4. Neurons (?) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 5. Oligodendrocytes (?) ----

FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN")
FeaturePlot(object = object_1, features = "PLP1")



#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) ----

DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)




## BT390 :: T,MG++,TC,OD of 2e tum-clone? CNV nodig ----


sid <- "BT390.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 3750,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 3750 &
                     nCount_RNA > 500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset
ElbowPlot(object_1, ndims = 45)


object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(1|4|9|8|6)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|2|10|12)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(11)$",paste0("\\1. T|OD?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("\\1. TC"),levels(object_1$seurat_clusters))

levels(object_1$seurat_clusters) <- gsub("^(3|5|7)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T",
  "2. T",
  "10. T",
  "11. T|OD?",
  "12. T",
  "3. OD",
  "5. OD",
  "7. OD",
  "1. TAM",
  "4. TAM",
  "6. TAM",
  "8. TAM",
  "9. TAM",
  "13. TC"  
))



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") +
  labs(subtitle=sid)




#### 1. Tumor (?) ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("GFAP","BTC","EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))


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


#### 2. Astrocyte (-) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")




#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 5. Oligodendrocytes + OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN")
FeaturePlot(object = object_1, features = "PLP1")




##### Figure S7e ----


tmp.c3 <- results.out |>
  dplyr::filter(!is.na(.data$C3.2022)) |> 
  dplyr::filter(.data$C3.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.opc <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.OPC)) |> 
  dplyr::filter(.data$neftel.meta.modules.OPC == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.c3.opc <- intersect(tmp.c3, tmp.opc)
tmp.c3 <- setdiff(tmp.c3, tmp.c3.opc)
tmp.opc <- setdiff(tmp.opc, tmp.c3.opc)


sid_print <- sid |> 
  stringr::str_replace(".filtered_gene_matrices","")


DotPlot(object = object_1, features =list('C3'=tmp.c3, 'OPC'=tmp.opc, 'C3+OPC'=tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC] in: ",sid_print, " (Couturier dataset)"))



ggsave(paste0("output/figures/2022_figure_S7e.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)





#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")

#### C3 (up) ----

f <- C3
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


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


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)



## BT397 [1+2/2] :: MG+++,TC+,T,OD,EN-,PE- :: mooie ----


rm(object_1)
gc()

sid <- "BT397_1of2.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

sid <- "BT397_2of2.filtered_gene_matrices"
object_1.tmp <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1.tmp <- CreateSeuratObject(counts = object_1.tmp, min.cells = 3, min.features = 200, project="Couturier")

object_1.m <- merge(object_1, y=object_1.tmp, add.cell.ids = c("1of2","2of2"), project="Couturier")

rm(object_1, object_1.tmp)
object_1 <- object_1.m
rm(object_1.m)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 650,col="red") +
  geom_hline(yintercept = 4500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 650 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- as.factor(paste( gsub("_.+$","",sid) , gsub("_.+$","",colnames(object_1)), sep='-'))
object_1$dataset <- as.character(object_1$state)


top10 <- head(VariableFeatures(object_1), 10)



print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))






# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0,1,11,7,14,2,5,8,6,10,15,3,9),"TAM", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(16),"OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(12),"TC", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(13,4,17) ,"T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(18) & object_1@reductions$umap@cell.embeddings[,2] >= 7.5,"EN", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(18) & object_1@reductions$umap@cell.embeddings[,2] < 7.5, "PE", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "4. T","13. T","17. T",
  "16. OD",
  "18. EN",
  "18. PE",
  "0. TAM","1. TAM","2. TAM","3. TAM","5. TAM","6. TAM","7. TAM","8. TAM","9. TAM" ,"10. TAM","11. TAM" ,"14. TAM","15. TAM",
  "12. TC"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_UMAP.png"),width=10,height=8)



# sum(object_1$seurat_clusters == "18. PE")



### Prepare for integration ----


object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_hline(yintercept=7.5, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0,1,11,7,14,2,5,8,6,10,15,3,9),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(16),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(12),"T-Cells", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(13,4,17) ,"Tumor", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(18) & object_1@reductions$umap@cell.embeddings[,2] >= 7.5,"Endothelial", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(18) & object_1@reductions$umap@cell.embeddings[,2] < 7.5, "Pericytes", object_1$youri_clusters)

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")


# Peri
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=15, linetype="dashed", color = "red") +
  geom_vline(xintercept=18, linetype="dashed", color = "red") +
  geom_hline(yintercept=4, linetype="dashed", color = "red") +
  geom_hline(yintercept=6, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 15 &
    object_1@reductions$umap@cell.embeddings[,1] <= 18 &
    object_1@reductions$umap@cell.embeddings[,2] >= 4 &
    object_1@reductions$umap@cell.embeddings[,2] <= 6 &
    object_1$youri_clusters != "Pericytes", "ambiguous", object_1$youri_clusters)


# T-cells
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=5, linetype="dashed", color = "red") +
  geom_vline(xintercept=10, linetype="dashed", color = "red") +
  geom_hline(yintercept=-15, linetype="dashed", color = "red") +
  geom_hline(yintercept=-7, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 5 &
    object_1@reductions$umap@cell.embeddings[,1] <= 10 &
    object_1@reductions$umap@cell.embeddings[,2] >= -15 &
    object_1@reductions$umap@cell.embeddings[,2] <= -7 &
    object_1$youri_clusters != "T-Cells", "ambiguous", object_1$youri_clusters)


# TAM/MG
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=5, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] < 5 &
    object_1$youri_clusters != "TAM/MG", "ambiguous", object_1$youri_clusters)


# OD
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=8, linetype="dashed", color = "red") +
  geom_vline(xintercept=12, linetype="dashed", color = "red") +
  geom_hline(yintercept=6, linetype="dashed", color = "red") +
  geom_hline(yintercept=11, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 8 &
    object_1@reductions$umap@cell.embeddings[,1] <= 12 &
    object_1@reductions$umap@cell.embeddings[,2] >= 6 &
    object_1@reductions$umap@cell.embeddings[,2] <= 11 &
    object_1$youri_clusters != "Oligodendrocytes", "ambiguous", object_1$youri_clusters)


# Tumor
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=10, linetype="dashed", color = "red") +
  geom_vline(xintercept=18, linetype="dashed", color = "red") +
  geom_hline(yintercept=-10, linetype="dashed", color = "red") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 10 &
    object_1@reductions$umap@cell.embeddings[,1] <= 18 &
    object_1@reductions$umap@cell.embeddings[,2] >= -10 &
    object_1@reductions$umap@cell.embeddings[,2] <= 0 &
    object_1$youri_clusters != "Tumor", "ambiguous", object_1$youri_clusters)




DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")
sum(object_1$youri_clusters == "ambiguous")

object_1.BT397 <- object_1



#### 1. Tumor (?) ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

# cluster 9:
# [1] "KIF2C"  "NUF2"   "ASPM"   "NEK2"   "CENPA"  "CKAP2L" "SGOL1"  "CENPE"  "CCNA2"  "PBK"    "MKI67"  "CDCA3"  "NUSAP1" "CCNB2"  "KIF23" 
# [16] "FAM64A" "AURKB"  "TOP2A"  "TPX2"   "CDC20" 
# mitosis genes

FeaturePlot(object = object_1, features = c("KIF2C", "NUF2", "OLIG1", "VIM", "S100B", "GFAP","NDRG1","RBFOX2"))


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



#### 2. Astrocyte (-) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")


#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----


FeaturePlot(object = object_1, features = "HBG1") # Tumor


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")



#### 5. Oligodendrocytes + OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN")
FeaturePlot(object = object_1, features = "PLP1")



##### Figure S7f ----


tmp.c3 <- results.out |>
  dplyr::filter(!is.na(.data$C3.2022)) |> 
  dplyr::filter(.data$C3.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.opc <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.OPC)) |> 
  dplyr::filter(.data$neftel.meta.modules.OPC == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.c3.opc <- intersect(tmp.c3, tmp.opc)
tmp.c3 <- setdiff(tmp.c3, tmp.c3.opc)
tmp.opc <- setdiff(tmp.opc, tmp.c3.opc)


sid_print <- sid |> 
  stringr::str_replace(".filtered_gene_matrices","") |> 
  stringr::str_replace("_2of2"," (1 & 2 of 2)")


DotPlot(object = object_1, features =list('C3'=tmp.c3, 'OPC'=tmp.opc, 'C3+OPC'=tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC] in: ",sid_print, " (Couturier dataset)"))



ggsave(paste0("output/figures/2022_figure_S7f.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)





#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### C3 (up) ----


endo <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_specificity') %>%
  dplyr::select(c('grand_mean', 'gene', 'Celltype')) %>%
  dplyr::filter(Celltype == 'end') %>% 
  dplyr::arrange(desc(grand_mean)) %>%
  dplyr::filter(gene %in% all.genes ) %>%
  dplyr::slice_head(n=25) %>%
  dplyr::mutate(grand_mean = NULL) %>% 
  dplyr::pull(gene)


C3.only <- setdiff(C3, endo)
C3.and.endo <- intersect(endo, C3)
endo.only <- setdiff(endo, C3)



DotPlot(object = object_1, features = list('C3'=C3.only, 'C3+endo'= C3.and.endo, 'endo'=endo.only,'pericyte'=c('PDGFRB','CD248','RGS5')), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C3 & top25 McKenzy endothelial cell markers] in: ",sid))


ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.pdf"),width=7.5, height=3,scale=2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C3.png"),width=7.5, height=3,scale=2)




#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)




FeaturePlot(object = object_1, features = C4A[1:4])
FeaturePlot(object = object_1, features = C4A[5:8])
FeaturePlot(object = object_1, features = C4A[9:12])
FeaturePlot(object = object_1, features = C4A[13:16])
FeaturePlot(object = object_1, features = C4A[17:19])

FeaturePlot(object = object_1, features = C4B[1:4])
FeaturePlot(object = object_1, features = C4B[5:9])



#### C5 (down) ----


DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)




FeaturePlot(object = object_1, features = C5[1:4])
FeaturePlot(object = object_1, features = C5[5:8])
FeaturePlot(object = object_1, features = C5[9:12])
FeaturePlot(object = object_1, features = C5[13:16])


#### C6 (up) :: corr met Endo+Pericytes ----


DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
#DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248","HEYL","CFH") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))

ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)




RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:8])
FeaturePlot(object = object_1, features = C6[9:16])
FeaturePlot(object = object_1, features = C6[17:24])
FeaturePlot(object = object_1, features = C6[25:33])

FeaturePlot(object = object_1, features = c("FBN1","COL1A1","RGS5", "CD248"),pt.size = 0.01 * 20)


#### CC-2022 (up) ----


DotPlot(object = object_1, features =list('C1'=
                                            results.out |> dplyr::filter(C1.2022) |> dplyr::pull(hugo_symbol)
                                            , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C1] in: ",sid))
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Couturier/2022-",sid,"_CC.png"),width=7.5, height=4,scale=1.2)




## BT400 [95-100% tumor, low res] ----

rm(sid, object_1)
gc()

sid <- "BT400.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 15000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 15000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")

#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor

#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = "CD14") # TAM/mg

#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")


## BT402 [95-100% tumor] ----


rm(sid, object_1)
gc()

sid <- "BT402.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4750,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4700 &
                     nCount_RNA > 500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




## BT407 [95-100% tumor] ----

rm(sid, object_1)
gc()


sid <- "BT407.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 

# 
# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 300,col="red") +
#   geom_hline(yintercept = 4500,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 20000,col="red") # + scale_y_log10()
# 
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 300 &
#                      nFeature_RNA < 4500 &
#                      nCount_RNA > 500 &
#                      nCount_RNA < 20000 &
#                      percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)


object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




## BT409 [95-100% tumor] ----

rm(sid, object_1)
gc()

sid <- "BT409.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 300,col="red") +
#   geom_hline(yintercept = 4500,col="red")
# 
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 20000,col="red") # + scale_y_log10()
# 
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 300 &
#                      nFeature_RNA < 4500 &
#                      nCount_RNA > 500 &
#                      nCount_RNA < 20000 &
#                      percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



#### 1. Tumor (?) ----

FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


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


#### 2. Astrocyte (-) ----


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))





#### 5. Oligodendrocytes (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN")
FeaturePlot(object = object_1, features = "PLP1")




## HFA567 CD133 :: No C6 ----

rm(sid, object_1)
gc()

sid <- "HFA567_cd133.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4500,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




#### 3A. TAM/mg/monocytes ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----

FeaturePlot(object = object_1, features = "HBG1") # Tumor



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


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)



## HFA567 Total :: One of the best HFA samples but no clear C6 ----


rm(sid, object_1)
gc()

sid <- "HFA567_total.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 6000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 25000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1000 &
                     nFeature_RNA < 6000 &
                     nCount_RNA > 1000 &
                     nCount_RNA < 25000 &
                     percent.mito < 0.15)

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

ElbowPlot(object_1, ndims = 45)

### cluster the cells

k = 40

object_1 <- FindNeighbors(object_1, dims = 1:k)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:k)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("OLIG/OPC"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("Immune cells (CD163-)"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(18)$",paste0("Hematopoietic stem cells?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(11)$",paste0("Big cluster w/ some endothelial & pericytes"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.13 <- FindMarkers(object_1, ident.1 = 13) # PDGFRA, OMG, OLIG1+2, S100B
# 9 and 10 show neuro-developmental genes
# tmp.9 <- FindMarkers(object_1, ident.1 = 9) # SCGN, DLX6, DLX6-AS, DLX5, DLX2, DLX1 (last two also relate to EGFRvIII)
# tmp.10 <- FindMarkers(object_1, ident.1 = 10) # SST, NXPH1, PLS3, LHX6 , DLX2, ERBB4, DLX5
# tmp.16 <- FindMarkers(object_1, ident.1 = 16) # PTGER3, KCNIP4, EBF1, CA8, LINGO2, FGF3, PRMT8, LHX1, PCP4

# FeaturePlot(object = object_1, features = "RBFOX2")




#### 1. Tumor (-) ----
#### 2. Astrocyte (+) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor



#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 3C. Hematopoietic stem cells? ----

FeaturePlot(object = object_1, features = "HBG1") # Tumor


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")

FeaturePlot(object = object_1, features = "VIM")
FeaturePlot(object = object_1, features = "S100B")
FeaturePlot(object = object_1, features = "CST3")

FeaturePlot(object = object_1, features = "RRM2")



#### 5. Oligodendrocytes (?) ----

FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN")
FeaturePlot(object = object_1, features = "PLP1")



#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



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


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)




## HFA570 CD133 :: Maybe C6 cluster? ----

rm(sid, object_1)
gc()

sid <- "HFA570_cd133.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1500,col="red") +
  geom_hline(yintercept = 5500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 4500,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1500 &
                     nFeature_RNA < 5500 &
                     nCount_RNA > 4500 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)

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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)


k <- 20
object_1 <- FindNeighbors(object_1, dims = 1:k)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
# head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:k)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


#levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.15 <- FindMarkers(object_1, ident.1 = 15) # SLN,CYP26B1,MS4A8,MCHR1,GDF10,FAM43A,CCDC140,LUM,MYLK,SLITRK6,IRX3,MSX2,LCAT,KCTD12,C5orf38,ASCL2,OPRK1,EMP1,ABCA8,HLA-DRB1,
# head(tmp.15,20)
# FeaturePlot(object = object_1, features = "SLN")


# some hemoglobins

FeaturePlot(object = object_1, features = "HBG1")
FeaturePlot(object = object_1, features = "HBG2")
FeaturePlot(object = object_1, features = "HBA1")
FeaturePlot(object = object_1, features = "HBA2")

#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")

FeaturePlot(object = object_1, features = "PROM1")



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



#### C6 ----

#f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6
# f <- c(C6, c("SLN","CYP26B1","MS4A8","MCHR1","GDF10","FAM43A","CCDC140")) - C15 markers
DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T) 
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])




## HFA570 Total :: No C6 ----


rm(object_1)
gc()

sid <- "HFA570_total.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 15000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 15000 &
                     percent.mito < 0.2)

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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)



object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.14 <- FindMarkers(object_1, ident.1 = 14)
# tmp.18 <- FindMarkers(object_1, ident.1 = 18)
# tmp.19 <- FindMarkers(object_1, ident.1 = 19) # NHLH2, MRLN, EBF3, SP5, LHX1, NGFR


FeaturePlot(object = object_1, features = c("HBG1", "HBG2")) # precursor red blood cells?


#### 1. Tumor (-) ----

#### 2. Astrocyte (-) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "AQP4")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "GJA1")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")



#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")

FeaturePlot(object = object_1, features = "VIM")
FeaturePlot(object = object_1, features = "S100B")
FeaturePlot(object = object_1, features = "CST3")

FeaturePlot(object = object_1, features = "RRM2")



#### 5. Oligodendrocytes (?) ----

FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


DotPlot(object = object_1, features=c("ABCB1","CD34","FLT4","TIE1","ITGA1","RGS5","PDGFRB","CD248"))



#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



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


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
#f <- C6

DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T) 
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])



## HFA571 CD133 :: No C6 ----

sid <- "HFA571_cd133.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1250,col="red") +
  geom_hline(yintercept = 4000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 2500,col="red") +
  geom_hline(yintercept = 12400,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1250 &
                     nFeature_RNA < 4000 &
                     nCount_RNA > 2500 &
                     nCount_RNA < 12400 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



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

#f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T) 
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])



## HFA571 Total :: No C6 ----

rm(sid, object_1)
gc()


sid <- "HFA571_total.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4000,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 10000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4000 &
                     nCount_RNA > 500 &
                     nCount_RNA < 10000 &
                     percent.mito < 0.2)


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


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)


object_1 <- FindNeighbors(object_1, dims = 1:20)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:20)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.9 <- FindMarkers(object_1, ident.1 = 9) # EOMES,PPP1R17,TMEM158,PENK,NEUROD4,NHLH1,CA12,RASGRP1,SSTR2,NEUROG1,XXbac-BPG32J3.19,NRN1,MFAP4,HES6,ZDHHC22,TTYH2,PRDX1,PAX6,SMOC1,CORO1C



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

#f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T) 
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])





## NSC1 CD133 :: C6 ! ----


rm(sid, object_1)
gc()

sid <- "NSC1_cd133.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 300,col="red") +
  geom_hline(yintercept = 4500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 15000,col="red") # + scale_y_log10()


object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 300 &
                     nFeature_RNA < 4500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 15000 &
                     percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- gsub("_cd","-CD",gsub("\\..+$","",sid))
object_1$dataset <- as.character(object_1$state)



top10 <- head(VariableFeatures(object_1), 10)


# plot variable features with and without labels

plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))     
plot1
plot2


#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)


object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 23
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



# tmp.10 <- FindMarkers(object_1, ident.1 = 10)
head(tmp.10,20) # CDC6,E2F1,TK1,MCM3,MCM6,DTL,KIAA0101,RMI2,DSN1,LOXL1,PCNA,ZNF367,GINS2,CHAF1A,RAD51AP1,MCM4,MCM5,GPX3,HES4,CDC45,
# tmp.11 <- FindMarkers(object_1, ident.1 = 11)
head(tmp.11,20) # BCAM


### Prepare for integration ----

#object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- 'NSC1 CD133+'


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=8, linetype="dashed", color = "red") +
  geom_vline(xintercept=11, linetype="dashed", color = "red") +
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  geom_hline(yintercept=4, linetype="dashed", color = "red") +
  geom_vline(xintercept=6, linetype="dashed", color = "blue") +
  geom_vline(xintercept=8, linetype="dashed", color = "blue") +
  geom_hline(yintercept=0, linetype="dashed", color = "blue") +
  geom_hline(yintercept=2, linetype="dashed", color = "blue")


object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 8 & 
    object_1@reductions$umap@cell.embeddings[,1] <= 11 & 
    object_1@reductions$umap@cell.embeddings[,2] >= 2 & 
    object_1@reductions$umap@cell.embeddings[,2] <= 4,
  "Pericytes", object_1$youri_clusters)

object_1$youri_clusters <- ifelse(
  object_1@reductions$umap@cell.embeddings[,1] >= 6 & 
    object_1@reductions$umap@cell.embeddings[,1] <= 8 & 
    object_1@reductions$umap@cell.embeddings[,2] >= 0 & 
    object_1@reductions$umap@cell.embeddings[,2] <= 2,
  "TAM/MG", object_1$youri_clusters)



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")

object_1.NSC1.CD133 <- object_1



#### 2. Astrocyte (+) ----


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")

DotPlot(object = object_1, features = c("GPR98","AQP4","BMPR1B","ETNPPL","GJB6","GJA1","FGFR3","SLC25A18","SLC1A2","SDC4","GFAP","EDNRB","RNF219-AS1","LINC00499","ALDH1L1","CHI3L1","CLDN10","AGT","SLCO1C1","SLC4A4","GPAM","SLC14A1")
        , group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#### 3A. TAM/mg/monocytes (?)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")


#### 4. Neurons (-) ----


FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
#FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")

FeaturePlot(object = object_1, features = "VIM")
FeaturePlot(object = object_1, features = "S100B")
FeaturePlot(object = object_1, features = "CST3")

FeaturePlot(object = object_1, features = "RRM2")



DotPlot(object = object_1, features = c("RBFOX3", "DLX6-AS1","SYNPR","RELN","CNR1","GAD2","OPRK1","GABRB2","RAB3C","SYT1","KCNC2","ZMAT4","RIMBP2","CHGB","GABRA1","MYT1L","GAD1","PTHLH","TAC3","OLFM3","RP11-679B19.1","ARL4C","TAC1","GABBR2","SCG2","ZCCHC12","GALNTL6"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




#### 5. Oligodendrocytes (?) ----

FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


DotPlot(object = object_1, features=c("ABCB1","CD34","FLT4","TIE1","ITGA1","RGS5","PDGFRB","CD248"),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### 6A. Endothelial (?) ----


FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")

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

#### C6 (up) :: in 10 en 11! ----

f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
#f <- c(C6 , c("ABCB1","CD34","FLT4","TIE1","ITGA1") )

f <- C6
DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])

## Combined ----


object_c <- merge(
  object_1.BT338,
  y=c(
  object_1.BT363,
  object_1.BT364,
  object_1.BT397
  #,object_1.NSC1.CD133
  ),
  add.cell.ids = c("BT338","BT363","BT364","BT397"
                   #,"NSC1_CD133"
                   ))
  

levels(as.factor(object_c$dataset))
object_c$dataset.short <- gsub("^([^-]+).+$","\\1",object_c$dataset)


object_c$seurat_clusters <- NULL
object_c <- NormalizeData(object = object_c, normalization.method = "LogNormalize", scale.factor = 1e4)
object_c <- FindVariableFeatures(object = object_c, selection.method = "vst", nfeatures = 2000)


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


d <- 26



object_c <- FindNeighbors(object_c, dims = 1:d)
object_c <- FindClusters(object_c, resolution = 1.2, algorithm=1)
object_c <- RunUMAP(object_c, dims = 1:d)
object_c@meta.data$pt = sapply(strsplit(rownames(object_c@meta.data), "[.]"), "[", 1)

#### name clusters ----

object_c <- FindClusters(object_c, resolution = 1.2, algorithm=1)
DimPlot(object_c, reduction = "umap", label = T, pt.size = .6, group.by = "seurat_clusters")


# levels(object_c$seurat_clusters) <- gsub("^(15)$","\\1. PE",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(30)$","\\1. EN",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(28)$","\\1. TC",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(12|21)$","\\1. OD",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(2|4|8|32|24|26|11|22|13)$","\\1. TAM/MG",levels(object_c$seurat_clusters))
# 
# levels(object_c$seurat_clusters) <- gsub("^(0|1|7|9|6|16|26|27)$","\\1. T BT364",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(20|29)$","\\1. T BT397",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(25|17|14)$","\\1. T BT338",levels(object_c$seurat_clusters))
# levels(object_c$seurat_clusters) <- gsub("^(3|5|10|18|19|23)$","\\1. T BT364",levels(object_c$seurat_clusters))

levels(object_c$seurat_clusters) <- gsub("^(25)$","\\1. TC",levels(object_c$seurat_clusters))

levels(object_c$seurat_clusters) <- gsub("^(15)$","\\1. PE",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(30)$","\\1. EN",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(12|21)$","\\1. OD",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(2|4|8|32|24|26|11|22|13)$","\\1. TAM/MG",levels(object_c$seurat_clusters))

levels(object_c$seurat_clusters) <- gsub("^(0|1|7|9|6|16|26|27)$","\\1. T BT364",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(20|29)$","\\1. T BT397",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(25|17|14)$","\\1. T BT338",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(3|5|10|18|19|23)$","\\1. T BT364",levels(object_c$seurat_clusters))

levels(object_c$seurat_clusters)

object_c$seurat_clusters <- factor(object_c$seurat_clusters, levels=c(
  "0. T BT364","1. T BT364","3. T BT364","5. T BT364","6. T BT364","7. T BT364","9. T BT364" ,"10. T BT364","18. T BT364","19. T BT364","16. T BT364","27. T BT364",
  "14. T BT338","17. T BT338","25. T BT338",
  "20. T BT397","29. T BT397",
  "23. T BT364",
  "12. OD","21. OD",
  "15. PE",
  "30. EN",
  "2. TAM/MG","4. TAM/MG","8. TAM/MG","11. TAM/MG","13. TAM/MG","22. TAM/MG","24. TAM/MG","26. TAM/MG",
  "28. TC"
))

DimPlot(object_c, reduction = "umap", label = F, pt.size = .8, group.by = "youri_clusters")


p1 <- DimPlot(object_c, reduction = "umap", label = F, pt.size = .6, group.by = "seurat_clusters",
              cols=c(
                rep("cornflowerblue",12),rep("darkblue",3),rep("#88CCEE",2),rep("#6699CC",1),
                rep("orange",2),rep("seagreen",1),rep("mediumvioletred",1),rep("tan3",8),rep("magenta",1))
              ) +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle="Yuan dataset")

p2 <- DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "dataset.short")

p1 + p2


# library(patchwork)



DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")  + DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "dataset.short")



#### 6A. Endothelial (?) ----


FeaturePlot(object = object_c, features = "ABCB1")
FeaturePlot(object = object_c, features = "CD34")
FeaturePlot(object = object_c, features = "FLT4")
FeaturePlot(object = object_c, features = "TIE1") # meh
FeaturePlot(object = object_c, features = "ITGA1") # endo + peri?


#### 6B. Pericytes (+) ----

FeaturePlot(object = object_c, features = "RGS5")
FeaturePlot(object = object_c, features = "PDGFRB")
FeaturePlot(object = object_c, features = "CD248")


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


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )


f <- C6
DotPlot(object_c, features=f, group.by = "seurat_clusters") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

VlnPlot(object = object_c, features = c(C6), group.by = "seurat_clusters",stack=T)

RidgePlot(object = object_c, features = c(C6), group.by = "seurat_clusters",stack=T)

VlnPlot(object = object_c, features = c(C6), group.by = "youri_clusters",stack=T)
RidgePlot(object = object_c, features = c(C6), group.by = "youri_clusters",stack=T)


FeaturePlot(object = object_c, features = C6)

FeaturePlot(object = object_c, features =  "CRABP2" )
FeaturePlot(object = object_c, features =  "CILP2" )
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


