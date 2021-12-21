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




# GSE138794_Diaz ----

## GSM4119521_SF10022 [HQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 65
# gender: Male
# chr7p+, chr8-, chr10-, chr12+


rm(sid, object_1)
gc()

sid <- 'GSM4119521_SF10022'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 2000,col="red") +
  geom_hline(yintercept = 6500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 4500,col="red") +
  geom_hline(yintercept = 22500,col="red") # + scale_y_log10()





object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 2000 &
                     nFeature_RNA < 6500 &
                     nCount_RNA > 4500 &
                     nCount_RNA < 22500 &
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

d <- 25
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(5)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9|11)$",paste0("TAM/microglia.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(15)$",paste0("Unknown\nEndothelial like?\nTBX3+"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(3|6|7|4|1|2|0)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("Tumor.G2"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("Tumor.G1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("Special endothelial\n or neurons?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("Neurons?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(10)$",paste0("Neurons type 2?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(8)$",paste0("Astrocytes?"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


#tmp.15 <- FindMarkers(object_1, ident.1 = 15)
head(tmp.15, 20) # A4GALT,CTD-2023N9.2,SNTG2,EFCC1,AC016768.1,  TWIST2  ,IL1R1,  KRT17  ,DLK1,MMP11,LGR6,FOXS1,CARMN,ESM1,FOXC2,FOXL1,BOK-AS1,ALDH1A3,SLC6A20,FOXA2,DMKN,IL18R1,SFRP5,LINC00261,  TBX3

# tmp.12 <- FindMarkers(object_1, ident.1 = 12) # AURKB, SKA1, BUB1B
# tmp.13 <- FindMarkers(object_1, ident.1 = 13) # DTL, MELK
# tmp.16 <- FindMarkers(object_1, ident.1 = 16) # SSTR2, SLC17A7, NHLH1, LHX1, CDH22
# tmp.8 <- FindMarkers(object_1, ident.1 = 8)
# tmp.10 <- FindMarkers(object_1, ident.1 = 10) # SUSD5, ETV1
# tmp.14 <- FindMarkers(object_1, ident.1 = 14)
# 
# 


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


#### 2. Astrocyte (+) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor



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


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
                                        ),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### 5. Oligodendrocytes (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


FeaturePlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"))


DotPlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2")) # other type of endothelial cells?
DotPlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2"), group.by = "seurat_clusters") # other type of endothelial cells?



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


#### C5 (down) ::T-cell & TAM/MG? ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) ----



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )


f <- C6
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])




## GSM4119522_SF10127 [LowQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 44
# gender: Female
# chr7+, chr10-


rm(sid, object_1)
gc()


sid <- 'GSM4119522_SF10127'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 200,col="red") +
  geom_hline(yintercept = 6500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 12500,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 200 &
                     nFeature_RNA < 6500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 12500 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)


d <- 25
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(5)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


## GSM4119523_SF11979sn [LowQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 76
# gender: Female
# chr7p+, chr10-, chr19q-


rm(sid, object_1)
gc()

sid <- 'GSM4119523_SF11979sn'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 200,col="red") +
  geom_hline(yintercept = 3750,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 9000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 200 &
                     nFeature_RNA < 3750 &
                     nCount_RNA > 500 &
                     nCount_RNA < 9000 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 10
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(5)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


## GSM4119524_SF12090 [T,OD,MG - LowQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 61
# gender: Male
# chr7+, chr10-


rm(sid, object_1)
gc()

sid <- 'GSM4119524_SF12090'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 200,col="red") +
  geom_hline(yintercept = 9000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 200,col="red") +
  geom_hline(yintercept = 25000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 200 &
                     nFeature_RNA < 9000 &
                     nCount_RNA > 200 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 10
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(5)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



## GSM4119525_SF12264 :: T,MG,TC+,OD [MHQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 51
# gender: Female
# chr6-, chr7+, chr14p-


rm(sid, object_1)
gc()

sid <- 'GSM4119525_SF12264'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 200,col="red") +
  geom_hline(yintercept = 12500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 50000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 200 &
                     nFeature_RNA < 12500 &
                     nCount_RNA > 500 &
                     nCount_RNA < 50000 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 25
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(2|9|12|8|11|4)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(5|0)$",paste0("Oligodendrocytes.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(1|10)$",paste0("TAM/MG\\.1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("T-Cells"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(7)$",paste0("Tumor [NPC]\n+ some endothelials\n+some TAM/MG"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(3|6)$",paste0("Tumor [NPC].\\1\n+ some T-Cells"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



# tmp.3.6.7 <- FindMarkers(object_1, ident.1 = c(3,6,7))
head(tmp.3.6.7, n=20)


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


#### 2. Astrocyte (+) ----


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


#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")


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


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


#### 6A. Endothelial (?) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


FeaturePlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                            "RGS5", "PDGFRB", "CD248"))


DotPlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2")) # other type of endothelial cells?
DotPlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2"), group.by = "seurat_clusters") # other type of endothelial cells?



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

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )


f <- C6
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])




## GSM4119526_SF4400 :: LowQ ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 51
# gender: Female
# chr7+, chr10-



rm(sid, object_1)
gc()

sid <- 'GSM4119526_SF4400'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 700,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 5000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 100 &
                     nFeature_RNA < 700 &
                     nCount_RNA > 100 &
                     nCount_RNA < 5000 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 15
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(5)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = 1.6, group.by = "seurat_clusters")


## GSM4119527_SF4297 :: T,OD,MG+,NE+,AC+ [MQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 39
# gender: Male
# chr3p-, chr7+, chr11q+, chr12-, chr14p+



rm(sid, object_1)
gc()

sid <- 'GSM4119527_SF4297'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 6500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 100 &
                     nFeature_RNA < 6500 &
                     nCount_RNA > 100 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 15
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(1|4|11)$",paste0("TAM/MG.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(8)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9|12)$",paste0("Neurons?\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(6|7)$",paste0("Astrocytes.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("Small Astrocyte subcluster?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|2|3|5|10)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.13 <- FindMarkers(object_1, ident.1 = 13) # LINC00943,AC007325.2,PRRX2,RP11-407A16.3,TPD52L1,GMPR,GPC5-AS1,RP3-395M20.12,NPNT,FXYD1,GJB6,EVC2,PLIN5,OTX1,ADIRF,CGNL1,RNF43,TMPRSS3,PLSCR4,SDC4,RP11-156K13.1,ALDH1A1,SDS,SHROOM3,LINC01515



#### 1. Tumor (+) ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "PRKG2") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


#### 2. Astrocyte (+) ----


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


#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")


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


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes (+) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")


#### 6A. Endothelial (-) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


FeaturePlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                            "RGS5", "PDGFRB", "CD248"))


DotPlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2")) # other type of endothelial cells?
DotPlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2"), group.by = "seurat_clusters") # other type of endothelial cells?



#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) :: (-) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) :: (-) ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )


f <- C6
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])






## GSM4119528_SF6996 ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 40
# gender: Female
# chr6q-, chr7+, chr10-, chr19+



rm(sid, object_1)
gc()

sid <- 'GSM4119528_SF6996'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1500,col="red") +
  geom_hline(yintercept = 10000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 3000,col="red") +
  geom_hline(yintercept = 40000,col="red") + ylim(0,65000)



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 1500 &
                     nFeature_RNA < 10000 &
                     nCount_RNA > 3000 &
                     nCount_RNA < 40000 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 20
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering :: T,OD+,MG [HQ]----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(0)$",paste0("TAM/MG?"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



# tmp.12 <- FindMarkers(object_1, ident.1 = 12) # PLP1, OPALIN, mitochondrial, dubbele kernen?

#### C5 (down) :: TAM? ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)






## GSM4119529_SF9259[S+R] :: T-,OD++,TM ----

#@todo pooling of the 2?

# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 73
# gender: Male
# chr7q+, chr10-, chr13-, chr14-

rm(sid, object_1)
gc()

sid <- 'GSM4119529_SF9259R'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

sid <- 'GSM4119530_SF9259S'
object_1.tmp <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
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
  geom_hline(yintercept = 200,col="red") +
  geom_hline(yintercept = 5000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 18000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 200 &
                     nFeature_RNA < 5000 &
                     nCount_RNA > 500 &
                     nCount_RNA < 18000 &
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

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1, ndims = 45)

d <- 20
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(0|1|7|3|4|5)$",paste0("Oligodendrocytes.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(6|11)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(10)$",paste0("Tumor OPC/BCAN.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(2|8|9)$",paste0("TAM/MG.\\1"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




# tmp.2.8.9 <- FindMarkers(object_1, ident.1 = c(2,8,9))
head(tmp.2.8.9, 25)

FeaturePlot(object = object_1, features = "BCAN") # Tumor



#### 1. Tumor (+) ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES
FeaturePlot(object = object_1, features = "PRKG2") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))



# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


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


#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")


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


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")


#### 6A. Endothelial (-) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


FeaturePlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                            "RGS5", "PDGFRB", "CD248"))


DotPlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2")) # other type of endothelial cells?
DotPlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2"), group.by = "seurat_clusters") # other type of endothelial cells?



#### 6B. Pericytes ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) :: (-) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) :: (-) ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )


f <- C6
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])






