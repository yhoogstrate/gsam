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

write.table(annotations.sub, file = "data/scRNA/GSE103224/annotation_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote=F)


# GSE103224 :: Yuan J. et. al ----

if(file.exists("tmp/GSE103224.scRNA.counts.Rds")) {
  GSE103224 <- readRDS("tmp/GSE103224.scRNA.counts.Rds")
} else {
  # diagnosis: Glioblastoma, WHO grade IV - idh1 status: R132H
  #a <- read.delim("data/scRNA/GSE103224/GSM2758471_PJ016.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  #`colnames<-`(c("ENS", "HGNC", paste0("GSM2758471.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-ampli
  b <- read.delim("data/scRNA/GSE103224/GSM2758472_PJ017.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758472_PJ017.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-non ampli
  c <- read.delim("data/scRNA/GSE103224/GSM2758473_PJ018.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758473_PJ018.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, idh1 status: wt, EGFR-ampli
  d <- read.delim("data/scRNA/GSE103224/GSM2758474_PJ025.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758474_PJ025.cell.",(1:(ncol(.)-2))+2) ))
  
  # Anaplastic Astrocytoma, WHO grade III, IDH-wt, EGFR-non ampli
  #e <- read.delim("data/scRNA/GSE103224/GSM2758475_PJ030.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  #`colnames<-`(c("ENS", "HGNC", paste0("GSM2758475.cell.",(1:(ncol(.)-2))+2) ))
  
  # GBM RECURRENT - IDH-wt, EGFR ampli
  f <- read.delim("data/scRNA/GSE103224/GSM2758476_PJ032.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758476_PJ032.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, recurrent, dh1 status: wt, egfr status: amplified in initial resection
  g <- read.delim("data/scRNA/GSE103224/GSM2758477_PJ035.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758477_PJ035.cell.",(1:(ncol(.)-2))+2) ))
  
  # Glioblastoma, WHO grade IV, IDH-wt, EGFR-non ampli
  h <- read.delim("data/scRNA/GSE103224/GSM2940098_PJ048.filtered.matrix.txt", stringsAsFactors = F,header=F) %>% 
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
  
  write10xCounts(path        = paste0("/home/youri/projects/gsam/data/scRNA/GSE103224/", sample),
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

sid <- 'GSM2758472_PJ017'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
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

ElbowPlot(object_1)

### cluster the cells ----

object_1 <- FindNeighbors(object_1, dims = 1:25)
object_1 <- FindClusters(object_1, resolution = 1)
head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:25)
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


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


dge <- FindMarkers(object_1, ident.1 = '9')
head(dge,15)

rm(dge)


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

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


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

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


FeaturePlot(object = object_1, features = C5) # LIF PI3 MMP7 tumor intrinsiek


#### C6 (up) ----


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

sid <- 'GSM2758473_PJ018'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
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

### scaling of data ----
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(object_1)
object_1 <- ScaleData(object_1, features = all.genes)

### PCA plot ----

object_1 <- RunPCA(object_1, features = VariableFeatures(object = object_1))
print(object_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object_1, dims = 1:2, reduction = "pca")
DimPlot(object_1, reduction = "pca")

#### estimation of the number of principle components in your dataset

ElbowPlot(object_1)

### cluster the cells ----

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

sid <- 'GSM2758474_PJ025'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
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

ElbowPlot(object_1)

### cluster the cells ----

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

tmp2 = FindMarkers(object = object_1, ident.1 = 10)


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

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


FeaturePlot(object = object_1, features = C6)

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
# GBM RECURRENT - IDH-wt, EGFR ampli

sid <- "GSM2758476_PJ032"
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
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

ElbowPlot(object_1)

### cluster the cells ----

object_1 <- FindNeighbors(object_1, dims = 1:25)
object_1 <- FindClusters(object_1, resolution = 1)
head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:25)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


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

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


FeaturePlot(object = object_1, features = C6)

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


sid <- 'GSM2758477_PJ035'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
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

ElbowPlot(object_1)

### cluster the cells ----


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


# ids <- read.table(paste0('data/scRNA/GSE103224/',sid,'/barcodes.tsv'))
# write.table(object_1$seurat_clusters %>%
#               data.frame(stringsAsFactors = F) %>%
#               tibble::rownames_to_column('cell') %>%
#               dplyr::rename(cluster = '.') %>% 
#               dplyr::mutate(cluster = paste0('SeuratCluster',cluster)) ,
#             file = paste0('data/scRNA/GSE103224/',sid,'/seurat_cluster_annotations.txt'),
#             row.names = F, col.names = F, sep = "\t",quote = FALSE)
# 


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

#### C3 - C6 ----


DotPlot(object = object_1, features = c(C6), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### C3 (down) :: endothelial ----

DotPlot(object = object_1, features = c(C3), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C3)


#### C4 (up) ----

DotPlot(object = object_1, features = c(C3), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(object = object_1, features = c(C4A,C4B), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C4A,C4B), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C4A,C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----

DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


DotPlot(object = object_1, features = c(C6), group.by = "seurat_clusters")
RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6)

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
#   file = paste0("data/scRNA/GSE103224/",sid,'/matrix.only-included-cells.txt'),
#   row.names = T, col.names = T, sep = "\t",quote = FALSE )


infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=paste0("data/scRNA/GSE103224/",sid,'/matrix.only-included-cells.txt'),
                                     annotations_file=paste0("data/scRNA/GSE103224/",sid,'/seurat_cluster_annotations.txt'),
                                     gene_order_file='data/scRNA/GSE103224/annotation_genes.txt',
                                     #annotations_file=system.file("extdata", "data/scRNA/GSE103224/GSM2758477_PJ035/seurat_cluster_annotations.txt", package = "infercnv"),
                                     #gene_order_file=system.file("extdata", "data/scRNA/GSE103224/annotation_genes.txt", package = "infercnv"),
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

sid <- 'GSM2940098_PJ048'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
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

ElbowPlot(object_1)

### cluster the cells ----

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


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")



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

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


DotPlot(object = object_1, features = c(C6), group.by = "seurat_clusters")
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




# scRNA HGG levi ----

## A :: Sample_Y ----

sid <- 'GSM2758472_PJ017'
object_1 <- Read10X(data.dir = "/home/youri/projects/gsam/data/scRNA_glim/Glioma_Y_and_O/Levi2_Glioma_Y/outs/raw_feature_bc_matrix")
object_1 <- CreateSeuratObject(counts = object_1,
                               min.cells = 3,
                               min.features = 200,
                               project = "glioma_glim")
mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 

    
ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1400,col="red") +
  geom_hline(yintercept = 4500,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 2200,col="red") +
  geom_hline(yintercept = 14000,col="red")
  #scale_y_log10()

object_1 <- subset(x = object_1, subset = 
                     nFeature_RNA > 1400 &
                     nFeature_RNA < 4500 & 
                     nCount_RNA > 2200 &
                     nCount_RNA < 14000 &
                     percent.mito < 0.025)

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
object_1 <- FindClusters(object_1, resolution = 0.8, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:40)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")

levels(object_1$seurat_clusters) <- gsub("^21$","Neuron",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(5|10|13|14)$","Immune + T-cells.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(4|6|12|19)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^16$","Endothelial",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^15$","Pericytes",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|9|8|11)$","Tumor.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^7$","Tumor/Dividing",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^17$","Tumor outlier",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^20$","Astrocyte",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^18$","Tumor/Apoptosis?",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^22$","Immune cell / OD hybrid?",levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


tmp.17 <- FindMarkers(object_1, ident.1 = 17)
head(tmp.17,20)

tmp.22 <- FindMarkers(object_1, ident.1 = 22)
head(tmp.22,20)






#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES


FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?

#### 2. Astrocyte (+) ----

FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor


# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


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

# NPC1
FeaturePlot(object = object_1, features = c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST")) # Tumor

# NPC2
FeaturePlot(object = object_1, features = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4")) # Tumor

# NPC1, NPC2, Neuron
FeaturePlot(object = object_1, features = c("BCAN", "NREP", "RBFOX3")) # Tumor
# BCAN


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


#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")


#### 4. Neurons (+) ----

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



#### 5. Oligodendrocytes (+) ----

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


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

#### C3 - C6 ----


DotPlot(object = object_1, features = c(C3, C4A, C4B, C5, C6), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#### C3 (down) :: endothelial ----

DotPlot(object = object_1, features = c(C3), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C3)


#### C4 (up) ----

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----

VlnPlot(object = object_1, features = C6, stack = T, sort = T)

DotPlot(object = object_1, features = C6) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(object = object_1, features = C6, stack = T, sort = T)

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




# GSE117891 ----

# GSE131928 :: Neftel ----

## convert to 10X/seurat ----

tmp <- read.table("data/scRNA/GSE131928_Neftel/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv",
                  header=T, 
                  stringsAsFactors=F) %>%
  tibble::column_to_rownames('GENE')

tmp <- tmp %>% dplyr::filter(rowSums(.) > 0)

tmp.genes  <- rownames(tmp)
tmp.cell.ids <- colnames(tmp)




# export each
for(sample in unique(gsub("^X([0-9]+).*$","\\1",tmp.cell.ids ))) {
  print(paste0("Exporting ",sample))
  
  tmp.subset <- tmp %>%
    as.data.frame() %>%
    dplyr::select(contains(sample))
  print(dim(tmp))
  
  tmp.cell.ids.subset <- colnames(tmp.subset)
  
  write10xCounts(path        = paste0("data/scRNA/GSE131928_Neftel/", sample),
                 x           = tmp.subset %>% as.matrix() %>% as("dgTMatrix"),
                 barcodes    = tmp.cell.ids.subset , 
#                 gene.id     = tmp.genes,
                 gene.symbol = tmp.genes)
  
  rm(GSE103224.subset, GSE103224.cell.ids.subset, sample)
  gc()
}



# load table and make 10X-like files

object_1 <- Read10X(data.dir = "/home/youri/projects/gsam/data/scRNA/GSE131928_Neftel/10x_obj_2/DropUtils")

## 102 :: ----
## 105 :: ----
## 114 :: ----
## 115 :: ----
## 118 :: ----
## 124 :: ----
## 125 :: ----
## 126 :: ----
## 143 :: ----



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
levels(object_1$seurat_clusters) <- gsub("^(15)$",paste0("Monocyte?"),levels(object_1$seurat_clusters))

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

rm(object_1)
gc()

sid <- "BT338_1of2.filtered_gene_matrices"
object_1 <- Read10X(data.dir = paste0("data/scRNA/EGAS00001004422_Couturier/filtered/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

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

object_1 <- FindNeighbors(object_1, dims = 1:35)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----


object_1 <- RunUMAP(object_1, dims = 1:35)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(1|4|11)$",paste0("Endothelial\n+Pericytes.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("TAM/microglia\n+T-cell"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(8)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))

levels(object_1$seurat_clusters) <- gsub("^(0|2|3|5|6|7|9|10|12)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



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
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M


#### 2. Astrocyte (-) ----


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor


#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (+) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")


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



#### 5. Oligodendrocytes (+) ----

FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


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

#### C3 - C6 ----


DotPlot(object = object_1, features = c(C3, C4A, C4B, C5, C6), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#### C3 (down) :: endothelial ----

DotPlot(object = object_1, features = c(C3), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C3)


#### C4 (up) ----

DotPlot(object = object_1, features = c(C4A,C4B))

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


DotPlot(object = object_1, features = C5)
FeaturePlot(object = object_1, features = C5)


#### C6 (up) ----


DotPlot(object = object_1, features = C6) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(object = object_1, features = C6, stack = T, sort = T)
VlnPlot(object = object_1, features = C6, stack = T, sort = T)

FeaturePlot(object = object_1, features = C6)


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
# levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("Pericytes.\\1"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


tmp.13 <- FindMarkers(object_1, ident.1 = 13)
head(tmp.13, 20) # LAMP3,IRF4,NCCRP1,CRIP1,SYNPO2,CCR7,EHF,CCL22,VTN,LSP1,CDX2,IDO1,TRAF3IP3,RP11-290F5.1,ANXA3,RP1-28O10.1,RP11-66B24.4,AP003774.1,NR4A3,CCL17
FeaturePlot(object = object_1, features = "LAMP3")


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


#### 3C. ? Leukocyte ?? ----


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


DotPlot(object = object_1, features = c(C4A, C4B))


VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)



#### C5 (down) ----


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters")

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)



#### C6 (up) ----


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)




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


levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|8|10|11|12|14)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(7|13)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(6|9)$","Tam/MG.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^15$","Endothelial",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^16$","Pericytes",levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")






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
FeaturePlot(object = object_1, features = "ELAVL4")# ? https://www.biorxiv.org/content/10.1101/775197v1



#### 5. Oligodendrocytes (+) ----

FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


#### 6A. Endothelial (-) ----

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



#### C6 (up) :: some cor w/ pericytes ----



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


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


## BT364 [1+2/2] :: T,MG,OD,AC(!) ----


rm(object_1)
gc()

object_1 <- Read10X(data.dir = "data/scRNA/EGAS00001004422_Couturier/filtered/BT364_1of2.filtered_gene_matrices/")
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

object_1.tmp <- Read10X(data.dir = "data/scRNA/EGAS00001004422_Couturier/filtered/BT364_2of2.filtered_gene_matrices/")
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


levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|9|12|13)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^15$","Tumor\nimmune cell recruiting?",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(10|14)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(6|11|8)$","TAM/MG\\1",levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters")



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


FeaturePlot(object = object_1, features = "STMN2") # Tumor
FeaturePlot(object = object_1, features = "ETNPPL") # Tumor



#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (1 of 2?) ----

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


#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
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



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6)



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

levels(object_1$seurat_clusters) <- gsub("^(1|4|9|8|6)$",paste0("TAM/microglia.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|11|2|10|12)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("T-Cells"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.3.5.7 <- FindMarkers(object_1, ident.1 = c(3,5,7))
head(tmp.3.5.7, 20) # STMN1


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


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (10?) ----

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

#### C3 (up) ----

f <- C3
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### C4 (up) ----

f <- c(C4A, C4B)

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)



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

levels(object_1$seurat_clusters) <- gsub("^(0|1|11|7|14|2|5|8|6|10|15|3|9)$",paste0("TAM/microglia.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13|4|17)$",paste0("Tumor.\\1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("T-cells"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("Oligodendrocytes"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(18)$",paste0("Endothelial\n+ Pericytes"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




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


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")



#### 3A. TAM/mg/monocytes (+)----

FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (10?) ----

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


#### C3 (up) ----

f <- C3
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### C4 (up) ----

f <- c(C4A, C4B)
DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----

DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: corr met Endo+Pericytes ----


f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = f, group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:8])
FeaturePlot(object = object_1, features = C6[9:16])
FeaturePlot(object = object_1, features = C6[17:24])
FeaturePlot(object = object_1, features = C6[25:33])

FeaturePlot(object = object_1, features = c("FBN1","COL1A1","RGS5", "CD248"),pt.size = 0.01 * 20)



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
tmp.9 <- FindMarkers(object_1, ident.1 = 9) # SCGN, DLX6, DLX6-AS, DLX5, DLX2, DLX1 (last two also relate to EGFRvIII)
tmp.10 <- FindMarkers(object_1, ident.1 = 10) # SST, NXPH1, PLS3, LHX6 , DLX2, ERBB4, DLX5

tmp.16 <- FindMarkers(object_1, ident.1 = 16) # PTGER3, KCNIP4, EBF1, CA8, LINGO2, FGF3, PRMT8, LHX1, PCP4


FeaturePlot(object = object_1, features = "RBFOX2")



#### hemoglobin ----

FeaturePlot(object = object_1, features = "HBG1") # Tumor


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


#### 3B. Til/T-cell (1 of 2?) ----

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

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)

DotPlot(object = object_1, features = c(C4A, C4B))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)



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


tmp.15 <- FindMarkers(object_1, ident.1 = 15) # SLN,CYP26B1,MS4A8,MCHR1,GDF10,FAM43A,CCDC140,LUM,MYLK,SLITRK6,IRX3,MSX2,LCAT,KCTD12,C5orf38,ASCL2,OPRK1,EMP1,ABCA8,HLA-DRB1,
head(tmp.15,20)

FeaturePlot(object = object_1, features = "SLN")


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


#### C6 ----

#f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
#f <- C6
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


#levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("TAM/microglia"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


# tmp.14 <- FindMarkers(object_1, ident.1 = 14)
# tmp.18 <- FindMarkers(object_1, ident.1 = 18)
# tmp.19 <- FindMarkers(object_1, ident.1 = 19) # NHLH2, MRLN, EBF3, SP5, LHX1, NGFR


FeaturePlot(object = object_1, features = c("HBG1", "HBG2")) # precursor red blood cells?


#### 1. Tumor (-) ----
#### 2. Astrocyte (+) ----


FeaturePlot(object = object_1, features = "STMN2")
FeaturePlot(object = object_1, features = "ETNPPL")



#### 3A. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (1 of 2?) ----

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

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)

DotPlot(object = object_1, features = c(C4A, C4B))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C5)



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


tmp.9 <- FindMarkers(object_1, ident.1 = 9) # EOMES,PPP1R17,TMEM158,PENK,NEUROD4,NHLH1,CA12,RASGRP1,SSTR2,NEUROG1,XXbac-BPG32J3.19,NRN1,MFAP4,HES6,ZDHHC22,TTYH2,PRDX1,PAX6,SMOC1,CORO1C


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


#### 3B. Til/T-cell (1 of 2?) ----

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



#### C5 (down) ----

f <- C5
DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### C6 (up) :: in 10 en 11! ----

f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
#f <- c(C6 , c("ABCB1","CD34","FLT4","TIE1","ITGA1") )

f <- C6
DotPlot(object = object_1, features = c(f), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


RidgePlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(f), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])



# GSE138794_Diaz ----

## GSM4119521_SF10022 [HQ] ----


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



DotPlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1", "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2")) # other type of endothelial cells?
DotPlot(object = object_1, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2"), group.by = "seurat_clusters") # other type of endothelial cells?



#### 6B. Pericytes (12) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B))
VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----

DotPlot(object = object_1, features = c(C5), group.by = "seurat_clusters")

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) ----



f <- c(C6 , c("RGS5", "PDGFRB", "CD248") )
f <- C6

DotPlot(object = object_1, features = f, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:6])


