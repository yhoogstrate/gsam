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


# GSE103224 ----

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


#### 2. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



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

FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh


#### 5B. Pericytes (3 cells ?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


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


#### 2. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))



#### 3. Neurons (-) ----

FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "GABRG2")

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

FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1")


#### 5B. Pericytes (+ heel dicht met endothelial?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### 6. Til/T-cell (?) ----

# FeaturePlot(object = object_1, features = "CD3D")
# FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "TRBC2")

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



## D :: GSM2758474_PJ025 :: T,MG,EN,PE ----
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


#### 2. TAM/mg/monocytes (+)----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


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

FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1")


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

FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1")


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

FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1")


#### 5B. Pericytes (+) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### 6. other cluster(s) ----

FeaturePlot(object = object_1, features = "COL1A1")
FeaturePlot(object = object_1, features = "VEGFA")
FeaturePlot(object = object_1, features = "NDRG1")




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

#levels(object_1$seurat_clusters) <- gsub("^0$","Tumor/CL~EGFR",levels(object_1$seurat_clusters))
#levels(object_1$seurat_clusters) <- gsub("^4$","Tumor/Neuronal",levels(object_1$seurat_clusters))
#levels(object_1$seurat_clusters) <- gsub("^1$","Tumor/OPC",levels(object_1$seurat_clusters))

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
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


#### 5A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1")


#### 5B. Pericytes (+) ----

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


# til/t-cell
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")


# scRNA HGG levi ----

## A :: H09-1815 ----

sid <- 'GSM2758472_PJ017'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")



# GSE117891 ----

# GSE131928 ----

# further processing ----