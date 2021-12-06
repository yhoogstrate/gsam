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


## make vis [ b: GSM2758472_PJ017 ] ----

object_1 <- Read10X(data.dir = "data/scRNA/GSE103224/GSM2758472_PJ017")
# sum(object_1[rownames(object_1) == 'RBFOX3',]) == 2
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")
# "EGFR" %in% rownames(object_1)

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & nFeature_RNA <10000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
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

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1)
head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "pt")

object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

#### Pericytes? ----
FeaturePlot(object = object_1, features = "RGS5", "PDGFRB", "CD248") # this endo gene defines the EM sig?!


#### TAM/mg ----
FeaturePlot(object = object_1, features = "CD163")
FeaturePlot(object = object_1, features = "CD163") # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg

#### tumor ----
FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor
FeaturePlot(object = object_1, features = "GFAP") # Tumor - ander gedeelte?

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # iedere cel?
FeaturePlot(object = object_1, features = "ANXA1") # weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")


#### Neurons ----
# counts tegen nul aan, te weinig cellen - niet tot nauwelijks, lijkt echt toename van neuronen
FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "DDN")
# FeaturePlot(object = object_1, features = "TNNT2")
# FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")


#### Oligodendrocytes ---
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN") # deze is niet weggefilterd, maar geen cellen hiermee?! - vraag levi
FeaturePlot(object = object_1, features = "PLP1") # wel/ook in tumor
FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "LGI3")
FeaturePlot(object = object_1, features = "PLPP2")

#### Astrocytes ----
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")
FeaturePlot(object = object_1, features = "GFAP")

# til/t-cell
FeaturePlot(object = object_1, features = "CD3E")
FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD8A")
FeaturePlot(object = object_1, features = "C1A1")

# C3
## endo ?
FeaturePlot(object = object_1, features = "TIE1")
FeaturePlot(object = object_1, features = "FLT1")

# C4
FeaturePlot(object = object_1, features = "NDRG1")
FeaturePlot(object = object_1, features = "MICAL2")

# C5 - lijkt niet-tumor
FeaturePlot(object = object_1, features = "PRF1")
FeaturePlot(object = object_1, features = "ARHGAP9") # mg?
# FeaturePlot(object = object_1, features = "FCMR")
FeaturePlot(object = object_1, features = "LXN") # top cluster?
FeaturePlot(object = object_1, features = "KCNE3") # top cluster??
FeaturePlot(object = object_1, features = "FPR2") # mg?
#FeaturePlot(object = object_1, features = "CCL13")
FeaturePlot(object = object_1, features = "MMP7") # tum/astr?
# FeaturePlot(object = object_1, features = "CALCR")
# FeaturePlot(object = object_1, features = "LRG1") mg?
FeaturePlot(object = object_1, features = "SAA2")
FeaturePlot(object = object_1, features = "PIE3")
FeaturePlot(object = object_1, features = "LIF")
FeaturePlot(object = object_1, features = "HSPA6")


# C6 genes
FeaturePlot(object = object_1, features = "COL1A1") # - mini cluster?
FeaturePlot(object = object_1, features = "COL1A2") # - mini cluster, little in tumor
FeaturePlot(object = object_1, features = "COL5A1") 
FeaturePlot(object = object_1, features = "ASS1") # tumor
FeaturePlot(object = object_1, features = "MFAP4") # tumor
FeaturePlot(object = object_1, features = "TWIST2")
#FeaturePlot(object = object_1, features = "TNNT3")
FeaturePlot(object = object_1, features = "COL5A1")
FeaturePlot(object = object_1, features = "COL12A1")
FeaturePlot(object = object_1, features = "MMP11")
FeaturePlot(object = object_1, features = "MFAP2")



## make vis [ c: GSM2758473_PJ018 ] ----

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

object_1 <- FindNeighbors(object_1, dims = 1:17)
object_1 <- FindClusters(object_1, resolution = 1) # 0.5 - 1 ?
#head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:17)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "pt")



#### TAM/mg ----

FeaturePlot(object = object_1, features = c("CD163", "CD14")) # TAM/mg


#### tumor ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

# GFAP en ANXA1 astro markers?
FeaturePlot(object = object_1, features = "GFAP") # Tumor - ander gedeelte?
FeaturePlot(object = object_1, features = "ANXA1") # weinig in GFAP-neg tumor cellen?

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")


#### Neurons (-) ----
# counts tegen nul aan, te weinig cellen - niet tot nauwelijks, lijkt echt toename van neuronen
FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")

#### Oligodendrocytes (+) ---
FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "OPALIN") # deze is niet weggefilterd, maar geen cellen hiermee?! - vraag levi
FeaturePlot(object = object_1, features = "LGI3")
FeaturePlot(object = object_1, features = "MOG")


#### OPC netfel ----
# lijkt wel in tumor te zitten(!)
FeaturePlot(object = object_1, features = "PLP1") # wel/ook in tumor [ is dit niet een Neftel OPC gene?]


#### C2 ~ neuronal ----
FeaturePlot(object = object_1, features = "PLPP2")



#### Astrocytes ----
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")
FeaturePlot(object = object_1, features = "GFAP")



## make vis [ d: GSM2758474_PJ025 ] ----

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

object_1 <- FindNeighbors(object_1, dims = 1:17)
object_1 <- FindClusters(object_1, resolution = 1) # 0.5 - 1 ?
#head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:17)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")


#### TAM/mg (+) ----


FeaturePlot(object = object_1, features = c("CD163", "CD14")) # TAM/mg


#### Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor

# GFAP en ANXA1 astro markers?
FeaturePlot(object = object_1, features = "GFAP") # Tumor - ander gedeelte?
FeaturePlot(object = object_1, features = "ANXA1") # weinig in GFAP-neg tumor cellen?

FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # niet alleen in tumor cellen (!)

FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")


#### Neurons (-) ----
# counts tegen nul aan, te weinig cellen - niet tot nauwelijks, lijkt echt toename van neuronen
FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")

#### Oligodendrocytes (-/?) ---
FeaturePlot(object = object_1, features = "TMEM144")
#FeaturePlot(object = object_1, features = "OPALIN") # deze is niet weggefilterd, maar geen cellen hiermee?! - vraag levi
FeaturePlot(object = object_1, features = "LGI3")
FeaturePlot(object = object_1, features = "MOG")

#### Pericytes? ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### Endothelial ----

FeaturePlot(object = object_1, features = "TIE1")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "CD34")

#### C6 ----

FeaturePlot(object = object_1, features = "DPT")


#### OPC netfel ----
# lijkt wel in tumor te zitten(!)
FeaturePlot(object = object_1, features = "PLP1") # wel/ook in tumor [ is dit niet een Neftel OPC gene?]


#### C2 ~ neuronal ----

FeaturePlot(object = object_1, features = "PLPP2")



#### Astrocytes ----
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")
FeaturePlot(object = object_1, features = "GFAP")






## make vis [ GSM2758476_PJ032 ] ----

object_1 <- Read10X(data.dir = "data/scRNA/GSE103224/GSM2758476_PJ032")
# sum(object_1[rownames(object_1) == 'RBFOX3',]) == 2
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")
# "EGFR" %in% rownames(object_1)

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & nFeature_RNA <10000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
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

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = 1)
head(Idents(object_1), 50)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "pt")

object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


#### TAM/mg ----


FeaturePlot(object = object_1, features = "CD163") # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor
FeaturePlot(object = object_1, features = "GFAP") # Tumor - ander gedeelte?
FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # iedere cel?
FeaturePlot(object = object_1, features = "ANXA1")
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

# C4
FeaturePlot(object = object_1, features = "NDRG1") 
FeaturePlot(object = object_1, features = "SLIT3") 



# neur - counts tegen nul aan, te weinig cellen
# FeaturePlot(object = object_1, features = "RBFOX3")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")

# oligod
# FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")

# astr
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")
FeaturePlot(object = object_1, features = "GFAP")

# til/t-cell
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")
FeaturePlot(object = object_1, features = "GFAP")


#C6 genes
FeaturePlot(object = object_1, features = "COL1A1") # - mini cluster?
FeaturePlot(object = object_1, features = "COL1A2") # - mini cluster?
FeaturePlot(object = object_1, features = "COL5A1") 
FeaturePlot(object = object_1, features = "ASS1") # tumor
FeaturePlot(object = object_1, features = "MFAP4") # tumor


# endo
FeaturePlot(object = object_1, features = "RGS5") # defines the EM sig?!
FeaturePlot(object = object_1, features = "TIE1")
FeaturePlot(object = object_1, features = "FLT1")


## make vis [ GSM2758477_PJ035 ] [ + EM ] ----


sid <- 'GSM2758477_PJ035'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE103224/",sid))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma")

mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.5, group.by = "orig.ident")
object_1 <- subset(x = object_1, subset = nFeature_RNA > 700 & nFeature_RNA <10000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
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

object_1 <- FindNeighbors(object_1, dims = 1:30)
object_1 <- FindClusters(object_1, resolution = .5)
head(Idents(object_1), 50)




### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:30)
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


#### TAM/mg ----


FeaturePlot(object = object_1, features = "CD163") # TAM/mg
FeaturePlot(object = object_1, features = "CD14") # TAM/mg



FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "EGFR") # Tumor
FeaturePlot(object = object_1, features = "GFAP") # Tumor - ander gedeelte?
FeaturePlot(object = object_1, features = "SOCS2") # Tumor?
FeaturePlot(object = object_1, features = "SOX2")
FeaturePlot(object = object_1, features = "VIM") # iedere cel?
FeaturePlot(object = object_1, features = "ANXA1") # weinig in GFAP-neg tumor cellen?
FeaturePlot(object = object_1, features = "PTEN")
FeaturePlot(object = object_1, features = "SETD5")
FeaturePlot(object = object_1, features = "SMAD5")
FeaturePlot(object = object_1, features = "TRIM24")

FeaturePlot(object = object_1, features = "PLP1")

#### Endothelial ----

#FeaturePlot(object = object_1, features = "CD31")
FeaturePlot(object = object_1, features = "CD34")



# neur - counts tegen nul aan, te weinig cellen - niet tot nauwelijks, lijkt echt toename van neuronen
FeaturePlot(object = object_1, features = "RBFOX3")
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "DDN")
# FeaturePlot(object = object_1, features = "TNNT2")
# FeaturePlot(object = object_1, features = "TMEM130")
# FeaturePlot(object = object_1, features = "GABRG2")
# FeaturePlot(object = object_1, features = "GABRA1")
# FeaturePlot(object = object_1, features = "GABRB2")


# oligod
FeaturePlot(object = object_1, features = "MOG")
#DefaultAssay(object = object_1 ) = 'RNA'
FeaturePlot(object = object_1, features = "OPALIN") # deze is niet weggefilterd, maar geen cellen hiermee?! - vraag levi
FeaturePlot(object = object_1, features = "PLP1") # wel/ook in tumor
FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "LGI3")
FeaturePlot(object = object_1, features = "PLPP2")

# astr
FeaturePlot(object = object_1, features = "GFAP")
FeaturePlot(object = object_1, features = "CACHD1")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "GPR37L1")

# til/t-cell
FeaturePlot(object = object_1, features = "CD3E")
FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD8A")
FeaturePlot(object = object_1, features = "C1A1")

# C3/endo
FeaturePlot(object = object_1, features = "TIE1")
FeaturePlot(object = object_1, features = "FLT1")

# C4
FeaturePlot(object = object_1, features = "NDRG1")
FeaturePlot(object = object_1, features = "MICAL2")

# C5 - lijkt niet-tumor
FeaturePlot(object = object_1, features = "PRF1")
FeaturePlot(object = object_1, features = "ARHGAP9") # mg?
# FeaturePlot(object = object_1, features = "FCMR")
FeaturePlot(object = object_1, features = "LXN") # top cluster?
FeaturePlot(object = object_1, features = "KCNE3") # top cluster??
FeaturePlot(object = object_1, features = "FPR2") # mg?
#FeaturePlot(object = object_1, features = "CCL13")
FeaturePlot(object = object_1, features = "MMP7") # tum/astr?
# FeaturePlot(object = object_1, features = "CALCR")
# FeaturePlot(object = object_1, features = "LRG1") mg?
FeaturePlot(object = object_1, features = "SAA2")
FeaturePlot(object = object_1, features = "PIE3")
FeaturePlot(object = object_1, features = "LIF")
FeaturePlot(object = object_1, features = "HSPA6")








# C6 genes
FeaturePlot(object = object_1, features = "COL1A1") # - mini cluster?
FeaturePlot(object = object_1, features = "COL1A2") # - mini cluster, little in tumor
FeaturePlot(object = object_1, features = "COL5A1") 
FeaturePlot(object = object_1, features = "ASS1") # tumor
FeaturePlot(object = object_1, features = "MFAP4") # tumor
FeaturePlot(object = object_1, features = "TWIST2")
#FeaturePlot(object = object_1, features = "TNNT3")
FeaturePlot(object = object_1, features = "COL5A1")
FeaturePlot(object = object_1, features = "COL12A1")
FeaturePlot(object = object_1, features = "MMP11")
FeaturePlot(object = object_1, features = "MFAP2")



### Find DE genes voor SC8 ----

dge <- FindMarkers(object_1, ident.1 = '8')

FeaturePlot(object = object_1, features = "CD248") # defines separate cluster - fibroblast marker
FeaturePlot(object = object_1, features = "RGS5") # defines separate cluster
FeaturePlot(object = object_1, features = "SPP1") # among most defining 



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





# read GSE117891, GSE131928 ----

# read GSE117891, GSE131928 ----

# further processing ----