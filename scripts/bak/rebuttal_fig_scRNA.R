f#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(AnnotationHub)
library(ensembldb)


# load data ----


source('scripts/R/subtype_genes.R')



# scRNA HGG levi Glimmunology ----

## A :: Sample_Y ----

rm(sid, object_1)
gc()


sid <- 'van_Hijfte_Sample_Y'
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


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))



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

# levels(object_1$seurat_clusters) <- gsub("^21$","\\1.Neuron",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(5|10|13|14)$","Immune + T-cells.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(4|6|12|19)$","Oligodendrocytes.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^16$","Endothelial",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^15$","Pericytes",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|9|8|11)$","Tumor.\\1",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^7$","Tumor/Dividing",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^17$","Tumor outlier",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^20$","Astrocyte",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^18$","Tumor/Apoptosis?",levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^22$","Immune cell / OD hybrid?",levels(object_1$seurat_clusters))

object_1 <- FindClusters(object_1, resolution = 0.8, algorithm=1)
object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(21), "NE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(5,10,13,14), "TAM/MG", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(4,6,12,19), "OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(16), "EN", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(15), "PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0,1,2,3,9,8,11), "T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(7), "T", object_1$class) # dividing
object_1$class <- ifelse(object_1$seurat_clusters %in% c(17), "T ?", object_1$class) #  Outlier?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(20), "AC", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(18), "T ?", object_1$class) # Apoptotic?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(22), "TAM/MG|OD", object_1$class)
object_1$class <- ifelse(object_1@reductions$umap@cell.embeddings[,1] >= 10 &
                           object_1@reductions$umap@cell.embeddings[,1] <= 11 &
                           object_1@reductions$umap@cell.embeddings[,2] >= 1.5 &
                           object_1@reductions$umap@cell.embeddings[,2] <= 3,
                         "TC", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters),". ",object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","1. T","2. T","3. T","7. T","8. T","9. T","11. T",
  "17. T ?","18. T ?",
  "20. AC",
  "21. NE",
  "4. OD","6. OD","12. OD","19. OD",
  "16. EN",
  "15. PE",
  "22. TAM/MG|OD" ,
  "5. TAM/MG","10. TAM/MG","13. TAM/MG","14. TAM/MG",
  "14. TC"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")  +
  labs(subtitle=sid) +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3)))






#### 1. Wang/MES ----

DotPlot(object = object_1, features = list('MES (Wang)'=subtype.mesenchymal$symbol,
                                           'pericyte'=c('PDGFRB','CD248','RGS5')
                                           
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


#### 2. Wang/CL ----

DotPlot(object = object_1, features = list('CL (Wang)'=subtype.classical$symbol,
                                           'pericyte'=c('PDGFRB','CD248','RGS5')
                                           
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


#### 3. Wang/PN ----

DotPlot(object = object_1, features = list('PN (Wang)'=subtype.proneural$symbol,
                                           'pericyte'=c('PDGFRB','CD248','RGS5')
                                           
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


#### 1. Neftel/MES1 ----

neftel.mes1 <- c('CHI3L1','ANXA2','ANXA1','CD44','VIM','MT2A','C1S','NAMPT','EFEMP1','C1R','SOD2','IFITM3','TIMP1','SPP1','A2M','S100A11','MT1X','S100A10','FN1','LGALS1','S100A16','CLIC1','MGST1','RCAN1','TAGLN2','NPC2','SERPING1','C8orf4','EMP1','APOE','CTSB','C3','LGALS3','MT1E','EMP3','SERPINA3','ACTN1','PRDX6','IGFBP7','SERPINE1','PLP2','MGP','CLIC4','GFPT2','GSN','NNMT','TUBA1C','GJA1','TNFRSF1A','WWTR1')

DotPlot(object = object_1, features = list('MES1 (Neftel)' = neftel.mes1,
                                           'PE'=c('PDGFRB','CD248','RGS5'),
                                           'TAM'=c('CD163','CD14'),
                                           'NE'=c('RBFOX3'),
                                           'OD'=c('TMEM144'),
                                           'TC'=c('CD2', 'TRAC')
                                           
                                           
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


#### 2. Neftel/MES2 ----

neftel.mes2 <- c('HILPDA','ADM','DDIT3','NDRG1','HERPUD1','DNAJB9','TRIB3','ENO2','AKAP12','SQSTM1','MT1X','ATF3','NAMPT','NRN1','SLC2A1','BNIP3','LGALS3','INSIG2','IGFBP3','PPP1R15A','VIM','PLOD2','GBE1','SLC2A3','FTL','WARS','ERO1L','XPOT','HSPA5','GDF15','ANXA2','EPAS1','LDHA','P4HA1','SERTAD1','PFKP','PGK1','EGLN3','SLC6A6','CA9','BNIP3L','RPL21','TRAM1','UFM1','ASNS','GOLT1B','ANGPTL4','SLC39A14','CDKN1A','HSPA9')

DotPlot(object = object_1, features = list('MES1 (Neftel)' = neftel.mes2,
                                           'PE'=c('PDGFRB','CD248','RGS5'),
                                           'TAM'=c('CD163','CD14'),
                                           'NE'=c('RBFOX3'),
                                           'OD'=c('TMEM144'),
                                           'TC'=c('CD2', 'TRAC')
                                           
                                           
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))




FeaturePlot(object = object_1, features =  "LOX" )
FeaturePlot(object = object_1, features =  "SERPINE1" )
FeaturePlot(object = object_1, features =  "ANXA2" )
FeaturePlot(object = object_1, features =  "TGFBI" )


# 


object_1 <- RunPCA(object_1, 
                   features = c(subtype.classical$symbol),
                   reduction.name = 'pca.subtype.cl',
                   reduction.key = 'PCcl')


object_1 <- RunPCA(object_1, 
                   features = c(subtype.mesenchymal$symbol),
                   reduction.name = 'pca.subtype.mes',
                   reduction.key = 'PCmes')


object_1 <- RunPCA(object_1, 
                   features = c(subtype.proneural$symbol),
                   reduction.name = 'pca.subtype.pn',
                   reduction.key = 'PCpn')

#DimPlot(object_2, reduction = "pca.subtype", label = TRUE, pt.size = .8, group.by = "seurat_clusters")



object_2 <- RunPCA(object_2, 
                   features = c(subtype.classical$symbol),
                   reduction.name = 'pca.subtype.cl',
                   reduction.key = 'PCcl')


object_2 <- RunPCA(object_2, 
                   features = c(subtype.mesenchymal$symbol),
                   reduction.name = 'pca.subtype.mes',
                   reduction.key = 'PCmes')


object_2 <- RunPCA(object_2, 
                   features = c(subtype.proneural$symbol),
                   reduction.name = 'pca.subtype.pn',
                   reduction.key = 'PCpn')

#DimPlot(object_2, reduction = "pca.subtype", label = TRUE, pt.size = .8, group.by = "seurat_clusters")






#### combined T/O ----

DotPlot(object = object_2, features = list('MES (Wang)'=subtype.mesenchymal$symbol,
                                           'CL (Wang)'=subtype.classical$symbol,
                                           'PN (Wang)'=subtype.proneural$symbol
                                           
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


#### 2. Wang/CL ----


FeaturePlot(object = object_1, features = c("PCcl_1","PCmes_1","PCpn_1"))

FeaturePlot(object = object_2, features = c("PCcl_1","PCmes_1","PCpn_1"))

