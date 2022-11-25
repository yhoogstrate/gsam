f#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)



# load data ----


source('scripts/load_results.out.R')


OPC <- c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4","LPPR1","PTPRZ1","VCAN","DBI","PMP2","CNP","TNS3","LIMA1","CA10","PCDHGC3","CNTN1","SCD5","P2RX7","CADM2","TTYH1","FGF12","TMEM206","NEU4","FXYD6","RNF13","RTKN","GPM6B","LMF1","ALCAM","PGRMC1","HRASLS","BCAS1","RAB31","PLLP","FABP5","NLGN3","SERINC5","EPB41L2","GPR37L1")
NPC1 <- c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST","ASCL1","BTG2","DCX","NXPH1","HN1","PFN2","SCG3","MYT1","CHD7","GPR56","TUBA1A","PCBP4","ETV1","SHD","TNR","AMOTL2","DBN1","HIP1","ABAT","ELAVL4","LMF1","GRIK2","SERINC5","TSPAN13","ELMO1","GLCCI1","SEZ6L","LRRN1","SEZ6","SOX11")
NPC2 <- c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4","DLX5","SOX4","MAP1B","RBFOX2","IGFBPL1","STMN1","HN1","TMEM161B-AS1","DPYSL3","SEPT3","PKIA","ATP1B1","DYNC1I1","CD200","SNAP25","PAK3","NDRG4","KIF5A","UCHL1","ENO2","KIF5C","DDAH2","TUBB2A","LBH","LOC150568","TCF4","GNG3","NFIB","DPYSL5","CRABP1","DBN1","NFIX","CEP170","BLCAP")



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
object_1@meta.data$pt <- sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)
DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")



object_1 <- FindClusters(object_1, resolution = 0.8, algorithm = 1)
object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(21), "NE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(5, 10, 13, 14), "TAM", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(4, 6, 12, 19), "OD", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(16), "EN", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(15), "PE", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(0, 1, 2, 3, 9, 8, 11), "T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(7), "T", object_1$class) # dividing
object_1$class <- ifelse(object_1$seurat_clusters %in% c(17), "T ?", object_1$class) #  Outlier?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(20), "AC", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(18), "T ?", object_1$class) # Apoptotic?
object_1$class <- ifelse(object_1$seurat_clusters %in% c(22), "TAM|OD", object_1$class)
object_1$class <- ifelse(object_1@reductions$umap@cell.embeddings[, 1] >= 10 &
  object_1@reductions$umap@cell.embeddings[, 1] <= 11 &
  object_1@reductions$umap@cell.embeddings[, 2] >= 1.5 &
  object_1@reductions$umap@cell.embeddings[, 2] <= 3,
"TC", object_1$class
)
object_1$seurat_clusters <- as.factor(paste0(as.character(object_1$seurat_clusters), ". ", object_1$class))

levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels = c(
  "0. T", "1. T", "2. T", "3. T", "7. T", "8. T", "9. T", "11. T",
  "17. T ?", "18. T ?",
  "20. AC",
  "21. NE",
  "4. OD", "6. OD", "12. OD", "19. OD",
  "16. EN",
  "15. PE",
  "22. TAM|OD",
  "5. TAM", "10. TAM", "13. TAM", "14. TAM",
  "14. TC"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters") +
  labs(subtitle = sid) +
  guides(col = guide_legend(ncol = 1, override.aes = list(size = 3)))



ggsave(paste0("output/figures/scRNA/Glimmunology/", sid, "_UMAP.pdf"), width = 10, height = 8)
ggsave(paste0("output/figures/scRNA/Glimmunology/", sid, "_UMAP.png"), width = 12, height = 10)

#od.markers <- FindMarkers(object_1, ident.1 = c(4,6,12,19,22))
#View(od.markers)



# 
# 
# tmp.17 <- FindMarkers(object_1, ident.1 = 17)
# head(tmp.17,20)
# 
# tmp.22 <- FindMarkers(object_1, ident.1 = 22)
# head(tmp.22,20)

tmp.15 <- FindMarkers(object_1, ident.1 = 15) # PE
View(tmp.15)




#### 1. Tumor (+) ----

FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR", max.cutoff = 4) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = "PDGFRB") # Tumor/MES


FeaturePlot(object = object_1, features = "EREG") # EREG
FeaturePlot(object = object_1, features = "EGF") # BTC
FeaturePlot(object = object_1, features = "BTC") # BTC

DotPlot(object = object_1, features =list( 'ligands'=c('EGF','EREG','AREG','BTC','EPGN','HBEGF') ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("EGFR and ligands in: ",sid))


FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?

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

FeaturePlot(object = object_1, features = "NODAL")


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


#FeaturePlot(object = object_1, features = "DCN") # DCN
#FeaturePlot(object = object_1, features = "COL1A2") # DCN
FeaturePlot(object = object_1, features = "ANPEP") # DCN



##### figure S6a ----


tmp.c4 <- results.out |>
  dplyr::filter(!is.na(.data$C4.2022)) |> 
  dplyr::filter(.data$C4.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.npc1 <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.NPC1)) |> 
  dplyr::filter(.data$neftel.meta.modules.NPC1 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.npc2 <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.NPC2)) |> 
  dplyr::filter(.data$neftel.meta.modules.NPC2 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.npc1.2 <- intersect(tmp.npc1, tmp.npc2)
tmp.npc1 <- setdiff(tmp.npc1, tmp.npc1.2)
tmp.npc2 <- setdiff(tmp.npc2, tmp.npc1.2)

tmp.c4.npc2 <- intersect(tmp.c4, tmp.npc2)
tmp.c4  <- setdiff(tmp.c4, tmp.c4.npc2)
tmp.npc2 <- setdiff(tmp.npc2, tmp.c4.npc2)


sid_print <- 'Samply Y (van Hijfte dataset - single nucleus RNA-seq)'



tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



ggsave(paste0("output/figures/2022_figure_S6a.pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)





#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "MOG")
#FeaturePlot(object = object_1, features = "OPALIN") # - clear small cluster
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "TMEM144")


FeaturePlot(object = object_1, features = "ST18") # OD?
FeaturePlot(object = object_1, features = "MBP") # OD?
FeaturePlot(object = object_1, features = "CTNNA3") # OD?
FeaturePlot(object = object_1, features = "SLC24A2") # OD?
FeaturePlot(object = object_1, features = "KIRREL3") # OD?
FeaturePlot(object = object_1, features = "NKAIN2") # OD?
FeaturePlot(object = object_1, features = "MAP7") # OD?
FeaturePlot(object = object_1, features = "RNF220") # OD?
FeaturePlot(object = object_1, features = "PEX5L") # OD?
FeaturePlot(object = object_1, features = "TMEM144") # OD?
FeaturePlot(object = object_1, features = "EDIL3") # OD?
FeaturePlot(object = object_1, features = "DOCK5") # OD?
FeaturePlot(object = object_1, features = "MOBP") # OD?
FeaturePlot(object = object_1, features = "UNC5C") # OD?
FeaturePlot(object = object_1, features = "CLDN11") # OD?
FeaturePlot(object = object_1, features = "SPOCK3") # OD?
FeaturePlot(object = object_1, features = "CNTNAP4") # OD?
FeaturePlot(object = object_1, features = "MAN2A1") # OD?
FeaturePlot(object = object_1, features = "PCSK6") # OD?
FeaturePlot(object = object_1, features = "TTLL7") # OD?

FeaturePlot(object = object_1, features = "OLIG2") # OD?



##### figure S6b ----


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



sid_print <- 'Samply Y (van Hijfte dataset - single nucleus RNA-seq)'



DotPlot(object = object_1, features =list('C3'=tmp.c3, 'OPC'=tmp.opc, 'C3+OPC'=tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC] in: ", sid_print))



ggsave(paste0("output/figures/2022_figure_S6b.pdf"),width=7.5*2, height=4,scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)




#### 6A. Endothelial (+) ----

FeaturePlot(object = object_1, features = "ABCB1")
FeaturePlot(object = object_1, features = "CD34")
FeaturePlot(object = object_1, features = "FLT4")
FeaturePlot(object = object_1, features = "TIE1") # meh
FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?



#### 6B. Pericytes (+) ----

FeaturePlot(object = object_1, features = c("RGS5","PDGFRB","CD248","PEAR1", "HEYL" , "CFH"))

FeaturePlot(object = object_1, features = c("RGS5"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))


#### F] Figure S14K ---

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  labs(subtitle = sid)
ggsave("output/figures/22022_Figure_S14K_labels.pdf", width = 6.5, height = 4, scale = 1.2) # to export cluster names


FeaturePlot(object = object_1, features = c("COL1A1", "COL1A2", "PDGFRB", "PECAM1"), min.cutoff = 1, order = T, pt.size = 0.15)
ggsave("output/figures/2022_Figure_S14K.pdf", width = 6.5, height = 4, scale = 1.2)




#### C0-2022 ----
##### F] Figure S12J - C0 ----


tmp.c0 <- results.out |>
  dplyr::filter(!is.na(.data$C0.2022)) |>
  dplyr::filter(.data$C0.2022 == T) |>
  dplyr::filter(!is.na(hugo_symbol)) |>
  dplyr::pull(hugo_symbol) |>
  unique()


sid_print <- "Samply Y (van Hijfte dataset - single nucleus RNA-seq)"


DotPlot(object = object_1, features = list("C0" = tmp.c0), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
  labs(x = paste0("Features [C0] in: ", sid_print))


ggsave(paste0("output/figures/2022_Figure_S12J.pdf"), width = 6.5, height = 4, scale = 1.2)
rm(tmp.c0, sid_print)



#### C1-2022 (up) ----
##### F] Figure S14F - C1 ----


tmp.c1 <- results.out |>
  dplyr::filter(!is.na(.data$C1.2022)) |> 
  dplyr::filter(.data$C1.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique() |> 
  sort()


sid_print <- 'Samply Y (van Hijfte dataset - single nucleus RNA-seq)'


DotPlot(object = object_1, features =list('C1'=tmp.c1, 'Peri'=c("RGS5", "PDGFRB", "CD248")), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C1] in: ",sid_print))



ggsave(paste0("output/figures/2022_Figure_S14F.pdf"),width=6.5, height=4, scale=1.2)
rm(tmp.c1, sid_print)




#### C2-2022 (Endo) (down) ----
##### figure S9c ----


tmp.c2 <- results.out |>
  dplyr::filter(.data$C2.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()

tmp.endo <- read_xlsx("data/McKenzie et al. Gene expression different cell types.xlsx", sheet='top_human_specificity') |>
  dplyr::select(c('grand_mean', 'gene', 'Celltype')) |>
  dplyr::filter(Celltype == 'end') |> 
  dplyr::arrange(desc(grand_mean)) |>
  dplyr::filter(gene %in% all.genes ) |>
  dplyr::slice_head(n=25) |>
  dplyr::mutate(grand_mean = NULL) |> 
  dplyr::pull(gene)

tmp.peri <- c('PDGFRB','CD248','RGS5')



tmp.c2 <- setdiff(tmp.c2, c(tmp.peri))
tmp.endo <- setdiff(tmp.endo, c(tmp.peri,tmp.c2))
tmp.peri <- setdiff(tmp.peri, c(tmp.c2, tmp.endo))


sid_print <- 'Samply Y (van Hijfte dataset - single nucleus RNA-seq)'


DotPlot(object = object_1, features = list('C2 (Endothelial)'=tmp.c2,
                                           'Endothelial'=tmp.endo,
                                           'Pericyte'=tmp.peri), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C2 & top25 McKenzie endothelial markers] in: ",sid_print))


ggsave(paste0("output/figures/2022_figure_S9c.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c2, tmp.peri, tmp.endo, sid_print)



#### MES ----


tmp <- results.out |> 
  # TCGA.subtype.marker == "TCGA-MES" | 
  dplyr::filter(neftel.meta.module.MES1 | neftel.meta.module.MES2) |> 
  dplyr::select(hugo_symbol,neftel.meta.module.MES1 , neftel.meta.module.MES2) |> 
  #dplyr::mutate(TCGA.subtype.marker = ifelse(is.na(TCGA.subtype.marker), F,T)) |> 
  tidyr::pivot_longer(cols=c( neftel.meta.module.MES1 , neftel.meta.module.MES2)) |> 
  dplyr::filter(value) |> 
  dplyr::mutate(value = NULL) |> 
  dplyr::distinct() |> 
  dplyr::group_by(hugo_symbol) |> 
  dplyr::summarise(str = paste0(name, collapse=",")) |> 
  dplyr::ungroup()



DotPlot(object = object_1, features = list(
  "MES1"                         = tmp |> dplyr::filter(str == "neftel.meta.module.MES1") |> dplyr::pull(hugo_symbol),
  #"TCGA.subtype.marker"                            = tmp |> dplyr::filter(str == "TCGA.subtype.marker") |> dplyr::pull(hugo_symbol),
  "MES2"                         = tmp |> dplyr::filter(str == "neftel.meta.module.MES2") |> dplyr::pull(hugo_symbol),
  "MES1+2" = tmp |> dplyr::filter(str == "neftel.meta.module.MES1,neftel.meta.module.MES2" ) |> dplyr::pull(hugo_symbol),
  'markers?' = c('DOCK7','DOCK10') # 'ZBTB20'
  #"TCGA.subtype.marker,neftel.meta.module.MES1"    = tmp |> dplyr::filter(str == "TCGA.subtype.marker,neftel.meta.module.MES1") |> dplyr::pull(hugo_symbol),
  #"TCGA.subtype.marker,neftel.meta.module.MES2"    = tmp |> dplyr::filter(str == "TCGA.subtype.marker,neftel.meta.module.MES2") |> dplyr::pull(hugo_symbol)
), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [MES] in: ",sid))



## de markers for MES [tumor intrinsic] vs [enviroment]
tmp.mes <- FindMarkers(object_1, 
                       ident.1 = c(1,9,18,20),
                       ident.2 = c(21,4,6,12,19,16,15,22,5,10,13,14,
                                   17,11,8,7,3,0
                                   ))

head(tmp.mes,n=15)


FeaturePlot(object = object_1, features = "ST6GALNAC3")
FeaturePlot(object = object_1, features = "RNF149")
FeaturePlot(object = object_1, features = "DOCK7")
FeaturePlot(object = object_1, features = "DOCK10")
FeaturePlot(object = object_1, features = "KCNN3")
FeaturePlot(object = object_1, features = "ZBTB20")





#' ## GLASS-NL purity ----

#' ### ACE ----


#' p1 <- DotPlot(object = object_1, features = list(
#'   'ACE.top.10'=c("NCBP3","AP002495.1","FAM106A","DEPDC5","RAD50","TRA2A","AL583810.2","BCS1L","LINC02804","CCDC32"),
#'   'ACE.bottom.10'=c("AL137139.3","RASL12","HTT-AS","CTXND1","AC104088.3","HKDC1","AL031056.1","SPOCK2","AC013553.3","AC104088.1")
#' ), group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#' ### ACE != 1.0 ----

#' p2 <- DotPlot(object = object_1, features = list(
#'   'ACE.below-1.top.10'=c("HSBP1L1","NAV3","PIWIL4","RGN","LGALS8","DNAJA4","SOWAHA","IGSF10","CES4A","SPATA6L"),
#'   'ACE.below-1.bottom.10'=c("MED20","HNRNPA1","PRR3","ZNF426","HOXD-AS2","GPR173","BAZ1A","ZNF121","ZNF333","ZNF747")),
#'   group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#' ### Erik curated fit ----

#' p3 <- DotPlot(object = object_1, features = list(
#'   'Erik.Manual.top.10'=c("KHNYN","TTC12","PLD1","SLC22A15","CYB5R2","SEMA4D","YPEL2","KCNJ2","GCA","NT5DC1"),
#'   'Erik.Manual.bottom.10'=c("BEST3","AL079301.1","VAX2","PRR3","RPAIN","PCGF2","AL645608.2","USP49","AC078846.1","AC006504.1")),
#'   group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#' p4 <- DotPlot(object = object_1, features = list(
#'   'VAF.top.10'=c("STPG1","HCG11","DLC1","NIPAL3","SYNJ2","UNC5C","CD55","EDIL3","SLC22A15","CBR1"),
#'   'VAF.bottom.10'=c("AC078846.1","USP49","BBS1","AL079301.1","ABCB6","BEST3","SOX2","RPAIN","PHF21B","PPOX")
#' ), group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#' p5 <- DotPlot(object = object_1, features = list(
#'   'meth.top.10'=c("RAB11FIP1","HCG11","SWAP70","TACC1","OSTF1","CD55","KHNYN","MCTP2","CCDC69","MYL12A"),
#'   'meth.bottom.10'=c("GPR173","C1orf61","SOX2","ATAT1","PCDHB9","ZNF462","PHF21B","SNCAIP","LINC00461","AC078846.1")
#' ), group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#' p1 + p2 + p3 + p4 + p5




#' DotPlot(object = object_1, features = list(
#'   'ACE.bottom.10'=c("NCBP3","AP002495.1","FAM106A","DEPDC5","RAD50","TRA2A","AL583810.2","BCS1L","LINC02804","CCDC32"),
#'   'ACE.top.10'=c("AL137139.3","RASL12","HTT-AS","CTXND1","AC104088.3","HKDC1","AL031056.1","SPOCK2","AC013553.3","AC104088.1"),
#'   
#'   'VAF.bottom.10'=c("STPG1","HCG11","DLC1","NIPAL3","SYNJ2","UNC5C","CD55","EDIL3","SLC22A15","CBR1"),
#'   'VAF.top.10'=c("AC078846.1","USP49","BBS1","AL079301.1","ABCB6","BEST3","SOX2","RPAIN","PHF21B","PPOX")
#' ), group.by = "seurat_clusters") + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




#' DotPlot(object = object_1, features = list(
#'   'ACE.top.10'=c("NCBP3","AP002495.1","FAM106A","DEPDC5","RAD50","TRA2A","AL583810.2","BCS1L","LINC02804","CCDC32"),
#'   'ACE.bottom.10'=c("AL137139.3","RASL12","HTT-AS","CTXND1","AC104088.3","HKDC1","AL031056.1","SPOCK2","AC013553.3","AC104088.1"),
#'   
#'   'meth.top.10'=c("RAB11FIP1","HCG11","SWAP70","TACC1","OSTF1","CD55","KHNYN","MCTP2","CCDC69","MYL12A"),
#'   'meth.bottom.10'=c("GPR173","C1orf61","SOX2","ATAT1","PCDHB9","ZNF462","PHF21B","SNCAIP","LINC00461","AC078846.1")
#' ), group.by = "seurat_clusters") + # T
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




#' FeaturePlot(object = object_1, features = "")



