f#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(infercnv)
library(AnnotationHub)
library(ensembldb)


# load data ----


source('scripts/load_results.out.R')



# cluster genes ----

# neuron.genes <- c("PTPRN","ICAM5","VWA7","PANX2","STRC","NPAS4","EGR4","OTOF","L1CAM","DOC2B","STXBP6","REPS2","DYNC1I1","AMPH",
#                   "PWWP3B","KCNQ5","ZNF98","KCNK12","BSPRY","YPEL4","PDZD7","ZNF365","SLC6A15","RBM11","CDKL2","KCTD4","RASGRF2",
#                   "IGFL2","TRIM54","TNNT2","ESYT3","RBP4","FSTL5","EPHB6","CAMKK1","SNCG","MATK","HS3ST2","FABP3","CBLN2","ATP8A2",
#                   "SULT4A1","ATP2B3","TMEM130","SSTR1","CNNM1","CABP1","SHISAL1","FGF13","CBLN1","NEFH","KCNS1","RYR2","NEFM",
#                   "ANKRD30B","KHDRBS2","SYCE1","SLC35F3","CDH18","CDH12","TRHDE","TCERG1L","FAM153CP","FAM153B","FAM153A","SOHLH1",
#                   "MTUS2","DLGAP2","FRMPD4","CCK","HOOK1","MAP7D2","GALNTL6","CDH9","NEGR1","LRFN5","GRM1","ACVR1C","PPP4R4","DRD1",
#                   "OCA2","CNGB1","RTN4RL1","GLP2R","MEPE","VSNL1","GABRB2","KRT222","HTR2A","SLITRK4","RXFP1","VSTM2L","SLC6A7",
#                   "ISLR2","CYP4X1","COL26A1","CABLES1","SYNGR3","ARHGDIG","PSD","NPM2","SGSM1","OLFM1","SRRM4","NEUROD2","UBE2QL1",
#                   "SVOP","CELF4","PCSK2","ATP1A3","SCN3B","CHGB","JPH3","DLGAP3","CKMT1B","KCNH3","DOC2A","KCNN1","UNC5A","SHANK1",
#                   "TMEM151B","VWA5B2","STX1B","FBXL16","ADAM11","TMEM63C","PDE2A","TMEM246","PRKCZ","SNAP91","SH3GL2","KCTD16",
#                   "NMNAT2","STXBP5L","TAC3","RIMBP2","CRHR2","CUX2","CLEC2L","SLC12A5","LRTM2","IQSEC3","WSCD2","HTR5A","GPR83",
#                   "NRGN","KCNS2","SLC30A3","DNAJC5G","SYT7","GABRD","TGFBR3L","KCNC2","GABRG2","GABRA1","HCN1","CREG2","CACNG3",
#                   "C1QL3","SLC6A17","SLC17A7","GABRG3","GABRA5","KCNV1","GPR26","GAD2","GPR22","CAMK1G","PDYN","SERTM1","HPCA",
#                   "CPLX2","SYT5","SLC8A2","CHRM1","RAB3A","GRIN1","CHGA","SNCB","FAM163B","SYN1","FSTL4","GPR61","SLC32A1","CPLX1",
#                   "SPTB","GLS2","SYT13","SNAP25","SYT1","SV2C","SYNPR","OLFM3","RBFOX1","UNC13C","RBFOX3","SYT2","PHYHIP","KCNA4",
#                   "HPCAL4","MPPED1","EMX1","PACSIN1","CALY","CACNA1B","TMEM132D","SV2B","CCKBR","RASAL1","DDN","C4orf50","HRH3",
#                   "CACNA1I","PHF24","MFSD4A","CAMK2A","SST","PRKCG","TBR1","SLC4A10","ARHGAP44","AJAP1","KCNT1","CHD5","GALNT9",
#                   "GRM4","NGEF","RTN4R","SCN2B","NAP1L2","DMTN","BRINP1","GRIN2A","GDA","SNCA","PRKCB","SERPINI1","NECAB1","KCNK1",
#                   "AK5","GABRA2","PPP2R2C","CPNE6","SOWAHA","ADARB2","NPY","KIRREL3","PNMA8B","FAIM2","DLG2","KIAA0319","SLC39A12",
#                   "ETNPPL","HPSE2")
# neuron.genes.ens <- c("ENSG00000187730","ENSG00000067606","ENSG00000196581","ENSG00000116254","ENSG00000121769","ENSG00000121905","ENSG00000116544","ENSG00000116983",
#                       "ENSG00000186377","ENSG00000134709","ENSG00000172260","ENSG00000154027","ENSG00000118733","ENSG00000156097","ENSG00000197106","ENSG00000157064",
#                       "ENSG00000118194","ENSG00000143858","ENSG00000174514","ENSG00000008118","ENSG00000135750","ENSG00000183780","ENSG00000198626","ENSG00000163032",
#                       "ENSG00000115155","ENSG00000115194","ENSG00000163793","ENSG00000138100","ENSG00000184261","ENSG00000135638","ENSG00000135625","ENSG00000175874",
#                       "ENSG00000123612","ENSG00000136535","ENSG00000144290","ENSG00000054356","ENSG00000066248","ENSG00000187094","ENSG00000163630","ENSG00000145087",
#                       "ENSG00000158220","ENSG00000163536","ENSG00000145198","ENSG00000157005","ENSG00000168993","ENSG00000181215","ENSG00000074211","ENSG00000151834",
#                       "ENSG00000138769","ENSG00000152595","ENSG00000145335","ENSG00000164089","ENSG00000171509","ENSG00000168843","ENSG00000174473","ENSG00000215218",
#                       "ENSG00000145526","ENSG00000154162","ENSG00000113100","ENSG00000164588","ENSG00000122012","ENSG00000113319","ENSG00000198944","ENSG00000053108",
#                       "ENSG00000183775","ENSG00000011083","ENSG00000070808","ENSG00000145864","ENSG00000022355","ENSG00000113327","ENSG00000184845","ENSG00000145920",
#                       "ENSG00000182230","ENSG00000074317","ENSG00000113763","ENSG00000170074","ENSG00000204677","ENSG00000137261","ENSG00000204396","ENSG00000124493",
#                       "ENSG00000124507","ENSG00000178233","ENSG00000112232","ENSG00000185760","ENSG00000065609","ENSG00000152822","ENSG00000122585","ENSG00000106113",
#                       "ENSG00000078053","ENSG00000158560","ENSG00000166448","ENSG00000160963","ENSG00000172209","ENSG00000236279","ENSG00000106123","ENSG00000157219",
#                       "ENSG00000198010","ENSG00000158806","ENSG00000158856","ENSG00000168490","ENSG00000104722","ENSG00000123119","ENSG00000156486","ENSG00000164794",
#                       "ENSG00000107295","ENSG00000122733","ENSG00000119125","ENSG00000165152","ENSG00000119411","ENSG00000078725","ENSG00000196990","ENSG00000130558",
#                       "ENSG00000165643","ENSG00000107147","ENSG00000176884","ENSG00000148408","ENSG00000185736","ENSG00000165985","ENSG00000148482","ENSG00000136750",
#                       "ENSG00000138311","ENSG00000173267","ENSG00000138207","ENSG00000172987","ENSG00000119946","ENSG00000186862","ENSG00000059915","ENSG00000154478",
#                       "ENSG00000176769","ENSG00000130643","ENSG00000171772","ENSG00000110148","ENSG00000182255","ENSG00000019505","ENSG00000166793","ENSG00000011347",
#                       "ENSG00000168539","ENSG00000174576","ENSG00000186642","ENSG00000150672","ENSG00000123901","ENSG00000149575","ENSG00000166257","ENSG00000154146",
#                       "ENSG00000149571","ENSG00000120645","ENSG00000166159","ENSG00000181418","ENSG00000135519","ENSG00000135472","ENSG00000135423","ENSG00000166863",
#                       "ENSG00000072657","ENSG00000166006","ENSG00000067715","ENSG00000072041","ENSG00000075035","ENSG00000166111","ENSG00000111249","ENSG00000111344",
#                       "ENSG00000139767","ENSG00000157782","ENSG00000151952","ENSG00000060709","ENSG00000182870","ENSG00000132932","ENSG00000132938","ENSG00000180440",
#                       "ENSG00000180332","ENSG00000102468","ENSG00000100884","ENSG00000168952","ENSG00000139874","ENSG00000165379","ENSG00000070182","ENSG00000165548",
#                       "ENSG00000100604","ENSG00000119698","ENSG00000186297","ENSG00000182256","ENSG00000104044","ENSG00000237289","ENSG00000242866","ENSG00000137766",
#                       "ENSG00000167178","ENSG00000185518","ENSG00000242173","ENSG00000127585","ENSG00000127561","ENSG00000078328","ENSG00000183454","ENSG00000122254",
#                       "ENSG00000166501","ENSG00000006116","ENSG00000149927","ENSG00000099365","ENSG00000102924","ENSG00000070729","ENSG00000154118","ENSG00000272636",
#                       "ENSG00000185924","ENSG00000004660","ENSG00000065325","ENSG00000006740","ENSG00000171532","ENSG00000213424","ENSG00000073670","ENSG00000167281",
#                       "ENSG00000180777","ENSG00000134508","ENSG00000101489","ENSG00000141668","ENSG00000007264","ENSG00000260001","ENSG00000105376","ENSG00000105642",
#                       "ENSG00000105649","ENSG00000197360","ENSG00000105409","ENSG00000204866","ENSG00000204851","ENSG00000118160","ENSG00000104888","ENSG00000161681",
#                       "ENSG00000126583","ENSG00000129990","ENSG00000101327","ENSG00000089199","ENSG00000132639","ENSG00000125851","ENSG00000132821","ENSG00000101438",
#                       "ENSG00000124134","ENSG00000124140","ENSG00000101180","ENSG00000185272","ENSG00000040608","ENSG00000167037","ENSG00000100285","ENSG00000100346",
#                       "ENSG00000186732","ENSG00000130540","ENSG00000138944","ENSG00000073150","ENSG00000169933","ENSG00000169891","ENSG00000184368","ENSG00000008056",
#                       "ENSG00000186462","ENSG00000157502","ENSG00000129682","ENSG00000179542","ENSG00000067842","ENSG00000198910")



# oligodendrocyte.genes <- c("RASGRF1","FBXO2","PPP1R16B","TPPP","SEC14L5","TMEM151A","LGI3","TMCC2","HHATL","RAB11FIP4","PDIA2",
#                            "HCN2","KCNJ9","DNAJC6","TUBB4A","ADCY5","CSDC2","AC118754.1","PLIN4","HSPA2","PI16","PTGDS","CDK18",
#                            "FA2H","AATK","NKX6-2","MAG","PLCH2","FAM131C","PPP1R14A","CHADL","TMEM88B","ABCA2","PLPP2","CERCAM",
#                            "BOK","ACP7","CYS1","ANXA3","TNFSF9","AVPI1","MYH11","ADH1B","TUBA4A","CORO6","IL12RB2","TESPA1",
#                            "MPP7","RSPO3","KCNJ12","OPN4","MKX","FRAS1","CPNE7")

OPC <- c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4","LPPR1","PTPRZ1","VCAN","DBI","PMP2","CNP","TNS3","LIMA1","CA10","PCDHGC3","CNTN1","SCD5","P2RX7","CADM2","TTYH1","FGF12","TMEM206","NEU4","FXYD6","RNF13","RTKN","GPM6B","LMF1","ALCAM","PGRMC1","HRASLS","BCAS1","RAB31","PLLP","FABP5","NLGN3","SERINC5","EPB41L2","GPR37L1")
NPC1 <- c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST","ASCL1","BTG2","DCX","NXPH1","HN1","PFN2","SCG3","MYT1","CHD7","GPR56","TUBA1A","PCBP4","ETV1","SHD","TNR","AMOTL2","DBN1","HIP1","ABAT","ELAVL4","LMF1","GRIK2","SERINC5","TSPAN13","ELMO1","GLCCI1","SEZ6L","LRRN1","SEZ6","SOX11")
NPC2 <- c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4","DLX5","SOX4","MAP1B","RBFOX2","IGFBPL1","STMN1","HN1","TMEM161B-AS1","DPYSL3","SEPT3","PKIA","ATP1B1","DYNC1I1","CD200","SNAP25","PAK3","NDRG4","KIF5A","UCHL1","ENO2","KIF5C","DDAH2","TUBB2A","LBH","LOC150568","TCF4","GNG3","NFIB","DPYSL5","CRABP1","DBN1","NFIX","CEP170","BLCAP")


# oligodendrocyte.genes.ens <- c("ENSG00000205116","ENSG00000149527","ENSG00000116661","ENSG00000185519","ENSG00000116675",
#                                "ENSG00000081985","ENSG00000162728","ENSG00000133069","ENSG00000117266","ENSG00000205795",
#                                "ENSG00000127824","ENSG00000176720","ENSG00000010282","ENSG00000173175","ENSG00000138759",
#                                "ENSG00000138772","ENSG00000196616","ENSG00000171368","ENSG00000164530","ENSG00000146374",
#                                "ENSG00000168481","ENSG00000167123","ENSG00000107317","ENSG00000107331","ENSG00000150051",
#                                "ENSG00000150054","ENSG00000122375","ENSG00000119986","ENSG00000148826","ENSG00000179292",
#                                "ENSG00000135426","ENSG00000126803","ENSG00000058335","ENSG00000185615","ENSG00000103184",
#                                "ENSG00000133392","ENSG00000103089","ENSG00000178773","ENSG00000183018","ENSG00000184185",
#                                "ENSG00000167549","ENSG00000131242","ENSG00000181409","ENSG00000141934","ENSG00000099822",
#                                "ENSG00000167676","ENSG00000104833","ENSG00000125657","ENSG00000105695","ENSG00000167641",
#                                "ENSG00000183760","ENSG00000101445","ENSG00000100399","ENSG00000172346")




# C3 <- c('VWF', 'TIE1', 'HIGD1B', 'MMRN1', 'CYSLTR2', 'MMP25','FLT4', 'BCL6B', 'GRAP', 'LAMC3', 'DPEP1', 'PXDNL', 'ANGPT2',
#         'PALD1', 'ADGRD1', 'GBP6', 'SLC52A3', 'CLDN5', 'VWA2', 'ABCB1', 'THSD7B', 'SPINK8', 'FOXQ1', 'ZIC3', 'NODAL')
# 
# C4A <- c('SOD3', "FSTL3", "FAM180A", "OSGIN1", "NDRG1", "AC010327.1","TRIM29", "HSPB7", "TNNT1", "CCN5", "MICAL2", "GLIS1", "SLIT3",
#         "CYP26B1", "NPR3", "FGF5", "CCBE1", "GPR68", "SH3RF2")
# C4B <- c("WNT11", "SCUBE3", "KRT17", "GPR78","CPZ","GLI1", "PRB2","MAFA","HAPLN1")
# 
# C5 <- c("PRF1", "ARHGAP9", "FCMR","LXN","KCNE3", "NR5A2","FPR2", "CCL13", "MMP7", "CALCR", "LRG1", "SAA2", "PI3", "LIF", "HSPA6")
# 
# C6 <- c('CRABP2', 'CILP2', 'DPT', 'FGF7', 'COL10A1', 'FBN1', 'GLT8D2',
#         'IRX3', 'MFAP5', 'MFAP4', "COL8A2", "FNDC1", "MMP11", "MFAP2",
#         "COL1A2", "COL1A1", "COL5A1", "ADAMTS2", "TPSB2", "KRT8", "OMD",
#         "OGN", "MME", "MLPH", "MRC1L1", "PTGFR", "TWIST2", "C5orf46",
#         "TNNT3", "ASS1", "PERP","KLHDC7B", "CCL8")





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



ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_UMAP.png"),width=12,height=10)


od.markers <- FindMarkers(object_1, ident.1 = c(4,6,12,19,22))
View(od.markers)



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
FeaturePlot(object = object_1, features = c("PDGFRB"))
FeaturePlot(object = object_1, features = c("CD248"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))




#### C3 - C6 ----



DotPlot(object = object_1, features = c(C3, C4A, C4B, C5, C6), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



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


ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C3.pdf"),width=7.5, height=3,scale=2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C3.png"),width=7.5, height=3,scale=2)



RidgePlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C3), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C3)



#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)




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

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


FeaturePlot(object = object_1, features = C5[1:4])
FeaturePlot(object = object_1, features = C5[5:8])
FeaturePlot(object = object_1, features = C5[9:12])
FeaturePlot(object = object_1, features = C5[13:16])


#### C6 (up) ----

# + "PEAR1", "HEYL" , "CFH"
DotPlot(object = object_1, features =list('C6'=C6 , 
                                          'Peri'=c("RGS5", "PDGFRB", "CD248",
                                         'FB1+FB2'=c('MMP2', 'ANPEP','ACTA2','GLI1')) ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))

ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Glimmunology/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)




VlnPlot(object = object_1, features = C6, stack = T, sort = T)
VlnPlot(object = object_1, features = C6, stack = T, sort = T)


FeaturePlot(object = object_1, features = C6)



VlnPlot(object = object_1, features = C6, stack = T, sort = T)

DotPlot(object = object_1, features = C6) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(object = object_1, features = C6, stack = T, sort = T)

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


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters") + NoLegend()
ggsave("output/figures/scRNA/Glimmunology/COL1A2_labels.svg", width=10, height=7)
FeaturePlot(object = object_1, features =  "COL1A2" )
ggsave("output/figures/scRNA/Glimmunology/COL1A2_expr.svg", width=10, height=7)




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


DotPlot(object = object_1, features = list('test'=c('GREM1')), group.by = "seurat_clusters") + # oligodendrocyte
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(object = object_1, features = list('test'=c('SOX10')), group.by = "seurat_clusters") + # oligodendrocyte
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(object = object_1, features = list('test'=c('WEE1','HOXA5','HOXA6','HOXD10','HOXD11','SHOX')), group.by = "seurat_clusters") + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



FeaturePlot(object = object_1, features =  c("TROAP","FOXM1","AURKB","TPX2","H3C3","H3C2","TOP2A","DTL"))



#### C0-2022 ----
##### figure S10j ----


tmp.c0 <- results.out |>
  dplyr::filter(!is.na(.data$C0.2022)) |> 
  dplyr::filter(.data$C0.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()



sid_print <- 'Samply Y (van Hijfte dataset - single nucleus RNA-seq)'



DotPlot(object = object_1, features =list('C0'=tmp.c0), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C0] in: ",sid_print))



ggsave(paste0("output/figures/2022_figure_S10j.pdf"),width=6.5, height=4,scale=1.2)
rm(tmp.c0, sid_print)



#### C1-2022 (up) ----
##### figure S12f ----


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



ggsave(paste0("output/figures/2022_figure_S12f.pdf"),width=6.5, height=4, scale=1.2)
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
  labs(x = paste0("Features [C2 & top25 McKenzy endothelial markers] in: ",sid_print))


ggsave(paste0("output/figures/2022_figure_S9c.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
rm(tmp.c2, tmp.peri, tmp.endo, sid_print)



## GLASS-NL purity ----

### ACE ----


p1 <- DotPlot(object = object_1, features = list(
  'ACE.top.10'=c("NCBP3","AP002495.1","FAM106A","DEPDC5","RAD50","TRA2A","AL583810.2","BCS1L","LINC02804","CCDC32"),
  'ACE.bottom.10'=c("AL137139.3","RASL12","HTT-AS","CTXND1","AC104088.3","HKDC1","AL031056.1","SPOCK2","AC013553.3","AC104088.1")
), group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



### ACE != 1.0 ----

p2 <- DotPlot(object = object_1, features = list(
  'ACE.below-1.top.10'=c("HSBP1L1","NAV3","PIWIL4","RGN","LGALS8","DNAJA4","SOWAHA","IGSF10","CES4A","SPATA6L"),
  'ACE.below-1.bottom.10'=c("MED20","HNRNPA1","PRR3","ZNF426","HOXD-AS2","GPR173","BAZ1A","ZNF121","ZNF333","ZNF747")),
  group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### Erik curated fit ----

p3 <- DotPlot(object = object_1, features = list(
  'Erik.Manual.top.10'=c("KHNYN","TTC12","PLD1","SLC22A15","CYB5R2","SEMA4D","YPEL2","KCNJ2","GCA","NT5DC1"),
  'Erik.Manual.bottom.10'=c("BEST3","AL079301.1","VAX2","PRR3","RPAIN","PCGF2","AL645608.2","USP49","AC078846.1","AC006504.1")),
  group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p4 <- DotPlot(object = object_1, features = list(
  'VAF.top.10'=c("STPG1","HCG11","DLC1","NIPAL3","SYNJ2","UNC5C","CD55","EDIL3","SLC22A15","CBR1"),
  'VAF.bottom.10'=c("AC078846.1","USP49","BBS1","AL079301.1","ABCB6","BEST3","SOX2","RPAIN","PHF21B","PPOX")
), group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p5 <- DotPlot(object = object_1, features = list(
  'meth.top.10'=c("RAB11FIP1","HCG11","SWAP70","TACC1","OSTF1","CD55","KHNYN","MCTP2","CCDC69","MYL12A"),
  'meth.bottom.10'=c("GPR173","C1orf61","SOX2","ATAT1","PCDHB9","ZNF462","PHF21B","SNCAIP","LINC00461","AC078846.1")
), group.by = "seurat_clusters", col.min = -2.8, col.max=2.8) + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p1 + p2 + p3 + p4 + p5




DotPlot(object = object_1, features = list(
  'ACE.bottom.10'=c("NCBP3","AP002495.1","FAM106A","DEPDC5","RAD50","TRA2A","AL583810.2","BCS1L","LINC02804","CCDC32"),
  'ACE.top.10'=c("AL137139.3","RASL12","HTT-AS","CTXND1","AC104088.3","HKDC1","AL031056.1","SPOCK2","AC013553.3","AC104088.1"),
  
  'VAF.bottom.10'=c("STPG1","HCG11","DLC1","NIPAL3","SYNJ2","UNC5C","CD55","EDIL3","SLC22A15","CBR1"),
  'VAF.top.10'=c("AC078846.1","USP49","BBS1","AL079301.1","ABCB6","BEST3","SOX2","RPAIN","PHF21B","PPOX")
), group.by = "seurat_clusters") + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




DotPlot(object = object_1, features = list(
  'ACE.top.10'=c("NCBP3","AP002495.1","FAM106A","DEPDC5","RAD50","TRA2A","AL583810.2","BCS1L","LINC02804","CCDC32"),
  'ACE.bottom.10'=c("AL137139.3","RASL12","HTT-AS","CTXND1","AC104088.3","HKDC1","AL031056.1","SPOCK2","AC013553.3","AC104088.1"),
  
  'meth.top.10'=c("RAB11FIP1","HCG11","SWAP70","TACC1","OSTF1","CD55","KHNYN","MCTP2","CCDC69","MYL12A"),
  'meth.bottom.10'=c("GPR173","C1orf61","SOX2","ATAT1","PCDHB9","ZNF462","PHF21B","SNCAIP","LINC00461","AC078846.1")
), group.by = "seurat_clusters") + # T
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




FeaturePlot(object = object_1, features = "")



