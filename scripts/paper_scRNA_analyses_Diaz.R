#!/usr/bin/env R

# load libs ----


library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(infercnv)
library(AnnotationHub)
library(ensembldb)
library(patchwork)



# cluster genes ----


neuron.genes <- c("PTPRN","ICAM5","VWA7","PANX2","STRC","NPAS4","EGR4","OTOF","L1CAM","DOC2B","STXBP6","REPS2","DYNC1I1","AMPH",
                  "PWWP3B","KCNQ5","ZNF98","KCNK12","BSPRY","YPEL4","PDZD7","ZNF365","SLC6A15","RBM11","CDKL2","KCTD4","RASGRF2",
                  "IGFL2","TRIM54","TNNT2","ESYT3","RBP4","FSTL5","EPHB6","CAMKK1","SNCG","MATK","HS3ST2","FABP3","CBLN2","ATP8A2",
                  "SULT4A1","ATP2B3","TMEM130","SSTR1","CNNM1","CABP1","SHISAL1","FGF13","CBLN1","NEFH","KCNS1","RYR2","NEFM",
                  "ANKRD30B","KHDRBS2","SYCE1","SLC35F3","CDH18","CDH12","TRHDE","TCERG1L","FAM153CP","FAM153B","FAM153A","SOHLH1",
                  "MTUS2","DLGAP2","FRMPD4","CCK","HOOK1","MAP7D2","GALNTL6","CDH9","NEGR1","LRFN5","GRM1","ACVR1C","PPP4R4","DRD1",
                  "OCA2","CNGB1","RTN4RL1","GLP2R","MEPE","VSNL1","GABRB2","KRT222","HTR2A","SLITRK4","RXFP1","VSTM2L","SLC6A7",
                  "ISLR2","CYP4X1","COL26A1","CABLES1","SYNGR3","ARHGDIG","PSD","NPM2","SGSM1","OLFM1","SRRM4","NEUROD2","UBE2QL1",
                  "SVOP","CELF4","PCSK2","ATP1A3","SCN3B","CHGB","JPH3","DLGAP3","CKMT1B","KCNH3","DOC2A","KCNN1","UNC5A","SHANK1",
                  "TMEM151B","VWA5B2","STX1B","FBXL16","ADAM11","TMEM63C","PDE2A","TMEM246","PRKCZ","SNAP91","SH3GL2","KCTD16",
                  "NMNAT2","STXBP5L","TAC3","RIMBP2","CRHR2","CUX2","CLEC2L","SLC12A5","LRTM2","IQSEC3","WSCD2","HTR5A","GPR83",
                  "NRGN","KCNS2","SLC30A3","DNAJC5G","SYT7","GABRD","TGFBR3L","KCNC2","GABRG2","GABRA1","HCN1","CREG2","CACNG3",
                  "C1QL3","SLC6A17","SLC17A7","GABRG3","GABRA5","KCNV1","GPR26","GAD2","GPR22","CAMK1G","PDYN","SERTM1","HPCA",
                  "CPLX2","SYT5","SLC8A2","CHRM1","RAB3A","GRIN1","CHGA","SNCB","FAM163B","SYN1","FSTL4","GPR61","SLC32A1","CPLX1",
                  "SPTB","GLS2","SYT13","SNAP25","SYT1","SV2C","SYNPR","OLFM3","RBFOX1","UNC13C","RBFOX3","SYT2","PHYHIP","KCNA4",
                  "HPCAL4","MPPED1","EMX1","PACSIN1","CALY","CACNA1B","TMEM132D","SV2B","CCKBR","RASAL1","DDN","C4orf50","HRH3",
                  "CACNA1I","PHF24","MFSD4A","CAMK2A","SST","PRKCG","TBR1","SLC4A10","ARHGAP44","AJAP1","KCNT1","CHD5","GALNT9",
                  "GRM4","NGEF","RTN4R","SCN2B","NAP1L2","DMTN","BRINP1","GRIN2A","GDA","SNCA","PRKCB","SERPINI1","NECAB1","KCNK1",
                  "AK5","GABRA2","PPP2R2C","CPNE6","SOWAHA","ADARB2","NPY","KIRREL3","PNMA8B","FAIM2","DLG2","KIAA0319","SLC39A12",
                  "ETNPPL","HPSE2")
neuron.genes.ens <- c("ENSG00000187730","ENSG00000067606","ENSG00000196581","ENSG00000116254","ENSG00000121769","ENSG00000121905","ENSG00000116544","ENSG00000116983",
                      "ENSG00000186377","ENSG00000134709","ENSG00000172260","ENSG00000154027","ENSG00000118733","ENSG00000156097","ENSG00000197106","ENSG00000157064",
                      "ENSG00000118194","ENSG00000143858","ENSG00000174514","ENSG00000008118","ENSG00000135750","ENSG00000183780","ENSG00000198626","ENSG00000163032",
                      "ENSG00000115155","ENSG00000115194","ENSG00000163793","ENSG00000138100","ENSG00000184261","ENSG00000135638","ENSG00000135625","ENSG00000175874",
                      "ENSG00000123612","ENSG00000136535","ENSG00000144290","ENSG00000054356","ENSG00000066248","ENSG00000187094","ENSG00000163630","ENSG00000145087",
                      "ENSG00000158220","ENSG00000163536","ENSG00000145198","ENSG00000157005","ENSG00000168993","ENSG00000181215","ENSG00000074211","ENSG00000151834",
                      "ENSG00000138769","ENSG00000152595","ENSG00000145335","ENSG00000164089","ENSG00000171509","ENSG00000168843","ENSG00000174473","ENSG00000215218",
                      "ENSG00000145526","ENSG00000154162","ENSG00000113100","ENSG00000164588","ENSG00000122012","ENSG00000113319","ENSG00000198944","ENSG00000053108",
                      "ENSG00000183775","ENSG00000011083","ENSG00000070808","ENSG00000145864","ENSG00000022355","ENSG00000113327","ENSG00000184845","ENSG00000145920",
                      "ENSG00000182230","ENSG00000074317","ENSG00000113763","ENSG00000170074","ENSG00000204677","ENSG00000137261","ENSG00000204396","ENSG00000124493",
                      "ENSG00000124507","ENSG00000178233","ENSG00000112232","ENSG00000185760","ENSG00000065609","ENSG00000152822","ENSG00000122585","ENSG00000106113",
                      "ENSG00000078053","ENSG00000158560","ENSG00000166448","ENSG00000160963","ENSG00000172209","ENSG00000236279","ENSG00000106123","ENSG00000157219",
                      "ENSG00000198010","ENSG00000158806","ENSG00000158856","ENSG00000168490","ENSG00000104722","ENSG00000123119","ENSG00000156486","ENSG00000164794",
                      "ENSG00000107295","ENSG00000122733","ENSG00000119125","ENSG00000165152","ENSG00000119411","ENSG00000078725","ENSG00000196990","ENSG00000130558",
                      "ENSG00000165643","ENSG00000107147","ENSG00000176884","ENSG00000148408","ENSG00000185736","ENSG00000165985","ENSG00000148482","ENSG00000136750",
                      "ENSG00000138311","ENSG00000173267","ENSG00000138207","ENSG00000172987","ENSG00000119946","ENSG00000186862","ENSG00000059915","ENSG00000154478",
                      "ENSG00000176769","ENSG00000130643","ENSG00000171772","ENSG00000110148","ENSG00000182255","ENSG00000019505","ENSG00000166793","ENSG00000011347",
                      "ENSG00000168539","ENSG00000174576","ENSG00000186642","ENSG00000150672","ENSG00000123901","ENSG00000149575","ENSG00000166257","ENSG00000154146",
                      "ENSG00000149571","ENSG00000120645","ENSG00000166159","ENSG00000181418","ENSG00000135519","ENSG00000135472","ENSG00000135423","ENSG00000166863",
                      "ENSG00000072657","ENSG00000166006","ENSG00000067715","ENSG00000072041","ENSG00000075035","ENSG00000166111","ENSG00000111249","ENSG00000111344",
                      "ENSG00000139767","ENSG00000157782","ENSG00000151952","ENSG00000060709","ENSG00000182870","ENSG00000132932","ENSG00000132938","ENSG00000180440",
                      "ENSG00000180332","ENSG00000102468","ENSG00000100884","ENSG00000168952","ENSG00000139874","ENSG00000165379","ENSG00000070182","ENSG00000165548",
                      "ENSG00000100604","ENSG00000119698","ENSG00000186297","ENSG00000182256","ENSG00000104044","ENSG00000237289","ENSG00000242866","ENSG00000137766",
                      "ENSG00000167178","ENSG00000185518","ENSG00000242173","ENSG00000127585","ENSG00000127561","ENSG00000078328","ENSG00000183454","ENSG00000122254",
                      "ENSG00000166501","ENSG00000006116","ENSG00000149927","ENSG00000099365","ENSG00000102924","ENSG00000070729","ENSG00000154118","ENSG00000272636",
                      "ENSG00000185924","ENSG00000004660","ENSG00000065325","ENSG00000006740","ENSG00000171532","ENSG00000213424","ENSG00000073670","ENSG00000167281",
                      "ENSG00000180777","ENSG00000134508","ENSG00000101489","ENSG00000141668","ENSG00000007264","ENSG00000260001","ENSG00000105376","ENSG00000105642",
                      "ENSG00000105649","ENSG00000197360","ENSG00000105409","ENSG00000204866","ENSG00000204851","ENSG00000118160","ENSG00000104888","ENSG00000161681",
                      "ENSG00000126583","ENSG00000129990","ENSG00000101327","ENSG00000089199","ENSG00000132639","ENSG00000125851","ENSG00000132821","ENSG00000101438",
                      "ENSG00000124134","ENSG00000124140","ENSG00000101180","ENSG00000185272","ENSG00000040608","ENSG00000167037","ENSG00000100285","ENSG00000100346",
                      "ENSG00000186732","ENSG00000130540","ENSG00000138944","ENSG00000073150","ENSG00000169933","ENSG00000169891","ENSG00000184368","ENSG00000008056",
                      "ENSG00000186462","ENSG00000157502","ENSG00000129682","ENSG00000179542","ENSG00000067842","ENSG00000198910")



oligodendrocyte.genes <- c("RASGRF1","FBXO2","PPP1R16B","TPPP","SEC14L5","TMEM151A","LGI3","TMCC2","HHATL","RAB11FIP4","PDIA2",
                           "HCN2","KCNJ9","DNAJC6","TUBB4A","ADCY5","CSDC2","AC118754.1","PLIN4","HSPA2","PI16","PTGDS","CDK18",
                           "FA2H","AATK","NKX6-2","MAG","PLCH2","FAM131C","PPP1R14A","CHADL","TMEM88B","ABCA2","PLPP2","CERCAM",
                           "BOK","ACP7","CYS1","ANXA3","TNFSF9","AVPI1","MYH11","ADH1B","TUBA4A","CORO6","IL12RB2","TESPA1",
                           "MPP7","RSPO3","KCNJ12","OPN4","MKX","FRAS1","CPNE7")

OPC <- c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4","LPPR1","PTPRZ1","VCAN","DBI","PMP2","CNP","TNS3","LIMA1","CA10","PCDHGC3","CNTN1","SCD5","P2RX7","CADM2","TTYH1","FGF12","TMEM206","NEU4","FXYD6","RNF13","RTKN","GPM6B","LMF1","ALCAM","PGRMC1","HRASLS","BCAS1","RAB31","PLLP","FABP5","NLGN3","SERINC5","EPB41L2","GPR37L1")
NPC1 <- c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST","ASCL1","BTG2","DCX","NXPH1","HN1","PFN2","SCG3","MYT1","CHD7","GPR56","TUBA1A","PCBP4","ETV1","SHD","TNR","AMOTL2","DBN1","HIP1","ABAT","ELAVL4","LMF1","GRIK2","SERINC5","TSPAN13","ELMO1","GLCCI1","SEZ6L","LRRN1","SEZ6","SOX11")
NPC2 <- c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4","DLX5","SOX4","MAP1B","RBFOX2","IGFBPL1","STMN1","HN1","TMEM161B-AS1","DPYSL3","SEPT3","PKIA","ATP1B1","DYNC1I1","CD200","SNAP25","PAK3","NDRG4","KIF5A","UCHL1","ENO2","KIF5C","DDAH2","TUBB2A","LBH","LOC150568","TCF4","GNG3","NFIB","DPYSL5","CRABP1","DBN1","NFIX","CEP170","BLCAP")


oligodendrocyte.genes.ens <- c("ENSG00000205116","ENSG00000149527","ENSG00000116661","ENSG00000185519","ENSG00000116675",
                               "ENSG00000081985","ENSG00000162728","ENSG00000133069","ENSG00000117266","ENSG00000205795",
                               "ENSG00000127824","ENSG00000176720","ENSG00000010282","ENSG00000173175","ENSG00000138759",
                               "ENSG00000138772","ENSG00000196616","ENSG00000171368","ENSG00000164530","ENSG00000146374",
                               "ENSG00000168481","ENSG00000167123","ENSG00000107317","ENSG00000107331","ENSG00000150051",
                               "ENSG00000150054","ENSG00000122375","ENSG00000119986","ENSG00000148826","ENSG00000179292",
                               "ENSG00000135426","ENSG00000126803","ENSG00000058335","ENSG00000185615","ENSG00000103184",
                               "ENSG00000133392","ENSG00000103089","ENSG00000178773","ENSG00000183018","ENSG00000184185",
                               "ENSG00000167549","ENSG00000131242","ENSG00000181409","ENSG00000141934","ENSG00000099822",
                               "ENSG00000167676","ENSG00000104833","ENSG00000125657","ENSG00000105695","ENSG00000167641",
                               "ENSG00000183760","ENSG00000101445","ENSG00000100399","ENSG00000172346")




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





# GSE138794_Diaz ----

## snRNA SF10022 [HQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 65
# gender: Male
# chr7p+, chr8-, chr10-, chr12+


rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119521_SF10022'
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
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))



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

d <- 25
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

levels(object_1$seurat_clusters) <- gsub("^(5)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9|11)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(15)$",paste0("\\1. EN&PE ?"),levels(object_1$seurat_clusters))# TBX3
levels(object_1$seurat_clusters) <- gsub("^(3|6|7|4|1|2|0)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. ?3"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(14)$",paste0("\\1. ?1"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(10)$",paste0("\\1. ?2"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(8)$",paste0("\\1. AC?"),levels(object_1$seurat_clusters))


levels(object_1$seurat_clusters)


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","1. T","2. T","3. T","4. T","6. T","7. T","12. T","13. T",
  "16. ?3",
  "15. EN&PE ?", # wss combi en lage depth
  "8. AC?", 
  "14. ?1", "10. ?2", # Neurons?
  "5. OD",
  "9. TAM", "11. TAM"
))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.png"),width=10,height=8)



#tmp.15 <- FindMarkers(object_1, ident.1 = 15)
head(tmp.15, 20)
# A4GALT,
# CTD-2023N9.2,
# SNTG2,EFCC1,AC016768.1,
# TWIST2  ,IL1R1,  KRT17,
# DLK1,MMP11,LGR6,FOXS1,CARMN,ESM1,FOXC2,FOXL1,
# BOK-AS1,ALDH1A3,SLC6A20,FOXA2,DMKN,IL18R1,SFRP5,LINC00261,
# TBX3

# tmp.12 <- FindMarkers(object_1, ident.1 = 12) # AURKB, SKA1, BUB1B
# tmp.13 <- FindMarkers(object_1, ident.1 = 13) # DTL, MELK
tmp.16 <- FindMarkers(object_1, ident.1 = 16)
# SSTR2, SLC17A7, NHLH1, LHX1, CDH22
View(tmp.16)

# tmp.8 <- FindMarkers(object_1, ident.1 = 8)
# tmp.10 <- FindMarkers(object_1, ident.1 = 10) # SUSD5, ETV1
# tmp.14 <- FindMarkers(object_1, ident.1 = 14)
#tmp.11 <- FindMarkers(object_1, ident.1 = 11)

tmp.8.10.14 <- FindMarkers(object_1, ident.1 = c(8,10,14))
tmp.c.8.10.14 <- FindConservedMarkers(object_1, ident.1 = c(8,10,14))

### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)


# DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
#   geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
#   geom_hline(yintercept=-0.75, linetype="dashed", color = "red")


object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(9,11),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(15),"Pericytes/Endothelial?", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(16),"Special neurons?", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(12,13),"Tumor [cycling?]", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(5),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(8),"Astrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(10,14),"Neurons", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(0:4,6:7),"Tumor", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")


#object_1.SF10022 <- object_1




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

FeaturePlot(object = object_1, features = c("BRINP3"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



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


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
                                        ),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### 5. Oligodendrocytes & OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")




DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)



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

FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))

FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


#### C2 (up) ----


#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)



VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ::T-cell & TAM/MG? ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) ----



DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)




RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])




## snRNA SF10127 [LowQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 44
# gender: Female
# chr7+, chr10-


rm(sid, object_1)
gc()


sid <- 'snRNA_GSM4119522_SF10127'
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


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))





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


## snRNA SF11979sn [LowQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 76
# gender: Female
# chr7p+, chr10-, chr19q-


rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119523_SF11979sn'
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


## snRNA SF12090 [T,OD,MG - LowQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 61
# gender: Male
# chr7+, chr10-


rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119524_SF12090'
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



## snRNA SF12264 :: T,MG,TC+,OD [MHQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 51
# gender: Female
# chr6-, chr7+, chr14p-


rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119525_SF12264'
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
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))


top10 <- head(VariableFeatures(object_1), 10)


print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))






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

# levels(object_1$seurat_clusters) <- gsub("^(2|9|12|8|11|4)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(5|0)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(1|10)$",paste0("\\1. TAM/MG"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("\\1. TC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(7)$",paste0("\\1. T + EN* + TAM/MG*"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(3|6)$",paste0("\\1. T + TC*"),levels(object_1$seurat_clusters))

object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
object_1$class <- as.character(object_1$seurat_clusters)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(2,12,8,9,11,4,3,6,7) & 
                           object_1@reductions$umap@cell.embeddings[,1] < 7
                         ,"T", object_1$class)
object_1$class <- ifelse(object_1$seurat_clusters %in% c(5,0),"OD", object_1$class)
object_1$class <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= 11.75 & 
                             object_1@reductions$umap@cell.embeddings[,1] <= 13 & 
                             object_1@reductions$umap@cell.embeddings[,2] >= 4 & 
                             object_1@reductions$umap@cell.embeddings[,2] <= 6
                           ,"TC", object_1$class)
object_1$class <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= 9 & 
                             object_1@reductions$umap@cell.embeddings[,1] <= 11 & 
                             object_1@reductions$umap@cell.embeddings[,2] >= 3 & 
                             object_1@reductions$umap@cell.embeddings[,2] <= 5
                           ,"BC", object_1$class) # other type?
object_1$class <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= 10 & 
                             object_1@reductions$umap@cell.embeddings[,1] <= 17 & 
                             object_1@reductions$umap@cell.embeddings[,2] >= -7 & 
                             object_1@reductions$umap@cell.embeddings[,2] <= 2
                           ,"TAM", object_1$class)
object_1$seurat_clusters <- as.factor(paste0(object_1$seurat_clusters,". ",object_1$class))



object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "2. T","3. T","4. T","6. T","7. T","8. T","9. T","11. T","12. T",
  "0. OD","5. OD",
  "7. TAM", "1. TAM", "10. TAM", "3. TAM",
  "1. TC","13. TC","2. TC","4. TC","7. TC",
  "3. BC", "7. BC"
))




DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.png"),width=10,height=8)


# tmp.t.cells <- FindMarkers(object_1, ident.1 = c("3. TC", "7. TC"), group.by = "seurat_clusters")
View(tmp.t.cells)




#tmp.t.cells <- FindMarkers(object_1, ident.1 = 13)
#head(tmp.t.cells,25)

# tmp.3.6.7 <- FindMarkers(object_1, ident.1 = c(3,6,7))
#head(tmp.3.6.7, n=20)


### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_vline(xintercept=10, linetype="dashed", color = "blue") +
  geom_vline(xintercept=17, linetype="dashed", color = "blue") +
  geom_hline(yintercept=-7, linetype="dashed", color = "blue") +
  geom_hline(yintercept=2, linetype="dashed", color = "blue")
  
  # geom_vline(xintercept=7, linetype="dashed", color = "red") +
  # 
  # geom_vline(xintercept=11.75, linetype="dashed", color = "blue") +
  # geom_vline(xintercept=13, linetype="dashed", color = "blue") +
  # geom_hline(yintercept=4, linetype="dashed", color = "blue") +
  # geom_hline(yintercept=6, linetype="dashed", color = "blue") +
  # 
  # geom_vline(xintercept=9, linetype="dashed", color = "darkgreen") +
  # geom_vline(xintercept=11, linetype="dashed", color = "darkgreen") +
  # geom_hline(yintercept=3, linetype="dashed", color = "darkgreen") +
  # geom_hline(yintercept=5, linetype="dashed", color = "darkgreen")




object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(2,12,8,9,11,4,3,6,7) & 
                                    object_1@reductions$umap@cell.embeddings[,1] < 7
                                  ,"Tumor", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(5,0),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= 11.75 & 
                                      object_1@reductions$umap@cell.embeddings[,1] <= 13 & 
                                      object_1@reductions$umap@cell.embeddings[,2] >= 4 & 
                                      object_1@reductions$umap@cell.embeddings[,2] <= 6
                                    ,"T-Cells", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= 9 & 
                                      object_1@reductions$umap@cell.embeddings[,1] <= 11 & 
                                      object_1@reductions$umap@cell.embeddings[,2] >= 3 & 
                                      object_1@reductions$umap@cell.embeddings[,2] <= 5
                                    ,"T-Cells", object_1$youri_clusters) # other type?
object_1$youri_clusters <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= 10 & 
                                      object_1@reductions$umap@cell.embeddings[,1] <= 17 & 
                                      object_1@reductions$umap@cell.embeddings[,2] >= -7 & 
                                      object_1@reductions$umap@cell.embeddings[,2] <= 2
                                    ,"TAM/MG", object_1$youri_clusters)



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")




object_1.SF12264 <- object_1




#tmp.u <- FindMarkers(object_1, ident.1 = "Unknown SF12264", group.by = "youri_clusters")
#head(tmp.u, 25)



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


FeaturePlot(object = object_1, features = c("CD163"),order=T) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12"),order=T) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14",order=T) # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"),order=T)
FeaturePlot(object = object_1, features = c("C1QC"),order=T)



#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "CD4",order=T)
FeaturePlot(object = object_1, features = "CD8A",order=T)
FeaturePlot(object = object_1, features = "CD8B",order=T)
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")

#### 3C. B-cells ----

FeaturePlot(object = object_1, features = c("IGLC3"),order=T)
FeaturePlot(object = object_1, features = c("CD19"),order=T)
FeaturePlot(object = object_1, features = c("CD79B"),order=T)


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

#### 5. Oligodendrocytes + OPC (+) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)




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



#### 6B. Pericytes (?) ----

FeaturePlot(object = object_1, features = "RGS5",order=T)
FeaturePlot(object = object_1, features = "PDGFRB",order=T)
FeaturePlot(object = object_1, features = "CD248",order=T)


#### C4 (up) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)



VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) ----


DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----



DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)



RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])




## snRNA SF4400 :: LowQ ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 51
# gender: Female
# chr7+, chr10-



rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119526_SF4400'
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


## snRNA SF4297 :: T,OD,MG+,NE+,AC+ [MQ] ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 39
# gender: Male
# chr3p-, chr7+, chr11q+, chr12-, chr14p+



rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119527_SF4297'
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
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))


top10 <- head(VariableFeatures(object_1), 10)



print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))





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

object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
object_1$seurat_clusters <- as.character(object_1$seurat_clusters)
object_1$seurat_clusters <- ifelse(  object_1@reductions$umap@cell.embeddings[,1] >= -10 & 
                             object_1@reductions$umap@cell.embeddings[,1] <= -9 & 
                             object_1@reductions$umap@cell.embeddings[,2] >= -6.2 & 
                             object_1@reductions$umap@cell.embeddings[,2] <= -5.2
                           , paste0(object_1$seurat_clusters,". PE|EN?"), object_1$seurat_clusters)
object_1$seurat_clusters <- as.factor(object_1$seurat_clusters)
levels(object_1$seurat_clusters) <- gsub("^(1|4|11)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(8)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9|12)$",paste0("\\1. NE"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(6|7)$",paste0("\\1. AC"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("\\1. ?"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|2|3|5|10)$",paste0("\\1. T"),levels(object_1$seurat_clusters))



object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","2. T","3. T","5. T","10. T",
  "6. AC","7. AC",
  "9. NE","12. NE",
  "8. OD",
  "6. PE|EN?",
  "13. ?",
  "1. TAM","4. TAM","11. TAM"
))
object_1$classes <- gsub("^[0-9\\. ]+","",object_1$seurat_clusters)



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


#  geom_vline(xintercept=-10, linetype="dashed", color = "blue") + 
#  geom_vline(xintercept=-9, linetype="dashed", color = "blue") + 
#  geom_hline(yintercept=-6.2, linetype="dashed", color = "blue") + 
#  geom_hline(yintercept=-5.2, linetype="dashed", color = "blue")

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.png"),width=10,height=8)




# tmp.13 <- FindMarkers(object_1, ident.1 = 13) #

# LINC00943,AC007325.2,PRRX2,RP11-407A16.3,TPD52L1,
# GMPR,GPC5-AS1,RP3-395M20.12,NPNT,FXYD1,GJB6,EVC2,
# PLIN5,OTX1,ADIRF,CGNL1,RNF43,TMPRSS3,PLSCR4,SDC4,
# RP11-156K13.1,ALDH1A1,SDS,SHROOM3,LINC01515

# tmp.9 <- FindMarkers(object_1, ident.1 = 9) # RBFOX3
# tmp.12 <- FindMarkers(object_1, ident.1 = 12) # 



### Prepare for integration ----

object_1 <- FindClusters(object_1, resolution = 1)
object_1$youri_clusters <- as.character(object_1$seurat_clusters)



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  geom_hline(yintercept=-5.2, linetype="dashed", color = "blue")





object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(1,4,11),"TAM/MG", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(9,12),"Neurons", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(8),"Oligodendrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(6,7,13) &
                                    object_1@reductions$umap@cell.embeddings[,2] >= -5.2
                                  ,"Astrocytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(6,7,13) &
                                    object_1@reductions$umap@cell.embeddings[,2] <= -5.2
                                  ,"Pericytes", object_1$youri_clusters)
object_1$youri_clusters <- ifelse(object_1$seurat_clusters %in% c(10,3,2,0,5),"Tumor", object_1$youri_clusters)


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters")





object_1.SF4297 <- object_1






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


#FeaturePlot(object = object_1, features = "STMN2") # Tumor
#FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

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
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")

FeaturePlot(object = object_1, features = "CD7")
FeaturePlot(object = object_1, features = "IL7R")


#### 4. Neurons (-) ----



FeaturePlot(object = object_1, features = "RBFOX3",order=T)
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_1, features = "DDN")
FeaturePlot(object = object_1, features = "TNNT2")
FeaturePlot(object = object_1, features = "TMEM130")
FeaturePlot(object = object_1, features = "GABRG2")
FeaturePlot(object = object_1, features = "GABRA1")
FeaturePlot(object = object_1, features = "GABRB2")
FeaturePlot(object = object_1, features = "SATB2") # 
FeaturePlot(object = object_1, features = "SLC17A7") # 
FeaturePlot(object = object_1, features = "GAD1") # 
FeaturePlot(object = object_1, features = "GAD2") # 


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")




tmp <- list('C1'=neuron.genes[neuron.genes %in% NPC2 == F][1:25] ,
            'NPC1'=NPC1[NPC1 %in% NPC2 == F] ,
            'NPC1+2' = intersect(NPC1, NPC2),
            'NPC2'=NPC2[NPC2 %in% c(NPC1, neuron.genes) == F],
            'NPC2 + C1' = intersect(neuron.genes, NPC2),
            'tumor'=c('EGFR','SOX2','GFAP'))

DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C1/NPC] in: ",sid))



# + "PEAR1", "HEYL" , "CFH"
DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))


##### figure S7a ----


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


tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)

DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid,  " (Diaz dataset)"))



ggsave(paste0("output/figures/2022_figure_S7a.pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2)




#### 5. Oligodendrocytes + OPC (+) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")



DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)





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

FeaturePlot(object = object_1, features = c("RGS5","PDGFRB","CD248","PEAR1", "HEYL" , "CFH"))

FeaturePlot(object = object_1, features = c("RGS5"))
FeaturePlot(object = object_1, features = c("PDGFRB"))
FeaturePlot(object = object_1, features = c("CD248"))
FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))



#### C4 (up) :: (-) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)


VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) :: (-) ----


#list('C5'=C5,'TAM/MG'=c("CD163","P2RY12","CD14","ITGB2","C1QC"))

DotPlot(object = object_1, features = C5,  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----
# for sure 2 cells are pericytes

DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248"
#                                                             ,"PEAR1", "HEYL" , "CFH"
                                                             ) ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)



RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])






## snRNA SF6996 ----
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 40
# gender: Female
# chr6q-, chr7+, chr10-, chr19+



rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119528_SF6996'
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

levels(object_1$seurat_clusters) <- gsub("^(0)$",paste0("TAM"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



# tmp.12 <- FindMarkers(object_1, ident.1 = 12) # PLP1, OPALIN, mitochondrial, dubbele kernen?



#### 5. Oligodendrocytes + OPC (+) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)



#### C5 (down) :: TAM? ----

DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)






## snRNA SF9259[S+R] :: T-,OD++,TM ----

#@todo pooling of the 2?

# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 73
# gender: Male
# chr7q+, chr10-, chr13-, chr14-

rm(sid, object_1)
gc()

sid <- 'snRNA_GSM4119529_SF9259R'
object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Couturier")

sid <- 'snRNA_GSM4119530_SF9259S'
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


levels(object_1$seurat_clusters) <- gsub("^(0|1|7|3|4|5)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(6|11)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(10)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(2|8|9)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))


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


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes + OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")





DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)




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



## scRNA SF11979 :: [LowQ] ----
# best wat cellen maar heel weinig reads?
# https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM4119531
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 76
# gender: Female


rm(sid, object_1)
gc()


sid <- 'scRNA_GSM4119531_SF11979'

object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"), gene.column = 1)
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 700,col="red") +
  geom_hline(yintercept = 6500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 30000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 700 &
                     nFeature_RNA < 6500 &
                     nCount_RNA > 100 &
                     nCount_RNA < 30000 &
                     percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))

# top10 <- head(VariableFeatures(object_1), 10)
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


d <- 18


object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)


levels(object_1$seurat_clusters) <- gsub("^(1)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(2|5|4|10|3|8|0|6)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. PE"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9)$",paste0("\\1. IL32 & CD69 +"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(7)$",paste0("\\1. EN?"),levels(object_1$seurat_clusters))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")



#tmp.7 <- FindMarkers(object_1, ident.1 = 7)# PKD1, GOLGA8A, GOLGA8B, SPRY4-IT1, COL16A1, MIRLET7BHG, XIST
#tmp.9 <- FindMarkers(object_1, ident.1 = 9) # IL32, CD69


#### 0. Cycling (+) ----

FeaturePlot(object = object_1, features = "TOP2A") # Tumor



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


#### 2. Astrocyte (?) ----


#FeaturePlot(object = object_1, features = "STMN2") # Tumor
#FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

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
FeaturePlot(object = object_1, features = "CD4",order=T)
FeaturePlot(object = object_1, features = "CD8A",order=T)
FeaturePlot(object = object_1, features = "CD8B",order=T)
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")

#### 3C. B-cells (?) ----

FeaturePlot(object = object_1, features = c("IGLC3"),order=T)
FeaturePlot(object = object_1, features = c("CD19"),order=T)
FeaturePlot(object = object_1, features = c("CD79B"),order=T)



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


FeaturePlot(object = object_1, features = "SYT1") 


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes + OPC (+) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")





DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)




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


FeaturePlot(object = object_1, features = C6)






## scRNA SF11977 :: OD [LowQ] ----
# https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM4119532
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 61
# gender: Female


rm(sid, object_1)
gc()


sid <- 'scRNA_GSM4119532_SF11977'


object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"), gene.column = 1)
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_violin() +
  geom_jitter(cex=1) +
  geom_hline(yintercept = 550,col="red") +
  geom_hline(yintercept = 2500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_violin() +
  geom_jitter(cex=1)  +
  geom_hline(yintercept = 1200,col="red") +
  geom_hline(yintercept = 20000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 550 &
                     nFeature_RNA < 2500 &
                     nCount_RNA > 1200 &
                     nCount_RNA < 20000 &
                     percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))

# top10 <- head(VariableFeatures(object_1), 10)
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

d <- 9
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)

### UMAP clustering ----

object_1 <- RunUMAP(object_1, dims = 1:d)
object_1@meta.data$pt = sapply(strsplit(rownames(object_1@meta.data), "[.]"), "[", 1)

# levels(object_1$seurat_clusters) <- gsub("^(1|4|11)$",paste0("TAM/MG.\\1"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")


tmp.9 <- FindMarkers(object_1, ident.1 = 9) # BCAN, OLIG1, PTPRZ1, SERPINE2 [BCAN]
tmp.4.8 <- FindMarkers(object_1, ident.1 = c(8,4))# [NPC2]
tmp.t <- FindMarkers(object_1, ident.1 = c(0,1,2,3,5,7))




## scRNA SF11956 :: [LowQ] ----
# https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM4119533
# tissue: glioma
# progression: Primary
# genotype/variation: IDHR132H WT GBM
# age: 63
# gender: Male


rm(sid, object_1)
gc()

sid <- 'scRNA_GSM4119533_SF11956'


object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"), gene.column = 1)
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



# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 100 &
#                      nFeature_RNA < 6500 &
#                      nCount_RNA > 100 &
#                      nCount_RNA < 20000 &
#                      percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))

# top10 <- head(VariableFeatures(object_1), 10)
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

# levels(object_1$seurat_clusters) <- gsub("^(1|4|11)$",paste0("TAM/MG.\\1"),levels(object_1$seurat_clusters))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters")




## scRNA SF11644 :: T,TC,OD,TM,EN/PE [MQ] ----
# https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM4119534
# tissue: glioma
# progression: Primary
# genotype/variation: GBM
# age: 57
# gender: Male


rm(sid, object_1)
gc()

sid <- 'scRNA_GSM4119534_SF11644'

object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"), gene.column = 1)
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1300,col="red") +
  geom_hline(yintercept = 7500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 45000,col="red") # + scale_y_log10()



# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 1300 &
#                      nFeature_RNA < 7500 &
#                      nCount_RNA > 100 &
#                      nCount_RNA < 45000 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))

# plot variable features with and without labels



print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))




top10 <- head(VariableFeatures(object_1), 10)
plot1 <- VariableFeaturePlot(object_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #CombinePlots(plots = list(plot1, plot2))
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


levels(object_1$seurat_clusters) <- gsub("^(18)$",paste0("\\1. TC"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(1|8|12|15|12|17)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|2|3|4|5|6|7|9|10|11|13)$",paste0("\\1. T"),levels(object_1$seurat_clusters))

object_1$seurat_clusters <- ifelse(grepl("^(14)$",as.character(object_1$seurat_clusters)) &  object_1@reductions$umap@cell.embeddings[,2] >= -11.4 ,"14. PE",as.character(object_1$seurat_clusters))
object_1$seurat_clusters <- ifelse(grepl("^(14)$",as.character(object_1$seurat_clusters)) &  object_1@reductions$umap@cell.embeddings[,2] <= -11.4 ,"14. EN",as.character(object_1$seurat_clusters))

object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","2. T","3. T","4. T","5. T","6. T","7. T","9. T","10. T","11. T","13. T",
  "16. OD",
  
  "1. TAM","8. TAM","12. TAM","15. TAM","17. TAM",
  "18. TC",
  
  "14. EN","14. PE"
  ))
object_1$cell.type <- as.factor(gsub("[0-9]+\\. (.+)$","\\1",object_1$seurat_clusters))



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  geom_hline(yintercept = -11.4,col="red", lty=2)  +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.png"),width=10,height=8)



# tmp.18 <- FindMarkers(object_1,  ident.1 = 18) CD7, IL7R



#### 1. Tumor (?) ----


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
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2", "GRAP", "", "AURKB")) # Tumor


#### 2. Astrocyte (?) ----


#FeaturePlot(object = object_1, features = "STMN2") # Tumor
#FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

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


#### 3A. TAM/mg/monocytes (?) ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (?) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")

FeaturePlot(object = object_1, features = "CD7")
FeaturePlot(object = object_1, features = "CD52")
FeaturePlot(object = object_1, features = "IL7R")



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


FeaturePlot(object = object_1, features = "SYT1") 


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes + OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")




DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)





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



#### 6B. Pericytes (1) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")


#### C3 :: endothelial (down) ----


C3 <- c('VWF', 'TIE1', 'HIGD1B', 'MMRN1', 'CYSLTR2', 'MMP25','FLT4', 'BCL6B', 'GRAP', 'LAMC3', 'DPEP1', 'PXDNL', 'ANGPT2',
        'PALD1', 'ADGRD1', 'GBP6', 'SLC52A3', 'CLDN5', 'VWA2', 'ABCB1', 'THSD7B', 'SPINK8', 'FOXQ1', 'ZIC3', 'NODAL')

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


ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C3.pdf"),width=7.5, height=3,scale=2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C3.png"),width=7.5, height=3,scale=2)




#### C4 (up) :: (-) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)



VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) :: (-) ----


DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)


RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----




DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)




RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6)




#### CC-2022 (up) ----


DotPlot(object = object_1, features =list('C1'=
                                            results.out |> dplyr::filter(C1.2022) |> dplyr::pull(hugo_symbol)
                                          , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C1] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/2022-",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/2022-",sid,"_C6.png"),width=7.5, height=4,scale=1.2)






## scRNA SF11681 ----
# https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM4658373
# progression: Recurrent
# strain: GBM, IDH1R132H WT
# age: 51
# gender: Male


rm(sid, object_1)
gc()

sid <- 'scRNA_GSM4658373_SF11681'


object_1 <- Read10X(data.dir = paste0("data/scRNA/GSE138794_Diaz/",sid,"/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 6500,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 100,col="red") +
  geom_hline(yintercept = 35000,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 500 &
                     nFeature_RNA < 6500 &
                     nCount_RNA > 100 &
                     nCount_RNA < 35000 &
                     percent.mito < 0.2)


object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "P1" 
object_1$dataset <- as.character(gsub("^[^_]+_","",sid))

# top10 <- head(VariableFeatures(object_1), 10)
# plot variable features with and without labels

print(paste0("Median(nCount_RNA) in ",sid, " = ",round(median(object_1$nCount_RNA))))
print(paste0("Median(nFeature_RNA) in ",sid, " = ",round(median(object_1$nFeature_RNA))))


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

levels(object_1$seurat_clusters) <- gsub("^(7)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(1|4|5|11)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|2|3|6|7|8|10|12)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(9)$",paste0("\\1. TC (+2PE?)"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("\\1. PE?"),levels(object_1$seurat_clusters))


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T","2. T","3. T","6. T","8. T","10. T","12. T",
  "7. OD",
  "1. TAM", "4. TAM", "5. TAM", "11. TAM",
  "9. TC (+2PE?)"))


DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.pdf"),width=10,height=8)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_UMAP.png"),width=10,height=8)




# tmp.9 <- FindMarkers(object_1,  ident.1 = 9) # IL7R
head(tmp.9, 25)



#### Tumor (?) ----


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
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2", "GFAP", "CDK4", "AURKB")) # Tumor


#### 2. Astrocyte (?) ----


#FeaturePlot(object = object_1, features = "STMN2") # Tumor
#FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

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


#### 3A. TAM/mg/monocytes (?) ----


FeaturePlot(object = object_1, features = c("CD163")) # TAM/mg
FeaturePlot(object = object_1, features = c("P2RY12")) # specifiek MG, niet Mac?
FeaturePlot(object = object_1, features = "CD14") # TAM/mg
FeaturePlot(object = object_1, features = c("ITGB2"))
FeaturePlot(object = object_1, features = c("C1QC"))


#### 3B. Til/T-cell (?) ----

FeaturePlot(object = object_1, features = "CD2")
FeaturePlot(object = object_1, features = "CD3D")
FeaturePlot(object = object_1, features = "TRBC2")
FeaturePlot(object = object_1, features = "TRAC")
FeaturePlot(object = object_1, features = "ICOS")
FeaturePlot(object = object_1, features = "GZMA")

FeaturePlot(object = object_1, features = "CD7")
FeaturePlot(object = object_1, features = "CD52")
FeaturePlot(object = object_1, features = "IL7R")


#### 3C. Hematopoietic stem cells? ----

FeaturePlot(object = object_1, features = "HBG1") # Tumor
FeaturePlot(object = object_1, features = "HBG2") # Tumor



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


FeaturePlot(object = object_1, features = "SYT1") 


DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters")

#### 5. Oligodendrocytes + OPC (?) ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")
FeaturePlot(object = object_1, features = "OPALIN")


DotPlot(object = object_1, features =list('C2'=oligodendrocyte.genes , 'OPC'=OPC ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C2/OPC] in: ",sid))

ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.pdf"),width=7.5*1.8, height=3.75,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C2_OPC.png"),width=7.5*1.8, height=3.75,scale=1.2)





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



#### 6B. Pericytes (?) ----

FeaturePlot(object = object_1, features = "RGS5")
FeaturePlot(object = object_1, features = "PDGFRB")
FeaturePlot(object = object_1, features = "CD248")



#### C4 (up) :: (-) ----


DotPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C4.png"),width=6.5, height=4,scale=1.2)




VlnPlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)
RidgePlot(object = object_1, features = c(C4A, C4B), group.by = "seurat_clusters",stack=T)

FeaturePlot(object = object_1, features = C4A)
FeaturePlot(object = object_1, features = C4B)


#### C5 (down) :: (-) ----


DotPlot(object = object_1, features = c(C5),  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C5] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.pdf"),width=6.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C5.png"),width=6.5, height=4,scale=1.2)




RidgePlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C5), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C5)



#### C6 (up) :: (-) ----


DotPlot(object = object_1, features =list('C6'=C6 , 'Peri'=c("RGS5", "PDGFRB", "CD248") ), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C6] in: ",sid))
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.pdf"),width=7.5, height=4,scale=1.2)
ggsave(paste0("output/figures/scRNA/Diaz/",sid,"_C6.png"),width=7.5, height=4,scale=1.2)



RidgePlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)
VlnPlot(object = object_1, features = c(C6), group.by = "seurat_clusters",stack=T)


FeaturePlot(object = object_1, features = C6[1:16])
FeaturePlot(object = object_1, features = C6[17:33])






## Combined ----



object_c <- merge(
  object_1.SF10022,
  y=c(
    object_1.SF12264,
    object_1.SF4297
  ),
  add.cell.ids = c("SF10022", "SF12264", "SF4297"))


levels(as.factor(object_c$dataset))


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


d <- 24



object_c <- FindNeighbors(object_c, dims = 1:d)
object_c <- FindClusters(object_c, resolution = 1.2, algorithm=1)
object_c <- RunUMAP(object_c, dims = 1:d)
object_c@meta.data$pt = sapply(strsplit(rownames(object_c@meta.data), "[.]"), "[", 1)


levels(object_c$seurat_clusters) <- gsub("^(26)$","Pericytes|Endothelial.\\1 SFF10022&SFF122",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(24)$","T-Cells.\\1 SF12264",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(19)$","Neurons.\\1 SF4297",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(16|8)$","Astrocytes\\1 SF4297",levels(object_c$seurat_clusters))
levels(object_c$seurat_clusters) <- gsub("^(3|6|4|27|10)$","Tumor.\\1 SF10022",levels(object_c$seurat_clusters))


DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "seurat_clusters")




DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "youri_clusters") +
  DimPlot(object_c, reduction = "umap", label = TRUE, pt.size = .8, group.by = "dataset")


#### 2. Astrocyte (+) ----


#FeaturePlot(object = object_1, features = "STMN2") # Tumor
#FeaturePlot(object = object_1, features = "ETNPPL") # Tumor

FeaturePlot(object = object_c, features = "GPR98")
FeaturePlot(object = object_c, features = "AQP4")
FeaturePlot(object = object_c, features = "BMPR1B")
FeaturePlot(object = object_c, features = "ETNPPL")
FeaturePlot(object = object_c, features = "GJB6")
FeaturePlot(object = object_c, features = "GJA1")
FeaturePlot(object = object_c, features = "FGFR3")
FeaturePlot(object = object_c, features = "SLC25A18")
FeaturePlot(object = object_c, features = "SLC1A2")
FeaturePlot(object = object_c, features = "SDC4")



#### 3B. Til/T-cell (-) ----

FeaturePlot(object = object_c, features = "CD2")
FeaturePlot(object = object_c, features = "CD3D")
FeaturePlot(object = object_c, features = "TRBC2")
FeaturePlot(object = object_c, features = "TRAC")
FeaturePlot(object = object_c, features = "ICOS")
FeaturePlot(object = object_c, features = "GZMA")



#### 4. Neurons ----


FeaturePlot(object = object_c, features = "RBFOX3")
FeaturePlot(object = object_c, features = "RBFOX1")
FeaturePlot(object = object_c, features = "RBFOX2") # NPC2 ~ Neftel
FeaturePlot(object = object_c, features = "DDN")
FeaturePlot(object = object_c, features = "TNNT2")
FeaturePlot(object = object_c, features = "TMEM130")
FeaturePlot(object = object_c, features = "GABRG2")
#FeaturePlot(object = object_c, features = "GABRA1")
FeaturePlot(object = object_c, features = "GABRB2")



DotPlot(object_c, features=c("RBFOX3","RELN","CNR1","GAD2","OPRK1","GABRB2","RAB3C",
                             "SYT1","KCNC2","ZMAT4","RIMBP2","CHGB","GABRA1","MYT1L","GAD1","TAC3","OLFM3","RP11-679B19.1","TAC1","GABBR2","SCG2","ZCCHC12"), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#### 6A. Endothelial (-) ----

FeaturePlot(object = object_c, features = "ABCB1")
FeaturePlot(object = object_c, features = "CD34")
FeaturePlot(object = object_c, features = "FLT4")
FeaturePlot(object = object_c, features = "TIE1") # meh
FeaturePlot(object = object_c, features = "ITGA1") # endo + peri?


FeaturePlot(object = object_c, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                            "RGS5", "PDGFRB", "CD248"))


DotPlot(object = object_c, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(object = object_c, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2")) # other type of endothelial cells?
DotPlot(object = object_c, features = c("SSTR2", "SST", "LHX1", "LHX6","NHLH1","NHLH2"), group.by = "seurat_clusters") # other type of endothelial cells?




#### 6B. Pericytes ----

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




# 
GSM4119529
