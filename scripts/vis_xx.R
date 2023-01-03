
reorder_levels <- function(seurat_obj, cell_type_ordered) {
  
  df1 <- data.frame(cell_type_ordered = cell_type_ordered) |> 
    dplyr::mutate(i = 1:n())
  
  df2 <- data.frame(seurat_clusters = object_1$seurat_clusters, 
                    cell_type = object_1$cell_type, 
                    annotated_clusters = object_1$annotated_clusters) |> 
    dplyr::distinct() |> 
    dplyr::mutate(seurat_clusters = as.numeric(seurat_clusters)) |> 
    dplyr::mutate(cell_type = as.character(cell_type)) |> 
    dplyr::mutate(annotated_clusters = as.character(annotated_clusters))
  
  p <- unique(df2$cell_type[df2$cell_type %in% df1$cell_type_ordered == F])
  if(length(p) > 0) {
    print(p)
  }
  stopifnot(df2$cell_type %in% df1$cell_type_ordered)
  
  lvls <- df2 |> 
    dplyr::left_join(df1, by=c('cell_type'='cell_type_ordered')) |> 
    dplyr::arrange(i, seurat_clusters)  |> 
    dplyr::pull(annotated_clusters) 

  seurat_obj$annotated_clusters <- factor(seurat_obj$annotated_clusters, levels=lvls)
  
  return(seurat_obj)
}


# sample dataset A ----

sid <- 'GSM4186981_MGH125'

object_1 <- Seurat::Read10X_h5(paste0("data/GSM4186981_Regev/MGH125/",sid,"_fresh_channel1_raw_feature_bc_matrix.h5"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 500,col="red") +
  geom_hline(yintercept = 6000,col="red")


ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 1000,col="red") +
  geom_hline(yintercept = 22500,col="red") # + scale_y_log10()



object_1 <- subset(x = object_1, subset =
                     nFeature_RNA > 500 &
                     nFeature_RNA < 6000 &
                     nCount_RNA > 1000 &
                     nCount_RNA < 22500 &
                     percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 30
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


# 1 = healthy AC
# 2,11 = OD
# 4,5,6 = TAM
# 7,10,9,0,3 = Tum
# 8 = T/T-cell


levels(object_1$seurat_clusters) <- gsub("^(1)$",paste0("\\1. AC"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(2|11)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(4|5|6)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))# TBX3
levels(object_1$seurat_clusters) <- gsub("^(0|3|7|9|10)$",paste0("\\1. T"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(8)$",paste0("\\1. TC"),levels(object_1$seurat_clusters))


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
  "0. T", "3. T", "7. T", "9. T", "10. T",
  "1. AC",
  "8. TC",
  "2. OD",  "11. OD",
  "4. TAM", "5. TAM", "6. TAM"
))




DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



#### 1. Tumor ----


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


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


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








DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")




##### C3/OPC -----


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
  stringr::str_replace("_", " / ")


DotPlot(object = object_1, features = list("C3" = tmp.c3, "OPC" = tmp.opc, "C3+OPC" = tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
  labs(x = paste0("Features [C3/OPC] in: ", sid_print, " (Regev dataset)"))


#ggsave(paste0("output/figures/2022_Figure_S8K.pdf"), width = 7.5 * 1.8, height = 3.75, scale = 1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)



#### 6A. Endothelial ----

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







# sample dataset C ----


sid <- 'H243_GBM'
object_1 <- Read10X(data.dir = paste0("data/",sid,"/H243_filtered_feature_bc_matrix/"))
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="GSE173278")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 30
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


# 1 = healthy AC
# 2,11 = OD
# 4,5,6 = TAM
# 7,10,9,0,3 = Tum
# 8 = T/T-cell


levels(object_1$seurat_clusters) <- gsub("^(2|3|6|14|23)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
levels(object_1$seurat_clusters) <- gsub("^(20)$",paste0("\\1. T&OD doublets"),levels(object_1$seurat_clusters)) # EGFR
levels(object_1$seurat_clusters) <- gsub("^(9)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # OPC - no EGFR or doublets with NE
levels(object_1$seurat_clusters) <- gsub("^(11)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(13|12|4|5|19|10|7|17|16|8|18|22|21)$",paste0("\\1. NE"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(15)$",paste0("\\1. AC"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(24)$",paste0("\\1. EN"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(0|1)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))


object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
"0. OD",
"1. OD",

"15. AC",

"24. EN",

"4. NE",
"5. NE",
"7. NE",
"8. NE",
"10. NE",
"12. NE",
"13. NE",
"16. NE",
"17. NE",
"18. NE",
"19. NE",
"21. NE",
"22. NE",

"2. T"   ,
"3. T",
"6. T",
"14. T",
"20. T",
"23. T",

"9. T" ,

"20. T&OD doublets",

"11. TAM"
 ))




DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters",
        dims=c(1,4),
        ) +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


#m.0.1 <- FindMarkers(object_1, ident.1 = c(0,1))
#m.9 <- FindMarkers(object_1, ident.1 = c(9))


#ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=seurat_clusters)) +
#  geom_jitter() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9))





#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel



DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




##### Figure Sxx - C4 ----


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


sid_print <- paste0("Bolleboom & Gao dataset - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



ggsave(paste0("output/figures/2022_Figure_S7Z.pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)




object_1 <- RunPCA(object_1, features = tmp.c4, reduction.name = 'pca.c4', reduction.key = 'C4',rev.pca=T)
object_1 <- RunPCA(object_1, features = tmp.npc, reduction.name = 'pca.npc', reduction.key = 'NPC')

DotPlot(object = object_1, features =c('C4_1','NPC_1'), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print, " (CPTAC-3 dataset)"))





#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")




##### C3/OPC -----


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


sid_print <- paste0("Bolleboom & Gao dataset - ", sid)



DotPlot(object = object_1, features = list("C3" = tmp.c3, "OPC" = tmp.opc, "C3+OPC" = tmp.c3.opc), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
  labs(x = paste0("Features [C3/OPC] in: ", sid_print, " (Diaz dataset)"))


#ggsave(paste0("output/figures/2022_Figure_S8K.pdf"), width = 7.5 * 1.8, height = 3.75, scale = 1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc, sid_print)



#### 6A. Endothelial ----

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

FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))

FeaturePlot(object = object_1, features = c("COL1A1"))
FeaturePlot(object = object_1, features = c("COL1A2"))

FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?


FeaturePlot(object = object_1, features = "VCAN")
FeaturePlot(object = object_1, features = "TNR")
FeaturePlot(object = object_1, features = "VIM")


# CPTAC-3 - A_CPT0167750015 ----
# https://portal.gdc.cancer.gov/repository?files_offset=20&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22files.cases.primary_site%22%2C%22value%22%3A%5B%22brain%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22transcriptome%20profiling%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Single%20Cell%20Analysis%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22scRNA-Seq%22%5D%7D%7D%5D%7D
# https://portal.gdc.cancer.gov/files/2558a6a4-c187-4dd1-a313-ce2ddf8a97a5
# https://www.sciencedirect.com/science/article/pii/S1535610821000507
# https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-022-01260-4#MOESM1


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'A_CPT0167750015'
origin_file <-  "data/CPTAC-3/A_CPT0167750015/dc8067dd-6147-4059-bae2-7d6d615496b5/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(12,0,2,6,11,14,10,5,4,3,20), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13,1,16,24,21,22), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(15,8), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(17), "NE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(7,9,19,18,23), "TC", object_1$cell_type)
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)



object_1 <- reorder_levels(object_1, c("TC", "TAM", "NE",  "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# m.x <- FindMarkers(object_1, ident.1 = c(22))


#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)



#### 0. Infer CNV ----


library(infercnv)


rm(infercnv_obj)
rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
stopifnot(colnames(rcm) == rownames(af))


fn <- paste0("cache/infercnv_",sid,".Rd")
if(file.exists(fn)) {
  infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")
} else {
  refgroup = data.frame(name = levels(object_1$seurat_clusters)) |> 
    dplyr::filter(
      grepl("^[0-9]+ (TAM|NE|TC|PE|EN|OD)$", name)
    )
  
# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
}
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

# infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)



#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


# CPTAC-3 - B_CPT0167860015 ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'B_CPT0167860015'
origin_file <-  "data/CPTAC-3/B_CPT0167860015/2558a6a4-c187-4dd1-a313-ce2ddf8a97a5/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(17,3,1,4,8,11,18,5,2,15,6,7,0,10,9), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(14,12), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(20,19,16), "NE", object_1$cell_type)
# object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(7,9,19,18,23), "TC", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(21), "TLL1+, non-T", object_1$cell_type) # some sort of strange RBFOX3 negative neuron?
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(22), "GABRG1+, non-T", object_1$cell_type) # some sort of strange RBFOX3 negative neuron?
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)



object_1 <- reorder_levels(object_1, c("", "TC", "TAM", "NE", "GABRG1+, non-T", "TLL1+, non-T", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# m.21 <- FindMarkers(object_1, ident.1 = c(21)) # CN neutral according to infercnv ...
# m.22 <- FindMarkers(object_1, ident.1 = c(22)) # CN neutral according to infercnv ...


#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
# dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
# dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)
DotPlot(object = object_1, features = "DoubletScore")


#### 0. Infer CNV ----



rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



# CPTAC-3 - C_CPT0205450015 [clean] ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'C_CPT0205450015'
origin_file <-  "data/CPTAC-3/C_CPT0205450015/ea055b54-8595-4ecf-a663-3a94871b291b/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(1,4,7,9,2,0,10), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(11,6, 5), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(3), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(12), "NE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(8), "TC", object_1$cell_type)
# object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(22), "GABRG1+, non-T", object_1$cell_type) # some sort of strange RBFOX3 negative neuron?
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)


object_1 <- reorder_levels(object_1, c("", "TC", "TAM", "NE", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# m.x <- FindMarkers(object_1, ident.1 = c(21))


#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
# dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
# dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)
DotPlot(object = object_1, features = "DoubletScore")



#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)



#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


#### 6A. Endothelial ----

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

FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))

FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?



# CPTAC-3 - D_CPT0205570014 ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'D_CPT0205570014'
origin_file <-  "data/CPTAC-3/D_CPT0205570014/c05cae5e-b505-45f9-af2b-b6801b4a21e4/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(10,6,5,3,2,1,9,0), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(8,4), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(12), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(11), "NE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13), "EN|PE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(7), "TC", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(22), "GABRG1+, non-T", object_1$cell_type) # some sort of strange RBFOX3 negative neuron?
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)



object_1 <- reorder_levels(object_1, c("", "EN|PE", "TC", "TAM", "NE", "GABRG1+, non-T", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


#m.l.s <- FindMarkers(object_1, ident.1 = c(10,6,5,3,2,1,9,0), ident.2 = c(8,4)) # CN neutral according to infercnv ...
#m.13 <- FindMarkers(object_1, ident.1 = c(13)) # en and pe genes


#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
# dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
# dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

DotPlot(object = object_1, features = "DoubletScore")
FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)


#### 0. Infer CNV ----


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

# m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


# CPTAC-3 - E_CPT0206000015 [nice sample, no neurons] ----
# infer CNV moet op deze worden gedaan

rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'E_CPT0206000015'
origin_file <-  "data/CPTAC-3/E_CPT0206000015/bb4be4b4-6678-4cfc-92d2-669f71024cab/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13,16,1,9,15), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(4,8,0,3,12,2,14,11), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(5,6,20), "OD", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(11), "NE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(23), "EN", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(21), "PE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(19,7,10,18), "TC", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(22), "BC", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(17), "NRXN1+", object_1$cell_type)
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)



object_1 <- reorder_levels(object_1, c("", "EN","PE", "BC", "TC", "TAM", "NE", "NRXN1+", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



#m.15 <- FindMarkers(object_1, ident.1 = c(15))  # T-cell CD7+
#m.17 <- FindMarkers(object_1, ident.1 = c(17)) # TNR, BCAN, NRNX1
#FeaturePlot(object = object_1, features = "BCAN")
#FeaturePlot(object = object_1, features = "TNR")
#FeaturePlot(object = object_1, features = "NRXN1")


#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
# dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
# dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

DotPlot(object = object_1, features = "DoubletScore")
FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)



#### 0. Infer CNV ----



rm(infercnv_obj)

path1 <- paste0("cache/infercnv_",sid,"_out_pdf")
path2 <- paste0("cache/infercnv_",sid,"_processed.Rds")
path3 <- paste0("cache/infercnv_",sid,".Rds")

if(!file.exists(path1)) {
  if(!file.exists(path2)) {
    if(!file.exists(path3)) {
      rcm <- object_1@assays$RNA@counts
      af <- data.frame(
        V1 = colnames(rcm)) |>
        dplyr::left_join(data.frame(V2 = object_1$annotated_clusters) |>
                           tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
        tibble::column_to_rownames('V1') |> 
        dplyr::mutate(V2 = as.character(V2)) 
      stopifnot(colnames(rcm) == rownames(af))
      
      refgroups = data.frame(annotated_clusters = as.character(object_1$annotated_clusters),
                             cell_type = as.character(object_1$cell_type)) |> 
        dplyr::filter(cell_type %in% c("TAM","NE","TC","BC", "PE","EN","OD")) |> 
        dplyr::pull(annotated_clusters) |> 
        unique()
      
      infercnv_obj = infercnv::CreateInfercnvObject(
        raw_counts_matrix= round(rcm),
        annotations_file = af ,
        gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
        ref_group_names=refgroups # group names for only the non-malignants
      )
      
      saveRDS(infercnv_obj, file=path3)
      rm(rcm, af, infercnv_obj)
      gc()
      
    } else {
      print(paste0("File: ", path3, " already present - skipping re-generation"))
    }
    
    infercnv_obj <- readRDS(file=path3)
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff=0.1,
                                  out_dir = path1,
                                  cluster_by_groups=TRUE,
                                  denoise=T,
                                  HMM=T,
                                  output_format="pdf")
    saveRDS(infercnv_obj, file=path2)
    rm(infercnv_obj)
    gc()
  }
}




#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)








FeaturePlot(object = object_1, features = "ETNPPL") # Tumor
FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 3C. B-cells ----


FeaturePlot(object = object_1, features = c("IGLC3"), order=T)
FeaturePlot(object = object_1, features = c("CD19"), order=T)
FeaturePlot(object = object_1, features = c("CD79B"), order=T)


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T, order=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F, order=T)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)


#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


#### 6A. Endothelial ----

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

FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))

FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?



# CPTAC-3 - F_CPT0206770002 - no neurons ----



rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'F_CPT0206770002'
origin_file <-  "data/CPTAC-3/F_CPT0206770002/1ef4a07d-2a1a-42b4-9d5f-fdf45c2345a0/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13,16,1,9,15), "T", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(4,8,0,3,12,2,14,11), "TAM", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(5,6,20), "OD", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(11), "NE", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(23), "EN", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(21), "PE", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(19,7,10,18), "TC", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(22), "BC", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(17), "NRXN1+", object_1$cell_type)
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)



object_1 <- reorder_levels(object_1, c("", "EN","PE", "BC", "TC", "TAM", "NE", "NRXN1+", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)




#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)


# CPTAC-3 - G_CPT0206880004 ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'G_CPT0206770002'
origin_file <-  "data/CPTAC-3/G_CPT0206880004/71494916-46c2-4235-9225-ffe218bbdf07/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(0,3,8,1,13), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(4,5,7), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(11,2), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(14,9,6,15,16,12), "NE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(18), "EN", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(19), "PE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(10), "TC", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(20), "AC", object_1$cell_type)# indeed CNV neutral
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(17), "TNR+, non-T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(21), "TACSTD2+, non-T", object_1$cell_type)
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)


object_1 <- reorder_levels(object_1, c("", "EN","PE", "BC", "TC", "TAM", "NE", "TNR+, non-T", "TACSTD2+, non-T", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# m.17 <- FindMarkers(object_1, ident.1 = c(17)) # non-T - AC124254.1, GALR1, COL9A1, SMOC1, CA10, AL450345.2, PCDH15, PDGFRA, PLPP4, BX284613.2, MEGF11, CHST9, FERMT1, TNR, AC023282.1, ITGA8, ADARB2, SCN9A, AFAP1L2, CRISPLD2, XYLT1, HIF3A, AL512308.1, MYRFL, LINC02283
# m.20 <- FindMarkers(object_1, ident.1 = c(20)) # AC? - ZNF98, MRVI1, AL137139.1, PRODH, SLC7A10, RNF219-AS1, AC018618.1, AC073941.1, SLC14A1, CD38, ETNPPL, LINC00499, RGS20, APLNR, AC012405.1, SLC39A12, FAM189A2, HPSE2, AC104574.2, RANBP3L, TPD52L1, FAM107A, AC097518.2, LINC00299, AC087482.1
# m.21 <- FindMarkers(object_1, ident.1 = c(21)) # ? TACSTD2, UGT2B7, PKHD1, ENSG00000272398, AOC1, LCN2, FXYD4, MMP7, AQP2, KRT18, KRT19, SLPI, WFDC2, CLDN8, FOLR1, SLC14A2, DEFB1, S100A11, PAX8, CDH16, EPCAM, AGR2, AC023421.1, TFPI2, KRT7

head(m.17, n=25)
head(m.20, n=25)
head(m.21, n=25)



#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
# dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
# dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

DotPlot(object = object_1, features = "DoubletScore")
FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)





#### 0. Infer CNV ----


rm(infercnv_obj)

path1 <- paste0("cache/infercnv_",sid,"_out_pdf")
path2 <- paste0("cache/infercnv_",sid,"_processed.Rds")
path3 <- paste0("cache/infercnv_",sid,".Rds")

if(!file.exists(path1)) {
  if(!file.exists(path2)) {
    if(!file.exists(path3)) {
      rcm <- object_1@assays$RNA@counts
      af <- data.frame(
        V1 = colnames(rcm)) |>
        dplyr::left_join(data.frame(V2 = object_1$annotated_clusters) |>
                           tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
        tibble::column_to_rownames('V1') |> 
        dplyr::mutate(V2 = as.character(V2)) 
      stopifnot(colnames(rcm) == rownames(af))
      
      refgroups = data.frame(annotated_clusters = as.character(object_1$annotated_clusters),
                             cell_type = as.character(object_1$cell_type)) |> 
        dplyr::filter(cell_type %in% c("TAM","NE","TC","BC", "PE","EN","OD")) |> 
        dplyr::pull(annotated_clusters) |> 
        unique()
      
      infercnv_obj = infercnv::CreateInfercnvObject(
        raw_counts_matrix= round(rcm),
        annotations_file = af ,
        gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
        ref_group_names=refgroups # group names for only the non-malignants
      )

      saveRDS(infercnv_obj, file=path3)
      rm(rcm, af, infercnv_obj)
      gc()
      
    } else {
      print(paste0("File: ", path3, " already present - skipping re-generation"))
    }
    
    infercnv_obj <- readRDS(file=path3)
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff=0.1,
                                  out_dir = path1,
                                  cluster_by_groups=TRUE,
                                  denoise=T,
                                  HMM=T,
                                  output_format="pdf")
    saveRDS(infercnv_obj, file=path2)
    rm(infercnv_obj)
    gc()
  }
}






#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=F, order=T)
FeaturePlot(object = object_1, features = "AQP4", label=F, order=T)
FeaturePlot(object = object_1, features = "TIMP3", label=F, order=T)
FeaturePlot(object = object_1, features = "NTRK2", label=F, order=T)
FeaturePlot(object = object_1, features = "KCNN3", label=F, order=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=F, order=T)




#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 3C. B-cells ----


FeaturePlot(object = object_1, features = c("IGLC3"), order=T)
FeaturePlot(object = object_1, features = c("CD19"), order=T)
FeaturePlot(object = object_1, features = c("CD79B"), order=T)


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T, order=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F, order=T)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144", order=T, label=F)
FeaturePlot(object = object_1, features = "TMEM144", order=T, label=T)
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



#### 6A. Endothelial ----


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

FeaturePlot(object = object_1, features = c("PEAR1"))
FeaturePlot(object = object_1, features = c("HEYL"))
FeaturePlot(object = object_1, features = c("CFH"))

FeaturePlot(object = object_1, features = "ITGA1") # endo + peri?

DotPlot(object = object_1, features = c("ABCB1", "CD34", "FLT4", "TIE1", "ITGA1",
                                        "RGS5", "PDGFRB", "CD248"), group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# CPTAC-3 - H_CPT0224600013 ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'H_CPT0224600013'
origin_file <-  "data/CPTAC-3/H_CPT0224600013/3351eae3-830b-4cdd-91aa-0f10b32ca2e8/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version


mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----


object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(4,6,2,11,7,5,8,1,0,9), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(10,3,14), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(12), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13), "NE", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(18), "EN", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(19), "PE", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(10), "TC", object_1$cell_type)
#object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(20), "AC", object_1$cell_type)# indeed CNV neutral
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)


object_1 <- reorder_levels(object_1, c("", "EN","PE", "BC", "TC", "TAM", "NE", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)









FeaturePlot(object = object_1, features = "ETNPPL") # Tumor
FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T, order=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F, order=T)

FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid_print), " (CPTAC-3 dataset)"))



ggsave(paste0("output/figures/2022_Figure_S7_ext_CPTAC-3_",sid,".pdf"),width=7.5*3, height=4,scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)


#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1",label=T)





# CPTAC-3 - I_CPT0167970014 - no neurons ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'I_CPT0167970014'
origin_file <-  "data/CPTAC-3/I_CPT0167970014/c59d4a7b-0a9c-4093-a801-b162d43f5027/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)








mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
# levels(object_1$seurat_clusters) <- gsub("^(13)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) # EGFR
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(22)$",paste0("\\1. GABRG1+"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# 
# m.22 <- FindMarkers(object_1, ident.1 = c(22))
# m.21 <- FindMarkers(object_1, ident.1 = c(21))
# 

#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)









FeaturePlot(object = object_1, features = "ETNPPL") # Tumor
FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----

FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel








DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





##### Figure Sxx - C4 ----

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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))





#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



# CPTAC-3 - J_CPT0189650015 [+] ----
# metadata cannot be trusted, data has same md5sum as M/CPT0207030018

rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'J_CPT0189650015'
origin_file <-  "data/CPTAC-3/J_CPT0189650015/62e997e9-0f83-4815-ab41-586d3af350a2/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)








mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(22)$",paste0("\\1. GABRG1+"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# 
# m.22 <- FindMarkers(object_1, ident.1 = c(22))
# m.21 <- FindMarkers(object_1, ident.1 = c(21))
# 

#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----

m.ac <- FindMarkers(object_1, ident.1 = c(1))


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)









FeaturePlot(object = object_1, features = "ETNPPL") # Tumor
FeaturePlot(object = object_1, features = "GPR98")
FeaturePlot(object = object_1, features = "BMPR1B")
FeaturePlot(object = object_1, features = "ETNPPL")
FeaturePlot(object = object_1, features = "GJB6")
FeaturePlot(object = object_1, features = "FGFR3")
FeaturePlot(object = object_1, features = "SLC25A18")
FeaturePlot(object = object_1, features = "SLC1A2")
FeaturePlot(object = object_1, features = "SDC4")


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----

FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel








DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





##### Figure Sxx - C4 ----

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


sid_print <- paste0("CPTAC-3 - ", sid)

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))





#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



# CPTAC-3 - K_CPT0125220004 [+] ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'K_CPT0125220004'
origin_file <-  "data/CPTAC-3/K_CPT0125220004/5b383196-0109-4452-8bce-459308d5b393/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)








mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
levels(object_1$seurat_clusters) <- gsub("^(7|12)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(22)$",paste0("\\1. GABRG1+"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


# 
# m.22 <- FindMarkers(object_1, ident.1 = c(22))
# m.21 <- FindMarkers(object_1, ident.1 = c(21))
# 

#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----

FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)
FeaturePlot(object = object_1, features = "RBFOX1")
FeaturePlot(object = object_1, features = "RBFOX2") # NPC2 ~ Neftel








DotPlot(object = object_1, features = c("RBFOX3",
                                        "CNR1","SYT1","SYNPR","GABRA1","RELN,","VIP",
                                        "CCT2","RUFY2","UBN2","ATP6V1H","HSPA4L","NASP","GNAO1","RAB6B","HLF","SLC25A36"
),group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





##### Figure Sxx - C4 ----

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


sid_print <- paste0("CPTAC-3 - ", gsub("^.+_","",sid))

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))





#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")




# CPTAC-3 - L_CPT0205890014 [+ but duplicates & noisy] ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'L_CPT0205890014'
origin_file <-  "data/CPTAC-3/L_CPT0205890014/2f54e97e-0961-410a-b429-f42d212c4755/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)








mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(22)$",paste0("\\1. GABRG1+"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



#### 0. Find Doublets ----


dbl <- Seurat::as.SingleCellExperiment(object_1)
dbl.out <- scDblFinder::findDoubletClusters(dbl, clusters=dbl$seurat_clusters)
dbl.out[1,]


set.seed(123456)
top.mam <- scran::getTopHVGs(dbl, prop=0.1)
dbl.dens <- scDblFinder::computeDoubletDensity(dbl,
                                               subset.row=top.mam, 
                                               d=ncol(reducedDim(dbl)))

dbl$DoubletScore <- dbl.dens
scater::plotUMAP(dbl, colour_by="DoubletScore")

stopifnot(colnames(dbl) == colnames(object_1))
object_1$DoubletScore <- dbl$DoubletScore
rm(dbl, dbl.out, dbl.dens, top.man)
gc()

FeaturePlot(object = object_1, features = "DoubletScore", label=T)
FeaturePlot(object = object_1, features = "DoubletScore", label=F)



#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)

FeaturePlot(object = object_1[,object_1$DoubletScore >= 5], features = "RBFOX3", label=F)
FeaturePlot(object = object_1[,object_1$DoubletScore >= 5], features = "CD163", label=F)
FeaturePlot(object = object_1[,object_1$DoubletScore >= 5], features = "TMEM144", label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", gsub("^.+_","",sid))

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")






# CPTAC-3 - M_CPT0207030018 [+ & clean] ----
# file has same MD5sum as J?

rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'M_CPT0207030018'
origin_file <-  "data/CPTAC-3/M_CPT0207030018/c868a478-c425-485b-9cfc-4ba7fe5365bc/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)








mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. GABRG1+, non-T"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


m.16 <- FindMarkers(object_1, ident.1 = c(16))
head(m.16, n =25)


#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=F) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)

FeaturePlot(object = object_1[,object_1$DoubletScore >= 5], features = "RBFOX3", label=F)
FeaturePlot(object = object_1[,object_1$DoubletScore >= 5], features = "CD163", label=F)
FeaturePlot(object = object_1[,object_1$DoubletScore >= 5], features = "TMEM144", label=F)


##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", gsub("^.+_","",sid))

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")



# CPTAC-3 - N_CPT0207030018 - insufficient neurons ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'N_CPT0168830014'
origin_file <-  "data/CPTAC-3/N_CPT0168830014/602a750e-a4ae-4089-b840-b8790d52179f/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)



mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 




object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


object_1 <- RunUMAP(object_1, dims = 1:d)



## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
#levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
#levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. GABRG1+, non-T"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


m.16 <- FindMarkers(object_1, ident.1 = c(16))
head(m.16, n =25)



#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)



# CPTAC-3 - O_CPT0168380014 - no neurons ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'O_CPT0168380014'
origin_file <-  "data/CPTAC-3/O_CPT0168380014/2ecf2dab-93f2-4b91-b6cd-acdd0a310b64/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)



mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)


object_1 <- RunUMAP(object_1, dims = 1:d)



## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
# levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. GABRG1+, non-T"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


m.16 <- FindMarkers(object_1, ident.1 = c(16))
head(m.16, n =25)




#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=F) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")





#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)


# CPTAC-3 - P_CPT0228220004 - no neurons ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'P_CPT0228220004'
origin_file <-  "data/CPTAC-3/P_CPT0228220004/052a06ca-73ea-4a40-aa6a-7118deee3078/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)




mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(6)$",paste0("\\1. T|AC"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. GABRG1+, non-T"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


m.16 <- FindMarkers(object_1, ident.1 = c(16))
head(m.16, n =25)


#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=F) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)



##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", gsub("^.+_","",sid))

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144", label=T)
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")

# CPTAC-3 - Q_CPT0087680014 ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'Q_CPT0087680014'
origin_file <-  "data/CPTAC-3/Q_CPT0087680014/299b7069-9be9-40ae-99a5-dc9d801cca0d/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)




mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
#levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(19)$",paste0("\\1. NE"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. GABRG1+, non-T"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


m.16 <- FindMarkers(object_1, ident.1 = c(16))
head(m.16, n =25)


#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=F) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)



##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", gsub("^.+_","",sid))

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144", label=T)
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")


# CPTAC-3 - R_CPT0168080014 -noisy nascent RNA sample, odd expression of VIM etc. ----


rhdf5::h5closeAll()
rm(object_1, lfile)
gc()


sid <- 'R_CPT0168080014'
origin_file <-  "data/CPTAC-3/R_CPT0168080014/9427ea4c-ac09-48b6-93ac-fcf628d0da71/seurat.loom"
target_file <- paste0("tmp/CPTAC-3_",sid,"_seurat.loom")
if(!file.exists(target_file)) {
  print("Copying loom file to rw cache")
  file.copy(origin_file, target_file)
  Sys.chmod(target_file,mode="0666")
}

lfile <- loomR::connect(filename = target_file, mode = "r+", skip.validate = TRUE)# skip.validate is needed because the provided loom file is in some old specification/version



mat <- t(lfile[['matrix']][,])

tt <- gencode.31 |> 
  dplyr::mutate(ENSG = gsub("\\..+$","",ENSG)) |> 
  dplyr::filter(!duplicated(ENSG))

gid <- lfile[["row_attrs/Gene"]][]
gid <- gsub("\\..+$","",gid)
gid <- data.frame(gid = gid)
gid <- gid |> 
  dplyr::left_join(tt, by=c('gid'='ENSG')) |> 
  dplyr::mutate(gid.new = ifelse(!is.na(GENE), GENE, gid))

gene.order <- gid |> 
  dplyr::select(gid.new, V1, V4, V5) |> 
  dplyr::rename(gene = gid.new) |> 
  dplyr::rename(chr = V1) |> 
  dplyr::rename(start = V4) |> 
  dplyr::rename(end = V5) |> 
  dplyr::filter(!is.na(chr)) |> 
  dplyr::arrange(order(gtools::mixedorder(chr)), start, end)

gid <- gid |> 
  dplyr::pull(gid.new)

rownames(mat) <- gid

barcode <- lfile[["col_attrs/CellID"]][]
colnames(mat) <- barcode



object_1 <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project="GSM4186981_MGH125")
rm(mat)




mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


# ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01) +
#   geom_hline(yintercept = 500,col="red") +
#   geom_hline(yintercept = 6000,col="red")
# 
# ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
#   geom_jitter(cex=0.01)  +
#   geom_hline(yintercept = 1000,col="red") +
#   geom_hline(yintercept = 22500,col="red") # + scale_y_log10()
# 
# object_1 <- subset(x = object_1, subset =
#                      nFeature_RNA > 500 &
#                      nFeature_RNA < 6000 &
#                      nCount_RNA > 1000 &
#                      nCount_RNA < 22500 &
#                      percent.mito < 0.2)



object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(object_1), 10)


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

d <- 35
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm=1)
head(Idents(object_1), 20)



object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
# 
# levels(object_1$seurat_clusters) <- gsub("^(0|1|2|3|4|5|7|8|9|10|11|17|18|21|15)$",paste0("\\1. T"),levels(object_1$seurat_clusters)) # EGFR
#levels(object_1$seurat_clusters) <- gsub("^(14|15)$",paste0("\\1. NE"),levels(object_1$seurat_clusters)) 
# levels(object_1$seurat_clusters) <- gsub("^(12)$",paste0("\\1. OD"),levels(object_1$seurat_clusters))
levels(object_1$seurat_clusters) <- gsub("^(19)$",paste0("\\1. NE"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(12|14)$",paste0("\\1. TAM"),levels(object_1$seurat_clusters))
# levels(object_1$seurat_clusters) <- gsub("^(16)$",paste0("\\1. GABRG1+, non-T"),levels(object_1$seurat_clusters)) # also PLA2G5+
# 
# 
# object_1$seurat_clusters <- factor(object_1$seurat_clusters, levels=c(
#   "0. T" , "1. T" , "2. T" , "3. T" , "4. T" ,
#   "5. T" , "7. T" , "8. T" , "9. T" ,
#   "10. T" , "11. T" , "15. T" , "17. T" , "18. T" , "21. T",
#   "6. T|AC" ,
#   "22. GABRG1+",
#   "13. OD" ,
#   "16. NE", "19. NE", "20. NE" ,
#   "12. TAM", "14. TAM" 
# ))
# 



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



DimPlot(object_1, reduction = "pca", label = TRUE, pt.size = .6, group.by = "seurat_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


m.16 <- FindMarkers(object_1, ident.1 = c(16))
head(m.16, n =25)


#### 0. Infer CNV ----


library(infercnv)


rcm <- object_1@assays$RNA@counts
af <- data.frame(
  V1 = colnames(rcm)) |>
  dplyr::left_join(data.frame(V2 = object_1$seurat_clusters) |>
                     tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
  tibble::column_to_rownames('V1') |> 
  dplyr::mutate(V2 = as.character(V2)) 
#dplyr::mutate(V2 = gsub("^[0-9]+\\. ","",V2))
stopifnot(colnames(rcm) == rownames(af))



# infercnv_obj = infercnv::CreateInfercnvObject(
#   raw_counts_matrix= round(rcm),
#   annotations_file = af ,
#   gene_order_file=gene.order |> dplyr::filter(!duplicated(gene)) |> tibble::column_to_rownames('gene'),
#   ref_group_names=c(
#     "14. TAM","12. TAM","20. NE","19. NE","16. NE","13. OD"
#     #"TAM","NE","OD"
#   ) # group names for only the non-malignants
#   )
# 
# saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015.Rds")

infercnv_obj <- readRDS(file="cache/infercnv_CPT0167860015.Rds")


infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir = "cache/infercnv_CPT0167860015_out_pdf",
                              cluster_by_groups=TRUE,
                              denoise=T,
                              HMM=T,
                              output_format="pdf"
)

#saveRDS(infercnv_obj , file="cache/infercnv_CPT0167860015_processed.Rds")



infercnv::plot_cnv(infercnv_obj,
                   out_dir = "tmp/infercnv_CPT0167860015_out/pdf/",
                   output_format = "pdf")


#### 1. Tumor ----


FeaturePlot(object = object_1, features = "ETV1") # Tumor
FeaturePlot(object = object_1, features = "CDK4") # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=T) # Tumor
FeaturePlot(object = object_1, features = "EGFR",label=F) # Tumor

FeaturePlot(object = object_1, features = "S100B") # Tumor/AC
FeaturePlot(object = object_1, features = "GFAP") # Tumor/AC
FeaturePlot(object = object_1, features = "OLIG1") # Tumor/OPC+NPC1
FeaturePlot(object = object_1, features = "VIM") # Tumor/MES

FeaturePlot(object = object_1, features = c("HSPA1A","HSPA1B","VEGFA")) # Apoptotic Tumor?
FeaturePlot(object = object_1, features = c("RRM2","PCNA","KIAA0101")) # G1/S
FeaturePlot(object = object_1, features = c("CCNB1","CDC20","CCNB2")) # G2/M
FeaturePlot(object = object_1, features = c("TMPO")) # G2/M
FeaturePlot(object = object_1, features = c("KIF2C","NUF2","ASPM","NEK2","CENPA","CKAP2L","SGOL1","CENPE","CCNA2","PBK","MKI67","CDCA3","NUSAP1","CCNB2","KIF23"))


FeaturePlot(object = object_1, features = c("EREG"))

# succes met vinden van een marker
FeaturePlot(object = object_1, features = c("EGFR","OLIG1","TMPO","VIM","STMN2",   "AURKB")) # Tumor


FeaturePlot(object = object_1, features = "CDH10")



#### 2. Astrocyte ----


FeaturePlot(object = object_1, features = "GJA1", label=T)
FeaturePlot(object = object_1, features = "AQP4", label=T)
FeaturePlot(object = object_1, features = "TIMP3", label=T)
FeaturePlot(object = object_1, features = "NTRK2", label=T)
FeaturePlot(object = object_1, features = "KCNN3", label=T)
FeaturePlot(object = object_1, features = "SLC14A1", label=T)


#### 3A. TAM ----


FeaturePlot(object = object_1, features = c("CD163"),label=T) # TAM/mg
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


#### 4. Neurons ----


FeaturePlot(object = object_1, features = "RBFOX3", label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)



##### Figure Sxx - C4 ----


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


sid_print <- paste0("CPTAC-3 - ", gsub("^.+_","",sid))

tmp <- list('C4'=tmp.c4,
            'NPC1'=tmp.npc1,
            'NPC1+2'=tmp.npc1.2,
            'NPC2'=tmp.npc2,
            'NPC2 + C4' = tmp.c4.npc2)


DotPlot(object = object_1, features = tmp, group.by = "seurat_clusters",
        cols = c("lightgrey", "purple")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ",sid_print))



#### 5. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144", label=T)
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")




