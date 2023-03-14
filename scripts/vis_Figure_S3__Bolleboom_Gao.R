#!/usr/bin/env R 


if(!exists('results.out')) {
  source('scripts/load_results.out.R')
}



# CX - OPC ----

ctx_G_CPT0206770004 = c("AC124254.1", "GALR1", "COL9A1", "SMOC1", "CA10", "AL450345.2", "PCDH15", "PDGFRA", "PLPP4", "BX284613.2", "MEGF11", "CHST9", "FERMT1", "TNR", "AC023282.1", "ITGA8", "ADARB2", "SCN9A", "AFAP1L2", "CRISPLD2", "XYLT1", "HIF3A", "AL512308.1", "MYRFL", "LINC02283")
ctx_K_CPT0125220004 = c("COL9A1", "PDGFRA", "GPR17", "ALDH1A3", "SMOC1", "PRKG2", "CA10", "PCDH15", "COL11A1", "CACNG5", "LINC02283", "AL512308.1", "TTLL6", "LUZP2", "LINC02223", "AC022034.3", "AC020584.1", "NOS1", "ADARB2", "NEU4", "SCN9A", "AC062021.1", "KCNJ16", "ATP2C2", "ADAMTS17")
all <- c(ctx_G_CPT0206770004, ctx_K_CPT0125220004)
dup <- all[duplicated(all)]

cell_type_opc <- list(
  'shared' = all[duplicated(all)],
  'G_CPT0206770004' = ctx_G_CPT0206770004[ctx_G_CPT0206770004 %in% dup == F][1:15],
  'K_CPT0125220004' = ctx_K_CPT0125220004[ctx_K_CPT0125220004 %in% dup == F][1:15]
)

rm(ctx_G_CPT0206770004, ctx_K_CPT0125220004, all, dup)



# CY ----

# astrocyte subtype? https://www.nature.com/articles/s41467-019-14198-8 fig3
cty_J_CPT0189650015 <- c("LRRC3B","GABRG1","ETNPPL","LINC00499","TPD52L1","AC002429.2","SLC39A12","EMX2","SPON1","WIF1","LINC00943","GJB6","ZNF98","FOXG1","SLC7A10","FOXG1-AS1","AL158058.1","AC023095.1","HPSE2","ABCC9","PRDM16","SIRPB3P","MRVI1","AC019270.1","AC106845.2")


all <- c(cty_J_CPT0189650015)
dup <- all[duplicated(all)]

cell_type_y <- list(
  'shared' = all[duplicated(all)],
  'J_CPT0189650015' = cty_J_CPT0189650015[cty_J_CPT0189650015 %in% dup == F][1:15]
)

rm(cty_J_CPT0189650015, all, dup)



# reorder func ----

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




# Bolleboom & GAO  H243-GBM ----
# chr7 gain, chr10 loss, chr9.q loss, chr13.p loss (InferCNV)
sid <- 'H243_GBM'
object_1 <- readRDS(file="tmp/seurat_bolleboom.2.Rds")

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
object_1 <- RunUMAP(object_1, dims = 1:d)


## clustering & annotation ----
### round 1 with doublets ----

object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1,random.seed=123456)


object_1$cell_type = ""
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(2,3,6,14,23), "T", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(9), "OPC", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(20), "Doublets", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(11), "TAM", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(0,1), "OD", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(13,12,4,5,19,10,7,17,16,8,18,22,21), "NE", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(15), "AC", object_1$cell_type)
object_1$cell_type = ifelse(object_1$seurat_clusters %in% c(24), "EN", object_1$cell_type)
object_1$annotated_clusters = paste0(object_1$seurat_clusters,". ",object_1$cell_type)



object_1 <- reorder_levels(object_1, c("", "Doublets" , "EN", "TC", "TAM", "NE", "OPC", "OD", "AC", "T"))

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_type") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)



### round 2 after doublet removal ----


object_1 <- object_1[,object_1[['cell_type']] != "Doublets"]

# object_1$seurat_clusters <- NULL
# object_1$cell_type <- NULL
# object_1$annotated_clusters <- NULL
object_1 <- FindNeighbors(object_1, dims = 1:d)
object_1 <- RunUMAP(object_1, dims = 1:d)
object_1 <- FindClusters(object_1, resolution = 1, algorithm = 1, random.seed=123456)



DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "annotated_clusters") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid)


#### FF] Figure S3-p06 - UMAP ----

DimPlot(object_1, reduction = "umap", label = TRUE, pt.size = .6, group.by = "cell_type") +
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) +
  labs(subtitle=sid) +
  scale_color_manual(values=c('T'='#F8766D', # neat reddish
                              'NE' = '#00B6EB', # neat blue
                              'OD'='#53B400', # neat green
                              'OPC'='gray65',
                              'EN'='gray65',
                              'TAM'='gray65',
                              'AC'='gray65',
                              'EN'='gray65'))


ggsave(paste0("output/figures/2022_Figure_S3-p06_", sid, "_UMAP_Bolleboom_Gao.pdf"), 
       width = 3.75 * 0.93 * (87/62) * 1.2 * 1.2,
       height = 3.75 * 0.93 * (87/62),scale=1.2)


#saveRDS(object_1, file="tmp/seurat_bolleboom.2.Rds")


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
        dplyr::filter(cell_type %in% c("TAM","NE","TC","BC", "PE","EN","OD", "OPC")) |> 
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


system(paste0("rm ",path1,"/*.txt"))
system(paste0("rm ",path1,"/*.dat"))
system(paste0("rm ",path1,"/*_obj"))



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


DotPlot(object = object_1, features = list('canonical'= c("GJA1","AQP4","TIMP3","NTRK2","KCNN3","SLC14A1"),
                                           'GABGRG1'= as.character(unlist(cell_type_y)))
        , group.by = "annotated_clusters",
        #cols = c("lightgrey", "purple")
) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("^.+_","",sid), " (CPTAC-3 dataset)"))


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

FeaturePlot(object = object_1, features = c("RBFOX3"), label=T)
FeaturePlot(object = object_1, features = "RBFOX3", label=F)
FeaturePlot(object = object_1, features = c("RBFOX3","DoubletScore"), label=F)





##### FF] Figure S3-p05 - C4 ----


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
            'NPC2'=tmp.npc2
            #'NPC2 + C4' = tmp.c4.npc2
            )


DotPlot(object = object_1, features = tmp, group.by = "annotated_clusters",
        cols = c("lightgrey", "purple")
        ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: ", gsub("_","-",sid), " (Bolleboom-Gao dataset)")) +
  guides(y = "none", y.sec = guide_axis())


ggsave(paste0("output/figures/2022_Figure_S3-p05_Bolleboom-Gao_",sid,"_C4.pdf"),
       width=7.5 * 1.8 * (612 / 744.4529) * (392.0837 / 412.4028) * 1.14,
       height=3.75 * 0.93 * (87/62), scale=1.2)
rm(tmp.c4, tmp.c4.npc2, tmp.npc1, tmp.npc1.2, tmp.npc2, sid_print)



#### 5A. Oligodendrocytes ----


FeaturePlot(object = object_1, features = "TMEM144")
FeaturePlot(object = object_1, features = "TMEM125")
FeaturePlot(object = object_1, features = "MOG")
FeaturePlot(object = object_1, features = "PLP1")




##### FF] Figure S3-p04 - C3 -----


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




DotPlot(object = object_1, features = list("C3" = tmp.c3, "OPC-like" = tmp.opc
                                           #, "C3 + OPC-like" = tmp.c3.opc
                                           ), group.by = "annotated_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
  labs(x = paste0("Features [C3/OPC-like] in: ", gsub("^.+_","",sid), " (Bolleboom-Gao dataset)"))


ggsave(paste0("output/figures/2022_Figure_S3-p04_Bolleboom-Gao_",sid,"_C3.pdf"),
       width=7.5 * 1.8 * (612 / 744.4529) * (392.0837 / 412.4028) * 1.14,
       height=3.75 * 0.93 * (87/62), scale=1.2)
rm(tmp.c3, tmp.opc, tmp.c3.opc)



#### 5B. OPC ----


DotPlot(object_1, features = cell_type_opc , group.by = "annotated_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
  labs(x = paste0("Features [CTX] in: ", gsub("^.+_","",sid), " (CPTAC-3 dataset)"))




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



