#!/usr/bin/env R 

# load libs ----


# load data ----

source('scripts/load_results.out.R')


# prepare data ----

## genes ----

metadata.genes <- results.out |> 
    dplyr::select(-contains("glass.tpc")) |>  # remove 2021/small batch
    
    dplyr::filter( !is.na(log2FoldChange.gsam.tpc.res)   ) %>%
    dplyr::filter( !is.na(`log2FoldChange.glass-2022.tpc.res`)  ) %>%
    dplyr::filter( !is.na(padj.gsam.tpc.res)   ) %>%
    dplyr::filter( !is.na(`padj.glass-2022.tpc.res`)  ) %>%
    
    dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down")) %>%
    dplyr::mutate(`direction.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` > 0 , "up", "down")) %>%
    
    dplyr::mutate(significant.2022 = 
                    padj.gsam.tpc.res < 0.01 &
                    abs(log2FoldChange.gsam.tpc.res) > 0.5 &
                    abs(`log2FoldChange.glass-2022.tpc.res`) > 0.5 & 
                    direction.gsam.tpc.res == `direction.glass-2022.tpc.res`) |> 
    
    dplyr::filter(significant.2022 == T)

print(paste0("Genes part of recursive clustering: ", nrow(metadata.genes)))




# plot ----


plt.labels <- metadata.genes %>%
  
  dplyr::mutate(direction_up = ifelse(log2FoldChange.gsam.tpc.res > 0 , T, F) ) %>% 
  dplyr::mutate(direction_down = ifelse(log2FoldChange.gsam.tpc.res < 0 , T, F) ) %>% 
  
  dplyr::mutate(TCGA.CL = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-CL" , T, F) ) %>% 
  dplyr::mutate(TCGA.MES = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-MES" , T, F) ) %>% 
  dplyr::mutate(TCGA.PN = ifelse(!is.na(TCGA.subtype.marker) & TCGA.subtype.marker == "TCGA-PN" , T, F) ) %>% 
  
  dplyr::mutate(Patel.Cell.cycle = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Cell cycle" , T, F) ) %>% 
  dplyr::mutate(Patel.Complete.Immune.response = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Complete/Immune response" , T, F) ) %>% 
  dplyr::mutate(Patel.Hypoxia = ifelse(!is.na(patel.scRNAseq.cluster) & patel.scRNAseq.cluster == "Hypoxia" , T, F) ) %>% 
  
  dplyr::rename(Neftel.AC = neftel.meta.module.AC) %>% 
  dplyr::rename(Neftel.NPC1 = neftel.meta.module.NPC1) %>% 
  dplyr::rename(Neftel.NPC2 = neftel.meta.module.NPC2) %>% 
  dplyr::rename(Neftel.OPC = neftel.meta.module.OPC) %>% 
  dplyr::rename(Neftel.MES1 = neftel.meta.module.MES1) %>% 
  dplyr::rename(Neftel.MES2 = neftel.meta.module.MES2) %>% 
  
  dplyr::rename(`Extracellular Matrix` = EM.struct.constituent) %>%
  
  dplyr::select(c('gid', 'hugo_symbol' , 'McKenzie_celltype_top_human_specificity', 'direction_up', 'direction_down',
                  TCGA.CL, TCGA.MES, TCGA.PN ,
                  Patel.Cell.cycle, Patel.Complete.Immune.response, Patel.Hypoxia, 
                  Neftel.AC, Neftel.OPC, Neftel.MES1, Neftel.MES2, Neftel.NPC1, Neftel.NPC2,
                  `Extracellular Matrix`
                  
                  ,significant.2022,
                  
                  direction.gsam.tpc.res ,log2FoldChange.gsam.tpc.res,  
                  `direction.glass-2022.tpc.res` ,`log2FoldChange.glass-2022.tpc.res`
  )) |> 
  reshape2::dcast(hugo_symbol +
                    direction_up + direction_down + 
                    TCGA.CL + TCGA.MES + TCGA.PN +
                    Neftel.AC + Neftel.NPC1 + Neftel.NPC2 + Neftel.OPC + Neftel.MES1 + Neftel.MES2 + 
                    Patel.Cell.cycle + Patel.Complete.Immune.response + Patel.Hypoxia +
                    `Extracellular Matrix` +
                    gid ~ `McKenzie_celltype_top_human_specificity`, fill = NA,fun.aggregate = as.logical) |> 
  dplyr::mutate('NA'=NULL) %>% 
  
  dplyr::mutate(direction_up = ifelse(direction_up , T, F) ) %>% # F functions as T, NA as F, for dcast2
  dplyr::mutate(direction_down = ifelse(direction_down , T, F) ) %>% # F functions as T, NA as F, for dcast2
  
  dplyr::mutate(TCGA.CL = ifelse(TCGA.CL , F, NA) ) %>%
  dplyr::mutate(TCGA.MES = ifelse(TCGA.MES , F, NA) ) %>%
  dplyr::mutate(TCGA.PN = ifelse(TCGA.PN , F, NA) ) %>%
  
  dplyr::mutate(Neftel.AC = ifelse(Neftel.AC , F, NA) ) %>%
  dplyr::mutate(Neftel.NPC1 = ifelse(Neftel.NPC1 , F, NA) ) %>%
  dplyr::mutate(Neftel.NPC2 = ifelse(Neftel.NPC2 , F, NA) ) %>%
  dplyr::mutate(Neftel.OPC = ifelse(Neftel.OPC , F, NA) ) %>%
  dplyr::mutate(Neftel.MES1 = ifelse(Neftel.MES1 , F, NA) ) %>%
  dplyr::mutate(Neftel.MES2 = ifelse(Neftel.MES2 , F, NA) ) %>%
  
  dplyr::mutate(`Extracellular Matrix` = ifelse(`Extracellular Matrix` , T, F) ) %>%
  
  dplyr::mutate(Patel.Cell.cycle = ifelse(Patel.Cell.cycle , F, NA) ) %>%
  dplyr::mutate(Patel.Complete.Immune.response = ifelse(Patel.Complete.Immune.response , F, NA) ) %>%
  dplyr::mutate(Patel.Hypoxia = ifelse(Patel.Hypoxia , F, NA) ) %>%
  
  dplyr::mutate(TCGA.CL = NULL , TCGA.MES = NULL , TCGA.PN = NULL ,
                Patel.Cell.cycle = NULL , Patel.Complete.Immune.response = NULL , Patel.Hypoxia = NULL)



gsam.gene.expression.all.vst <- gsam.gene.expression.all %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() 



sel = gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F & pat.with.IDH == F) |> 
  dplyr::mutate(sid = paste0("GSAM-",sid)) |> 
  dplyr::filter(sid %in% colnames(gsam.gene.expression.all.vst)) |>  # low res
  dplyr::pull(sid)


gsam.gene.expression.all.vst |>
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(sel) 



# substitute gene names
plt.data <- plt.labels %>%
  dplyr::select('gid','hugo_symbol') %>%
  dplyr::left_join(gsam.gene.expression.all.vst %>%
                     as.data.frame %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid'='gid')) %>% 
  tibble::column_to_rownames('hugo_symbol') %>%
  dplyr::mutate(gid = NULL)


# check if data and 
stopifnot(rownames(plt.data) == plt.labels$hugo_symbol)
sum(duplicated(plt.labels$hugo_symbol)) == 0
sum(duplicated(rownames(plt.data))) == 0

recursiveCorPlot(plt.data, plt.labels |> 
                   dplyr::mutate(EM.GO.0031012 = hugo_symbol %in% go.0031012) |> 
                   tibble::column_to_rownames('hugo_symbol') |> 
                   dplyr::mutate(gid = NULL) |> 
                   dplyr::mutate(astrocyte = ifelse(is.na( astrocyte ), F, astrocyte)) |> 
                   dplyr::mutate(endothelial = ifelse(is.na( endothelial ), F, endothelial)) |> 
                   dplyr::mutate(neuron = ifelse(is.na( neuron ), F, neuron)) |> 
                   dplyr::mutate(oligodendrocyte = ifelse(is.na( oligodendrocyte ), F, oligodendrocyte)) |> 
                   dplyr::mutate(Neftel.AC = NULL,Neftel.NPC1 = NULL,Neftel.NPC2 = NULL,Neftel.OPC = NULL,Neftel.MES1 = NULL,Neftel.MES2 = NULL )
                 , 3.6, 3)



ggsave("output/figures/recursiveCorPlot.2022.svg", height=30, width=30) # then average, then ward.D1?


hclust <- recursiveCorPlot(plt.data |> 
                             tibble::rownames_to_column('hugo_symbol') |> 
                             dplyr::left_join(plt.labels |> dplyr::select(hugo_symbol, gid),by=c('hugo_symbol'='hugo_symbol')) |> 
                             tibble::column_to_rownames('gid') |> 
                             dplyr::mutate(hugo_symbol=NULL)
                             , plt.labels |> 
                   dplyr::mutate(EM.GO.0031012 = hugo_symbol %in% go.0031012) |> 
                   tibble::column_to_rownames('gid') |> 
                   dplyr::mutate(hugo_symbol = NULL) |> 
                   dplyr::mutate(astrocyte = ifelse(is.na( astrocyte ), F, astrocyte)) |> 
                   dplyr::mutate(endothelial = ifelse(is.na( endothelial ), F, endothelial)) |> 
                   dplyr::mutate(neuron = ifelse(is.na( neuron ), F, neuron)) |> 
                   dplyr::mutate(oligodendrocyte = ifelse(is.na( oligodendrocyte ), F, oligodendrocyte)) |> 
                   dplyr::mutate(Neftel.AC = NULL,Neftel.NPC1 = NULL,Neftel.NPC2 = NULL,Neftel.OPC = NULL,Neftel.MES1 = NULL,Neftel.MES2 = NULL )
                 , 3.6, 3,"ward.D2", TRUE)

#saveRDS(hclust, "cache/h.2022.Rds")



#h <- cor_cor_plot(plt, labels, 3,6, method="ward.D2")
#ggsave("output/figures/cor_ward.D2x.pdf", height=30, width=30) # then average, then ward.D1?
#ggsave("output/figures/cor_ward.D2x.svg", height=30, width=30) # then average, then ward.D1?


#recursiveCorPlot(G.SAM.corrected.DE.genes.VST, G.SAM.corrected.DE.labels, 3 , 3)
# 
# sum(go.0031012 %in% de)
# go.0031012[go.0031012 %in% de & (go.0031012 %in% CC.2022 == F)]
# go.0031012[go.0031012 %in% de & (go.0031012 %in% CC.2022 == T)]
# 
# de = plt.labels$hugo_symbol



### COL1A2 ----


plt <- glass.gbm.rnaseq.expression.all.samples.vst |> 
  tibble::rownames_to_column('ensembl_id') |> 
  dplyr::filter(ensembl_id == "ENSG00000164692") |> 
  tibble::column_to_rownames('ensembl_id') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('aliquot_barcode') |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples, by=c('aliquot_barcode'='aliquot_barcode')
  )


ggplot(plt, aes(x= ifelse(is.primary,"prim","recur"), y=ENSG00000164692, group=case_barcode, label = sid.label, col=aliquot_batch_synapse)) +
  facet_grid(cols = vars(aliquot_batch_synapse), scales = "free")  + 
  geom_line() +
  geom_point() +
  ggrepel::geom_text_repel()



## obtain h-clust ----





