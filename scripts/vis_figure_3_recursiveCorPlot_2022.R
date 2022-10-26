#!/usr/bin/env R 

# load libs ----


library(recursiveCorPlot)


# load data ----


source('scripts/load_results.out.R')
source('scripts/load_G-SAM_expression_data.R')


# prepare data ----


## genes ----


metadata.genes <- results.out |> 
    dplyr::select(-contains("glass.tpc")) |>  # remove 2021/small batch
    
    dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) |> 
    dplyr::filter(!is.na(`log2FoldChange.glass-2022.tpc.res`)) |> 
    dplyr::filter(!is.na(padj.gsam.tpc.res)) |> 
    dplyr::filter(!is.na(`padj.glass-2022.tpc.res`)) |> 
    
    dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down")) |> 
    dplyr::mutate(`direction.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` > 0 , "up", "down")) |> 
    
    dplyr::mutate(significant.2022 = 
                    padj.gsam.tpc.res < 0.01 &
                    abs(log2FoldChange.gsam.tpc.res) > 0.5 &
                    abs(`log2FoldChange.glass-2022.tpc.res`) > 0.5 & 
                    direction.gsam.tpc.res == `direction.glass-2022.tpc.res`)


# stats - DE Genes G-SAM:
metadata.genes |> 
  dplyr::filter(
    padj.gsam.tpc.res < 0.01 &
      abs(log2FoldChange.gsam.tpc.res) > 0.5
  ) |>  dim()



metadata.genes <- metadata.genes |> 
  dplyr::filter(significant.2022 == T)



k.genes <- nrow(metadata.genes)
print(paste0("Genes part of recursive clustering: ", k.genes))





# plot ----


plt.labels <- metadata.genes %>%
  dplyr::mutate(`McKenzie_celltype_top_human_specificity` = as.factor(.data$`McKenzie_celltype_top_human_specificity`)) |> 
  dplyr::rename(`Extracellular Matrix` = EM.struct.constituent) %>%
  dplyr::select(c('gid', 'hugo_symbol' , `McKenzie_celltype_top_human_specificity`, 
                  `Extracellular Matrix`, `EM.GO.0031012`,
                  `direction.glass-2022.tpc.res`
  )) |> 
  
  dplyr::mutate(decoy = T) |> 
  tidyr::pivot_wider(names_from = `direction.glass-2022.tpc.res`, values_from = decoy, values_fill = FALSE) |> 
  dplyr::mutate(`NA` = NULL) |> 
  as.data.frame() |> 

  dplyr::mutate(decoy = T) |> 
  tidyr::pivot_wider(names_from = `McKenzie_celltype_top_human_specificity`, values_from = decoy, values_fill = FALSE) |> 
  dplyr::mutate(`NA` = NULL) |> 
  as.data.frame()

if(!'TAM' %in% colnames(plt.labels)) {
  plt.labels$TAM <- FALSE
}


### export Table S2:

tmp.order <- readRDS('tmp/h.2022.Rds')
tmp.order <- data.frame(hugo_symbol = tmp.order$labels, order = tmp.order$order)

exp.ts2 <- plt.labels |> 
  dplyr::rename(`increases over time`=up) |> 
  dplyr::rename(`decreases over time`=down) |> 
  
  dplyr::mutate(oligodendrocyte = NULL) |> 
  dplyr::mutate(endothelial = NULL) |> 
  dplyr::mutate(neuron = NULL) |> 
  dplyr::mutate(astrocyte = NULL) |> 
  dplyr::mutate(TAM = NULL) |> 
  
  dplyr::left_join(tmp.order, by=c('hugo_symbol'='hugo_symbol'), suffix=c('','')) |> 
  
  dplyr::left_join(
    results.out |> 
      dplyr::select(gid, ensembl_id, cluster.2022 ,`McKenzie_celltype_top_human_specificity`),
    by=c('gid'='gid'),suffix=c('','')
  ) |> 
  dplyr::arrange(desc(order(order))) |> 
  dplyr::mutate(order=  NULL) |> 
  dplyr::mutate(gid =  NULL)

stopifnot(nrow(tmp.order) == 484)



openxlsx::createWorkbook(
  creator = "Dr. Youri Hoogstrate",
  title = "Gene cluster information G-SAM study",
  subject = "Gene clusters",
  category = "G-SAM study"
) -> wb
openxlsx::addWorksheet(wb, "Sheet1 - Gene clusters G-SAM")
openxlsx::writeDataTable(wb, sheet = "Sheet1 - Gene clusters G-SAM", x = exp.ts2)

openxlsx::saveWorkbook(wb, file = "output/tables/tab_S2_DGE_gene_clusters.xlsx", overwrite = T)




# :: ::: 



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

recursiveCorPlot::recursiveCorPlot(plt.data, plt.labels |> 
                   tibble::column_to_rownames('hugo_symbol') |> 
                   dplyr::mutate(gid = NULL),
                 font_scale = 3.6,
                 legend_scale = 3,
                 method = "ward.D2",
                 caption = paste0("Genes: n=",       k.genes, ", ",
                                  "G-SAM samples: n=",ncol(gsam.gene.expression.all.vst) ," ",
                                  "(", ncol(gsam.gene.expression.all.vst)," >= 15%, 0 < 15%)"),
                 return_h_object = FALSE
                 )





ggsave("output/figures/2022_figure_3_recursiveCorPlot.pdf", height=30, width=30) # then average, then ward.D1?
ggsave("output/figures/2022_figure_3_recursiveCorPlot.svg", height=30, width=30) # then average, then ward.D1?




h <- recursiveCorPlot(plt.data, plt.labels |> 
                   tibble::column_to_rownames('hugo_symbol') |> 
                   dplyr::mutate(gid = NULL),
                 font_scale = 3.6,
                 legend_scale = 3,
                 method = "ward.D2", 
                 return_h_object = TRUE)

#saveRDS(h, "cache/h.2022.Rds")



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





