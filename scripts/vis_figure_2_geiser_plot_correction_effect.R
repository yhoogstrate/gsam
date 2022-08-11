#!/usr/bin/env R



## Figure 3 B,C,D,E re-styled double ----


# show CD14, CD74, CD163,
# CD33, CD84, FCER1G, GPR34, TMEM119



plt <- results.out %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) |>  # TPC correlation G-SAM
  dplyr::filter(!is.na(`statistic.t.glass-2022.cor.tpc`)) |> # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) |> 
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) |> 
  dplyr::filter(!is.na(`log2FoldChange.glass-2022.res`)) |> 
  dplyr::filter(!is.na(`log2FoldChange.glass-2022.tpc.res`))  |> 
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > 2, 2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -2, -2 , log2FoldChange.gsam.res)) %>%
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 2, 2 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -2, -2 , log2FoldChange.gsam.tpc.res)) %>%
  dplyr::mutate(is.limited.glass.res = as.character(abs(`log2FoldChange.glass-2022.res`) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(`log2FoldChange.glass-2022.res` = ifelse(`log2FoldChange.glass-2022.res` > 2, 2 , `log2FoldChange.glass-2022.res`)) %>%
  dplyr::mutate(`log2FoldChange.glass-2022.res` = ifelse(`log2FoldChange.glass-2022.res` < -2, -2 , `log2FoldChange.glass-2022.res`)) %>%
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(`log2FoldChange.glass-2022.tpc.res`) > 2)) %>% # change pch to something that is limited
  dplyr::mutate(`log2FoldChange.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` > 2, 2 , `log2FoldChange.glass-2022.tpc.res`)) %>%
  dplyr::mutate(`log2FoldChange.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` < -2, -2 , `log2FoldChange.glass-2022.tpc.res`)) %>% 
  dplyr::mutate(is.mg.or.tam = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == "microglia/TAM")



# pivot longer is complicated because both the log2fc & cor.stat need to be pivotted

plt.expanded <- rbind(
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.res) %>%
    dplyr::rename(padj = padj.gsam.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(padj = padj.gsam.tpc.res) %>%
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) %>%
    dplyr::rename(is.limited = is.limited.gsam.tpc.res) %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::mutate(DESeq2.tcp.corrected = T) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.res`) %>%
    dplyr::rename(padj = `padj.glass-2022.res`) %>%
    dplyr::rename(cor.stat = `statistic.t.glass-2022.cor.tpc`) %>%
    dplyr::rename(is.limited = is.limited.glass.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = F) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt %>%
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.tpc.res`) %>%
    dplyr::rename(padj = `padj.glass-2022.tpc.res`) %>%
    dplyr::rename(cor.stat = `statistic.t.glass-2022.cor.tpc`) %>%
    dplyr::rename(is.limited = is.limited.glass.tpc.res) %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::mutate(DESeq2.tcp.corrected = T) %>%
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
) %>%
  dplyr::mutate(DESeq2.tcp.corrected = ifelse(DESeq2.tcp.corrected == T, "Tumor cell-% corrected" , "Tumor cell-% uncorrected")) %>%
  dplyr::mutate(DESeq2.tcp.corrected = factor(DESeq2.tcp.corrected, levels=
                                                c("Tumor cell-% uncorrected" , "Tumor cell-% corrected"))) %>%
  dplyr::mutate(is.limited = is.limited == 'TRUE') %>% 
  dplyr::mutate(col = case_when(
    is.mg.or.tam == F ~ "no-mg-or-tam",
    is.mg.or.tam == T ~ dataset
  )) %>%
  dplyr::mutate(show.label = hugo_symbol %in% c("CD33", "CD84", "FCER1G", "GPR34", "CD14","CD74", "CD163")) %>% 
  dplyr::rename(`Purity corrected` = DESeq2.tcp.corrected) 




ggplot(plt.expanded, aes(x = log2FoldChange ,
                         y =  cor.stat,
                         shape = is.limited ,
                         col =  is.limited ,
                         label = hugo_symbol,
                         #size = is.limited  ,
                         fill = col  ) ) +
  facet_grid(rows = vars(dataset), cols=vars(`Purity corrected`), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.mg.or.tam` == F)) +
  geom_point(size=2.2, data = subset(plt.expanded, `is.mg.or.tam` == T)) +
  geom_smooth(data = subset(plt.expanded, abs(cor.stat) > 7.5 &  is.limited == FALSE),
              aes(group=1),
              method="lm",
              se = FALSE,  formula=y ~ x, orientation="y", col="#ff2929", show.legend=F ) +
  ggrepel::geom_text_repel(data=subset(plt.expanded, show.label  & log2FoldChange > 0), size=2.5 , 
                           segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1, 
                           direction = "y", hjust = "left", col="black", nudge_y = -4.8 ) +
  ggrepel::geom_text_repel(data=subset(plt.expanded, show.label  & log2FoldChange < 0), size=2.5 ,
                           segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1, 
                           direction = "y", hjust = "right", col="black", nudge_y = -4.8 ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e69356', 'GLASS' = '#69a4d5', 'no-mg-or-tam'='white') ) + 
  youri_gg_theme + 
  labs(x = "log2FC R1 vs. R2",
       y="Correlation t-statistic with tumour percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.5,2.5)



##ggsave('output/figures/geiser_plot_mckenzy_immune_tam_double.pdf', height=5.75 * 1.1 * 1.1,width=12 * 1.1)



