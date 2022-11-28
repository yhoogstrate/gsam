#!/usr/bin/env R


# load data -----


if(!exists('results.out')) {
  source('scripts/load_results.out.R')
}


## F] Figure M2C, M2D, M2E & M2F ----

# show CD14, CD74, CD163,
# CD33, CD84, FCER1G, GPR34, TMEM119


l.gsam <- 2.0
l.glass <- 2.5

plt <- results.out |>
  dplyr::filter(!is.na(statistic.gsam.cor.tpc)) |>  # TPC correlation G-SAM
  dplyr::filter(!is.na(`statistic.t.glass-2022.cor.tpc`)) |> # TPC correlation GLASS
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) |> 
  dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)) |> 
  dplyr::filter(!is.na(`log2FoldChange.glass-2022.res`)) |> 
  dplyr::filter(!is.na(`log2FoldChange.glass-2022.tpc.res`)) |> 
  
  dplyr::mutate(is.limited.gsam.res = as.character(abs(log2FoldChange.gsam.res) > l.gsam)) |> # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res > l.gsam, l.gsam, log2FoldChange.gsam.res)) |>
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(log2FoldChange.gsam.res < -l.gsam, -l.gsam , log2FoldChange.gsam.res)) |>
  dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > l.gsam)) |> # change pch to something that is limited
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > l.gsam, l.gsam , log2FoldChange.gsam.tpc.res)) |>
  dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -l.gsam, -l.gsam , log2FoldChange.gsam.tpc.res)) |>
  
  dplyr::mutate(is.limited.glass.res = as.character(abs(`log2FoldChange.glass-2022.res`) > l.glass)) |> # change pch to something that is limited
  dplyr::mutate(`log2FoldChange.glass-2022.res` = ifelse(`log2FoldChange.glass-2022.res` > l.glass, l.glass , `log2FoldChange.glass-2022.res`)) |>
  dplyr::mutate(`log2FoldChange.glass-2022.res` = ifelse(`log2FoldChange.glass-2022.res` < -l.glass, -l.glass , `log2FoldChange.glass-2022.res`)) |>
  dplyr::mutate(is.limited.glass.tpc.res = as.character(abs(`log2FoldChange.glass-2022.tpc.res`) > l.glass)) |> # change pch to something that is limited
  dplyr::mutate(`log2FoldChange.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` > l.glass, l.glass , `log2FoldChange.glass-2022.tpc.res`)) |>
  dplyr::mutate(`log2FoldChange.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` < -l.glass, -l.glass , `log2FoldChange.glass-2022.tpc.res`)) |> 
  dplyr::mutate(is.mg.or.tam = !is.na(McKenzie_celltype_top_human_specificity) & McKenzie_celltype_top_human_specificity == "microglia/TAM")



# pivot longer is complicated because both the log2fc & cor.stat need to be pivotted

plt.expanded.gsam <- rbind(
  plt |>
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.res) |>
    dplyr::rename(padj = padj.gsam.res) |>
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) |>
    dplyr::rename(is.limited = is.limited.gsam.res) |>
    dplyr::mutate(dataset = "G-SAM") |>
    dplyr::mutate(DESeq2.tcp.corrected = F) |>
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt |>
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) |>
    dplyr::rename(padj = padj.gsam.tpc.res) |>
    dplyr::rename(cor.stat = statistic.gsam.cor.tpc) |>
    dplyr::rename(is.limited = is.limited.gsam.tpc.res) |>
    dplyr::mutate(dataset = "G-SAM") |>
    dplyr::mutate(DESeq2.tcp.corrected = T) |>
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
) |>
  dplyr::mutate(DESeq2.tcp.corrected = ifelse(DESeq2.tcp.corrected == T, "Tumor cell-% corrected" , "Tumor cell-% uncorrected")) |>
  dplyr::mutate(DESeq2.tcp.corrected = factor(DESeq2.tcp.corrected, levels=
                                                c("Tumor cell-% uncorrected" , "Tumor cell-% corrected"))) |>
  dplyr::mutate(is.limited = is.limited == 'TRUE') |> 
  dplyr::mutate(col = case_when(
    is.mg.or.tam == F ~ "no-mg-or-tam",
    is.mg.or.tam == T ~ dataset
  )) |>
  dplyr::mutate(show.label = hugo_symbol %in% c("CD33", "CD84", "FCER1G", "GPR34", "CD14","CD74", "CD163"
                                                #, "TMEM144","RBFOX3"
  )) |> 
  dplyr::rename(`Purity corrected` = DESeq2.tcp.corrected) 



plt.expanded.glass <- rbind(
  plt |>
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.res`) |>
    dplyr::rename(padj = `padj.glass-2022.res`) |>
    dplyr::rename(cor.stat = `statistic.t.glass-2022.cor.tpc`) |>
    dplyr::rename(is.limited = is.limited.glass.res) |>
    dplyr::mutate(dataset = "GLASS") |>
    dplyr::mutate(DESeq2.tcp.corrected = F) |>
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
  ,
  plt |>
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.tpc.res`) |>
    dplyr::rename(padj = `padj.glass-2022.tpc.res`) |>
    dplyr::rename(cor.stat = `statistic.t.glass-2022.cor.tpc`) |>
    dplyr::rename(is.limited = is.limited.glass.tpc.res) |>
    dplyr::mutate(dataset = "GLASS") |>
    dplyr::mutate(DESeq2.tcp.corrected = T) |>
    dplyr::select(hugo_symbol, log2FoldChange, cor.stat, is.limited, dataset, DESeq2.tcp.corrected, padj, is.mg.or.tam)
) |>
  dplyr::mutate(DESeq2.tcp.corrected = ifelse(DESeq2.tcp.corrected == T, "Tumor cell-% corrected" , "Tumor cell-% uncorrected")) |>
  dplyr::mutate(DESeq2.tcp.corrected = factor(DESeq2.tcp.corrected, levels=
                                                c("Tumor cell-% uncorrected" , "Tumor cell-% corrected"))) |>
  dplyr::mutate(is.limited = is.limited == 'TRUE') |> 
  dplyr::mutate(col = case_when(
    is.mg.or.tam == F ~ "no-mg-or-tam",
    is.mg.or.tam == T ~ dataset
  )) |>
  dplyr::mutate(show.label = hugo_symbol %in% c("CD33", "CD84", "FCER1G", "GPR34", "CD14","CD74", "CD163"
                                                #, "TMEM144","RBFOX3"
  )) |> 
  dplyr::rename(`Purity corrected` = DESeq2.tcp.corrected) 




p1 <- ggplot(plt.expanded.gsam, aes(x = log2FoldChange ,
                                    y =  cor.stat,
                                    shape = is.limited ,
                                    col =  is.limited ,
                                    label = hugo_symbol,
                                    #size = is.limited  ,
                                    fill = col  ) ) +
  facet_grid(rows = vars(dataset), cols=vars(`Purity corrected`), scales = "free") +
  #facet_grid(rows = vars(dataset), cols=vars(`Purity corrected`), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded.gsam, `is.mg.or.tam` == F)) +
  geom_point(size=2.2, data = subset(plt.expanded.gsam, `is.mg.or.tam` == T)) +
  geom_smooth(data = subset(plt.expanded.gsam, abs(cor.stat) > 6.5 &  is.limited == FALSE),
              aes(group=1),
              method="lm",
              se = FALSE,  formula=y ~ x, orientation="y",col="#ff2929", show.legend=F ) +
  ggrepel::geom_text_repel(data=subset(plt.expanded.gsam, show.label  & log2FoldChange > 0), size=2.5 ,
                           segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1,
                           direction = "y", hjust = "left", col="black", nudge_y = -4.8 ) +
  ggrepel::geom_text_repel(data=subset(plt.expanded.gsam, show.label  & log2FoldChange < 0), size=2.5 ,
                           segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1,
                           direction = "y", hjust = "right", col="black", nudge_y = -4.8 ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e69356', 'GLASS' = '#69a4d5', 'no-mg-or-tam'='white') ) + 
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x = "log2FC primary vs. recurrence",
       y="Correlation t-statistic with tumor percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis") +
  xlim(-2.25,2.25) +
  ggpubr::stat_cor(
    aes(col=NULL,fill=NULL,shape=NULL,
        label = ..r.label..),
    method="spearman") +
  labs(caption="287 - 216")


p2 <- ggplot(plt.expanded.glass, aes(x = log2FoldChange ,
                                     y =  cor.stat,
                                     shape = is.limited ,
                                     col =  is.limited ,
                                     label = hugo_symbol,
                                     #size = is.limited  ,
                                     fill = col  ) ) +
  facet_grid(rows = vars(dataset), cols=vars(`Purity corrected`), scales = "free") +
  #facet_grid(rows = vars(dataset), cols=vars(`Purity corrected`), scales = "free") +
  geom_hline(yintercept = 0, col="gray60",lwd=0.5) +
  geom_vline(xintercept = 0, col="gray60",lwd=0.5) +
  geom_point(size=2.2, data = subset(plt.expanded.glass, `is.mg.or.tam` == F)) +
  geom_point(size=2.2, data = subset(plt.expanded.glass, `is.mg.or.tam` == T)) +
  geom_smooth(data = subset(plt.expanded.glass, abs(cor.stat) > 6.5 &  is.limited == FALSE),
              aes(group=1),
              method="lm",
              se = FALSE,  formula=y ~ x, orientation="y",col="#ff2929", show.legend=F ) +
  ggrepel::geom_text_repel(data=subset(plt.expanded.glass, show.label  & log2FoldChange > 0), size=2.5 ,
                           segment.size = 0.25, segment.linetype = 1, nudge_x = 3.1,
                           direction = "y", hjust = "left", col="black", nudge_y = -4.8 ) +
  ggrepel::geom_text_repel(data=subset(plt.expanded.glass, show.label  & log2FoldChange < 0), size=2.5 ,
                           segment.size = 0.25, segment.linetype = 1, nudge_x = -3.1,
                           direction = "y", hjust = "right", col="black", nudge_y = -4.8 ) +
  scale_shape_manual(values = c('TRUE'= 23, 'FALSE' = 21)   )  +
  scale_color_manual(values = c('FALSE' = '#00000040', 'TRUE' = '#000000aa')) +
  scale_fill_manual(values = dataset_colors <- c('G-SAM' = '#e69356', 'GLASS' = '#69a4d5', 'no-mg-or-tam'='white') ) + 
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x = "log2FC primary vs. recurrence",
       y="Correlation t-statistic with tumor percentage",
       shape = "Truncated at x-axis",
       size="Truncated at x-axis",
       col="Truncated at x-axis"
       ) +
  xlim(-2.75,2.75) +
  ggpubr::stat_cor(
    aes(col=NULL,fill=NULL,shape=NULL,
        label = ..r.label.. # nice hack to only obtain the cor [https://github.com/kassambara/ggpubr/issues/188#issuecomment-501976246]
        ),
    method="spearman") +
  labs(caption="G-SAM: n=287  -  GLASS: n=216  samples")



p1 / p2


ggsave('output/figures/2022_Figure_M2C_M2D_M2E_M2F.pdf', width=8.3 / 2,height=8.3/2.25, scale=2.1)


rm(l.gsam, l.glass, plt, plt.expanded.gsam, plt.expanded.glass, p1, p2)



