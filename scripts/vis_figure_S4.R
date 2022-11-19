#!/usr/bin/env R

# load data ----


source('scripts/load_results.out.R')


# Figure S4ab chr7, chr10, PTEN + EGFR ----



plt <- results.out |> 
  dplyr::mutate(mark.chr7 = chr == "chr7") |> 
  dplyr::mutate(mark.chr10 = chr == "chr10") |> 
  dplyr::mutate(show.label.chr7 = grepl("EGFR",hugo_symbol)) |> 
  dplyr::mutate(show.label.chr10 = grepl("PTEN|MGMT",hugo_symbol)) # TET1 really interesting




plt.expanded <- rbind(
  plt %>% 
    dplyr::mutate(status = case_when(show.label.chr7 == T ~ "chr7 label", mark.chr7 == T ~ "chr7", TRUE ~ "other")) %>% 
    dplyr::mutate(panel = "chr7 genes") %>%
    dplyr::select(hugo_symbol, statistic.gsam.cor.tpc, log2FoldChange.gsam.res, status, panel, chr)
  ,
  plt %>% 
    dplyr::mutate(status = case_when(show.label.chr10 == T ~ "chr10 label", mark.chr10 == T ~ "chr10", TRUE ~ "other")) %>% 
    dplyr::mutate(panel = "chr10 genes") %>%
    dplyr::select(hugo_symbol, statistic.gsam.cor.tpc, log2FoldChange.gsam.res, status, panel, chr)
) %>%
  dplyr::filter(!is.na(statistic.gsam.cor.tpc) & !is.na(log2FoldChange.gsam.res)) %>%
  dplyr::mutate(limited = abs(log2FoldChange.gsam.res) > 2.5) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(limited & log2FoldChange.gsam.res < 0, -2, log2FoldChange.gsam.res)) %>%
  dplyr::mutate(log2FoldChange.gsam.res = ifelse(limited & log2FoldChange.gsam.res > 0, 2, log2FoldChange.gsam.res)) %>%
  dplyr::mutate(panel = factor(panel, levels = c('chr7 genes', 'chr10 genes') )) |> 
  dplyr::mutate()



ggplot(plt.expanded, aes(x = log2FoldChange.gsam.res, y = statistic.gsam.cor.tpc, shape=limited,
                         fill = status, col = status, label=hugo_symbol)) +
  facet_grid(cols = vars(panel), scales = "free") +
  geom_point(data = subset(plt.expanded, status == "other"), size=1.8, col="#00000044") +
  geom_point(data = subset(plt.expanded, status %in% c("chr7", "chr10")), size=2.5, col="black") +
  geom_point(data = subset(plt.expanded, grepl("label", status)), size=2.75, fill="red", col="black") +
  ggrepel::geom_text_repel(data=subset(plt.expanded, status == "chr7 label"),  size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, show.legend  = FALSE) + 
  ggrepel::geom_text_repel(data=subset(plt.expanded, status == "chr10 label"),  size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, show.legend  = FALSE) + 
  scale_shape_manual(values = c('TRUE' = 23, 'FALSE' = 21)) +
  scale_fill_manual(values=c('chr7'='#6ba6e5AA','chr10'='#eab509AA','other'='white' , "chr10 label"='red',"chr7 label"='red'), guide="none") + 
  scale_color_manual(values=c('chr7'='#6ba6e5AA','chr10'='#eab509AA','other'='white' , "chr10 label"='red',"chr7 label"='red'), guide="none") + 
  theme_bw()  +
  scale_x_continuous(breaks = c(0), limits=c(-2.25,2.25)) +
  scale_y_continuous(breaks = c(0)) +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    
    panel.grid.major.x = element_line(colour = 'grey20', size=0.5),
    panel.grid.major.y = element_line(colour = 'grey20', size=0.5),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),

    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=8)) +
  labs(caption = "G-SAM: n=287 samples", 
       x="log2FC DESeq2 primary vs. recurrent (G-SAM; naive)",
       y="t-stat correlation tumor purity (G-SAM)")



ggsave("output/figures/2022_figure_S4ab_geiser_chr7_chr10.pdf", width=8.3 / 2,height=8.3/5, scale=2.25)
ggsave("output/figures/2022_figure_S4ab_geiser_chr7_chr10.svg", width=8.3 / 2,height=8.3/5, scale=2.25)








# Figure S4c [stats chr7&chr10] ----


plt <-
  results.out |> 
  dplyr::filter(chr %in% c('chrX','chrY') == F) |> 
  dplyr::mutate(pos = as.numeric(gsub(".+\\-([0-9]+)...$","\\1",gid))) |> 
  dplyr::mutate(chr = as.character(chr)) |> 
  dplyr::filter(!is.na(log2FoldChange.gsam.res)) |> 
  dplyr::group_by(chr) |> 
  dplyr::mutate(order = median(log2FoldChange.gsam.res)) |> 
  dplyr::ungroup()



## stats ----



tmp <- data.frame(chr = unique(plt$chr)) |> 
  dplyr::mutate(pval = lapply(chr , function(x, data){
    return(
      t.test(
        data |> dplyr::filter(chr != x) |> dplyr::pull(log2FoldChange.gsam.res),
        data |> dplyr::filter(chr == x) |> dplyr::pull(log2FoldChange.gsam.res),
        var.equal=T )$p.value
    )
  }, data=plt)) |>
  dplyr::mutate(padj.fdr = p.adjust(pval, method="fdr"))


tmp |> 
  dplyr::filter(padj.fdr < 0.01)



## add stats to plt ----


plt <- plt |> 
  dplyr::left_join(
    tmp |>
                     dplyr::select(chr, padj.fdr) |> 
                     dplyr::mutate(padj.fdr = format.pval(padj.fdr,digits=1))
                     , by=c('chr'='chr'),suffix=c('','')) |> 
  dplyr::mutate(padj.fdr = ifelse(duplicated(chr),NA,padj.fdr)) |> 
  dplyr::arrange(order) |> 
  # add newlines for or after label:
  dplyr::mutate(tmp = !is.na(padj.fdr)) |> 
  dplyr::group_by(tmp, ) |> 
  dplyr::mutate(prefix = ifelse(1:n() %% 2 == 0, "before", "after")) |> 
  dplyr::mutate(padj.fdr = case_when(
    tmp == T & prefix == "before" ~ paste0("\n",padj.fdr),
    tmp == T & prefix == "after"  ~ paste0(padj.fdr,"\n"),
    T ~ as.character(NA)
  )) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(tmp = NULL, prefix=NULL)




## plot ----



ggplot(plt, aes(x=reorder(chr, order),y=log2FoldChange.gsam.res, label=padj.fdr)) +
  ggbeeswarm::geom_quasirandom(cex=0.5,width=0.29, alpha=0.5) +
  geom_boxplot(outlier.shape = NA, col=alpha("white",1),fill=NA,lwd=1.1, width=0.70,coef=0) +
  geom_boxplot(outlier.shape = NA, col=alpha("red",1),fill=NA,alpha=0.75, width=0.70,coef=0) +
  geom_label(y = -0.72,col="blue",alpha=NA, fill=NA, label.size = NA, cex=3.2) +
  theme_bw() +
  labs(caption = "G-SAM: n=287 samples",
       x = NULL,
       y = "log2FC DESeq2 primary vs.\nrecurrent (G-SAM; naive)") +
  scale_x_discrete(expand=expansion(add = 0.85)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25), limits=c(-0.82,0.72)) + # leave space for labels
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(0.8)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    
    #panel.grid.major.x = element_line(colour = 'grey20', size=0.5),
    panel.grid.major.y = element_line(colour = 'blue', size=0.5),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


ggsave("output/figures/2022_figure_S4c_LFC_per_chr.pdf", width=8.3 / 2,height=8.3/6, scale=2)
ggsave("output/figures/2022_figure_S4c_LFC_per_chr.svg", width=8.3 / 2,height=8.3/6, scale=2)




# Figure S4def-ghi ----


plt <- results.out |> 
  
  #dplyr::filter(!is.na(log2FoldChange.gsam.tpc.res)  & !is.na(statistic.gsam.cor.tpc) ) |> 
  
  #dplyr::mutate(is.limited.gsam.tpc.res = as.character(abs(log2FoldChange.gsam.tpc.res) > 3)) |> 
  #dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 3, 3 , log2FoldChange.gsam.tpc.res)) |> 
  #dplyr::mutate(log2FoldChange.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res < -3, -3 , log2FoldChange.gsam.tpc.res)) |> 
  
  dplyr::mutate(direction.gsam.tpc.res = ifelse(log2FoldChange.gsam.tpc.res > 0 , "up", "down") ) |> 
  dplyr::mutate(`direction.glass-2022.tpc.res` = ifelse(`log2FoldChange.glass-2022.tpc.res` > 0 , "up", "down") ) |> 
  dplyr::mutate(significant = 
                  padj.gsam.tpc.res < 0.01 &
                  abs(log2FoldChange.gsam.tpc.res) > 0.5 & 
                  abs(`log2FoldChange.glass-2022.tpc.res`) > 0.5 & 
                  `direction.gsam.tpc.res` == `direction.glass-2022.tpc.res` 
  ) %>%
  dplyr::mutate(show.label.gains = hugo_symbol %in% c(
    "AKT1","AKT3","EGFR", "PDGFRA","MET", "PIK3C2B", "MDM2","MDM4", "CDK4",
    "CDK6","SOX2","FGFR3","MYCN","MYC","CCND1","CCND2","BMI1" # gains
  )) %>%
  dplyr::mutate(show.label.other = hugo_symbol %in% c( 
    "TP53", "PIK3CA","PIK3R1","NF1","SPTA1","GABRA6","ABCC6","CXCL12",
    "LTBP4", "TGFB1","PREX1","MSH6", "MSH2", "MLH1","VEGFA", "PTPN11",
    "STAT3","DCC", "MGMT", "PTCH1","VEGFA","GLI1",
    "NODAL",
    "PTPN11", # mut [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5627776/]
    "LTBP4", "MSH6", "PRDM2", "IGF1R"
  )) %>%
  dplyr::mutate(show.label.losses = hugo_symbol %in% c(
    "CDKN2A","CDKN2B","PTEN","RB1","NF1","QKI","CDKN2C","TP53","MTAP","ELAVL2","PARK2"
  ))


#PTPN11, LTBP4, MSH6, PRDM2, IGF1R



plt.expanded <- rbind(
  plt %>% 
    dplyr::rename(cor.tpc = statistic.gsam.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(show.label = show.label.gains) %>% 
    dplyr::mutate(marker.genes = "Commonly gained") %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.gsam.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(show.label = show.label.other) %>% 
    dplyr::mutate(marker.genes = "Other") %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = statistic.gsam.cor.tpc) %>%
    dplyr::rename(log2FoldChange = log2FoldChange.gsam.tpc.res) %>%
    dplyr::rename(show.label = show.label.losses) %>% 
    dplyr::mutate(marker.genes = "Commonly lost") %>%
    dplyr::mutate(dataset = "G-SAM") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = `statistic.t.glass-2022.cor.tpc`) %>%
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.tpc.res`) %>%
    dplyr::rename(show.label = show.label.gains) %>% 
    dplyr::mutate(marker.genes = "Commonly gained") %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = `statistic.t.glass-2022.cor.tpc`) %>%
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.tpc.res`) %>%
    dplyr::rename(show.label = show.label.other) %>% 
    dplyr::mutate(marker.genes = "Other") %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
  ,
  plt %>% 
    dplyr::rename(cor.tpc = `statistic.t.glass-2022.cor.tpc`) %>%
    dplyr::rename(log2FoldChange = `log2FoldChange.glass-2022.tpc.res`) %>%
    dplyr::rename(show.label = show.label.losses) %>% 
    dplyr::mutate(marker.genes = "Commonly lost") %>%
    dplyr::mutate(dataset = "GLASS") %>%
    dplyr::select(hugo_symbol, cor.tpc, log2FoldChange, show.label, dataset, significant, marker.genes)
) %>%
  dplyr::mutate(status = case_when(show.label == F & significant == F ~ "other",
                                   show.label == F & significant == T ~ "significant",
                                   show.label == T & significant == F ~ "show label",
                                   show.label == T & significant == T ~ "label & show significant")) %>%
  dplyr::filter(!is.na(cor.tpc) & !is.na(log2FoldChange)) %>%
  dplyr::mutate(marker.genes = factor(marker.genes, levels = c("Commonly gained", "Other", "Commonly lost"))) %>%
  dplyr::mutate(limited = abs(log2FoldChange) > 2.5) %>%
  dplyr::mutate(log2FoldChange = ifelse(limited & log2FoldChange < 0, -2.5, log2FoldChange)) %>%
  dplyr::mutate(log2FoldChange = ifelse(limited & log2FoldChange > 0, 2.5, log2FoldChange)) %>% 
  dplyr::mutate(limited = as.character(limited))



ggplot(plt.expanded, aes(x = log2FoldChange, y = cor.tpc, shape=limited,
                         fill = status, col = status, label=hugo_symbol)) +
  facet_grid(cols = vars(marker.genes), rows = vars(dataset), scales = "free") +
  geom_point(data = subset(plt.expanded, status == "other"), size=1.8, col="#00000044") +
  geom_point(data = subset(plt.expanded, status == "significant"), size=2, col="black") +
  geom_vline(xintercept = c(-0.5,0.5), lty=1, lwd=2, col="#FFFFFF88") +
  geom_vline(xintercept = c(-0.5,0.5), lty=2, col="red") +
  ggrepel::geom_text_repel(data=subset(plt.expanded, grepl("label",status) & log2FoldChange > 0), col="blue", size=2.5 ,nudge_x = 3.1, direction = "y", hjust = "left", segment.size=0.35, show.legend  = FALSE) + 
  ggrepel::geom_text_repel(data=subset(plt.expanded, grepl("label",status) & log2FoldChange < 0), col="blue", size=2.5 , nudge_x = -3.1, direction = "y", hjust = "right", segment.size=0.35, show.legend  = FALSE) + 
  geom_point(data = subset(plt.expanded, status == "show label"), size=2.75, col="black") +
  geom_point(data = subset(plt.expanded, status == "label & show significant"), size=2.75, col="black") +
  scale_shape_manual(values = c('TRUE' = 23, 'FALSE' = 21)) +
  scale_fill_manual(values = c("other" = "white", "significant" = "#eab509DD", "show label" = "#6ba6e5DD", "label & show significant" = "red")) + 
  scale_color_manual(values = c("other" = "white", "significant" = "#eab509DD", "show label" = "#6ba6e5DD", "label & show significant" = "red"),guide="none") + 
  theme_bw() +
  xlim(-3, 3) +
  #scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) +
  #scale_y_continuous(breaks = c(10,7.5,5,2.5,0,-2.5,-5,-7.5,-10)) +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    
    #panel.grid.major.x = element_line(colour = 'grey20', size=0.5),
    #panel.grid.major.y = element_line(colour = 'grey20', size=0.5),
    
    #panel.grid.minor.x = element_blank(),
    #panel.grid.minor.y = element_blank(),
    
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
    ) +
  labs(x="log2FC primary vs. recurrent (purity corrected)",
       y="t-stat correlation tumor purity",
       caption = "G-SAM: n=287  -  GLASS: n=216  samples")


ggsave("output/figures/2022_figure_S4def-ghi_geiser_gains_losses_other.pdf",  width=8.3 / 2,height=8.3/3, scale=2)
ggsave("output/figures/2022_figure_S4def-ghi_geiser_gains_losses_other.svg",  width=8.3 / 2,height=8.3/3, scale=2)





