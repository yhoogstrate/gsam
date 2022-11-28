

stat.gsam <- results.out |>
  dplyr::filter(TCGA.subtype.marker == "TCGA-MES") |>
  dplyr::select(gid,statistic.gsam.cor.tpc,log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res ) 

median(stat.gsam$log2FoldChange.gsam.res)
median(stat.gsam$log2FoldChange.gsam.tpc.res )

plot(sort(stat.gsam$log2FoldChange.gsam.tpc.res - stat.gsam$log2FoldChange.gsam.res))
sum((stat.gsam$log2FoldChange.gsam.tpc.res - stat.gsam$log2FoldChange.gsam.res) > 0)
sum((stat.gsam$log2FoldChange.gsam.tpc.res - stat.gsam$log2FoldChange.gsam.res) <= 0)
abline(h=0)



# GSAM ----
## MES ----

plt.gsam <- results.out |>
  #dplyr::filter(TCGA.subtype.marker == "TCGA-MES") |>
  dplyr::select(gid,hugo_symbol,statistic.gsam.cor.tpc,log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res , TCGA.subtype.marker )  |> 
  dplyr::mutate(rank = rank(log2FoldChange.gsam.tpc.res)) |> 
  tidyr::pivot_longer(c(log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res)) |> 
  dplyr::rename(log2FC = value) |> 
  dplyr::rename(col = name ) |> 
  dplyr::mutate(col = case_when(
    TCGA.subtype.marker == "TCGA-MES" & col == "log2FoldChange.gsam.res" ~ "TCGA & naive",
    TCGA.subtype.marker == "TCGA-MES" & col == "log2FoldChange.gsam.tpc.res" ~ "TCGA & purity corrected",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-MES") & col == "log2FoldChange.gsam.res" ~ "remaining & naive",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-MES") & col == "log2FoldChange.gsam.tpc.res" ~ "remaining & purity corrected",
  )) |> 
  dplyr::filter(col != "remaining & naive")


col.cor <- '#A15E5E'
col.naive <- '#55864E'

median.corrected <- plt.gsam |>
  dplyr::filter(col == "TCGA & purity corrected") |>
  dplyr::pull(log2FC) |>
  median()
median.naive <- plt.gsam |>
  dplyr::filter(col == "TCGA & naive") |>
  dplyr::pull(log2FC) |>
  median()

labels <- data.frame(log2FC = c(median.naive, median.corrected),
                     col=c('TCGA & naive','TCGA & purity corrected')) |> 
  dplyr::mutate(label = paste0("LFC median: ", round(log2FC,2))) |> 
  dplyr::mutate(hugo_symbol = label) |>  
  dplyr::mutate(gid = col) |>  
  dplyr::mutate(fill = col) |> 
  dplyr::mutate(statistic.gsam.cor.tpc=7.5)



ggplot(plt.gsam, aes(x=log2FC, y=statistic.gsam.cor.tpc,fill=col,label=hugo_symbol,group=gid )) +
  geom_vline(xintercept=0, lwd=0.5, color = "gray20") +  
  geom_point(data = plt.gsam |>  dplyr::filter(col == "remaining & purity corrected"),cex=0.1,col="gray") +
  

  geom_vline(xintercept=median.naive, lwd=0.75, color = col.naive,lty=2) +
  geom_vline(xintercept=median.corrected, lwd=0.75, color = col.cor,lty=2) +
  ggrepel::geom_text_repel(data=labels) +
  #geom_point(data = plt.gsam |>  dplyr::filter(col == "TCGA & naive"),pch=21,size=2) +
  ggplot2::geom_path(
    data = plt.gsam |>  dplyr::filter(grepl("TCGA",col)),
    lineend = "butt",
    linejoin = "mitre",
    lwd=0.5,
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125*0.65 , "inches")))  +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits=c(-2.25,2.25)) +
  scale_color_manual(values=c('TCGA & purity corrected'=col.cor,
                              'TCGA & naive'=col.naive))+
  scale_fill_manual(values=c('TCGA & purity corrected'=col.cor,
                              'TCGA & naive'=col.naive)) +
  geom_rug(data = plt.gsam |> dplyr::filter(col == "TCGA & purity corrected"),sides="b",aes(col=col),show.legend=F)+
  geom_rug(data = plt.gsam |> dplyr::filter(col == "TCGA & naive"),sides="t",aes(col=col),show.legend=F) +
  labs(x = "log2FC primary vs. recurrence", y="Correlation t-statistic with tumour purity")




ggsave('output/figures/2022_figure_SX_MES.pdf', width=8.3 / 2 ,height=8.3/2.9 , scale=2.1)



## PN ----


plt.gsam <- results.out |>
  dplyr::select(gid,hugo_symbol,statistic.gsam.cor.tpc,log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res , TCGA.subtype.marker )  |> 
  dplyr::mutate(rank = rank(log2FoldChange.gsam.tpc.res)) |> 
  tidyr::pivot_longer(c(log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res)) |> 
  dplyr::rename(log2FC = value) |> 
  dplyr::rename(col = name ) |> 
  dplyr::mutate(col = case_when(
    TCGA.subtype.marker == "TCGA-PN" & col == "log2FoldChange.gsam.res" ~ "TCGA & naive",
    TCGA.subtype.marker == "TCGA-PN" & col == "log2FoldChange.gsam.tpc.res" ~ "TCGA & purity corrected",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-PN") & col == "log2FoldChange.gsam.res" ~ "remaining & naive",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-PN") & col == "log2FoldChange.gsam.tpc.res" ~ "remaining & purity corrected",
  )) |> 
  dplyr::filter(col != "remaining & naive")


col.cor <- '#A15E5E'
col.naive <- '#55864E'

median.corrected <- plt.gsam |>
  dplyr::filter(col == "TCGA & purity corrected") |>
  dplyr::pull(log2FC) |>
  median()
median.naive <- plt.gsam |>
  dplyr::filter(col == "TCGA & naive") |>
  dplyr::pull(log2FC) |>
  median()

labels <- data.frame(log2FC = c(median.naive, median.corrected),
                     col=c('TCGA & naive','TCGA & purity corrected')) |> 
  dplyr::mutate(label = paste0("LFC median: ", round(log2FC,2))) |> 
  dplyr::mutate(hugo_symbol = label) |>  
  dplyr::mutate(gid = col) |>  
  dplyr::mutate(fill = col) |> 
  dplyr::mutate(statistic.gsam.cor.tpc=7.5)



ggplot(plt.gsam, aes(x=log2FC, y=statistic.gsam.cor.tpc,fill=col,label=hugo_symbol,group=gid )) +
  geom_vline(xintercept=0, lwd=0.5, color = "gray20") +  
  geom_point(data = plt.gsam |>  dplyr::filter(col == "remaining & purity corrected"),cex=0.1,col="gray") +
  
  
  geom_vline(xintercept=median.naive, lwd=0.75, color = col.naive,lty=2) +
  geom_vline(xintercept=median.corrected, lwd=0.75, color = col.cor,lty=2) +
  ggrepel::geom_text_repel(data=labels) +
  geom_text(data = plt.gsam |>  dplyr::filter(grepl("TCGA & naive",col) & log2FC > 0.3 & statistic.gsam.cor.tpc < 0),hjust = 0) +
  #geom_point(data = plt.gsam |>  dplyr::filter(col == "TCGA & naive"),pch=21,size=2) +
  ggplot2::geom_path(
    data = plt.gsam |>  dplyr::filter(grepl("TCGA",col)),
    lineend = "butt",
    linejoin = "mitre",
    lwd=0.5,
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125*0.65 , "inches")))  +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits=c(-2.25,2.25)) +
  scale_color_manual(values=c('TCGA & purity corrected'=col.cor,
                              'TCGA & naive'=col.naive))+
  scale_fill_manual(values=c('TCGA & purity corrected'=col.cor,
                             'TCGA & naive'=col.naive)) +
  geom_rug(data = plt.gsam |> dplyr::filter(col == "TCGA & purity corrected"),sides="b",aes(col=col),show.legend=F)+
  geom_rug(data = plt.gsam |> dplyr::filter(col == "TCGA & naive"),sides="t",aes(col=col),show.legend=F) +
  labs(x = "log2FC primary vs. recurrence", y="Correlation t-statistic with tumour purity")




## CL ----


plt.gsam <- results.out |>
  dplyr::select(gid,hugo_symbol,statistic.gsam.cor.tpc,log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res , TCGA.subtype.marker )  |> 
  dplyr::mutate(rank = rank(log2FoldChange.gsam.tpc.res)) |> 
  tidyr::pivot_longer(c(log2FoldChange.gsam.res,log2FoldChange.gsam.tpc.res)) |> 
  dplyr::rename(log2FC = value) |> 
  dplyr::rename(col = name ) |> 
  dplyr::mutate(col = case_when(
    TCGA.subtype.marker == "TCGA-CL" & col == "log2FoldChange.gsam.res" ~ "TCGA & naive",
    TCGA.subtype.marker == "TCGA-CL" & col == "log2FoldChange.gsam.tpc.res" ~ "TCGA & purity corrected",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-CL") & col == "log2FoldChange.gsam.res" ~ "remaining & naive",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-CL") & col == "log2FoldChange.gsam.tpc.res" ~ "remaining & purity corrected",
  )) |> 
  dplyr::filter(col != "remaining & naive")


col.cor <- '#A15E5E'
col.naive <- '#55864E'

median.corrected <- plt.gsam |>
  dplyr::filter(col == "TCGA & purity corrected") |>
  dplyr::pull(log2FC) |>
  median()
median.naive <- plt.gsam |>
  dplyr::filter(col == "TCGA & naive") |>
  dplyr::pull(log2FC) |>
  median()

labels <- data.frame(log2FC = c(median.naive, median.corrected),
                     col=c('TCGA & naive','TCGA & purity corrected')) |> 
  dplyr::mutate(label = paste0("LFC median: ", round(log2FC,2))) |> 
  dplyr::mutate(hugo_symbol = label) |>  
  dplyr::mutate(gid = col) |>  
  dplyr::mutate(fill = col) |> 
  dplyr::mutate(statistic.gsam.cor.tpc=7.5)



ggplot(plt.gsam, aes(x=log2FC, y=statistic.gsam.cor.tpc,fill=col,label=hugo_symbol,group=gid )) +
  geom_vline(xintercept=0, lwd=0.5, color = "gray20") +  
  geom_point(data = plt.gsam |>  dplyr::filter(col == "remaining & purity corrected"),cex=0.1,col="gray") +
  
  
  geom_vline(xintercept=median.naive, lwd=0.75, color = col.naive,lty=2) +
  geom_vline(xintercept=median.corrected, lwd=0.75, color = col.cor,lty=2) +
  ggrepel::geom_text_repel(data=labels) +
  geom_text(data = plt.gsam |>  dplyr::filter(grepl("TCGA & naive",col) & log2FC > 0.3 & statistic.gsam.cor.tpc < 0),hjust = 0) +
  #geom_point(data = plt.gsam |>  dplyr::filter(col == "TCGA & naive"),pch=21,size=2) +
  ggplot2::geom_path(
    data = plt.gsam |>  dplyr::filter(grepl("TCGA",col)),
    lineend = "butt",
    linejoin = "mitre",
    lwd=0.5,
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125*0.65 , "inches")))  +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits=c(-2.25,2.25)) +
  scale_color_manual(values=c('TCGA & purity corrected'=col.cor,
                              'TCGA & naive'=col.naive))+
  scale_fill_manual(values=c('TCGA & purity corrected'=col.cor,
                             'TCGA & naive'=col.naive)) +
  geom_rug(data = plt.gsam |> dplyr::filter(col == "TCGA & purity corrected"),sides="b",aes(col=col),show.legend=F)+
  geom_rug(data = plt.gsam |> dplyr::filter(col == "TCGA & naive"),sides="t",aes(col=col),show.legend=F) +
  labs(x = "log2FC primary vs. recurrence", y="Correlation t-statistic with tumour purity")





# GLASS ----




plt.glass <- results.out |>
  dplyr::select(gid,hugo_symbol,`statistic.t.glass-2022.cor.tpc`,`log2FoldChange.glass-2022.res`,`log2FoldChange.glass-2022.tpc.res` , TCGA.subtype.marker )  |> 
  dplyr::mutate(rank = rank(`log2FoldChange.glass-2022.tpc.res`)) |> 
  tidyr::pivot_longer(c(`log2FoldChange.glass-2022.res`,`log2FoldChange.glass-2022.tpc.res`)) |> 
  dplyr::rename(log2FC = value) |> 
  dplyr::rename(col = name ) |> 
  dplyr::mutate(col = case_when(
    TCGA.subtype.marker == "TCGA-MES" & col == "log2FoldChange.glass-2022.res" ~ "TCGA & naive",
    TCGA.subtype.marker == "TCGA-MES" & col == "log2FoldChange.glass-2022.tpc.res" ~ "TCGA & purity corrected",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-MES") & col == "log2FoldChange.glass-2022.res" ~ "remaining & naive",
    (is.na(TCGA.subtype.marker) | TCGA.subtype.marker != "TCGA-MES") & col == "log2FoldChange.glass-2022.tpc.res" ~ "remaining & purity corrected",
  )) |> 
  dplyr::filter(col != "remaining & naive")


col.cor <- '#A15E5E'
col.naive <- '#55864E'

median.corrected <- plt.glass |>
  dplyr::filter(col == "TCGA & purity corrected") |>
  dplyr::pull(log2FC) |>
  median()
median.naive <- plt.glass |>
  dplyr::filter(col == "TCGA & naive") |>
  dplyr::pull(log2FC) |>
  median()

labels <- data.frame(log2FC = c(median.naive, median.corrected),
                     col=c('TCGA & naive','TCGA & purity corrected')) |> 
  dplyr::mutate(label = paste0("LFC median: ", round(log2FC,3))) |> 
  dplyr::mutate(hugo_symbol = label) |>  
  dplyr::mutate(gid = col) |>  
  dplyr::mutate(fill = col) |> 
  dplyr::mutate(`statistic.t.glass-2022.cor.tpc`=7.5)



ggplot(plt.glass, aes(x=log2FC, y=`statistic.t.glass-2022.cor.tpc`,fill=col,label=hugo_symbol,group=gid )) +
  geom_point(data = plt.glass |>  dplyr::filter(col == "remaining & purity corrected"),cex=0.1,col="gray") +
  geom_vline(xintercept=median.naive, lwd=0.75, color = col.naive,lty=2) +
  geom_vline(xintercept=median.corrected, lwd=0.75, color = col.cor,lty=2) +
  ggrepel::geom_text_repel(data=labels) +
  #geom_point(data = plt.glass |>  dplyr::filter(col == "TCGA & naive"),pch=21,size=2) +
  ggplot2::geom_path(
    data = plt.glass |>  dplyr::filter(grepl("TCGA",col)),
    lineend = "butt",
    linejoin = "mitre",
    lwd=0.5,
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125*0.65 , "inches")))  +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    panel.grid.major.x =  element_line(colour = 'grey20', size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_x_continuous(breaks = c(0), limits=c(-2.25,2.25)) +
  scale_color_manual(values=c('TCGA & purity corrected'=col.cor,
                              'TCGA & naive'=col.naive))+
  scale_fill_manual(values=c('TCGA & purity corrected'=col.cor,
                             'TCGA & naive'=col.naive)) +
  geom_rug(data = plt.glass |> dplyr::filter(col == "TCGA & purity corrected"),sides="b",aes(col=col),show.legend=F)+
  geom_rug(data = plt.glass |> dplyr::filter(col == "TCGA & naive"),sides="t",aes(col=col),show.legend=F)



# delta NMF ----

plt <- gsam.rna.metadata |> 
  dplyr::filter(.data$blacklist.pca == F) |>
  dplyr::filter(.data$pat.with.IDH == F) |>
  dplyr::filter(.data$sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
  dplyr::filter(.data$batch != "old") |>
  dplyr::filter(.data$tumour.percentage.dna >= 15) |>
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::select(sid,pid,resection,`NMF:150:1`, `NMF:150:2`, `NMF:150:3`) |> 
  tidyr::pivot_wider(id_cols = pid, names_from = resection, values_from = c(sid, `NMF:150:1`, `NMF:150:2`, `NMF:150:3`)) |> 
  dplyr::mutate(delta.NMF.1 = `NMF:150:1_r2` - `NMF:150:1_r1`) |> 
  dplyr::mutate(delta.NMF.2 = `NMF:150:2_r2` - `NMF:150:2_r1`) |> 
  dplyr::mutate(delta.NMF.3 = `NMF:150:3_r2` - `NMF:150:3_r1`) |> 
  
  dplyr::mutate(rank = order(order(delta.NMF.3))) |> 
  
  dplyr::mutate( sid_r1 = NULL ) |> 
  dplyr::mutate( sid_r2 = NULL ) |> 
  dplyr::mutate( `NMF:150:1_r1` = NULL ) |> 
  dplyr::mutate( `NMF:150:1_r2` = NULL ) |> 
  dplyr::mutate( `NMF:150:2_r1` = NULL ) |> 
  dplyr::mutate( `NMF:150:2_r2` = NULL ) |> 
  dplyr::mutate( `NMF:150:3_r1` = NULL ) |> 
  #dplyr::mutate( `NMF:150:3_r2` = NULL ) |> 
  
  tidyr::pivot_longer(cols = c(delta.NMF.1, delta.NMF.2, delta.NMF.3), names_to = "NMF.metafeature", values_to = "NMF.value")


ggplot(plt, aes(x=reorder(pid, rank), y=NMF.value)) +
  facet_grid(rows = vars(NMF.metafeature), scales = "free", space="free_y") +
  geom_point()



