#!/usr/bin/env R 


# load libs ----

# get data ----

# make plot ----


plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    #dplyr::filter(tumour.percentage.dna >= 15) |> - has no effect on RNA levels here, so keep included
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      is.primary,
      tumour.percentage.dna
    ) |> 
    dplyr::rename(purity = tumour.percentage.dna) |>
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    #dplyr::filter(tumour.percentage.2022 >= 15) |> - has no effect on RNA levels here, so keep included
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      is.primary,
      tumour.percentage.2022,
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(purity = tumour.percentage.2022) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::mutate(is.primary = dplyr::recode(as.character(is.primary), 'TRUE'='Primary', 'FALSE'='Recurrence'))



tmp.n.gsam <- plt |> 
  dplyr::filter(dataset == "G-SAM") |> 
  nrow()
tmp.n.glass <- plt |> 
  dplyr::filter(dataset == "GLASS") |> 
  nrow()



ggplot(plt, aes(x=is.primary,y=purity, fill = is.primary)) +
  facet_grid(cols=vars(dataset), scales = "free", space="free_y") + 
  geom_hline(yintercept = 15, color = "gray60") +
  geom_violin(draw_quantiles = c(), col=NA, alpha=0.2) +
  #geom_jitter( position=position_jitter(0.2), size=2.5, pch=21, col="black") +
  ggbeeswarm::geom_quasirandom(pch=21,size=3,col="black",alpha=0.85) +
  geom_violin(draw_quantiles = c(0.5), col="black", fill=alpha("white",0)) +
  ylim(0, 100) +
  ggsignif::geom_signif(
    comparisons = list(c("Primary" ,  "Recurrence")),
    test="wilcox.test",
    col="black"
  ) + 
  scale_fill_manual(values =  c('Primary'=resection_colors[['R1']],'Recurrence'=resection_colors[['R2']]), name="Resection", guide="none" ) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(
    x = NULL, col = NA, y = "Tumor purity",
    caption=paste0("G-SAM: n=",tmp.n.gsam,"  -  GLASS: n=",tmp.n.glass)  ) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1)  )



ggsave("output/figures/2022_figure_1de_purity_scatter.pdf", width=8.3 / 4,height=8.3/4, scale=2)
ggsave("output/figures/2022_figure_1de_purity_scatter.svg", width=8.3 / 4,height=8.3/4, scale=2)




