#!/usr/bin/env R 


# load libs ----


if(!exists('gsam.rna.metadata')) {
  source('scripts/load_G-SAM_metadata.R')
}

if(!exists('glass.gbm.rnaseq.metadata.all.samples')) {
  source('scripts/load_GLASS_data.R')
}

source('scripts/R/palette.R')



# F] Figure M1E, M1F - resection ~ purity scatters ----


plt <- rbind(
  gsam.rna.metadata |>
    dplyr::filter(blacklist.pca == F) |>
    dplyr::filter(pat.with.IDH == F) |>
    dplyr::filter(sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
    # dplyr::filter(tumour.percentage.dna >= 15) |> - has no effect on RNA levels here, so keep included
    dplyr::mutate(is.primary = resection == "r1") |>
    dplyr::select(
      sid,
      is.primary,
      tumour.percentage.dna
    ) |>
    dplyr::rename(purity = tumour.percentage.dna) |>
    dplyr::mutate(dataset = "G-SAM"),
  glass.gbm.rnaseq.metadata.all.samples |>
    # dplyr::filter(tumour.percentage.2022 >= 15) |> - has no effect on RNA levels here, so keep included
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
  dplyr::mutate(is.primary = dplyr::recode(as.character(is.primary), "TRUE" = "Primary", "FALSE" = "Recurrence"))




tmp.n.gsam <- plt |>
  dplyr::filter(dataset == "G-SAM") |>
  nrow()
tmp.n.glass <- plt |>
  dplyr::filter(dataset == "GLASS") |>
  nrow()



# add median stats
plt <- plt |> 
  rbind(
    plt |>
      dplyr::filter(dataset == "G-SAM") |>
      dplyr::filter(is.primary == "Primary" & !is.na(purity)) |>
      dplyr::mutate(purity = median(purity)) |>
      head(n = 1) |>
      dplyr::mutate(sid = "median"),
    plt |>
      dplyr::filter(dataset == "G-SAM") |>
      dplyr::filter(is.primary != "Primary" & !is.na(purity)) |>
      dplyr::mutate(purity = median(purity)) |>
      head(n = 1) |>
      dplyr::mutate(sid = "median"),
    plt |>
      dplyr::filter(dataset == "GLASS") |>
      dplyr::filter(is.primary == "Primary" & !is.na(purity)) |>
      dplyr::mutate(purity = median(purity)) |>
      head(n = 1) |>
      dplyr::mutate(sid = "median"),
    plt |>
      dplyr::filter(dataset == "GLASS") |>
      dplyr::filter(is.primary != "Primary" & !is.na(purity)) |>
      dplyr::mutate(purity = median(purity)) |>
      head(n = 1) |>
      dplyr::mutate(sid = "median")
  ) |> 
  dplyr::mutate(purity = round(purity,1))




ggplot(plt |> dplyr::filter(sid != 'median'), aes(x = is.primary, y = purity, fill = is.primary)) +
  facet_grid(cols = vars(dataset), scales = "free", space = "free_y") +
  geom_hline(yintercept = 15, color = "gray60") +
  geom_violin(draw_quantiles = c(), col = NA, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom(pch = 21, size = 3, col = "black", alpha = 0.85) +
  ylim(0, 100) +
  ggsignif::geom_signif(
    comparisons = list(c("Primary", "Recurrence")),
    test = "wilcox.test",
    col = "black",
    tip_length = 0
  ) +
  geom_hline(data = plt |> dplyr::filter(sid == 'median'), 
             aes(yintercept = purity,col = is.primary)) +
  scale_color_manual(values = c("Primary" = resection_colors[["R1"]], 
                               "Recurrence" = resection_colors[["R2"]]), name = "Resection", guide = "none") +
  scale_fill_manual(values = c("Primary" = resection_colors[["R1"]], 
                               "Recurrence" = resection_colors[["R2"]]), name = "Resection", guide = "none") +
  geom_text(data = plt |> dplyr::filter(sid == 'median'), aes(label = purity),nudge_y=3,nudge_x=-0.25) +
  theme_bw() +
  theme(
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(
    x = NULL, 
    y = "Tumor purity",
    col = NA, 
    caption = paste0("G-SAM: n=", tmp.n.gsam, "  -  GLASS: n=", tmp.n.glass)
  ) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.1))

ggsave("output/figures/2022_Figure_M1E_M1F.pdf", width=8.3 / 4,height=8.3/4, scale=2)
rm(plt)



# F] Figure S2D- subtype ~ purity ----


plt <- rbind(
  gsam.rna.metadata |>
    dplyr::filter(blacklist.pca == F) |>
    dplyr::filter(pat.with.IDH == F) |>
    dplyr::filter(sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |> 
    dplyr::filter(tumour.percentage.dna >= 15) |> 
    dplyr::mutate(is.primary = resection == "r1") |>
    dplyr::select(
      sid,
      is.primary,
      tumour.percentage.dna,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(purity = tumour.percentage.dna) |>
    dplyr::mutate(dataset = "G-SAM"),
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |> 
    dplyr::mutate(is.primary = resection == "TP") |>
    dplyr::select(
      aliquot_barcode,
      is.primary,
      tumour.percentage.2022,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(purity = tumour.percentage.2022) |>
    dplyr::mutate(dataset = "GLASS")
) |>
  dplyr::mutate(is.primary = dplyr::recode(as.character(is.primary), "TRUE" = "Primary", "FALSE" = "Recurrence")) |> 
  dplyr::mutate(GITS.150.svm.2022.subtype = dplyr::recode(GITS.150.svm.2022.subtype, 'Mesenchymal'='MES', 'Classical'='CL', 'Proneural'='PN'))


tmp.n.gsam <- plt |>
  dplyr::filter(dataset == "G-SAM") |>
  nrow()
tmp.n.glass <- plt |>
  dplyr::filter(dataset == "GLASS") |>
  nrow()


stats <- ggpubr::compare_means(`purity` ~ `GITS.150.svm.2022.subtype`, 
                               data = plt,
                               method="wilcox.test",
                               p.adjust.method ="holm",
                               ) |> # , group.by = "GITS.150.svm.2022.subtype_primary" geen facet
  dplyr::mutate(y_pos = c(108,96,104), p.adj = format.pval(p.adj, digits = 1)) |> 
  dplyr::mutate(GITS.150.svm.2022.subtype = NA)


ggplot(plt, aes(x = GITS.150.svm.2022.subtype, y = purity, fill = GITS.150.svm.2022.subtype)) +
  ggplot2::geom_violin(draw_quantiles = c(), col = NA, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom(pch = 21, size = 3, col = "black", alpha = 0.85) +
  ggplot2::geom_violin(draw_quantiles = c(0.5), col = "black", fill = alpha("white", 0)) +
  ggplot2::scale_y_continuous(limits = c(2, 111),breaks=c(0,25,50,75,100)) +
  ggsignif::geom_signif(
    data = stats,
    aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y_pos),
    manual = TRUE,
    tip_length = 0
  ) +
  scale_fill_manual(values = c(
    "MES" = as.character(subtype_colors['Mesenchymal']),
    "CL" = as.character(subtype_colors['Classical']),
    "PN" = as.character(subtype_colors['Proneural'])
                         ), name = "Subtype", guide = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.1)
  ) +
  labs(
    x = NULL, col = NA, y = "Tumor purity",
    caption = paste0("G-SAM: n=", tmp.n.gsam, "  -  GLASS: n=", tmp.n.glass)
  )

ggsave("output/figures/2022_Figure_S2D.pdf", width=8.3 / ((3/2) * 4),height=1.5, scale=2)
rm(plt, tmp.stats, tmp.n.glass, tmp.n.gsam)




