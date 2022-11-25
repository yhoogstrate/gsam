#!/bin/env R


source('scripts/R/palette.R')


# Figure S8C ----

plt <- gsam.rna.metadata |>
  dplyr::filter(.data$blacklist.pca == F) |>
  dplyr::filter(.data$pat.with.IDH == F) |>
  dplyr::filter(.data$sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
  dplyr::filter(.data$batch != "old") |>
  dplyr::filter(.data$tumour.percentage.dna >= 15) |> 
  dplyr::mutate(resection = dplyr::recode(resection,'r1'='Primary','r2'='Recurrence'))


ggplot(plt, aes(x=tumour.percentage.dna, y=rna.signature.C4.neuron.2022, fill=resection))+
  geom_point(pch = 21, size = 2, col = "black", alpha = 0.85) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL, label = ..r.label..), method = "pearson",label.x=70) +
  #stat_smooth(method="loess",se=FALSE,aes(fill=NULL),col="red") +
  scale_fill_manual(values = c("Primary" = resection_colors[["R1"]],
                               "Recurrence" = resection_colors[["R2"]]), name = "Resection") +
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
    y = "C4 (neuron) signature", 
    x = "Tumor purity",
    col = NA, 
    caption = paste0("G-SAM: n=", length(unique(plt$sid)))
  ) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.1))


ggsave("output/figures/2022_figure_S8C.pdf", width=8.3 / 4,height=8.3 / 4, scale=2)




# Figure S8D ----

plt <- gsam.rna.metadata |>
  dplyr::filter(.data$blacklist.pca == F) |>
  dplyr::filter(.data$pat.with.IDH == F) |>
  dplyr::filter(.data$sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
  dplyr::filter(.data$batch != "old") |>
  dplyr::filter(.data$tumour.percentage.dna >= 15) |> 
  dplyr::mutate(resection = dplyr::recode(resection,'r1'='Primary','r2'='Recurrence')) |> 
  dplyr::left_join(
    gsam.gene.expression.all.vst |>
      as.data.frame() |> 
      tibble::rownames_to_column('gid') |> 
      dplyr::filter(grepl("RBFOX3",gid)) |> 
      tibble::column_to_rownames('gid') |> 
      t() |> 
      as.data.frame() |> 
      tibble::rownames_to_column('sid'),
    by=c('sid'='sid'),suffix=c('','')  )



ggplot(plt, aes(x=tumour.percentage.dna, y=`ENSG00000167281.19_4|RBFOX3|chr17:77085427-77512230(-)`, fill=resection))+
  geom_point(pch = 21, size = 2, col = "black", alpha = 0.85) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL, label = ..r.label..), method = "pearson",label.x=70) +
  scale_fill_manual(values = c("Primary" = resection_colors[["R1"]],
                               "Recurrence" = resection_colors[["R2"]]), name = "Resection") +
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
    y = "RBFOX3", 
    x = "Tumor purity",
    col = NA, 
    caption = paste0("G-SAM: n=", length(unique(plt$sid)))
  ) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.1))



ggsave("output/figures/2022_figure_S8D.pdf", width=8.3 / 4,height=8.3 / 4, scale=2)



# Figure S8 CD cor ----
# 
# 
# tmp <- gsam.gene.expression.all.vst |>
#   as.data.frame() |> 
#   tibble::rownames_to_column('gid') |> 
#   dplyr::filter(gid %in% (results.out |> dplyr::filter(C4.2022) |>  dplyr::pull('gid'))) |> 
#   dplyr::left_join(results.out |>  dplyr::select(gid, hugo_symbol), by=c('gid'='gid'),suffix=c('','')) |> 
#   tibble::column_to_rownames('hugo_symbol') |> 
#   dplyr::mutate(gid=NULL) |> 
#   t() |> 
#   as.data.frame() |> 
#   tibble::rownames_to_column('sid')
# 
# 
# plt <- gsam.rna.metadata |>
#   dplyr::filter(.data$blacklist.pca == F) |>
#   dplyr::filter(.data$pat.with.IDH == F) |>
#   dplyr::filter(.data$sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
#   dplyr::filter(.data$batch != "old") |>
#   dplyr::filter(.data$tumour.percentage.dna >= 15) |> 
#   #dplyr::mutate(resection = dplyr::recode(resection,'r1'='Primary','r2'='Recurrence')) |> 
#   dplyr::select(sid, tumour.percentage.dna) |> 
#   dplyr::left_join(tmp , by=c('sid'='sid'),suffix=c('',''))  |> 
#   tibble::column_to_rownames('sid')
# 
# corrplot::corrplot(cor(plt))
# 
