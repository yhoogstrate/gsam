#!/usr/bin/env R

# settings ----


options(warnPartialMatchDollar = TRUE) # https://stackoverflow.com/questions/32854683/data-frames-in-r-name-autocompletion


# load libs ----


# load data ----


source('scripts/load_results.out.R')


## figure S8a (NE / NPC) ----
# a. eerst losse PCA bepaling, dan correlatie daar tussen?


labels <- results.out |>
  dplyr::filter(!is.na(gid)) |>
  dplyr::filter(.data$C4.2022 | .data$neftel.meta.modules.NPC1 | .data$neftel.meta.modules.NPC2) |>
  dplyr::select(gid, ensembl_id, hugo_symbol, `C4.2022`, neftel.meta.module.NPC1, neftel.meta.module.NPC2) |>
  dplyr::mutate(col = case_when(
    `C4.2022` & !`neftel.meta.module.NPC1` & !`neftel.meta.module.NPC2` ~ "C4 (Neuron)",
    !`C4.2022` &  `neftel.meta.module.NPC1` & !`neftel.meta.module.NPC2` ~ "NPC1",
    !`C4.2022` & !`neftel.meta.module.NPC1` &  `neftel.meta.module.NPC2` ~ "NPC2",
    !`C4.2022` &  `neftel.meta.module.NPC1` &  `neftel.meta.module.NPC2` ~ "NPC1 & 2",
    T ~ "ambiguous"
  )) |>
  dplyr::filter(col != "ambiguous") |>
  dplyr::mutate(C4.2022 = NULL) |>
  dplyr::mutate(neftel.meta.modules.NPC1 = NULL) |>
  dplyr::mutate(neftel.meta.modules.NPC2 = NULL)

labels <- labels |> 
  dplyr::mutate(col = ifelse(col == "C4 (Neuron)",paste0(col," [n=",sum(labels$col == "C4 (Neuron)"),"]"),col)) |> 
  dplyr::mutate(col = ifelse(col == "NPC1",paste0(col," [n=",sum(labels$col == "NPC1"),"]"),col)) |> 
  dplyr::mutate(col = ifelse(col == "NPC2",paste0(col," [n=",sum(labels$col == "NPC2"),"]"),col)) |> 
  dplyr::mutate(col = ifelse(col == "NPC1 & 2",paste0(col," [n=",sum(labels$col == "NPC1 & 2"),"]"),col))


# labels |> 
#   dplyr::pull(col) |> 
#   unique()



plt <- labels |>
  dplyr::select("gid", "hugo_symbol") |>
  dplyr::left_join(gsam.gene.expression.all.vst |>
                     as.data.frame() |>
                     tibble::rownames_to_column("gid"),
                   by = c("gid" = "gid")
  ) |>
  dplyr::mutate(gid = NULL) |>
  tibble::column_to_rownames("hugo_symbol")

res.pca <- prcomp(t(plt), scale = TRUE)

factoextra::fviz_pca_var(res.pca, col.var = labels$col, repel = TRUE,geom=c('arrow')) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.1)
  ) +
  scale_color_manual(values = c(
    `C4 (Neuron) [n=254]` = "#ff5f68",
    `NPC1 [n=39]` = "#6ba6e5",
    `NPC2 [n=33]` = "#6be5d9",
    `NPC1 & 2 [n=10]` = "#6bcee5"
  )) +
  labs(caption = paste0("G-SAM: n=", ncol(plt), " samples"))


ggsave("output/figures/2022_figure_S8a.pdf", width = 8.3 / 4, height = 8.3 / 4, scale = 2)




## figure S8b (OD / OPC) ----



labels <- results.out |> 
  dplyr::filter(!is.na(gid)) |> 
  dplyr::filter(.data$C3.2022 | .data$neftel.meta.modules.OPC) |>
  dplyr::select(gid, ensembl_id, hugo_symbol, `C3.2022`, neftel.meta.module.OPC) |>
  dplyr::mutate(col = case_when(
    `C3.2022` & !`neftel.meta.module.OPC` ~ "C3 (Oligodendrocyte)",
    !`C3.2022` &  `neftel.meta.module.OPC` ~ "OPC",
    T ~ "ambiguous"
  )) |>
  dplyr::filter(col != "ambiguous") |>
  dplyr::mutate(C3.2022 = NULL) |> 
  dplyr::mutate(neftel.meta.modules.OPC = NULL)


labels <- labels |> 
  dplyr::mutate(col = ifelse(col == "C3 (Oligodendrocyte)",paste0(col," [n=",sum(labels$col == "C3 (Oligodendrocyte)"),"]"),col)) |> 
  dplyr::mutate(col = ifelse(col == "OPC",paste0(col," [n=",sum(labels$col == "OPC"),"]"),col))


labels |>
  dplyr::pull(col) |>
  unique()



plt <- labels |>
  dplyr::select("gid", "hugo_symbol") |>
  dplyr::left_join(gsam.gene.expression.all.vst |>
                     as.data.frame() |>
                     tibble::rownames_to_column("gid"),
                   by = c("gid" = "gid")
  ) |>
  dplyr::mutate(gid = NULL) |>
  tibble::column_to_rownames("hugo_symbol")

res.pca <- prcomp(t(plt), scale = TRUE)

factoextra::fviz_pca_var(res.pca, col.var = labels$col, repel = T,geom=c('arrow')) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.1)
  ) +
  scale_color_manual(values = c(
    `C3 (Oligodendrocyte) [n=125]` = "#ff5f68",
    `OPC [n=48]` = "#6ba6e5"
  )) +
  labs(caption = paste0("G-SAM: n=", ncol(plt), " samples"))


ggsave("output/figures/2022_figure_S8b.pdf", width = 8.3 / 4, height = 8.3 / 4, scale = 2)


