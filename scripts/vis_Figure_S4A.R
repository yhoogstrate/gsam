#!/usr/bin/env R 


if(!exists('results.out')) {
  source('scripts/load_results.out.R')
}


## FF] Figure S4-p01 - C4 ----


tmp.c4 <- results.out |>
  dplyr::filter(!is.na(.data$C4.2022)) |> 
  dplyr::filter(.data$C4.2022 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.npc1 <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.NPC1)) |> 
  dplyr::filter(.data$neftel.meta.modules.NPC1 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.npc2 <- results.out |> 
  dplyr::filter(!is.na(.data$neftel.meta.modules.NPC2)) |> 
  dplyr::filter(.data$neftel.meta.modules.NPC2 == T) |> 
  dplyr::filter(!is.na(hugo_symbol)) |> 
  dplyr::pull(hugo_symbol) |> 
  unique()
tmp.npc1.2 <- intersect(tmp.npc1, tmp.npc2)
tmp.npc1 <- setdiff(tmp.npc1, tmp.npc1.2)
tmp.npc2 <- setdiff(tmp.npc2, tmp.npc1.2)

tmp.c4.npc2 <- intersect(tmp.c4, tmp.npc2)
tmp.c4  <- setdiff(tmp.c4, tmp.c4.npc2)
tmp.npc2 <- setdiff(tmp.npc2, tmp.c4.npc2)


plt <- readxl::read_xlsx("data/1-s2.0-S0092867422005360-mmc3.xlsx",skip=1) |> 
  dplyr::rename(gene_symbol = NAME) |> 
  tidyr::pivot_longer(cols = -gene_symbol) |> 
  dplyr::filter(gene_symbol %in% c(tmp.c4, tmp.npc1, tmp.npc1.2, tmp.npc2, tmp.c4.npc2)) |> 
  #dplyr::filter(gene_symbol %in% c( tmp.c4.npc2)) |> 
  dplyr::mutate(facet = factor(case_when(
    gene_symbol %in% tmp.c4 ~ "C4",
    gene_symbol %in% tmp.npc1 ~ "NPC1",
    gene_symbol %in% tmp.npc1.2 ~ "NPC1+2",
    gene_symbol %in% tmp.npc2 ~ "NPC2",
    gene_symbol %in% tmp.c4.npc2 ~ "C4 + NPC2"
  ), levels=c("C4","NPC1","NPC1+2","NPC2", "C4 + NPC2"))) |> 
  dplyr::mutate(facet_y = factor(ifelse(grepl("tumor", name), "tumor","non-tumor"), levels=c("tumor","non-tumor"))) |> 
  dplyr::mutate(value = log(value))


ggplot(plt, aes(x = gene_symbol, y=name, size=value, col=value)) +
  facet_grid(cols = vars(facet),rows=vars(facet_y), scales = "free", space = "free") +
  geom_point() +
  scale_colour_gradient(low="gray", high="purple") +
  theme(line = element_blank()) +
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C4/NPC] in: log GLASS CYBERSORTx signatures")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  )


ggsave(paste0("output/figures/2022_Figure_S4-p01.pdf"),width=7.5*1.8, height=3.75,scale=1.2)


## FF] Figure S4-p02 - C3/OD & OPC-L -----


tmp.c3 <- results.out |>
  dplyr::filter(!is.na(.data$C3.2022)) |>
  dplyr::filter(.data$C3.2022 == T) |>
  dplyr::filter(!is.na(hugo_symbol)) |>
  dplyr::pull(hugo_symbol) |>
  unique()
tmp.opc <- results.out |>
  dplyr::filter(!is.na(.data$neftel.meta.modules.OPC)) |>
  dplyr::filter(.data$neftel.meta.modules.OPC == T) |>
  dplyr::filter(!is.na(hugo_symbol)) |>
  dplyr::pull(hugo_symbol) |>
  unique()
tmp.c3.opc <- intersect(tmp.c3, tmp.opc)
tmp.c3 <- setdiff(tmp.c3, tmp.c3.opc)
tmp.opc <- setdiff(tmp.opc, tmp.c3.opc)



plt <- readxl::read_xlsx("data/1-s2.0-S0092867422005360-mmc3.xlsx",skip=1) |> 
  dplyr::rename(gene_symbol = NAME) |> 
  tidyr::pivot_longer(cols = -gene_symbol) |> 
  dplyr::filter(gene_symbol %in% c(tmp.c3, tmp.opc, tmp.c3.opc)) |> 
  dplyr::mutate(facet = factor(case_when(
    gene_symbol %in% tmp.c3 ~ "C3",
    gene_symbol %in% tmp.opc ~ "OPC-like",
    gene_symbol %in% tmp.c3.opc ~ "C3 + OPC-like"
  ), levels=c("C3","OPC-like","C3 + OPC-like"))) |> 
  dplyr::mutate(facet_y = factor(ifelse(grepl("tumor", name), "tumor","non-tumor"), levels=c("tumor","non-tumor"))) |> 
  dplyr::mutate(value = log(value))


ggplot(plt, aes(x = gene_symbol, y=name, size=value, col=value)) +
  facet_grid(cols = vars(facet),rows=vars(facet_y), scales = "free", space = "free") +
  geom_point() +
  scale_colour_gradient(low="gray", high="blue") +
  theme(line = element_blank()) +
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5)) +
  labs(x = paste0("Features [C3/OPC-like] in: log GLASS CYBERSORTx signatures")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  )

ggsave(paste0("output/figures/2022_Figure_S4-p02.pdf"),width=7.5*1.8, height=3.75,scale=1.2)



