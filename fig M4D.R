#!/usr/bin/env R


library(tidyverse)
library(ggplot2)



resection_colors <- c(
  "primary" = "#55864e",
  "Primary Res." = "#55864e",
  "R1" = "#55864e" ,
  
  "recurrence" = "#a15e5e",
  "Recurrent Res." = "#a15e5e",
  "R2" = "#a15e5e"  # #a1805e
)




data.per.tile <- readRDS(file="/tmp/data.per.tile.Rds")



## F] Figure M4D ----


plt <- data.per.tile |> 
  dplyr::filter(!is.na(n.cells)) |>  #
  dplyr::mutate(facet_x = grepl("control", sid)) |>
  dplyr::mutate(facet_x = ifelse(facet_x, "control", "tumor")) |>
  dplyr::mutate(facet_x = factor(facet_x, levels = c("tumor", "control")))


median.lines <- plt |>
  dplyr::filter(resection %in% c("primary", "recurrence")) |>
  dplyr::group_by(sid) |>
  dplyr::summarise(
    log.ratio.neurons = log2(sum(n.neuron) / (sum(n.tumor) + sum(n.negative))),
    frac.neurons = sum(n.neuron) / sum(n.cells) * 100,
    Image = unique(Image),
    sid = unique(sid),
    pid = unique(pid),
    resection = unique(resection),
    n.tumor.low = sum(tumor.low.high == "tumor low"),
    n.tumor.high = sum(tumor.low.high == "tumor high"),
    median.n.cells = median(n.cells)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(facet_x = 'tumor') |> 
  dplyr::mutate(proportion.tumor.low.tiles = n.tumor.low / (n.tumor.low + n.tumor.high))



# stats should be paired
stats <- median.lines |>
  tidyr::pivot_wider(id_cols = pid, names_from = resection, values_from = c(log.ratio.neurons, frac.neurons))

wilcox.test(stats$frac.neurons_primary, stats$frac.neurons_recurrence, paired=T)
p = wilcox.test(stats$frac.neurons_primary, stats$frac.neurons_recurrence, paired=T)$p.value



n.pairs <- plt |>
  dplyr::filter(resection != "control") |> 
  dplyr::pull(pid) |> 
  unique() |>
  length()
n.controls <- plt |>
  dplyr::filter(resection == "control") |> 
  dplyr::pull(sid) |> 
  unique() |>
  length()





ggplot(plt, aes(x = sid, y = frac.neurons, col = resection, group = pid)) +
  facet_grid(cols = vars(facet_x), scales = "free", space = "free") +
  ggbeeswarm::geom_quasirandom(cex = 0.3) +
  geom_line(data = median.lines, col = "black", lwd = 0.5) +
  scale_color_manual(values = c(resection_colors[c("primary", "recurrence")], "control" = "gray40")) +
  theme_bw() +
  geom_text(data = data.frame(sid = "AAV1", pid = "-", facet_x = 'tumor', frac.neurons = 75, resection = "primary", label = paste0("P = ", format.pval(p, digits = 2))), aes(label = label), hjust = 0, col = "black") +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.text.x = element_text(angle = 90, vjust = 0.45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(
    caption = paste0("tumor pairs=", n.pairs, " normals: n=", n.controls),
    x = NULL,
    y = "fraction neurons in tile"
  )


