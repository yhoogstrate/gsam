#!/usr/bin/env R

# load data ----


if(!exists('gsam.rna.metadata')) {
  source("scripts/load_G-SAM_metadata.R")
}


# overview ----


plt.ids <- gsam.rna.metadata |> 
  dplyr::filter(.data$blacklist.pca == F) |>
  dplyr::filter(.data$pat.with.IDH == F) |>
  dplyr::filter(.data$sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
  dplyr::filter(.data$batch != "old") |>
  dplyr::filter(.data$tumour.percentage.dna >= 15) |>
  dplyr::select(c("pid", "sid", "blacklist.pca")) |>
  dplyr::mutate(resection = gsub("^...(.).*$", "R\\1", sid)) |>
  dplyr::mutate(pid = as.factor(as.character(pid))) # refactor


plt.single <- plt.ids |>
  dplyr::left_join(
    gsam.rna.metadata |>
      dplyr::select(c("sid", "tumour.percentage.dna", "wt.reads.v3", "vIII.reads.v3", "vIII.percentage")),
    by = c("sid" = "sid")
  ) |>
  dplyr::filter(!is.na(vIII.percentage)) |>
  dplyr::filter(!is.na(tumour.percentage.dna))


# only complete pairs
filter1 <- plt.single |>
  dplyr::group_by(pid) |>
  dplyr::filter(n() == 2) |>
  dplyr::pull(pid)

# of which EGFRvIII percentage >= 1.0
filter2 <- plt.single |>
  dplyr::filter(vIII.percentage >= 1.0) |>
  dplyr::pull(pid)

filter <- intersect(filter1, filter2)


n.paired.v3 <- length(filter)



## F] Figure S6C, S6D - arrow plot ----


plt <- plt.single |>
  dplyr::filter(pid %in% filter) |>
  dplyr::arrange(pid, sid) |>
  dplyr::mutate(col = ifelse(vIII.percentage < 1.0, "EGFRvIII percentage < 1.0", "EGFRvIII percentage >= 1.0"))

status <- dplyr::full_join(
  plt |>  dplyr::filter(resection == "R1") |>  dplyr::select(pid, vIII.percentage) |> dplyr::rename_with(~ paste0(.x, ".R1")),
  plt |>  dplyr::filter(resection == "R2") |>  dplyr::select(pid, vIII.percentage) |> dplyr::rename_with(~ paste0(.x, ".R2")),
  by = c("pid.R1" = "pid.R2")
) |>
  dplyr::mutate(status.v3 = ifelse(vIII.percentage.R1 > vIII.percentage.R2, "EGFRvIII percentage decreases", "EGFRvIII percentage increases")) |>
  dplyr::rename(pid = pid.R1)

plt <- plt |>
  dplyr::left_join(status |> plyr::select(pid, status.v3), by = c("pid" = "pid"))


ggplot(plt, aes(x = vIII.percentage, y = tumour.percentage.dna, group = pid)) +
  facet_grid(cols = vars(status.v3), scales = "free") +
  geom_point(data = subset(plt, resection == "R2"), pch = 22, cex = 3, aes(fill = col)) +
  geom_point(data = subset(plt, resection == "R2"), pch = 19, cex = 0.5, col = "black") +
  geom_path(
    lineend = "butt",
    linejoin = "mitre",
    arrow = arrow(
      ends = "last",
      type = "closed",
      angle = 15,
      length = unit(0.15, "inches")
    ),
    alpha = 0.7
  ) +
  geom_point(data = subset(plt, resection == "R1"), pch = 21, cex = 3, aes(fill = col)) +
  geom_point(data = subset(plt, resection == "R1"), pch = 19, cex = 0.5, col = "black") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 100)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = "EGFRvIII percentage [ EGFRvIII / (EGFRwt+EGFRvIII) ]",
    y = "Tumor cell percentage",
    caption = paste0("G-SAM: n=", n.paired.v3, " pairs (EGFRvIII positive)")
  ) +
  scale_fill_manual(values = c("EGFRvIII percentage < 1.0" = "red", "EGFRvIII percentage >= 1.0" = "white"), name = NULL) +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    # panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  )


ggsave("output/figures/2022_Figure_S6C_S6D.pdf", width = 8.3 / 2, height = 8.3 / 3.4, scale = 2)



## F] Figure S6A, S6B  - vertical plot ----


delta.v3.percentage <-
  dplyr::full_join(
    plt.single %>% dplyr::filter(resection == "R1") %>% `colnames<-`(paste0(colnames(.), ".R1")),
    plt.single %>% dplyr::filter(resection == "R2") %>% `colnames<-`(paste0(colnames(.), ".R2")),
    by = c("pid.R1" = "pid.R2")
  ) %>%
  dplyr::rename(pid = pid.R1) %>%
  dplyr::filter(pid %in% filter) %>%
  dplyr::mutate(delta.v3.percentage = vIII.percentage.R2 - vIII.percentage.R1) %>%
  dplyr::select(pid, delta.v3.percentage)


plt <- plt.single %>%
  dplyr::filter(pid %in% filter) %>%
  dplyr::left_join(delta.v3.percentage, by = c("pid" = "pid")) %>%
  dplyr::rename(order = delta.v3.percentage) %>%
  dplyr::mutate(tumour.percentage.dna = ifelse(resection == "R1", tumour.percentage.dna + 0.001, tumour.percentage.dna))


plt <- rbind(
  plt %>% dplyr::select(pid, resection, order, vIII.percentage) %>%
    dplyr::rename(percentage = vIII.percentage) %>%
    dplyr::mutate(status = "EGFRvIII"),
  plt %>% dplyr::select(pid, resection, order, tumour.percentage.dna) %>%
    dplyr::rename(percentage = tumour.percentage.dna) %>%
    dplyr::mutate(status = "Tumor cell")
) %>%
  dplyr::mutate(status = factor(status, levels = c("Tumor cell", "EGFRvIII"))) %>%
  dplyr::mutate(col = ifelse(percentage < 1, "-", "+"))


ggplot(plt, aes(x = reorder(pid, order), y = percentage)) +
  facet_grid(rows = vars(status), scales = "free") +
  # geom_point(data = subset(plt, resection == "R2"), pch=22, cex=2.5, aes(fill=col)) +
  # geom_line(lwd=1.2) +
  geom_path(
    lineend = "butt",
    linejoin = "mitre",
    arrow = arrow(ends = "last", type = "closed", angle = 15, length = unit(0.1, "inches")), lwd = 1.2
  ) +
  geom_point(data = subset(plt, resection == "R1"), pch = 19, cex = 2.5, col = "black") +
  geom_point(data = subset(plt, resection == "R1"), pch = 19, cex = 1.25, col = "white") +


  # geom_point(data = subset(plt, resection == "R2"),  pch=19, cex=2.5, col="red") +

  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_manual(values = c("-" = "red", "+" = "white")) +
  labs(
    y = "Percentage",
    x = "Patient",
    caption = paste0("G-SAM: n=", length(unique(as.character(plt$pid))), " pairs (EGFRvIII positive)")
  )


ggsave("output/figures/2022_Figure_S6A_S6B.pdf", width = 8.3 / 2, height = 8.3 / 3.3, scale = 2)


