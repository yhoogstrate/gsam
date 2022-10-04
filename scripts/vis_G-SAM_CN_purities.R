#!/usr/bin/env R 

# load libs ----


# load data ----


dat <- read.delim("data/gsam/output/tables/cnv_copynumber-ratio.cns_all.txt", stringsAsFactors = F) |>
  dplyr::mutate(gene = NULL) |>
  dplyr::mutate(segment.length = end - start)


purities <- read.table("output/tables/cnv/tumor.percentage.estimate.txt")


stopifnot(purities$sample %in% dat$patient.id)



# analysis ----


render_gsam <- function(bc) {
  # bc = unique(purities$sample)[1]
  # print(bc)


  dat.pat <- dat |>
    dplyr::filter(.data$patient.id == bc) |>
    dplyr::mutate(chromosome = paste0("chr", chromosome)) |>
    dplyr::rename(segment.chr = chromosome) |>
    dplyr::filter(segment.chr %in% c("chrX", "chrY") == F) |>
    dplyr::filter(segment.length >= 5000000) |>
    dplyr::filter(log2 < 1.1)



  purity <- purities |>
    dplyr::filter(.data$sample == bc) |>
    dplyr::pull(pct) / 100






  plt <- dat.pat |>
    dplyr::mutate(id = paste0("id.", 1:n())) |>
    tidyr::pivot_longer(cols = c(`start`, `end`), names_to = "type", values_to = "pos") |>
    dplyr::mutate(segment.chr = factor(segment.chr, levels = gtools::mixedsort(unique(as.character(segment.chr))))) |>
    dplyr::mutate(pos = pos / 1000000)





  ggplot(plt, aes(x = pos, y = log2, group = id, col = segment.chr)) +
    facet_grid(cols = vars(segment.chr), scales = "free", space = "free") +
    geom_hline(yintercept = 0, lwd = 0.5, lty = 3, col = "black", alpha = 0.5) +
    geom_line(lwd = 2) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 4) / 2), lwd = 0.7, lty = 2, col = "black", alpha = 0.35) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 3) / 2), lwd = 0.7, lty = 2, col = "black", alpha = 0.35) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 1) / 2), lwd = 0.7, lty = 2, col = "black", alpha = 0.35) +
    ylim(-2, 2) +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold", size = rel(1)),
      # axis.text.x = element_blank(),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.45, hjust = 1)
    ) +
    scale_color_discrete(guide = "none") +
    labs(
      x = NULL, y = "CNVKit log2 copy ratio",
      caption = paste0("", dat.pat$aliquot_barcode[1], "  -  purity estimate: ", purity)
    ) +
    scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300))

  ggsave(paste0("output/figures/cnv/g-sam/2022/", bc, "_estimate_gg.pdf"), width = 8.3, height = 3.5, scale = 1.75)
}



pbapply::pblapply(unique(purities$sample), render_gsam)
