#!/usr/bin/env R 

# library ----


source('scripts/R/palette.R')

library(ggplot2)
library(dplyr)



# load data ----

raw.per.cell       = as.data.frame(Reduce(rbind, readRDS("data/gsam/Stainings/IF (old)/Experiments revisions/GFAP_NeuN/data_analysis/data_correct_order.RDS")))
data.per.tile      = readRDS("data/gsam/Stainings/IF (old)/Experiments revisions/GFAP_NeuN/data_analysis/data_correct_order_ann.RDS")
thresholds         = readxl::read_xlsx(path = "data/gsam/Stainings/IF (old)/Experiments revisions/GFAP_NeuN/data_analysis/sample_annotation_NeuN.xlsx")
data.per.resection = readRDS("data/gsam/Stainings/IF (old)/Experiments revisions/GFAP_NeuN/data_analysis/gsam-signatures.Rds") |> 
  tidyr::drop_na(rna.signature.C4.neuron.2022) |>
  dplyr::mutate(sid = substr(sid, 1,4))



# processing ----


#check distribution individual ROIs
raw.per.cell |>
  dplyr::filter(Image == "batch 1 Levi GFAP Neun.czi - Scene #1") |>
  ggplot(aes(as.numeric(Centroid.X.Âµm), as.numeric(Centroid.Y.Âµm))) +
  geom_point(cex=0.1) +
  coord_equal()


tmp <- raw.per.cell |> 
  dplyr::select(c("Image", "Name", "Parent", "Centroid.X.Âµm", "Centroid.Y.Âµm", "Nucleus..FITC.mean", "Cell..Red610.mean"))  |> 
  dplyr::rename(Tile = Parent) |> 
  dplyr::mutate(`Image.filename` = Image) |> 
  dplyr::mutate(`Tile.QuPath.Name` = Tile) |> 
  dplyr::rename(`Centroid.X.µm` = "Centroid.X.Âµm") |> 
  dplyr::rename(`Centroid.Y.µm` = "Centroid.Y.Âµm") |> 
  dplyr::rename(`NeuN_nucleus` = 'Nucleus..FITC.mean') |> 
  dplyr::rename(`GFAP_cell` = 'Cell..Red610.mean')

tmp$Image     = gsub("batch ",                                 "b",     tmp$Image)
tmp$Image     = gsub(" Levi GFAP Neun.czi - Scene #",          "_s",    tmp$Image)
tmp$Image     = gsub(" levi GFAP Neun.czi - Scene #",          "_s",    tmp$Image)
tmp$Image     = gsub("GFAP-NE_control_rescan.czi - Scene #",   "bc1_s", tmp$Image)
tmp$Image     = gsub("GFAP-NE_control_rescan_2.czi - Scene #", "bc2_s", tmp$Image)
tmp$Tile      = gsub("Tile ",                                  "t",     tmp$Tile)
tmp$Name      = paste(tmp$Image, tmp$Tile, sep = "_")



data.per.tile$Image     = gsub("batch ",                                 "b",     data.per.tile$Image)
data.per.tile$Image     = gsub(" Levi GFAP Neun.czi - Scene #",          "_s",    data.per.tile$Image)
data.per.tile$Image     = gsub(" levi GFAP Neun.czi - Scene #",          "_s",    data.per.tile$Image)
data.per.tile$Image     = gsub("GFAP-NE_control_rescan.czi - Scene #",   "bc1_s", data.per.tile$Image)
data.per.tile$Image     = gsub("GFAP-NE_control_rescan_2.czi - Scene #", "bc2_s", data.per.tile$Image)
data.per.tile$Name2     = gsub("Tile ",                                  "t",     data.per.tile$Name2)
data.per.tile$Name      = paste(data.per.tile$Image, data.per.tile$Name2, sep = "_")
data.per.tile <- data.per.tile |> 
  dplyr::left_join(
    tmp |>
      dplyr::select(Image.filename , Image) |> 
      dplyr::distinct(),
    by=c('Image'='Image'),suffix=c('','')
  )  |> # Remove tiles full with eri's
  dplyr::filter(
    (
      (Image == "b3_s5") &
        (Name2 %in% c(
          "t138","t151","t159","t160","t164","t165","t169","t209","t212","t215","t225","t227",
          "t228","t229","t241","t248","t254","t26","t261","t263","t275","t298","t318","t345",
          "t350","t370","t390","t393","t399","t401","t406","t415","t430","t443","t445","t450",
          "t457","t464","t468","t494","t5","t501","t517","t520","t521","t53","t530","t538",
          "t550","t552","t556","t557","t566","t63","t66","t91","t180","t272","t528","t118",
          "t437","t485","t403","t540","t238","t373","t326","t346"))
    )
    == F)



tmp3 = tmp |>
  dplyr::group_by(Name) |>
  dplyr::summarize(xran = range(Centroid.X.µm)[2] - range(Centroid.X.µm)[1],
            yran = range(Centroid.Y.µm)[2] - range(Centroid.Y.µm)[1]) |>
  dplyr::filter(xran < 200, yran <200) |>
  dplyr::pull(Name)


tmp4 = data.per.tile |>
  dplyr::filter(Perimeter.Âµm < 801 & Perimeter.Âµm > 799) |>
  dplyr::pull(Name)


tmp = tmp |>
  dplyr::filter(Name %in% tmp3) |>
  dplyr::filter(Name %in% tmp4) |>
  dplyr::select(Image, Image.filename, Name, Tile, Tile.QuPath.Name, NeuN_nucleus, GFAP_cell, `Centroid.X.µm` , `Centroid.Y.µm` )

thr = thresholds |>
  dplyr::mutate(batch = paste0("b", batch),
         scene = paste0("s", scene),
         Image = paste(batch, scene, sep = "_"),
         resection = substr(sample_correct, 4,4)) |>
  dplyr::select(Image, sample_correct, resection, Neun_Threshold, GFAP_Threshold) |>
  dplyr::rename(sample = sample_correct)

data.per.cell = tmp |>
  dplyr::left_join(thr, by = "Image") |>
  dplyr::mutate(type = dplyr::case_when(NeuN_nucleus >= Neun_Threshold & GFAP_cell >= GFAP_Threshold ~ "Neuron",
                          NeuN_nucleus >= Neun_Threshold & GFAP_cell <= GFAP_Threshold ~ "Neuron",
                          NeuN_nucleus <= Neun_Threshold & GFAP_cell >= GFAP_Threshold ~ "Tumor",
                          NeuN_nucleus <= Neun_Threshold & GFAP_cell <= GFAP_Threshold ~ "Neg")) |>
  dplyr::select(!c(Neun_Threshold, GFAP_Threshold))



# generic loadings ----


tmp <- data.per.cell |>
  dplyr::group_by(Name) |>
  dplyr::summarise(
    n.tumor = sum(`type` == "Tumor"),
    n.neuron = sum(`type` == "Neuron"),
    n.negative = sum(`type` == "Neg"),
    #Image = unique(Image),
    Tile = unique(Tile),
    sid = unique(sample),
    resection = unique(resection)
  ) |>
  dplyr::mutate(sid = dplyr::recode(sid, # eye candy
    "control1" = "control-A",
    "control2" = "control-B",
    "control4" = "control-C",
    "control8" = "control-D",
    "control9" = "control-E",
  )) |>
  dplyr::mutate(resection = dplyr::case_when(
    resection == "1" ~ "primary",
    resection == "2" ~ "recurrence",
    resection == "t" ~ "control"
  )) |>
  dplyr::mutate(resection = factor(resection, levels=c('primary','recurrence','control'))) |> 
  dplyr::mutate(pid = gsub("^(...).+$", "\\1", sid)) |>
  dplyr::mutate(pid = ifelse(pid == "con", "control", pid)) |>
  dplyr::mutate(n.cells = n.neuron + n.negative + n.tumor) |>
  dplyr::mutate(frac.neurons = n.neuron / n.cells * 100) |>
  dplyr::mutate(log.ratio.neurons = log2((n.neuron + 1) / (n.tumor + n.negative + 1))) |>
  dplyr::mutate(tumor.low.high = dplyr::case_when(
    resection == "control" ~ "tumor low",
    n.cells > 125 ~ "tumor high",
    n.cells <= 125 ~ "tumor low"
  )) 


## medians:
tmp |> 
  dplyr::filter(resection == "control") |> 
  dplyr::pull(n.cells) |> 
  median()

tmp |> 
  dplyr::filter(resection != "control" & n.cells > 125) |> 
  dplyr::pull(n.cells) |> 
  median()

tmp |> 
  dplyr::filter(resection != "control" & n.cells < 125) |> 
  dplyr::pull(n.cells) |> 
  median()



## means:
tmp |> 
  dplyr::filter(resection == "control") |> 
  dplyr::pull(n.cells) |> 
  mean()

tmp |> 
  dplyr::filter(resection != "control" & n.cells > 125) |> 
  dplyr::pull(n.cells) |> 
  mean()

tmp |> 
  dplyr::filter(resection != "control" & n.cells < 125) |> 
  dplyr::pull(n.cells) |> 
  mean()



#stopifnot(length(intersect(colnames(data.per.tile), colnames(tmp))) == 1) # Name, identifier

data.per.tile <- data.per.tile |> 
  dplyr::left_join(tmp, by=c('Name'='Name'),suffix=c('',''))




# 1: fraction neurons increases ----
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


ggsave("output/figures/2022_Figure_M4C.pdf", width=8.3 / 2,height=8.3/5, scale=2)


saveRDS(data.per.tile, file="/tmp/data.per.tile.Rds") # export for sharing




# 2: separate low high tumor regions ----

## FF] Figure S5C - cut-off tumor low/high ----

plt <- data.per.tile |>
  dplyr::filter(!is.na(n.cells)) |>
  dplyr::arrange(n.cells) |>
  dplyr::mutate(x = 1:n()) |>
  dplyr::mutate(colr = ifelse(resection == "control", "control", "tumor"))

p1 <- ggplot(plt, aes(y = n.cells, x = x, col = colr)) +
  geom_hline(yintercept = 125, lwd = 0.65, color = "gray40", lty = "dotted") +
  geom_point(data = (plt |> dplyr::filter(colr != "control")), cex = 0.1, alpha = 0.5) +
  geom_point(data = (plt |> dplyr::filter(colr == "control")), cex = 0.1, alpha = 0.5) +
  scale_y_continuous(trans = "log2", limits = c(1, 512)) +
  scale_color_manual(values = c("tumor" = "red", "control" = "darkgreen")) +
  labs(
    x = "tiles ranked by cell count", y = "per tile cell count", color = NULL,
    caption = paste0(
      "n-samples: ", formatC(length(unique(plt$sid)), format = "f", big.mark = ",", digits = 0),
      "  -  n-tiles: ", formatC(nrow(plt), format = "f", big.mark = ",", digits = 0),
      "  -  n-cells: ", formatC(sum(plt$n.cells), format = "f", big.mark = ",", digits = 0)
    )
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  annotate(geom = "text", x = 1000, y = 320, label = "\nCut-off: 125", size = 3.5, hjust = 0)



plt <- data.per.tile |>
  dplyr::mutate(n.cells = n.neuron + n.negative + n.tumor) |>
  dplyr::arrange(n.cells) |>
  dplyr::mutate(x = 1:n()) |>
  dplyr::mutate(colr = ifelse(resection == "control", "control", "tumor")) |>
  dplyr::filter(colr == "tumor")

p2 <- ggplot(plt, aes(y = n.cells, x = x, col = colr)) +
  geom_hline(yintercept = 125, lwd = 0.65, color = "gray40", lty = "dotted") +
  geom_point(data = (plt |> dplyr::filter(colr != "control")), cex = 0.1, alpha = 0.5) +
  geom_point(data = (plt |> dplyr::filter(colr == "control")), cex = 0.1, alpha = 0.5) +
  scale_y_continuous(trans = "log2", limits = c(1, 512)) +
  scale_color_manual(values = c("tumor" = "red", "control" = "darkgreen")) +
  labs(
    x = "tiles ranked by cell count", y = "per tile cell count", color = NULL,
    caption = paste0(
      "n-samples: ", formatC(length(unique(plt$sid)), format = "f", big.mark = ",", digits = 0),
      "  -  n-tiles: ", formatC(nrow(plt), format = "f", big.mark = ",", digits = 0),
      "  -  n-cells: ", formatC(sum(plt$n.cells), format = "f", big.mark = ",", digits = 0)
    )
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  annotate(geom = "text", x = 1000, y = 320, label = "\nCut-off: 125", size = 3.5, hjust = 0)


p1 + p2



ggsave("output/figures/2022_Figure_S5C.pdf", width = 8.3 / 2, height = 8.3 / 4, scale = 2)




## FF] Figure S5E - high vs. low ~ diff neurons ----


plt <- data.per.tile |>
  dplyr::filter(!is.na(n.cells)) |> #
  dplyr::mutate(facet_x = grepl("control", sid)) |>
  dplyr::mutate(facet_x = ifelse(facet_x, "control", "tumor")) |>
  dplyr::mutate(facet_x = factor(facet_x, levels = c("tumor", "control")))


per.tile.stats <- data.per.tile |>
  dplyr::filter(resection != "control") |>
  dplyr::filter((pid == "CDD" & tumor.low.high == "tumor high") == F) |> # all neurons are in tumor high, cannot estimate median
  dplyr::filter((pid == "AAV" & tumor.low.high == "tumor low") == F) |> # too few data points
  dplyr::filter((pid == "AZG" & tumor.low.high == "tumor low") == F) |> # too few data points
  dplyr::filter((pid == "AZH" & tumor.low.high == "tumor low") == F) |> # too few data points
  dplyr::filter((pid == "BAT" & tumor.low.high == "tumor low") == F) |> # too few data points
  dplyr::filter((pid == "JAG" & tumor.low.high == "tumor low") == F) |> # too few data points
  dplyr::filter((pid == "JAK" & tumor.low.high == "tumor low") == F) |> # too few data points
  dplyr::group_by(sid, tumor.low.high) |>
  dplyr::summarise(
    # log.ratio.neurons = median(log.ratio.neurons), # this is wrong
    log.ratio.neurons = log2(sum(n.neuron) / (sum(n.tumor) + sum(n.negative))),
    frac.neurons = sum(n.neuron) / sum(n.cells) * 100, # fraction, human intuitive and extremely skewed, for non-parametric tests
    sid = unique(sid),
    pid = unique(pid),
    resection = unique(resection),
    tumor.low.high = unique(tumor.low.high),
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(facet_x = "tumor")


per.pair.stats <- per.tile.stats |>
  tidyr::pivot_wider(id_cols = pid, names_from = c(resection, tumor.low.high), values_from = c(log.ratio.neurons, frac.neurons))


p.tumor.high <- wilcox.test(
  per.pair.stats |>
    dplyr::filter(!is.na(`frac.neurons_primary_tumor high`) & !is.na(`frac.neurons_recurrence_tumor high`)) |>
    dplyr::pull(`frac.neurons_primary_tumor high`),
  per.pair.stats |>
    dplyr::filter(!is.na(`frac.neurons_primary_tumor high`) & !is.na(`frac.neurons_recurrence_tumor high`)) |>
    dplyr::pull(`frac.neurons_recurrence_tumor high`),
  paired = T
)$p.value

p.tumor.low <- wilcox.test(
  per.pair.stats |>
    dplyr::filter(!is.na(`frac.neurons_primary_tumor low`) & !is.na(`frac.neurons_recurrence_tumor low`)) |>
    dplyr::pull(`frac.neurons_primary_tumor high`),
  per.pair.stats |>
    dplyr::filter(!is.na(`frac.neurons_primary_tumor low`) & !is.na(`frac.neurons_recurrence_tumor low`)) |>
    dplyr::pull(`frac.neurons_recurrence_tumor high`),
  paired = T
)$p.value


df.pvals <- data.frame(
  sid = rep("AAV1", 2), # align above first sample
  pid = rep("-", 2),
  facet_x = rep("tumor", 2),
  frac.neurons = c(65, 90),
  resection = rep("primary", 2),
  tumor.low.high = c("tumor high", "tumor low"),
  label = c(
    paste0("P = ", format.pval(p.tumor.high, digits = 2)),
    paste0("P = ", format.pval(p.tumor.low, digits = 2))
  )
)


ggplot(
  plt,
  aes(x = sid, y = frac.neurons, col = resection, group = pid)
) +
  facet_grid(cols = vars(facet_x), rows = vars(tumor.low.high), scales = "free", space = "free_x") +
  ggbeeswarm::geom_quasirandom(cex = 0.3) +
  geom_line(data = per.tile.stats, col = "black", lwd = 0.5) +
  scale_color_manual(values = c(resection_colors[c("primary", "recurrence")], "control" = "gray40")) +
  theme_bw() +
  geom_text(data = df.pvals, aes(label = label), hjust = 0, col = "black") +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.text.x = element_text(angle = 90, vjust = 0.45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(
    caption = paste0("tumor pairs=", n.pairs, " normals: n=", n.controls),
    x = NULL,
    y = "fraction neurons"
  )


ggsave("output/figures/2022_Figure_S5E.pdf", width=8.3 / 2,height=8.3/4, scale=2)



# 3: per sample tile plots ----
## F] Figure M4C -- aggregate counts/tile ann ----


ssids <- data.per.tile |> 
  dplyr::filter(!is.na(sid) & !grepl("ontrol", sid)) |>
  dplyr::pull(Image.filename) |> 
  unique()

for(ssid in ssids) {
  #ssid = 'batch 1 Levi GFAP Neun.czi - Scene #1'
  ssid = "batch 3 Levi GFAP Neun.czi - Scene #5"
  print(ssid)
  
  plt.per.tile <- data.per.tile |> 
    dplyr::filter(!grepl("ontrol", sid)) |> 
    dplyr::filter(`Image.filename` == ssid) |> 
    dplyr::filter(!is.na(n.cells)) |> 
    dplyr::mutate(`Centroid.X.µm` = xmin) |> 
    dplyr::mutate(`Centroid.Y.µm` = ymin) |> 
    dplyr::mutate(neuron.fraction = cut(frac.neurons,breaks=10)) # nice trick - select ALL frac.neurons values and cut it in 10 slices / quantiles - use these as green overlay in the tiles
  
  # make quantiles numerical proportions
  plt.per.tile <- plt.per.tile |> 
    dplyr::mutate(nnf = 0) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[1], 1/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[2], 2/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[3], 3/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[4], 4/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[5], 5/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[6], 6/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[7], 7/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[8], 8/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[9], 9/11,nnf)) |> 
    dplyr::mutate(nnf = ifelse(neuron.fraction == levels(neuron.fraction)[10], 10/11,nnf)) |> 
    dplyr::mutate(filly = ymin + ((ymax-ymin) * nnf))
  
  
  plt.per.cell <- data.per.cell |> 
    dplyr::filter(Image.filename == ssid) |> 
    dplyr::mutate(col = "Nucleus")
    #dplyr::mutate(col = ifelse(type == "Neuron", "Neuron", "Not neuron"))
  
  ggplot(plt.per.cell, aes(x=`Centroid.X.µm`, y=`Centroid.Y.µm`)) +
    geom_point(cex=0.005,pch=21, col="gray80") +
    geom_rect(data=plt.per.tile |> dplyr::mutate(fill='proportion neurons (interval)'), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=filly,fill=fill), col=NA) +
    geom_rect(data=plt.per.tile, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,col=tumor.low.high), fill=NA,lwd=0.15) +
    coord_equal() +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold", size = rel(1)),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.45, hjust = 1),
      axis.ticks.x = element_blank(),
      
      axis.text.y = element_blank(), #element_text(angle = 90, vjust = 0.45, hjust = 1),
      axis.ticks.y = element_blank(),
      
      panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
    ) +
    labs(
      x = "location x-axis IF slide",
      y = "location y-axis IF slide",
      col=NULL,
      fill=NULL,
      caption = paste0(unique(plt.per.tile$sid), ": n-tiles=",format(length(unique(plt.per.tile$Name)),nsmall=1, big.mark=","),"  n-cells=",format(nrow(plt.per.cell),nsmall=1, big.mark=","))
    ) +
    scale_color_manual(values=c('tumor low'='blue',
                                'tumor high'='red'
    ))+
    scale_fill_manual(values=c('proportion neurons (interval)'=alpha("darkgreen",0.5)
    ))
  
  ggsave(paste0("output/figures/2022_Figure_M4C__",unique(plt.per.tile$sid),".pdf"), width=8.3 / 2,height=(8.3/ 2) * 1.25, scale=3)
}



