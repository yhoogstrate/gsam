#!/usr/bin/env R

# load data ----

source('scripts/load_G-SAM_metadata.R')
source('scripts/load_GLASS_data.R')



# cached shuffled distances
v2 <- readRDS('tmp/PCA.eucledian.distances.shuffled.Rds') |>
  dplyr::pull(.data$`NMF:150:PCA:eucledian.dist`)



# F] Figure S2A - Distance real vs shuffled ----

v1.gsam <- gsam.rna.metadata |> 
  dplyr::filter(.data$blacklist.pca == F) |> 
  dplyr::filter(.data$pat.with.IDH == F) |> 
  dplyr::filter(.data$sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F) |> 
  dplyr::filter(.data$tumour.percentage.dna >= 15) |>
  dplyr::filter(!is.na(.data$`NMF:150:PCA:eucledian.dist`)) |> 
  dplyr::filter(.data$resection == "r1") |>  # skip duplicate
  dplyr::pull(.data$`NMF:150:PCA:eucledian.dist`)
stopifnot(length(v1.gsam) == 122)


v1.glass <- glass.gbm.rnaseq.metadata.all.samples |>
  dplyr::filter(.data$resection == "TP") |> 
  dplyr::filter(!is.na(.data$`NMF:150:PCA:eucledian.dist`)) |> 
  dplyr::pull(.data$`NMF:150:PCA:eucledian.dist`)
stopifnot(length(v1.glass) == 87)


v1 <- c(v1.gsam, v1.glass)
k <- length(v1)


l <- k * (k-1)
stopifnot(length(v2) == l)


plt <- rbind(
  data.frame(dist = v1) |> 
    plyr::mutate(dataset = "Actual pairs"),
  data.frame(dist = v2) |> 
    plyr::mutate(dataset = "Shuffled pairs")
)

ggplot(plt, aes(x=dataset, y=dist, size=dataset)) +
  ggbeeswarm::geom_quasirandom(width=0.29) +
  geom_boxplot(outlier.shape = NA, col=alpha("red",1),fill=NA,lwd=0.4, width=0.70,coef=0) +
  ggsignif::geom_signif(comparisons = list(c("Actual pairs" ,  "Shuffled pairs")), test="wilcox.test", col="black", tip_length = 0,y_position=5.4) +
  scale_size_manual(values=c("Actual pairs"=1,"Shuffled pairs"=0.4 ), guide="none") +
  theme_bw() +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x=element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25),
    strip.text = element_text(size = 8) # facet size
  ) +
  labs(y = "Distance GITS", 
       x=NULL, caption=paste0("G-SAM: n=",length(v1.gsam),"  -  GLASS: n=",length(v1.glass)," pairs\nshuffled: n * (n-1) = ",l)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,.15)))

ggsave("output/figures/2022_Figure_S2A.pdf", width=8.3 / ((3/2) * 4),height=8.3 / ((3/2) * 4), scale=2)
ggsave("output/figures/2022_Figure_S2A.png", width=8.3 / ((3/2) * 4),height=1, scale=2)


print(wilcox.test(v1, v2))



