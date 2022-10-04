#!/usr/bin/env R 

# load libs ----


# load data ----

source("scripts/R/chrom_sizes.R")

if(!exists('glass.gbm.rnaseq.metadata.all.samples')) {
  source('scripts/load_glass')
}


# analysis ----



dat <- 'data/gsam/data/GLASS_GBM_R1-R2/variants_gatk_seg.syn31121137.all_samlpes.tsv' |> 
  read.table(sep="\t",header=T,stringsAsFactors = F) |> 
  dplyr::mutate(ROW_ID=NULL, ROW_VERSION=NULL) |> 
  dplyr::mutate(portion_barcode = gsub("^([^\\-]+)-([^\\-]+)-([^\\-]+)-([^\\-]+)-([0-9]+).+$$","\\1-\\2-\\3-\\4-\\5",aliquot_barcode)) |> 
  dplyr::filter(portion_barcode %in% sel$portion_barcode) |> 
  dplyr::mutate(chrom = paste0("chr",chrom) )  |> 
  #dplyr::filter(grepl('WXS', aliquot_barcode)) %>% # WGS are not as clean, and only present for TCGA  
  dplyr::filter(chrom %in% c('chr23', 'chr24') == F) |>  # chr23 is most likely chrX
  dplyr::filter(aliquot_barcode != "GLSS-CU-R017-R1-01D-WXS-C0CBCW") |>  # too noisy data
  dplyr::filter(aliquot_barcode != "GLSS-HF-3118-R1-01D-WXS-QRF6VZ") |>  # too noisy data
  dplyr::filter(aliquot_barcode != "GLSS-HF-3162-R1-01D-WXS-3HVQJ6") |> 
  dplyr::filter(aliquot_barcode != "GLSS-HF-3162-TP-01D-WXS-HTQZ6B") |> 
  dplyr::filter(aliquot_barcode != "GLSS-HF-EE74-R1-01D-WXS-5AS51O") |> 
  dplyr::mutate(aliquot_barcode = as.factor(aliquot_barcode)) 


purities <-   purity <- read.table('output/tables/cnv/tumor.percentage.estimate_glass.2022.all_samples.txt') 


render <- function(bc) {
  #print(bc)
  
  
  dat.pat <- dat %>%
    dplyr::filter(aliquot_barcode == bc) 
  
  
  purity <- purities |> 
    dplyr::filter(.data$aliquot_barcode.wgs == dat.pat$aliquot_barcode[1]) |> 
    dplyr::pull(tumor.purity.cnv.pct.2022) / 100
  
  
  
  center <-  dat.pat
  center <- median(rep(center$log2_copy_ratio, center$num_points))
  
  
  dat.pat <- dat.pat %>%
    dplyr::mutate(outlier = num_points <= 60 ) %>%
    dplyr::filter(outlier == F) %>%
    dplyr::mutate(log2_copy_ratio = log2_copy_ratio - center) %>%
    dplyr::filter(log2_copy_ratio < 1.1) %>%
    dplyr::filter(log2_copy_ratio > -2.0)
  
  
  
  plt <- dat.pat |> 
    dplyr::mutate(id = paste0("id.",1:n())) |> 
    tidyr::pivot_longer(cols=c(`start`,`end`),names_to = 'type', values_to='pos') |>
    dplyr::mutate(chrom = factor(chrom, levels=gtools::mixedsort(unique(as.character(chrom))) )) |> 
    dplyr::mutate(pos = pos / 1000000)
  
  
  
  
  
  ggplot(plt, aes(x=pos, y=log2_copy_ratio, group=id, col=chrom)) +
    facet_grid(cols = vars(chrom), scales = "free", space="free") +
    geom_hline(yintercept = 0, lwd=0.5, lty=3, col="black",alpha=0.5) +
    geom_line(lwd=2) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 4) / 2), lwd=0.7, lty=2, col="black",alpha=0.35) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 3) / 2), lwd=0.7, lty=2, col="black",alpha=0.35) +
    geom_hline(yintercept = log2(((1 - purity) * 2 + purity * 1) / 2), lwd=0.7, lty=2, col="black",alpha=0.35) +
    ylim(-2,2) +
    theme_bw()  +
    theme(
      axis.title = element_text(face = "bold",size = rel(1)),
      #axis.text.x = element_blank(),
      legend.position = 'bottom',
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.45, hjust = 1)
    ) +
    scale_color_discrete(guide="none") +
    labs(x=NULL, y="log2 copy ratio Synapse portal",
         caption = paste0("aliquot_barcode W[X/G]S: ", dat.pat$aliquot_barcode[1], "  -  purity estimate: ", purity )) +
    scale_x_continuous(breaks = c(0,50,100,150,200,250,300))
  
  ggsave(paste0("output/figures/cnv/glass/2022/",dat.pat$aliquot_barcode[1], "_estimate_gg.pdf"), width=8.3, height=3.5, scale=1.75)
}



pbapply::pblapply(unique(dat$aliquot_barcode), render)

