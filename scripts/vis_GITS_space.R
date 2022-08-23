#!/usr/bin/env R

# load libs ----


# load data ----


source('scripts/load_G-SAM_metadata.R')
source('scripts/load_GLASS_data.R')

source('scripts/R/palette.R')


# plots ----

## fig 1c ----



plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
)




# Add initial subtype
tmp.initial.subtype <- plt |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::filter(is.primary) |> 
  dplyr::select(pid, GITS.150.svm.2022.subtype) |> 
  dplyr::rename(subtype.primary = GITS.150.svm.2022.subtype)

plt <- plt |> 
  dplyr::left_join(tmp.initial.subtype, by=c('pid'='pid'), suffix=c('',''))
rm(tmp.initial.subtype)


# Add primary's as xend and yend to the pairs
tmp.segment <- plt |>
  dplyr::mutate(subtype.primary = NULL) |>
  dplyr::mutate(is.primary = ifelse(is.primary,"primary","recurrence")) |> 
  dplyr::group_by(pid) |>
  dplyr::filter(n() == 2) |>
  dplyr::ungroup() |>
  tidyr::pivot_wider(names_from = is.primary,
                     values_from = c(sid, `NMF:150:PC1`, `NMF:150:PC2`, `GITS.150.svm.2022.subtype`)) |>
  dplyr::mutate(
    `NMF:150:PC1_primary`  = NULL,
    `NMF:150:PC2_primary` = NULL,
    `GITS.150.svm.2022.subtype_primary` = NULL,
    `sid_recurrence` = NULL,
    `pid` = NULL
  )

plt <- plt |> 
  dplyr::left_join(tmp.segment, by=c('sid'='sid_primary'), suffix=c('',''))
rm(tmp.segment)


# plt <- plt |> 
#   dplyr::filter(dataset == "GLASS")

  
# make facets
plt.expanded <- rbind(
  plt |> dplyr::mutate(facet = "Classical") |>
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Mesenchymal")  |> 
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Proneural")  |> 
    dplyr::mutate(highlight = subtype.primary == facet)
)



plt.contours <- readRDS("cache/analysis_GITS_space_GITS_contours.Rds") |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype` = class, `pid` = NA)



ggplot(plt.expanded, aes(x=-`NMF:150:PC1`,y=-`NMF:150:PC2`, group=pid, col =`GITS.150.svm.2022.subtype`,fill =`GITS.150.svm.2022.subtype`)) +
  geom_raster(data = plt.contours, alpha=0.05) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)), 
               colour="gray40", 
               size=0.25, 
               lty=2,
               breaks=c(1.75,2.25)) +
  facet_grid(cols = vars(facet)) +
  geom_point(data = plt.expanded |> dplyr::filter(highlight == F) , size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Classical') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Classical'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Mesenchymal') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Mesenchymal'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Proneural') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Proneural'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt.expanded |>  dplyr::filter(highlight == T & is.primary), aes(fill = `GITS.150.svm.2022.subtype`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors) +
  labs(x = "PC1 on NMF meta-features", y = "PC2 on NMF meta-features", fill = "Subtype") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size = rel(1))) +
  theme(legend.position = 'bottom') +
  theme(text = element_text(family = 'Arial'))

rm(plt, plt.expanded)



### MES ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM") |> 
    dplyr::mutate(batch = dataset)
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      GITS.150.svm.2022.subtype,
      aliquot_batch_synapse
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS") |> 
    dplyr::rename(batch = aliquot_batch_synapse)
)




# Add initial subtype
tmp.initial.subtype <- plt |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::filter(is.primary) |> 
  dplyr::select(pid, GITS.150.svm.2022.subtype) |> 
  dplyr::rename(subtype.primary = GITS.150.svm.2022.subtype)

plt <- plt |> 
  dplyr::left_join(tmp.initial.subtype, by=c('pid'='pid'), suffix=c('',''))
rm(tmp.initial.subtype)


# Add primary's as xend and yend to the pairs
tmp.segment <- plt |>
  dplyr::mutate(batch = NULL) |> 
  dplyr::mutate(subtype.primary = NULL) |>
  dplyr::mutate(is.primary = ifelse(is.primary,"primary","recurrence")) |> 
  dplyr::group_by(pid) |>
  dplyr::filter(n() == 2) |>
  dplyr::ungroup() |>
  tidyr::pivot_wider(names_from = is.primary,
                     values_from = c(sid, `NMF:150:PC1`, `NMF:150:PC2`, `GITS.150.svm.2022.subtype`)) |>
  dplyr::mutate(
    `NMF:150:PC1_primary`  = NULL,
    `NMF:150:PC2_primary` = NULL,
    `GITS.150.svm.2022.subtype_primary` = NULL,
    `sid_recurrence` = NULL,
    `pid` = NULL
  )

plt <- plt |> 
  dplyr::left_join(tmp.segment, by=c('sid'='sid_primary'), suffix=c('',''))
rm(tmp.segment)


# plt <- plt |> 
#   dplyr::filter(dataset == "GLASS")


# make facets
plt.expanded <- rbind(
  plt |> dplyr::mutate(facet = "Classical") |>
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Mesenchymal")  |> 
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Proneural")  |> 
    dplyr::mutate(highlight = subtype.primary == facet)
) |> 
  dplyr::filter(facet ==  "Mesenchymal")



plt.contours <- readRDS("cache/analysis_GITS_space_GITS_contours.Rds") |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype` = class, `pid` = NA, batch=NA)



ggplot(plt.expanded, aes(x=-`NMF:150:PC1`,y=-`NMF:150:PC2`, group=pid, col =`batch`,
                         fill =`batch`,
                         label = batch)) +
  geom_raster(data = plt.contours, alpha=0.05) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)), 
               colour="gray40", 
               size=0.25, 
               lty=2,
               breaks=c(1.75,2.25)) +
  facet_grid(cols = vars(facet)) +
  geom_point(data = plt.expanded |> dplyr::filter(highlight == F) , size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Classical') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Classical'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  #geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Mesenchymal') , # fill.arrow is no aesthetic (yet?)
  #             aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Mesenchymal'], col=rgb(0,0,0,0.6),
  #             lwd = 0.35,
  #             arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Proneural') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Proneural'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt.expanded |>  dplyr::filter(highlight == T & is.primary), aes(fill = `GITS.150.svm.2022.subtype`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  geom_text(data = plt.expanded |> 
                             dplyr::filter(is.primary == F &
                                           subtype.primary == "Mesenchymal" &
                                             GITS.150.svm.2022.subtype != "Mesenchymal"
                                             )
                             ) +
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors) +
  labs(x = "PC1 on NMF meta-features", y = "PC2 on NMF meta-features", fill = "Subtype") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size = rel(1))) +
  theme(legend.position = 'bottom')

rm(plt, plt.expanded)


### 3050 ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::filter(grepl("3050",pid))




# Add initial subtype
tmp.initial.subtype <- plt |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::filter(is.primary) |> 
  dplyr::select(pid, GITS.150.svm.2022.subtype) |> 
  dplyr::rename(subtype.primary = GITS.150.svm.2022.subtype)

plt <- plt |> 
  dplyr::left_join(tmp.initial.subtype, by=c('pid'='pid'), suffix=c('',''))
rm(tmp.initial.subtype)


# Add primary's as xend and yend to the pairs
tmp.segment <- plt |>
  dplyr::mutate(subtype.primary = NULL) |>
  dplyr::mutate(is.primary = ifelse(is.primary,"primary","recurrence")) |> 
  dplyr::group_by(pid) |>
  dplyr::filter(n() == 2) |>
  dplyr::ungroup() |>
  tidyr::pivot_wider(names_from = is.primary,
                     values_from = c(sid, `NMF:150:PC1`, `NMF:150:PC2`, `GITS.150.svm.2022.subtype`)) |>
  dplyr::mutate(
    `NMF:150:PC1_primary`  = NULL,
    `NMF:150:PC2_primary` = NULL,
    `GITS.150.svm.2022.subtype_primary` = NULL,
    `sid_recurrence` = NULL,
    `pid` = NULL
  )

plt <- plt |> 
  dplyr::left_join(tmp.segment, by=c('sid'='sid_primary'), suffix=c('',''))
rm(tmp.segment)


# plt <- plt |> 
#   dplyr::filter(dataset == "GLASS")


# make facets
plt.expanded <- rbind(
  plt |> dplyr::mutate(facet = "Classical") |>
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Mesenchymal")  |> 
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Proneural")  |> 
    dplyr::mutate(highlight = subtype.primary == facet)
)



plt.contours <- readRDS("cache/analysis_GITS_space_GITS_contours.Rds") |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype` = class, `pid` = NA)



ggplot(plt.expanded, aes(x=-`NMF:150:PC1`,y=-`NMF:150:PC2`, group=pid, col =`GITS.150.svm.2022.subtype`,fill =`GITS.150.svm.2022.subtype`)) +
  geom_raster(data = plt.contours, alpha=0.05) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)), 
               colour="gray40", 
               size=0.25, 
               lty=2,
               breaks=c(1.75,2.25)) +
  facet_grid(cols = vars(facet)) +
  geom_point(data = plt.expanded |> dplyr::filter(highlight == F) , size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Classical') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Classical'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Mesenchymal') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Mesenchymal'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data= plt.expanded |> dplyr::filter(highlight == T & `GITS.150.svm.2022.subtype_recurrence` == 'Proneural') , # fill.arrow is no aesthetic (yet?)
               aes(xend = -`NMF:150:PC1_recurrence`, yend=-`NMF:150:PC2_recurrence`), arrow.fill = subtype_colors['Proneural'], col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt.expanded |>  dplyr::filter(highlight == T & is.primary), aes(fill = `GITS.150.svm.2022.subtype`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors) +
  labs(x = "PC1 on NMF meta-features", y = "PC2 on NMF meta-features", fill = "Subtype") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size = rel(1))) +
  theme(legend.position = 'bottom') +
  theme(text = element_text(family = 'Arial'))

rm(plt, plt.expanded)





## fig s1a ----



plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      `sid`,
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      `NMF:150:membership`,
      `GITS.150.svm.2022.subtype`
    )
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      `aliquot_barcode`,
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      `NMF:150:membership`,
      `GITS.150.svm.2022.subtype`
    ) |>
    dplyr::rename(sid = aliquot_barcode)) |> 
  
  dplyr::mutate(max.1 = pmax(`NMF:150:2`, `NMF:150:3`) )  |> 
  dplyr::mutate(max.2 = pmax(`NMF:150:1`, `NMF:150:3`) ) |> 
  dplyr::mutate(max.3 = pmax(`NMF:150:1`, `NMF:150:2`) ) |> 
  
  dplyr::mutate(diff.1 = ifelse(`NMF:150:1` > max.1, `NMF:150:1` , 0) ) |>
  dplyr::mutate(diff.2 = ifelse(`NMF:150:2` > max.2, `NMF:150:2`, 0) ) |>
  dplyr::mutate(diff.3 = ifelse(`NMF:150:3` > max.3, `NMF:150:3` , 0) ) |>
  
  dplyr::mutate(order =  order(order(-diff.1, -diff.2, -diff.3))) |>
  dplyr::mutate(max.1 = NULL, max.2 = NULL, max.3 = NULL, diff.1 = NULL, diff.2 = NULL, diff.3 = NULL) |>
  dplyr::arrange(order) |> 
  
  tidyr::pivot_longer(cols=c(`NMF:150:1`, `NMF:150:2`, `NMF:150:3`),names_to = 'NMF H-matrix meta-feature') |> 
  dplyr::mutate(`NMF H-matrix meta-feature` = gsub("NMF:150:","NMF meta-feature ",`NMF H-matrix meta-feature`)) 



plt.line.based <- rbind(
  plt |>  dplyr::mutate(value=0),
  plt
)



ggplot(plt.line.based, aes(x=reorder(sid, order), xend=reorder(sid, order), y =value, yend = 0,
                           #fill=`GITS.150.svm.2022.subtype`,
                           col=`GITS.150.svm.2022.subtype`
                           )) +
  #geom_point()  +
  #geom_segment() +
  #geom_bar(stat="identity",lwd=0.2) +
  geom_line(aes(group=sid),lwd=0.425) +
  facet_grid(rows = vars(`NMF H-matrix meta-feature`)) +
  scale_color_manual(values = mixcol(subtype_colors, c("black","black","black"),0.15) ) + 
  labs(col = 'ssGSEA subtype', y = "NMF H-matrix values", x = NULL) + 
  theme_bw() +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x=element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25),
    strip.text = element_text(size = 8) # facet size
  ) 

ggsave("output/figures/2022_figure_S1A.pdf", width=8.3 / 2,height=8.3/4, scale=2)




## fig s1b ----


tmp.pca <- readRDS('tmp/analysis_GITS_space_GITS_PCA_150.Rds')


plt <- data.frame(PC = colnames(tmp.pca$x),
                  var_explained=(tmp.pca$sdev)^2/sum((tmp.pca$sdev)^2)) %>%
  dplyr::mutate(label = paste0(round(var_explained * 100,1),"%"))

ggplot(plt, aes(x=PC,y=var_explained, group=1, label=label)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_col(fill="gray60", col="black", lwd=0.4) +
  geom_line()+
  geom_label(y = 0.35, size=4) +
  geom_point(size=4)+
  labs(title="Scree plot: PCA on NMF H-matrix", x= NULL, y="Variance by PC explained") +
  theme_bw()  +
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
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) 

ggsave("output/figures/2022_figure_S1B.pdf", width=8.3 / 4,height=8.3/4, scale=2)


rm(tmp.pca)


## fig s1c ----


tmp.pca <- readRDS('tmp/analysis_GITS_space_GITS_PCA_150.Rds')


# PC being a prcomp object
data <- data.frame(obsnames=row.names(tmp.pca$x), tmp.pca$x) |> 
  dplyr::mutate(PC1 = -PC1, PC2 = -PC2) # same orientation as first rivision
datapc <- data.frame(varnames=rownames(tmp.pca$rotation), tmp.pca$rotation) |> 
  dplyr::mutate(PC1 = -PC1, PC2 = -PC2) |>  # same orientation as first rivision
  dplyr::mutate(varnames = gsub('NMF:150:','NMF meta-feature ', varnames))
  

x = "PC1"
y = "PC2"
mult <- min((max(data[, y]) - min(data[, y]) / (max(datapc[, y]) - min(datapc[, y]))),
            (max(data[, x]) - min(data[, x]) / (max(datapc[, x]) - min(datapc[, x]))))
datapc <- transform(datapc,
                    v1 = .7 * mult * (get(x)),
                    v2 = .7 * mult * (get(y)))


ggplot(data, aes(x=PC1, y=PC2, col=varnames)) +
  geom_point(size=3,  col='black', fill='gray60', pch=21, alpha=0.85) +
  coord_equal() +
  geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.4,"cm")), lwd=1.75, alpha=0.9) +
  ggrepel::geom_label_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3,
                            vjust=0.6, show.legend=F, col="black") +
  labs(col = "Associated with") +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill = "Sub-type (GlioVis)") +
  scale_color_manual(values = mixcol(subtype_colors_nmf2022, rep("black",length(subtype_colors_nmf2022)),0.15),
                     labels = c('NMF meta-feature 2'='PN',
                                                              'NMF meta-feature 1'='CL',
                                                              'NMF meta-feature 3'='MES')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) 

ggsave("output/figures/2022_figure_S1C.pdf", width=8.3 / 4,height=8.3/4, scale=2)






## fig s1d ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
)




ggplot(plt, aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=ssGSEA.2022.subtype)) +
  geom_point(size=3,  col='black', pch=21, alpha=0.85) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill = "Subtype (ssGSEA)") +
  scale_fill_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.15),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2))

ggsave("output/figures/2022_figure_S1D.pdf", width=8.3 / 4,height=8.3/4, scale=2)







## fig s1e ----



plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
)


plt.contours <- readRDS("cache/analysis_GITS_space_GITS_contours.Rds") |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = class, `pid` = NA)


ggplot(plt, aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=ssGSEA.2022.subtype, group=pid)) +
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)), 
               colour="gray40", 
               size=0.25, 
               lty=2,
               breaks=c(1.5,2.5)) +
  geom_point(size=3,  col='black', pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill = "Subtype (ssGSEA)") +
  scale_fill_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.1),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2))

ggsave("output/figures/2022_figure_S1E.pdf", width=8.3 / 4,height=8.3/4, scale=2)






## fig s1f ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      is.primary,
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::mutate(
    misclass = case_when(
      grepl("|",ssGSEA.2022.subtype,fixed=T) ~ "?",
      ssGSEA.2022.subtype != GITS.150.svm.2022.subtype ~ ".",
      T ~ " "
    )
  )
  #dplyr::filter(misclass != " ")


plt.contours <- readRDS("cache/analysis_GITS_space_GITS_contours.Rds") |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = class, `pid` = NA, misclass=".")


ggplot(plt, aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=ssGSEA.2022.subtype, group=pid, shape=misclass)) +
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)), 
               colour="gray40", 
               size=0.25, 
               lty=2,
               breaks=c(1.5,2.5)) +
  geom_point(data = plt |>  dplyr::filter(misclass == " ") ,size=3,  col='black', pch=21, alpha=0.7) +
  geom_point(data = plt |>  dplyr::filter(misclass == "?") ,size=3,  col='black', pch=21, alpha=0.7) +
  geom_point(data = plt |>  dplyr::filter(misclass == ".") ,size=3,  col='black', pch=21, alpha=0.7) +
  geom_point(data = plt |>  dplyr::filter(misclass == "?"),
             col = 'black',  size=1.2) + # pch=4,
  geom_point(data = plt |>  dplyr::filter(misclass == "."), col = 'black',fill="white",pch=21, size=1.4,stroke=0.65) + # pch=19, 
  coord_equal() +
  labs(x="PC1 on NMF meta-features",
       y="PC2 on NMF meta-features", 
       fill = "Subtype (ssGSEA)", 
       shape="") +
  scale_fill_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.1),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  scale_shape_manual(values=c("?"=4,"."=19),
                     label=c('?'='ssGSEA undecisive','.'='GITS ~ ssGSEA discordant')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill="none") #guide_legend(ncol=2)

ggsave("output/figures/2022_figure_S1F.pdf", width=8.3 / 4,height=8.3/4, scale=2)




## fig s1g ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,

      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      
      ssGSEA.2022.Proneural.enrichment_score,
      ssGSEA.2022.Mesenchymal.enrichment_score,
      ssGSEA.2022.Classical.enrichment_score,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,
      
      
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      
      ssGSEA.2022.Proneural.enrichment_score,
      ssGSEA.2022.Mesenchymal.enrichment_score,
      ssGSEA.2022.Classical.enrichment_score,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
  tidyr::pivot_longer(cols=c(`NMF:150:1`, `NMF:150:2`, `NMF:150:3`),names_to="NMF meta-feature",values_to="NMF contribution") |> 
  tidyr::pivot_longer(cols=c(`ssGSEA.2022.Proneural.enrichment_score`,`ssGSEA.2022.Mesenchymal.enrichment_score`,`ssGSEA.2022.Classical.enrichment_score`),
                      names_to="ssGSEA subtype", values_to="ssGSEA enrichment score") |> 
  dplyr::mutate(`ssGSEA subtype` = case_when(
    `ssGSEA subtype` == "ssGSEA.2022.Proneural.enrichment_score" ~ "PN",
    `ssGSEA subtype` == "ssGSEA.2022.Mesenchymal.enrichment_score" ~ "MES",
    `ssGSEA subtype` == "ssGSEA.2022.Classical.enrichment_score" ~ "CL"
  )) |> 
  dplyr::filter((`NMF meta-feature` == "NMF:150:1" &  `ssGSEA subtype` == "CL") |
                (  `NMF meta-feature` == "NMF:150:2" &  `ssGSEA subtype` == "PN") |
                 ( `NMF meta-feature` == "NMF:150:3" &  `ssGSEA subtype` == "MES")) |> 
  dplyr::mutate(facet = paste0(`NMF meta-feature`," ~ " , `ssGSEA subtype`)) |> 
  dplyr::mutate(needle = case_when(`ssGSEA subtype` == "CL" ~ "Classical",
                                     `ssGSEA subtype` == "MES" ~ "Mesenchymal",
                                     `ssGSEA subtype` == "PN" ~ "Proneural"
                                     )) |> 
  dplyr::mutate(facet.target = stringr::str_detect(ssGSEA.2022.subtype, `needle`)) |> 
  dplyr::mutate(`needle` = NULL) |> 
  dplyr::mutate(
    misclass = case_when(
      facet.target & grepl("|",ssGSEA.2022.subtype,fixed=T) ~ " ",
      facet.target & ssGSEA.2022.subtype != GITS.150.svm.2022.subtype ~ ".",
      T ~ " "
    )
  ) 




ggplot(plt, aes(x=`ssGSEA enrichment score`, y=`NMF contribution`, fill=`ssGSEA.2022.subtype`, shape=misclass)) +
  facet_grid(cols = vars(`facet`), scales = "free") +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(size=3, data=plt |>  dplyr::filter(misclass != "."), col='black', pch=21, alpha=0.7) +
  geom_point(size=3, data=plt |>   dplyr::filter(misclass == ".") ,col='black', pch=21, alpha=0.7) +
  #geom_point(data = plt |>  dplyr::filter(misclass == "?"), col = 'white',  size=1.2) + # pch=4,
  geom_point(data = plt |>  dplyr::filter(misclass == "."), col = 'black',fill="white",pch=21, size=0.9,stroke=0.65) +
  labs(y="NMF meta-feature score", 
       fill = "Subtype (ssGSEA)",
       shape="") +
  scale_fill_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.1),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  scale_shape_manual(values=c("?"=4,"."=19),
                     label=c('?'='ssGSEA undecisive','.'='GITS ~ ssGSEA discordant')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides() #guide_legend(ncol=2)

ggsave("output/figures/2022_figure_S1G.pdf", width=8.3 / 2,height=8.3/5, scale=2)



## fig s1h ----



plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::select(
      sid,
      pid,
      
      `NMF:150:PC1`,
      `NMF:150:PC2`,
      
      tumour.percentage.dna,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    )  |> 
    dplyr::rename(purity = tumour.percentage.dna) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
 glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      
      `NMF:150:PC1`,
      `NMF:150:PC2`,

      tumour.percentage.2022,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::rename(purity = tumour.percentage.2022) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::mutate(purity.binned = case_when(purity < 33 ~ "<33",
                                          purity >= 33 & purity < 66 ~ "33 - 66",
                                          T ~ ">= 66" )) |> 
  dplyr::mutate(col = purity.binned )



plt.contours <- readRDS("cache/analysis_GITS_space_GITS_contours.Rds") |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = class, `pid` = NA, misclass=".", col = class)


ggplot(plt, aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=col, group=pid)) +
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)), 
               colour="gray40", 
               size=0.25, 
               lty=2,
               breaks=c(1.5,2.5))  +
  geom_point(size=2,  col='black', pch=21, alpha=0.7,stroke=0.4) +
  coord_equal() +
  
  scale_fill_manual(values = c(mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.1),
                               "<33" = mixcol("#55864e","white",0.66),
                               "33 - 66" = mixcol("#55864e","white",0.33),
                               ">= 66"   = mixcol("#55864e","white",0)
                               ),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  labs(x="PC1 on NMF meta-features",
       y="PC2 on NMF meta-features", 
       fill = "", 
       shape="") +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )

ggsave("output/figures/2022_figure_S1H.pdf", width=8.3 / 4,height=8.3/4, scale=2)


## fig s1k ----


tmp <- rbind(
    gsam.rna.metadata |>
      
      dplyr::filter(blacklist.pca == F) %>%
      dplyr::filter(pat.with.IDH == F) %>%
      dplyr::filter(
        sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
      ) %>%
      dplyr::filter(tumour.percentage.dna >= 15) |>
      
      dplyr::select(
        sid,
        pid,
        
        ssGSEA.2022.Proneural.enrichment_score,
        ssGSEA.2022.Mesenchymal.enrichment_score,
        ssGSEA.2022.Classical.enrichment_score,
        
        ssGSEA.2022.subtype,
        GITS.150.svm.2022.subtype
      )  |> 
      dplyr::mutate(dataset = "G-SAM")
    ,
    glass.gbm.rnaseq.metadata.all.samples |>
      dplyr::select(
        aliquot_barcode,
        case_barcode,
        
        ssGSEA.2022.Proneural.enrichment_score,
        ssGSEA.2022.Mesenchymal.enrichment_score,
        ssGSEA.2022.Classical.enrichment_score,
        
        ssGSEA.2022.subtype,
        GITS.150.svm.2022.subtype
      ) |>
      dplyr::rename(sid = aliquot_barcode) |>
      dplyr::rename(pid = case_barcode) |> 
      dplyr::mutate(dataset = "GLASS")
  )

tmp.pca <- tmp |> 
  tibble::column_to_rownames('sid') |>
  dplyr::select(`ssGSEA.2022.Proneural.enrichment_score`, `ssGSEA.2022.Mesenchymal.enrichment_score`, `ssGSEA.2022.Classical.enrichment_score`) |>
  dplyr::rename(`ssGSEA PN score`= `ssGSEA.2022.Proneural.enrichment_score`) |> 
  dplyr::rename(`ssGSEA MES score`= `ssGSEA.2022.Mesenchymal.enrichment_score`) |> 
  dplyr::rename(`ssGSEA CL score` = `ssGSEA.2022.Classical.enrichment_score`)  |> 
  prcomp()  


# PC being a prcomp object
data <- data.frame(obsnames=row.names(tmp.pca$x), tmp.pca$x) |> 
  dplyr::mutate(PC2 = -PC2) |>  # same orientation as first rivision
  dplyr::left_join(tmp |> dplyr::select(sid, ssGSEA.2022.subtype), by=c('obsnames'='sid'),suffix=c('',''))
datapc <- data.frame(varnames=rownames(tmp.pca$rotation), tmp.pca$rotation) |> 
  dplyr::mutate(PC2 = -PC2) |>  # same orientation as first rivision
  dplyr::mutate(varnames = gsub('NMF:150:','NMF meta-feature ', varnames))


x = "PC1"
y = "PC2"
mult <- min((max(data[, y]) - min(data[, y]) / (max(datapc[, y]) - min(datapc[, y]))),
            (max(data[, x]) - min(data[, x]) / (max(datapc[, x]) - min(datapc[, x]))))
datapc <- transform(datapc,
                    v1 = .7 * mult * (get(x)),
                    v2 = .7 * mult * (get(y)))


ggplot(data, aes(x=PC1, y=PC2, col=varnames)) +
  geom_point(size=3,  col='black', fill='gray60', pch=21, alpha=0.85) +
  coord_equal() +
  geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.4,"cm")), lwd=1.75, alpha=0.9) +
  ggrepel::geom_label_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3,
                            vjust=0.6, show.legend=F, col="black") +
  labs(col = "Associated with") +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill = "Sub-type (GlioVis)") +
  scale_color_manual(values = mixcol(subtype_colors_ssGSEA, rep("black",length(subtype_colors_ssGSEA)),0.15),
                     labels = c('ssGSEA PN score'='PN',
                                'ssGSEA CL score'='CL',
                                'ssGSEA MES score'='MES')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) 

ggsave("output/figures/2022_figure_S1K.pdf", width=8.3 / 4,height=8.3/4, scale=2)


## fig s1l ----

# probs x class & misclass



plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::select(
      sid,
      
      ssGSEA.2022.Classical_pval,
      ssGSEA.2022.Mesenchymal_pval,
      ssGSEA.2022.Proneural_pval,

      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      aliquot_barcode,

      ssGSEA.2022.Classical_pval,
      ssGSEA.2022.Mesenchymal_pval,
      ssGSEA.2022.Proneural_pval,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::mutate(ssGSEA.err1 = case_when(
    ssGSEA.2022.subtype == "Classical" ~ ssGSEA.2022.Classical_pval^2 ,
    ssGSEA.2022.subtype == "Proneural" ~ ssGSEA.2022.Proneural_pval^2 ,
    ssGSEA.2022.subtype == "Mesenchymal" ~ ssGSEA.2022.Mesenchymal_pval^2,
    ssGSEA.2022.subtype == "Proneural|Classical" ~ ssGSEA.2022.Classical_pval * ssGSEA.2022.Proneural_pval
  )) |> 
  dplyr::mutate(ssGSEA.err2 = case_when(
    ssGSEA.2022.subtype == "Classical" ~ sqrt((1 - ssGSEA.2022.Proneural_pval)^2 + (1 - ssGSEA.2022.Mesenchymal_pval)^2) ,
    ssGSEA.2022.subtype == "Proneural" ~ sqrt((1 - ssGSEA.2022.Classical_pval)^2 + (1 - ssGSEA.2022.Mesenchymal_pval)^2) ,
    ssGSEA.2022.subtype == "Mesenchymal" ~ sqrt((1 - ssGSEA.2022.Classical_pval)^2 + (1 - ssGSEA.2022.Proneural_pval)^2) ,
    ssGSEA.2022.subtype == "Proneural|Classical" ~ ssGSEA.2022.Mesenchymal_pval
  )) |> 
  tidyr::pivot_longer(cols=c( ssGSEA.2022.Classical_pval,
                              ssGSEA.2022.Mesenchymal_pval,
                              ssGSEA.2022.Proneural_pval), names_to = "ssGSEA subtype", values_to = "ssGSEA pval") |> 
  dplyr::mutate(`ssGSEA subtype` = gsub("ssGSEA.2022.","",`ssGSEA subtype`)) |> 
  dplyr::mutate(`ssGSEA subtype` = gsub("_pval","",`ssGSEA subtype`)) |> 
  dplyr::mutate(order = order(order(ssGSEA.err1, ssGSEA.err2)))



plt.line.based <- rbind(
  plt |>  dplyr::mutate(`ssGSEA pval`=0),
  plt |>  dplyr::mutate(`ssGSEA pval` = -log(`ssGSEA pval`))
)


plt.line.based |> dplyr::filter(grepl("|",`GITS.150.svm.2022.subtype`))


ggplot(plt.line.based, aes(x=reorder(sid, order), y=`ssGSEA pval`, group=sid, col=`GITS.150.svm.2022.subtype`)) +
  #geom_point(size=3,  col='black', pch=21, alpha=0.85) +
  geom_line(lwd=0.8) +
  facet_grid(cols = vars(ssGSEA.2022.subtype), rows=vars(`ssGSEA subtype`), scale="free_x", space="free_x") +
  scale_color_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.15),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2))

ggsave("output/figures/2022_figure_S1L.pdf", width=8.3 / 4,height=8.3/4, scale=2)



## fig s2 ----


plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::mutate(is.primary = ifelse(resection == "r1","primary","recurrent")) |> 
    
    dplyr::select(
      `sid`,
      `pid`,
      `is.primary`,
      
      `NMF:150:PCA:eucledian.dist`,
      
      `GITS.150.svm.2022.subtype`
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = ifelse(resection == "TP","primary","recurrent")) |> 
    dplyr::select(
      `aliquot_barcode`,
      `case_barcode`,
      `is.primary`,
      
      `NMF:150:PCA:eucledian.dist`,
      
      `GITS.150.svm.2022.subtype`
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::group_by(pid) |> # filter for matching pairs
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |>
  dplyr::mutate(`GITS.150.svm.2022.subtype` = case_when(
    `GITS.150.svm.2022.subtype` == "Proneural" ~ "PN",
    `GITS.150.svm.2022.subtype` == "Classical" ~ "CL",
    `GITS.150.svm.2022.subtype` == "Mesenchymal" ~ "MES"
  )) |> 
  tidyr::pivot_wider(id_cols = pid, 
                     names_from = c(is.primary),
                     values_from = c(sid,`GITS.150.svm.2022.subtype`, `NMF:150:PCA:eucledian.dist`, `dataset`)) |> 
  dplyr::mutate(stable = ifelse(`GITS.150.svm.2022.subtype_primary` == `GITS.150.svm.2022.subtype_recurrent`, "stable","transition"))

stopifnot(plt$`NMF:150:PCA:eucledian.dist_primary` == plt$`NMF:150:PCA:eucledian.dist_recurrent`)


plt <- plt |>
  dplyr::rename(`NMF:150:PCA:eucledian.dist` = `NMF:150:PCA:eucledian.dist_primary`) |> 
  dplyr::mutate(`NMF:150:PCA:eucledian.dist_recurrent` = NULL) |> 
  dplyr::rename(`dataset` = `dataset_primary`) |> 
  dplyr::mutate(`dataset_recurrent` = NULL) |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype_primary` = paste0("from: ", `GITS.150.svm.2022.subtype_primary`)) |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype_recurrent` = paste0("to: ", `GITS.150.svm.2022.subtype_recurrent`))


#, group = type
ggplot(plt, aes(x=GITS.150.svm.2022.subtype_recurrent, y=`NMF:150:PCA:eucledian.dist`)) +
  facet_grid(rows = vars(GITS.150.svm.2022.subtype_primary), space="free_x") +
  #ylim(0, 4.3) +
  geom_violin(width=1.05,aes(fill=stable)) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  ggbeeswarm::geom_quasirandom(aes(fill=dataset), pch=21,size=3, cex=4) +
  labs(x = NULL, y = "Euclidean distance GITS space", fill = "") +
  scale_fill_manual(values = c('stable'='gray90', 'transition'='white',dataset_colors['G-SAM'], dataset_colors['GLASS'] )) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    #strip.text = element_text(size = 7),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2))



ggsave("output/figures/2022_figure_S2.pdf", width=8.3 / 4,height=8.3/4, scale=2)


rm(plt)



## fig s3a ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::mutate(is.primary = ifelse(resection == "r1","primary","recurrent")) |> 
    
    dplyr::select(
      `sid`,
      `pid`,
      `is.primary`,
      
      `NMF:150:PC1.n`,
      `NMF:150:PC2.n`,
      
      `GITS.150.svm.2022.subtype`
    ) 
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::mutate(is.primary = ifelse(resection == "TP","primary","recurrent")) |> 
    dplyr::select(
      `aliquot_barcode`,
      `case_barcode`,
      `is.primary`,
      
      `NMF:150:PC1.n`,
      `NMF:150:PC2.n`,
      
      `GITS.150.svm.2022.subtype`
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) 
) |> 
  dplyr::group_by(pid) |> # filter for matching pairs
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |>
  tidyr::pivot_wider(id_cols = pid, 
                     names_from = c(is.primary),
                     values_from = c(sid,`GITS.150.svm.2022.subtype`,  `NMF:150:PC1.n`, `NMF:150:PC2.n`)) |> 
  dplyr::mutate(`NMF:150:PC1.n_recurrent` = `NMF:150:PC1.n_recurrent` - `NMF:150:PC1.n_primary`)|> 
  dplyr::mutate(`NMF:150:PC2.n_recurrent` = `NMF:150:PC2.n_recurrent` - `NMF:150:PC2.n_primary`)|> 
  dplyr::mutate(`NMF:150:PC1.n_primary` = `NMF:150:PC1.n_primary` - `NMF:150:PC1.n_primary`)|> 
  dplyr::mutate(`NMF:150:PC2.n_primary` = `NMF:150:PC2.n_primary` - `NMF:150:PC2.n_primary`) |> 
  dplyr::mutate(subtype.primary = `GITS.150.svm.2022.subtype_primary`) |> 
  tidyr::pivot_longer(
    cols = c(
      `sid_primary`,
      `sid_recurrent`,
      
      `GITS.150.svm.2022.subtype_primary`,
      `GITS.150.svm.2022.subtype_recurrent`,
      
      `NMF:150:PC1.n_primary`,
      `NMF:150:PC1.n_recurrent`,
      
      `NMF:150:PC2.n_primary`,
      `NMF:150:PC2.n_recurrent`
    ),
    names_to = c('.value','resection'),
    names_pattern = "(.+)_(.+)"
  )


m.r <- plt %>%
  dplyr::filter(resection == "recurrent") |> 
  dplyr::group_by(subtype.primary) %>%
  dplyr::summarise(`NMF:150:PC1.n` = mean(`NMF:150:PC1.n`), 
                   `NMF:150:PC2.n` = mean(`NMF:150:PC2.n`)) |> 
  dplyr::ungroup() %>%
  dplyr::mutate(pid = "Mean",
            sid = "Mean",
            `GITS.150.svm.2022.subtype` = "Mean",
            resection = "recurrent")

m.p <- m.r |> dplyr::mutate(`NMF:150:PC1.n` = 0 , `NMF:150:PC2.n` = 0 , resection = "primary")


plt <- rbind(plt, m.p, m.r)
rm(m.p, m.r)



ggplot(plt, aes(x = -`NMF:150:PC1.n`, y= -`NMF:150:PC2.n`, group=pid, col=subtype.primary, fill=GITS.150.svm.2022.subtype)) +
  facet_grid(cols = vars(subtype.primary)) +
  geom_line(data = plt |> dplyr::filter(pid != "Mean") ,alpha=0.75, lwd=1) +
  geom_point(data = plt |> dplyr::filter(resection == "recurrent" & sid != "Mean"),size=3, pch=21,col="black",stroke=0.8)  +
  geom_line(data = plt |> dplyr::filter(pid == "Mean") , lwd=1.2, col="black") +
  geom_point(data = plt |> dplyr::filter(resection == "recurrent" & sid == "Mean"),size=3, pch=21,col="black",fill="white",stroke=0.8)  +
  scale_color_manual(values = subtype_colors, 
                     guide = "none") +
  scale_fill_manual(values = subtype_colors) +
  scale_x_continuous(breaks = c(0)) +
  scale_y_continuous(breaks = c(0)) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill="Subtype") +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    legend.position = 'bottom',
    #panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #strip.text = element_text(size = 7),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=3))



ggsave("output/figures/2022_figure_S3a.pdf", width=8.3 / 2,height=8.3/4, scale=2)


rm(plt)



## fig s3b ----


plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::select(
      `pid`,
      `NMF:150:PCA:eucledian.dist`,
    ) |> 
    dplyr::left_join(gsam.patient.metadata |> 
                       dplyr::select(`studyID`, `HM`),
                     by=c('pid'='studyID'),suffix=c('','')) |> 
    dplyr::mutate(HM = case_when(
      HM == "No" ~ F,
      HM == "Yes" ~ T, 
      T ~ as.logical(NA)
    ))
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::select(
      `case_barcode`,
      `NMF:150:PCA:eucledian.dist`,
      `HM`
    ) |>
    dplyr::rename(pid = case_barcode) 
) |> 
  dplyr::filter(!is.na(`NMF:150:PCA:eucledian.dist`)) |> 
  dplyr::distinct() |> 
  dplyr::mutate(pid = as.character(pid)) |> 
  dplyr::mutate(facet.HM.determinated = is.na(HM))


# load transition segments 
plt.expanded <- readRDS('tmp/analysis_GITS_space.transition.segments.Rds') |> 
  dplyr::left_join(plt, by=c('pid'='pid')) |> 
  dplyr::mutate(x = pct_transition * `NMF:150:PCA:eucledian.dist` ) |> 
  dplyr::group_by(pid) |> 
  dplyr::mutate(d = max(x)) |> 
  dplyr::mutate(init = x == min(x)) |> 
  dplyr::mutate(second = order(x) == 2) |> 
  dplyr::mutate(end = x == max(x)) |> 
  dplyr::mutate(n = n()) |> 
  dplyr::ungroup() |> 
  dplyr::rename(subtype = pred)


# add initial subtype to each segment, for facetting
plt.expanded <- plt.expanded |> 
  dplyr::left_join(
    plt.expanded |> 
      dplyr::filter(init) |> 
      dplyr::select(pid, subtype) |> 
      dplyr::rename(subtype.init = subtype),
    by=c('pid'='pid'), suffix=c('','')
  )

# add position of second segment, to normalize for right facet
plt.expanded <- plt.expanded |> 
  dplyr::left_join(
    plt.expanded |> 
      dplyr::filter(second) |> 
      dplyr::select(pid, x) |> 
      dplyr::rename(x.init = x),
    by=c('pid'='pid'), suffix=c('','')
  )

# plt.expanded <- plt.expanded |>
#  dplyr::filter(pid == 'GLSS-HF-DE05')
# print(plt.expanded)


# add HM status
plt.expanded <- rbind(
  plt.expanded,
  plt.expanded |>
    dplyr::filter(init == T) |>
    dplyr::filter(!is.na(HM) & HM == T) |>
    dplyr::mutate(x = - 0.025) |>
    dplyr::mutate(subtype = "Hyper-mutant") |>
    dplyr::mutate(segment = paste0("segment.",pid, ".HM")),
  
  plt.expanded |>
    dplyr::filter(init == T) |>
    dplyr::filter(!is.na(HM) & HM == T) |>
    dplyr::mutate(x = - 0.05) |>
    dplyr::mutate(subtype = "Hyper-mutant") |>
    dplyr::mutate(segment = paste0("segment.",pid, ".HM"))
)
print(plt.expanded)



# make facettes
plt.expanded <- rbind(
  plt.expanded |>
    dplyr::mutate(type ="a")
  ,
  plt.expanded |>
    dplyr::mutate(type ="b") |>
    dplyr::mutate(x = case_when(
      subtype == "Hyper-mutant" & x == -0.025 ~ -2.75,
      subtype == "Hyper-mutant" & x == -0.05 ~ -2.75 - 0.04,
      T ~ x - x.init)) |>
    dplyr::mutate(d = ifelse(n == 2, -x.init,d - x.init))
) |>
  dplyr::mutate(dataset = grepl("TCGA|GLSS",pid))




p1 <- ggplot(plt.expanded |> dplyr::filter(type=="a"), aes(x = x , y = reorder(pid, d), col = subtype, group=segment)) +
  geom_line(lwd = 0.8) + 
  facet_grid(rows = vars(facet.HM.determinated, subtype.init), scales = "free", space="free_y") + 
  labs(x = "Distance GITS space",
       y=NULL,
       col = "Distance within Subtype space") +
  theme_bw()  +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  scale_color_manual(values = c(subtype_colors, 'Hyper-mutant' = 'gray40'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


p2 <- ggplot(plt.expanded |> dplyr::filter(type=="b"), aes(x = x , y = reorder(pid, d), col = subtype, group=segment)) +
  #geom_point(cex=1) +
  geom_line(lwd = 0.8) + 
  facet_grid(rows = vars(facet.HM.determinated, subtype.init), scales = "free", space="free_y") + 
  labs(x = "Distance GITS space",
       y=NULL,
       col = "Distance within Subtype space") +
  theme_bw()  +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  scale_color_manual(values = c(subtype_colors, 'Hyper-mutant' = 'gray40'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


p1 + p2
ggsave("output/figures/2022_figure_S3b.pdf", width=8.3 / 2,height=8.3/2, scale=2)




## NMF 7k 3 x purity ----



plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::select(
      `sid`,
      `tumour.percentage.dna`,
      `NMF:150:error.v.per.M`,
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      
      `NMF:7k:error.v.per.M`,
      `NMF:7k:1`,
      `NMF:7k:2`,
      `NMF:7k:3`,
      `NMF:7k:4`
    )  |> 
    dplyr::rename(purity = tumour.percentage.dna) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(`tumour.percentage.2022` >= 15) |>
    dplyr::select(
      `aliquot_barcode`,
      `tumour.percentage.2022`,
      
      `NMF:150:error.v.per.M`,
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      
      `NMF:7k:error.v.per.M`,
      `NMF:7k:1`,
      `NMF:7k:2`,
      `NMF:7k:3`,
      `NMF:7k:4`
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(purity = tumour.percentage.2022) |> 
    dplyr::mutate(dataset = "GLASS")
)
 

ggplot(plt, aes(x=`NMF:7k:3` , y =`purity` , col = dataset, label=sid)) +
  facet_grid(rows = vars(dataset), scales = "free", space="free_y") + 
  geom_point() +
  theme_bw()  +
  #ggrepel::geom_label_repel() +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) + 
  stat_smooth(method="lm") +
  ggpubr::stat_cor( method = "pearson") 



# cor.test(plt$`NMF:150:error.v.per.M`,plt$purity, method="pearson")



# 〰 © Dr. Youri Hoogstrate 〰 ----

