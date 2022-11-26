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
    dplyr::filter(tumour.percentage.2022 >= 15) |>
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

tmp.n.gsam <- plt |> 
  dplyr::filter(dataset == "G-SAM") |> 
  nrow()
tmp.n.glass <- plt |> 
  dplyr::filter(dataset == "GLASS") |> 
  nrow()

tmp.p.gsam <- plt |> 
  dplyr::filter(dataset == "G-SAM") |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() > 1) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()
tmp.p.glass <- plt |> 
  dplyr::filter(dataset == "GLASS") |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() > 1) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()




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
  facet_grid(rows = vars(facet)) +
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
  theme(axis.title = element_text(face = "bold", size = rel(1))) +
  theme(legend.position = 'bottom') +
  labs(caption=paste0("G-SAM: n=",tmp.n.gsam, ", pairs=", tmp.p.gsam, " GLASS: n=",tmp.n.glass, ", pairs=", tmp.p.glass)  ) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.1)  )


ggsave("output/figures/2022_figure_1c.pdf", width=8.3 / 4,height=8.3/4 * 2.8, scale=2)
ggsave("output/figures/2022_figure_1c.svg", width=8.3 / 4,height=8.3/4 * 2.8, scale=2)



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
    #dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      `sid`,
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      `NMF:150:membership`,
      `GITS.150.svm.2022.subtype`
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    #dplyr::filter(tumour.percentage.2022 >= 15) |> 
    dplyr::mutate(is.primary = resection == "TP") |> 
    dplyr::select(
      `aliquot_barcode`,
      `NMF:150:1`,
      `NMF:150:2`,
      `NMF:150:3`,
      `NMF:150:membership`,
      `GITS.150.svm.2022.subtype`
    ) |>
    dplyr::rename(sid = aliquot_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
  ) |> 
  
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
  dplyr::mutate(`NMF H-matrix meta-feature` = gsub("NMF:150:","NMF meta-feature ",`NMF H-matrix meta-feature`)) |> 
  dplyr::filter(!is.na(value))


n.glass <- plt |>
  dplyr::filter(!is.na(value)) |> 
  dplyr::filter(!duplicated(.data$sid)) |>
  dplyr::pull(.data$dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |>
  dplyr::filter(!is.na(value)) |> 
  dplyr::filter(!duplicated(.data$sid)) |>
  dplyr::pull(.data$dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')



plt.line.based <- rbind(
  plt |>dplyr::mutate(value=0),
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
  labs(col = 'ssGSEA subtype',
       y = "NMF H-matrix values",
       x = NULL,
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples" )) + 
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


rm(n.gsam, n.glass, plt.line.based, plt)




## fig s1b ----


tmp.pca <- readRDS('tmp/analysis_GITS_space_GITS_PCA_150.Rds')
stopifnot(nrow(tmp.pca$x) == 287+216)


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


sel <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::mutate(dataset = "G-SAM") |> 
    dplyr::select(sid, dataset)
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |> 
    dplyr::rename(sid = aliquot_barcode) |> 
    dplyr::mutate(dataset = "GLASS") |> 
    dplyr::select(sid, dataset)
)

n.glass <- table(sel |> dplyr::pull(.data$dataset))['GLASS']
n.gsam <- table(sel |> dplyr::pull(.data$dataset))['G-SAM']


tmp.pca <- readRDS('tmp/analysis_GITS_space_GITS_PCA_150.Rds')
stopifnot(sel$sid %in% rownames(tmp.pca$x))

tmp.pca$x <- tmp.pca$x |> 
  as.data.frame(stringsAsFactors=F) |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::filter(.data$sid %in% sel$sid) |> 
  tibble::column_to_rownames('sid')


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
  labs(col = "Associated with",
       x="PC1 on NMF meta-features",
       y="PC2 on NMF meta-features",
       fill = "Subtype (GlioVis)",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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

rm(n.glass, n.gsam, mult, tmp.pca, data, x, y, datapc, sel)





## fig s1d ----



plt <- rbind(
  gsam.rna.metadata |>

    dplyr::filter(blacklist.pca == F) |> 
    dplyr::filter(pat.with.IDH == F) |> 
    dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F) |> 
    dplyr::filter(tumour.percentage.dna >= 15) |> # avoid NA values
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
    dplyr::filter(tumour.percentage.2022 >= 15) |> # avoid NA values
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


n.glass <- plt |> dplyr::pull(.data$dataset) |> table() |> purrr::pluck('GLASS')
n.gsam <- plt |> dplyr::pull(.data$dataset) |> table() |> purrr::pluck('G-SAM')


ggplot(plt, aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=ssGSEA.2022.subtype)) +
  geom_point(size=3,  col='black', pch=21, alpha=0.85) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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


rm(n.glass, n.gsam, plt)




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
    dplyr::filter(tumour.percentage.2022 >= 15) |> # avoid NA values
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


n.glass <- plt |> dplyr::pull(.data$dataset) |> table() |> purrr::pluck('GLASS')
n.gsam <- plt |> dplyr::pull(.data$dataset) |> table() |> purrr::pluck('G-SAM')


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
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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

rm(plt, n.glass, n.gsam, plt.contours)




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
    dplyr::filter(tumour.percentage.2022 >= 15) |> # avoid NA values
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


n.glass <- plt |> dplyr::pull(.data$dataset) |> table() |> purrr::pluck('GLASS')
n.gsam <- plt |> dplyr::pull(.data$dataset) |> table() |> purrr::pluck('G-SAM')


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
       shape="",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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


rm(plt, plt.contours, n.glass, n.gsam)




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
    dplyr::filter(tumour.percentage.2022 >= 15) |> # avoid NA values
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


n.glass <- plt |> 
  dplyr::filter(!is.na(`ssGSEA enrichment score`)) |> 
  dplyr::filter(!is.na(`NMF contribution`)) |> 
  dplyr::filter(!duplicated(sid)) |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |> 
  dplyr::filter(!is.na(`ssGSEA enrichment score`)) |> 
  dplyr::filter(!is.na(`NMF contribution`)) |> 
  dplyr::filter(!duplicated(sid)) |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')



ggplot(plt, aes(x=`ssGSEA enrichment score`, y=`NMF contribution`, fill=`ssGSEA.2022.subtype`, shape=misclass)) +
  #facet_grid(cols = vars(`facet`), scales = "free") +
  facet_wrap(~facet, scales = "free") +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL, label = ..r.label..), method = "pearson") +
  geom_point(size=3, data=plt |>  dplyr::filter(misclass != "."), col='black', pch=21, alpha=0.7) +
  geom_point(size=3, data=plt |>   dplyr::filter(misclass == ".") ,col='black', pch=21, alpha=0.7) +
  #geom_point(data = plt |>  dplyr::filter(misclass == "?"), col = 'white',  size=1.2) + # pch=4,
  geom_point(data = plt |>  dplyr::filter(misclass == "."), col = 'black',fill="white",pch=21, size=0.9,stroke=0.65) +
  labs(y="NMF meta-feature score", 
       fill = "Subtype (ssGSEA)",
       shape="",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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

rm(plt, n.glass, n.gsam)



## fig s1g equivalent for NMF:7k ----

# 1 = CL, 2 = PN, 4 = MES
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
      
      `NMF:7k:1`,
      `NMF:7k:2`,
      `NMF:7k:3`,
      `NMF:7k:4`,
      
      ssGSEA.2022.Classical.enrichment_score,
      ssGSEA.2022.Proneural.enrichment_score,
      tumour.percentage.dna,
      ssGSEA.2022.Mesenchymal.enrichment_score,
      
      ssGSEA.2022.subtype
    ) |> 
    dplyr::rename(purity = tumour.percentage.dna) |>
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    #dplyr::filter(tumour.percentage.2022 >= 15) |> 
    dplyr::select(
      aliquot_barcode,
      
      `NMF:7k:1`,
      `NMF:7k:2`,
      `NMF:7k:3`,
      `NMF:7k:4`,
      
      ssGSEA.2022.Classical.enrichment_score,
      ssGSEA.2022.Proneural.enrichment_score,
      tumour.percentage.2022,
      ssGSEA.2022.Mesenchymal.enrichment_score,
      
      ssGSEA.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(purity = tumour.percentage.2022) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
  tidyr::pivot_longer(cols=c(`NMF:7k:1`, `NMF:7k:2`,`NMF:7k:3`,`NMF:7k:4`),names_to="NMF meta-feature",values_to="NMF contribution") |> 
  tidyr::pivot_longer(cols=c(`ssGSEA.2022.Proneural.enrichment_score`,`ssGSEA.2022.Mesenchymal.enrichment_score`,purity,`ssGSEA.2022.Classical.enrichment_score`),
                      names_to="ssGSEA subtype", values_to="ssGSEA enrichment score") |> 
  dplyr::mutate(`ssGSEA subtype` = case_when(
    `ssGSEA subtype` == "ssGSEA.2022.Proneural.enrichment_score" ~ "PN",
    `ssGSEA subtype` == "ssGSEA.2022.Mesenchymal.enrichment_score" ~ "MES",
    `ssGSEA subtype` == "ssGSEA.2022.Classical.enrichment_score" ~ "CL",
    T ~ "purity"
  )) |> 
  dplyr::filter((`NMF meta-feature` == "NMF:7k:1" &  `ssGSEA subtype` == "CL") |
                  (  `NMF meta-feature` == "NMF:7k:2" &  `ssGSEA subtype` == "PN") |
                  (  `NMF meta-feature` == "NMF:7k:3" &  `ssGSEA subtype` == "purity") |
                  ( `NMF meta-feature` == "NMF:7k:4" &  `ssGSEA subtype` == "MES")) |> 
  dplyr::mutate(facet = paste0(`NMF meta-feature`," ~ " , `ssGSEA subtype`)) |> 
  dplyr::mutate(needle = case_when(`ssGSEA subtype` == "CL" ~ "Classical",
                                   `ssGSEA subtype` == "MES" ~ "Mesenchymal",
                                   `ssGSEA subtype` == "PN" ~ "Proneural",
                                   `ssGSEA subtype` == "purity" ~ "Purity"
  )) |> 
  dplyr::mutate(facet.target = stringr::str_detect(ssGSEA.2022.subtype, `needle`)) |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = as.character(`ssGSEA.2022.subtype`)) |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = ifelse(facet == "NMF:7k:3 ~ purity","-",`ssGSEA.2022.subtype`)) |> 
  dplyr::mutate(`needle` = NULL)


n.glass <- plt |> 
  dplyr::filter(`NMF meta-feature` == "NMF:7k:1") |> 
  dplyr::filter(!is.na(`NMF contribution`)) |> 
  dplyr::pull(.data$dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |> 
  dplyr::filter(`NMF meta-feature` == "NMF:7k:1") |> 
  dplyr::filter(!is.na(`NMF contribution`)) |> 
  dplyr::pull(.data$dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')


ggplot(plt, aes(x=`ssGSEA enrichment score`, y=`NMF contribution`, fill=`ssGSEA.2022.subtype`)) +
  facet_wrap(~facet, ncol=4, scales = "free") +
  geom_point(size=3, col='black', pch=21, alpha=0.7) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL, label = ..r.label..), method = "pearson") +
  labs(x = "ssGSEA enrichment score or tumor purity",
        y="NMF meta-feature score", 
       fill = "Subtype (ssGSEA)",
       shape="",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
  scale_fill_manual(values = c(mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.1),"-"="gray60"),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL',"-"="")) +
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


ggsave("output/figures/2022_figure_S1G_equivalent_for_NMF_7K.pdf", width=8.3 / 2,height=8.3/5, scale=2)


rm(plt, n.glass, n.gsam)



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
   dplyr::filter(tumour.percentage.2022 >= 15) |>  # erase NA's
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


n.glass <- plt |> 
  dplyr::filter(!is.na(`NMF:150:PC1`)) |> 
  dplyr::filter(!is.na(`NMF:150:PC2`)) |> 
  dplyr::filter(!duplicated(sid)) |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |> 
  dplyr::filter(!is.na(`NMF:150:PC1`)) |> 
  dplyr::filter(!is.na(`NMF:150:PC2`)) |> 
  dplyr::filter(!duplicated(sid)) |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')



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
       shape = "",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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


rm(plt, plt.contours, n.glass, n.gsam)




## fig s1k - ssGSEA PCA ----


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
    dplyr::filter(tumour.percentage.2022 >= 15) |>  # erase NA's
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


n.glass <- tmp |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- tmp |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')


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
  labs(x="PC1 on ssGSEA enrichment scores",
       y="PC2 on ssGSEA enrichment scores",
       fill = "Subtype (GlioVis)",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
  ) +
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

## fig s1k - 7k PCA equivalent ----


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
      
      `NMF:7k:1`,
      `NMF:7k:2`,
      `NMF:7k:3`,
      `NMF:7k:4`,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    )  |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |>  # erase NA's
    dplyr::select(
      aliquot_barcode,
      case_barcode,
      
      `NMF:7k:1`,
      `NMF:7k:2`,
      `NMF:7k:3`,
      `NMF:7k:4`,
      
      ssGSEA.2022.subtype,
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
)


n.glass <- tmp |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- tmp |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')


tmp.pca <- tmp |> 
  tibble::column_to_rownames('sid') |>
  dplyr::select(
    `NMF:7k:1`,
    `NMF:7k:2`,
    `NMF:7k:4`
  ) |>
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


ggplot(data, aes(x = PC1, y = PC2, col = varnames)) +
  geom_point(size = 3, col = "black", fill = "gray60", pch = 21, alpha = 0.85) +
  coord_equal() +
  geom_segment(data = datapc, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.4, "cm")), lwd = 1.75, alpha = 0.9) +
  ggrepel::geom_label_repel(
    data = datapc, aes(x = v1, y = v2, label = varnames), size = 3,
    vjust = 0.6, show.legend = F, col = "black"
  ) +
  labs(col = "Associated with") +
  labs(
    x = "PC1 on NMF [7k+ genes]",
    y = "PC2 on NMF [7k+ genes]",
    fill = "Subtype (GlioVis)",
    caption = paste0("G-SAM: n=", n.gsam, "  -  GLASS: n=", n.glass, " samples")
  ) +
  scale_color_manual(
    values = mixcol(
      c(
        "NMF:7k:1" = as.character(subtype_colors["Classical"]),
        "NMF:7k:2" = as.character(subtype_colors["Proneural"]),
        "NMF:7k:4" = as.character(subtype_colors["Mesenchymal"])
      ),
      rep("black", length(subtype_colors_ssGSEA)), 0.15
    ),
    labels = c("NMF:7k:1" = "CL", "NMF:7k:2" = "PN", "NMF:7k:4" = "MES")
  ) +
  theme_bw() +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  )


ggsave("output/figures/2022_figure_S1K_NMF-7k_eq.pdf", width=8.3 / 4,height=8.3/4, scale=2)


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
  dplyr::mutate(`ssGSEA subtype` = dplyr::recode(`ssGSEA subtype`,
                                                 'Classical' = 'CL',
                                                 'Mesenchymal' = 'MES',
                                                 'Proneural' = 'PN',
                                                 'Proneural|Classical'='PN|CL'
  )) |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = dplyr::recode(`ssGSEA.2022.subtype`,
                                                 'Classical' = 'CL',
                                                 'Mesenchymal' = 'MES',
                                                 'Proneural' = 'PN',
                                                 'Proneural|Classical'='PN|CL'
  )) |> 
  dplyr::mutate(`ssGSEA subtype` = paste0("ssGSEA ", `ssGSEA subtype`, " pval")) |> 
  dplyr::mutate(`ssGSEA.2022.subtype` = paste0("ssGSEA call: ", `ssGSEA.2022.subtype`, "")) |> 
  dplyr::mutate(order = order(order(ssGSEA.err1, ssGSEA.err2)))


n.glass <- plt |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype)) |> 
  dplyr::filter(!duplicated(sid)) |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype)) |> 
  dplyr::filter(!duplicated(sid)) |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')


plt.line.based <- rbind(
  plt |>  dplyr::mutate(`ssGSEA pval`=0),
  plt |>  dplyr::mutate(`ssGSEA pval` = -log(`ssGSEA pval`))
)



ggplot(plt.line.based, aes(x=reorder(sid, order), y=`ssGSEA pval`,
                           group=sid,
                           col=`GITS.150.svm.2022.subtype`)) +
  #geom_point(size=3,  col='black', pch=21, alpha=0.85) +
  geom_line(lwd=0.8) +
  facet_grid(cols = vars(ssGSEA.2022.subtype), 
             rows=vars(`ssGSEA subtype`),
             scale="free_x",
             space="free_x") +
  scale_color_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.15),
                    label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
  labs(x = "Patients", 
       y="-log10(ssGSEA pval)", 
       col="GITS subtype",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
       ) +
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
    axis.text.x=element_blank(),
    
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2))


ggsave("output/figures/2022_figure_S1L.pdf", width=8.3 / 4,height=8.3/4, scale=2)

rm(plt, plt.line.based, n.glass, n.gsam)



## figure S2a ----


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
    dplyr::filter(tumour.percentage.2022 >= 15) |> 
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


n.glass <- table(plt$dataset_primary)['GLASS']
n.gsam <- table(plt$dataset_primary)['G-SAM']



plt <- plt |>
  dplyr::rename(`NMF:150:PCA:eucledian.dist` = `NMF:150:PCA:eucledian.dist_primary`) |> 
  dplyr::mutate(`NMF:150:PCA:eucledian.dist_recurrent` = NULL) |> 
  dplyr::rename(`dataset` = `dataset_primary`) |> 
  dplyr::mutate(`dataset_recurrent` = NULL) |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype_primary` = paste0("from: ", `GITS.150.svm.2022.subtype_primary`)) |> 
  dplyr::mutate(`GITS.150.svm.2022.subtype_recurrent` = paste0("to: ", `GITS.150.svm.2022.subtype_recurrent`))



# FDR + geom_signif -> https://github.com/kassambara/ggpubr/issues/65#issuecomment-400918671
stats <- ggpubr::compare_means(`NMF:150:PCA:eucledian.dist` ~ `GITS.150.svm.2022.subtype_recurrent`,
                               group.by = "GITS.150.svm.2022.subtype_primary", 
                               data = plt
                               ) |>
  dplyr::mutate(y_pos = 4.25 + ((1:n() %% 3) * 0.4)) |> 
  dplyr::mutate(p.adj = format.pval(p.adj, digits = 1))


ggplot(plt, aes(x = GITS.150.svm.2022.subtype_recurrent, y = `NMF:150:PCA:eucledian.dist`)) +
  facet_grid(rows = vars(GITS.150.svm.2022.subtype_primary), space = "free_x") +
  ggsignif::geom_signif(
    data = stats,
    aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y_pos),
    manual = TRUE,
    textsize = 2.5,
    size = 0.375
  ) +
  geom_violin(width = 1.05, aes(fill = stable)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5, fill = NA, col = "black") +
  ggbeeswarm::geom_quasirandom(aes(fill = dataset), pch = 21, size = 3, cex = 4) +
  labs(
    x = NULL,
    y = "Euclidean distance GITS space",
    fill = "",
    caption = paste0("G-SAM: n=", n.gsam, "  -  GLASS: n=", n.glass, " pairs")
  ) +
  scale_fill_manual(values = c("stable" = "gray90", "transition" = "white", dataset_colors["G-SAM"], dataset_colors["GLASS"])) +
  scale_y_continuous(limits = c(0, 5.25)) + # for the signif
  theme_bw() +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    # strip.text = element_text(size = 7),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.1)
  ) +
  guides(fill = guide_legend(ncol = 4))



ggsave("output/figures/2022_figure_S2a.pdf", width=8.3 / ((3/2) * 2),height=8.3/3.2, scale=2)
ggsave("output/figures/2022_figure_S2a.svg", width=8.3 / ((3/2) * 2),height=8.3/3.2, scale=2)


rm(n.glass, n.gsam)
rm(plt)



## figure 1d ----




plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    
    dplyr::mutate(resection = ifelse(resection == "r1","primary","recurrent")) |> 
    dplyr::select(
      `resection`,
      `pid`,
      `GITS.150.svm.2022.subtype`
    ) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |> 
    dplyr::mutate(resection = ifelse(resection == "TP","primary","recurrent")) |>
    dplyr::select(
      `resection`,
      `case_barcode`,
      `GITS.150.svm.2022.subtype`
    ) |>
    #dplyr::rename(sid = aliquot_barcode) |> 
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(GITS.150.svm.2022.subtype = recode(GITS.150.svm.2022.subtype,
                                                   'Mesenchymal'='MES',
                                                   'Proneural'='PN',
                                                   'Classical'='CL'
                                                   )) |> 
  
  tidyr::pivot_wider(id_cols = pid, 
                     names_from = c(resection),
                     values_from = c(GITS.150.svm.2022.subtype, dataset)) |> 
  dplyr::mutate(dataset_primary = NULL) |> 
  dplyr::rename(dataset = dataset_recurrent) |> 
  dplyr::mutate(label = paste0(GITS.150.svm.2022.subtype_primary, " -> ", GITS.150.svm.2022.subtype_recurrent))

stopifnot(sum(duplicated(plt$pid)) == 0)



stats <- data.frame(table(plt$label)) |> 
  dplyr::left_join(data.frame(table(plt$label) / nrow(plt) * 100) |> dplyr::rename(Percentage = Freq), by=c('Var1'='Var1')) |> 
  dplyr::mutate(stats = paste0(Freq," (",round(Percentage,1),"%)" ) )




# plot 3 dots
l <- 10
plt.dots <- data.frame(x = c(0, l, sinpi(60/360) * l),
                  y= c(0,0, cospi(60/360) * l),
                  label = c("MES","CL","PN"),
                  type="dot"
                  )


plt.labels <- data.table::CJ(plt.dots$label, plt.dots$label) |>
  dplyr::rename(from = V1, to = V2) |>
  dplyr::left_join(plt.dots, by = c("from" = "label"), suffix = c("", "")) |>
  dplyr::mutate(type = NULL) |> 
  dplyr::rename(x_from = x) |> 
  dplyr::rename(y_from = y)  |>
  dplyr::left_join(plt.dots, by = c("to" = "label"), suffix = c("", "")) |>
  dplyr::mutate(type = NULL) |> 
  dplyr::rename(x_to = x) |> 
  dplyr::rename(y_to = y) |> 
  dplyr::mutate(id = paste0(from , " -> ", to)) |> 
  dplyr::left_join(stats, by=c('id'='Var1'),suffix=c('','')) |> 
  dplyr::mutate(label = paste0(id, ": " , stats)) |> 
  dplyr::mutate(x = (0.6 * x_from) + (0.4 * x_to)) |> 
  dplyr::mutate(y = (0.6 * y_from) + (0.4 * y_to)) |> 
  dplyr::mutate(x_from = NULL,   y_from = NULL, x_to = NULL,     y_to = NULL)



ggplot(plt.dots, aes(x=x,y=y,label=label,fill=label)) +
  geom_point(cex=30, pch=21) +
  geom_text(data=rbind(plt.dots |> dplyr::select(x,y,label), plt.labels|> dplyr::select(x,y,label)), cex=5, aes(fill=NULL)) +
  scale_x_continuous(expand = expansion(mult = .18*1.5)) +
  scale_y_continuous(expand = expansion(mult = .18*1.5)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )  + 
  coord_equal() +
  labs(x=NULL, y=NULL,fill=NULL) +
  scale_fill_manual(values = 
                      mixcol(c('MES' = as.character(subtype_colors['Mesenchymal']), 'PN' = as.character(subtype_colors['Proneural']), 'CL' = as.character(subtype_colors['Classical'])),rep("white",3),0.0)
    )


ggsave("output/figures/2022_figure_1d.pdf", width=8.3 / 2,height=8.3/2, scale=2)



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
    ) |> 
    dplyr::filter(!is.na(.data$`NMF:150:PC1.n`)) |> 
    dplyr::filter(!is.na(.data$`NMF:150:PC2.n`)) |> 
    
    dplyr::mutate(dataset = "G-SAM")
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
    dplyr::filter(!is.na(.data$`NMF:150:PC1.n`)) |> 
    dplyr::filter(!is.na(.data$`NMF:150:PC2.n`)) |> 
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::group_by(pid) |> # filter for matching pairs
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup()


n.glass <- plt |> 
  dplyr::filter(!is.na(.data$`NMF:150:PC1.n`)) |> 
  dplyr::filter(!is.na(.data$`NMF:150:PC2.n`)) |> 
  dplyr::filter(!duplicated(.data$pid)) |> 
  dplyr::pull(.data$dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |> 
  dplyr::filter(!is.na(.data$`NMF:150:PC1.n`)) |> 
  dplyr::filter(!is.na(.data$`NMF:150:PC2.n`)) |> 
  dplyr::filter(!duplicated(.data$pid)) |> 
  dplyr::pull(.data$dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')



plt <- plt |> 
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
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features", 
       fill="Subtype",
       caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," pairs" )
       ) +
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


rm(plt, n.glass, n.gsam)



## F] Figure S3A - all & merged ----


plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) |> 
    dplyr::filter(pat.with.IDH == F) |> 
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
    )) |> 
    dplyr::mutate(dataset = "G-SAM")
  ,
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |>
    dplyr::select(
      `case_barcode`,
      `NMF:150:PCA:eucledian.dist`,
      `HM`
    ) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::filter(!is.na(`NMF:150:PCA:eucledian.dist`)) |> 
  dplyr::group_by(pid) |> # filter for matching pairs
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |>
  dplyr::distinct() |> 
  dplyr::mutate(pid = as.character(pid)) |> 
  dplyr::mutate(facet.HM.determinated = is.na(HM))


n.glass <- table(plt$dataset)['GLASS']
n.gsam <- table(plt$dataset)['G-SAM']



stopifnot(readRDS('tmp/analysis_GITS_space.transition.segments.Rds') |>
            dplyr::filter(grepl("^GLSS|TCGA",pid)) |> 
            dplyr::filter(!duplicated(pid)) |> 
            nrow() == 87)
stopifnot(readRDS('tmp/analysis_GITS_space.transition.segments.Rds') |>
            dplyr::filter(grepl("^GLSS|TCGA",pid) == F) |> 
            dplyr::filter(!duplicated(pid)) |> 
            nrow() == 122)


# load transition segments 
plt.expanded <- readRDS('tmp/analysis_GITS_space.transition.segments.Rds') |> 
  dplyr::inner_join(plt, by=c('pid'='pid')) |> 
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
    by=c('pid'='pid'), suffix=c('',''))


# add position of second segment, to normalize for right facet
plt.expanded <- plt.expanded |> 
  dplyr::left_join(
    plt.expanded |> 
      dplyr::filter(second) |> 
      dplyr::select(pid, x) |> 
      dplyr::rename(x.init = x),
    by=c('pid'='pid'), suffix=c('','')
  )




# make facettes
plt.expanded <- rbind(
  plt.expanded |>
    dplyr::mutate(type ="a")
  ,
  plt.expanded |>
    dplyr::mutate(type ="b") |>
    dplyr::mutate(x = x - x.init) |>
    dplyr::mutate(d = ifelse(n == 2, -x.init,d - x.init))
) |>
  dplyr::mutate(dataset = grepl("TCGA|GLSS",pid))




p1 <- ggplot(plt.expanded |> dplyr::filter(type=="a"), aes(x = x , y = reorder(pid, d), col = subtype, group=segment)) +
  geom_line(lwd = 0.8) + 
  facet_grid(rows = vars(subtype.init), scales = "free", space="free_y") + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_color_manual(values = c(subtype_colors), 
                     label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL'),
                     guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_y_discrete(expand=expansion(add = 0.8))


p2 <- ggplot(plt.expanded |> dplyr::filter(type=="b"), aes(x = x , y = reorder(pid, d), col = subtype, group=segment)) +
  #geom_point(cex=1) +
  geom_line(lwd = 0.8) + 
  facet_grid(rows = vars(subtype.init), scales = "free", space="free_y") + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_color_manual(values = c(subtype_colors), 
                     label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL'),
                     guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_y_discrete(expand=expansion(add = 0.8))


p1 + p2 + patchwork::plot_annotation(caption =  paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," pairs" ))


ggsave("output/figures/2022_Figure_S3A.pdf", width=8.3 / 2,height=8.3/2.5, scale=2)



rm(n.gsam, n.glass, p1 , p2 , plt.expanded, plt)





## F] Figure S3B - HM status samples only ----


plt <- rbind(
  gsam.rna.metadata |>
    dplyr::filter(blacklist.pca == F) |>
    dplyr::filter(pat.with.IDH == F) |>
    dplyr::filter(sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |> # replicates
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::select(`pid`, `NMF:150:PCA:eucledian.dist`, ) |>
    dplyr::left_join(gsam.patient.metadata |>
      dplyr::select(`studyID`, `HM`),
    by = c("pid" = "studyID"), suffix = c("", "")
    ) |>
    dplyr::mutate(HM = case_when(
      HM == "No" ~ F,
      HM == "Yes" ~ T,
      T ~ as.logical(NA)
    )) |>
    dplyr::mutate(dataset = "G-SAM"),
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |>
    dplyr::select(
      `case_barcode`,
      `NMF:150:PCA:eucledian.dist`,
      `HM`
    ) |>
    dplyr::rename(pid = case_barcode) |>
    dplyr::mutate(dataset = "GLASS")
) |>
  dplyr::filter(!is.na(`NMF:150:PCA:eucledian.dist`)) |>
  dplyr::distinct() |>
  dplyr::mutate(pid = as.character(pid)) |>
  dplyr::filter(!is.na(HM))


n.glass <- table(plt$dataset)["GLASS"]
n.gsam <- table(plt$dataset)["G-SAM"]



# load transition segments
plt.expanded <- readRDS("tmp/analysis_GITS_space.transition.segments.Rds") |>
  dplyr::inner_join(plt, by = c("pid" = "pid")) |>
  dplyr::mutate(x = pct_transition * `NMF:150:PCA:eucledian.dist`) |>
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
    by = c("pid" = "pid"), suffix = c("", "")
  )

# add position of second segment, to normalize for right facet
plt.expanded <- plt.expanded |>
  dplyr::left_join(
    plt.expanded |>
      dplyr::filter(second) |>
      dplyr::select(pid, x) |>
      dplyr::rename(x.init = x),
    by = c("pid" = "pid"), suffix = c("", "")
  )


# add HM status
plt.expanded <- rbind(
  plt.expanded,
  plt.expanded |>
    dplyr::filter(init == T) |>
    dplyr::filter(!is.na(HM) & HM == T) |>
    dplyr::mutate(x = -0.025) |>
    dplyr::mutate(subtype = "Hyper-mutant") |>
    dplyr::mutate(segment = paste0("segment.", pid, ".HM")),
  plt.expanded |>
    dplyr::filter(init == T) |>
    dplyr::filter(!is.na(HM) & HM == T) |>
    dplyr::mutate(x = -0.05) |>
    dplyr::mutate(subtype = "Hyper-mutant") |>
    dplyr::mutate(segment = paste0("segment.", pid, ".HM"))
)
print(plt.expanded)



# make facettes
plt.expanded <- rbind(
  plt.expanded |>
    dplyr::mutate(type = "a"),
  plt.expanded |>
    dplyr::mutate(type = "b") |>
    dplyr::mutate(x = case_when(
      subtype == "Hyper-mutant" & x == -0.025 ~ -2.75,
      subtype == "Hyper-mutant" & x == -0.05 ~ -2.75 - 0.04,
      T ~ x - x.init
    )) |>
    dplyr::mutate(d = ifelse(n == 2, -x.init, d - x.init))
) |>
  dplyr::mutate(dataset = grepl("TCGA|GLSS", pid))




p1 <- ggplot(plt.expanded |> dplyr::filter(type=="a"),
             aes(x = x , y = reorder(pid, d), col = subtype, group=segment, size=subtype)) +
  geom_line(lwd = 0.8) + 
  facet_grid(rows = vars( subtype.init), scales = "free", space="free_y") + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_size_manual(values = c('Mesenchymal'= 0.8,'Proneural'= 0.8,'Classical'= 0.8,'Hyper-mutant' = 1.3 )  ) +
  scale_color_manual(values = c(mixcol(subtype_colors, rep("gray75",3), 0.85), 'Hyper-mutant' = 'gray40'),
                     label = c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL'),
                     guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_y_discrete(expand=expansion(add = 1.5))


p2 <- ggplot(plt.expanded |> dplyr::filter(type=="b"), aes(x = x , y = reorder(pid, d), col = subtype, group=segment)) +
  #geom_point(cex=1) +
  geom_line(lwd = 0.8) + 
  facet_grid(rows = vars( subtype.init), scales = "free", space="free_y") + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  scale_size_manual(values = c('Mesenchymal'= 0.8,'Proneural'= 0.8,'Classical'= 0.8,'Hyper-mutant' = 1.3 )  ) +
  scale_color_manual(values = c(mixcol(subtype_colors, rep("gray75",3), 0.7), 'Hyper-mutant' = 'gray40'),
                     label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL'),
                     guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_y_discrete(expand=expansion(add = 1.5))


p1 + p2 + patchwork::plot_annotation(caption =  paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," pairs [with HM status]" ))


ggsave("output/figures/2022_Figure_S3B.pdf", width=8.3 / 2,height=8.3/4.1, scale=2)
rm(n.gsam, n.glass, p1 , p2 , plt.expanded, plt)



## NMF 7k 3 x purity ----

#@comment redudendt fig because of s1g for NMF7k


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


n.glass <- plt |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('GLASS')
n.gsam <- plt |> 
  dplyr::pull(dataset) |> 
  table() |> 
  purrr::pluck('G-SAM')


 

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
  ggpubr::stat_cor( method = "pearson", col="black") 



# cor.test(plt$`NMF:150:error.v.per.M`,plt$purity, method="pearson")



# 〰 © Dr. Youri Hoogstrate 〰 ----

