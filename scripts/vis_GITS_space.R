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
    )
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
    dplyr::rename(pid = case_barcode)
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


  
# make facets
plt.expanded <- rbind(
  plt |> dplyr::mutate(facet = "Classical") |>
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Mesenchymal")  |> 
    dplyr::mutate(highlight = subtype.primary == facet),
  plt |> dplyr::mutate(facet = "Proneural")  |> 
    dplyr::mutate(highlight = subtype.primary == facet)
)



ggplot(plt.expanded, aes(x=-`NMF:150:PC1`,y=-`NMF:150:PC2`, group=pid, col =`GITS.150.svm.2022.subtype`,fill =`GITS.150.svm.2022.subtype`)) +
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
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors)



ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.classical | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.classical & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  #geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  #geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title=NULL) +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors) +
  facet_grid(cols = vars(facet)) +
  #  geom_text_repel(data = subset(plt, pid %in% c('CDA','CDD','CBP','CBR','CCZ') )) + 
  theme(strip.background = element_blank(), strip.text = element_blank())









## fig s1a ----

## fig s1b ----

## fig s1c ----

## fig s1d ----

## fig s1e ----

## fig s1f ----

## fig s1g ----

## fig s2 ----

## fig s3a ----

## fig s3b ----




#'@todo error.m x purity



plt <- data.frame(PC= paste0("PC",1:3),
                  var_explained=(tmp.pca$sdev)^2/sum((tmp.pca$sdev)^2)) %>%
  dplyr::mutate(label = paste0(round(var_explained * 100,1),"%"))

ggplot(plt, aes(x=PC,y=var_explained, group=1, label=label))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_col(fill="gray60", col="black", lwd=0.4) +
  geom_line()+
  geom_label(y = 0.35, size=4) +
  geom_point(size=4)+
  labs(title="Scree plot: PCA on NMF (for p=3)", x= NULL, y="Variance in NMF V-matrix explained") +
  theme_bw()




### ssGSEA x 150 ----

# combined 
plt <- tmp.out |>
  dplyr::filter(!is.na(ssGSEA.2022.subtype)) 

p1 = ggplot(plt , aes(x=`NMF:150:1`, y=ssGSEA.2022.Classical.enrichment_score,
                 shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p2 = ggplot(plt , aes(x=`NMF:150:2`, y=`ssGSEA.2022.Proneural.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p3 = ggplot(plt , aes(x=`NMF:150:3`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw()

p1 + p2 + p3


#ggsave("output/figures/")



# # GLASS
# plt <- tmp.out |>
#   dplyr::filter(!is.na(ssGSEA.Synapse.subtype.2022))
# ggplot(plt , aes(x=`NMF:150:1`, y=ssGSEA.Synapse.subtype.2022.Classical.enrichment_score, col=ssGSEA.Synapse.subtype.2022)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:2`, y=ssGSEA.Synapse.subtype.2022.Proneural.enrichment_score, col=ssGSEA.Synapse.subtype.2022)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:3`, y=ssGSEA.Synapse.subtype.2022.Mesenchymal.enrichment_score, col=ssGSEA.Synapse.subtype.2022)) +
#   geom_point() +
#   theme_bw()
# 
# 
# 
# # GSAM
# plt <- tmp.out |>
#   dplyr::filter(!is.na(gliovis.ssGSEA.Mesenchymal.score))
# ggplot(plt , aes(x=`NMF:150:1`, y=gliovis.ssGSEA.Classical.score, col=gliovis.gsea_call)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:2`, y=gliovis.ssGSEA.Proneural.score, col=gliovis.gsea_call)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:3`, y=gliovis.ssGSEA.Mesenchymal.score, col=gliovis.gsea_call)) +
#   geom_point() +
#   theme_bw()


### ssGSEA x 7k ----



### ssGSEA x 150 ----

# combined 
plt <- tmp.out |>
  dplyr::filter(!is.na(ssGSEA.2022.subtype)) 

p1 = ggplot(plt , aes(x=`NMF:7k:1`, y=ssGSEA.2022.Classical.enrichment_score,
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p2 = ggplot(plt , aes(x=`NMF:7k:2`, y=`ssGSEA.2022.Proneural.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p3 = ggplot(plt , aes(x=`NMF:7k:4`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw()

p1 + p2 + p3



# GLASS - match NMF's w/ subtypes
plt <- tmp.out |>
  dplyr::filter(!is.na(ssGSEA.2022.subtype))

# GLASS - match NMF's w/ subtypes
ggplot(plt , aes(x=`NMF:7k:1`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
  geom_point() +
  theme_bw()
# ggplot(plt , aes(x=`NMF:7k:2`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:3`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:4`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()

ggplot(plt , aes(x=`NMF:7k:2`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
  geom_point() +
  theme_bw()
# ggplot(plt , aes(x=`NMF:7k:1`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:3`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:4`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()

ggplot(plt , aes(x=`NMF:7k:4`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
  geom_point() +
  theme_bw()
# ggplot(plt , aes(x=`NMF:7k:1`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:2`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:3`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()


### ssGSEA [glass] ----

### PCA ----


tmp.pca <- tmp.out |>
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022)) |> 
  dplyr::select(sid, 
                GBM.transcriptional.subtype.Synapse.Proneural.enrichment_score,
                GBM.transcriptional.subtype.Synapse.Classical.enrichment_score,
                GBM.transcriptional.subtype.Synapse.Mesenchymal.enrichment_score
  ) |> 
  tibble::column_to_rownames('sid') |> 
  prcomp()


tmp.out <- tmp.out |> 
  dplyr::left_join(
    tmp.pca |> 
      purrr::pluck('x') |> 
      as.data.frame() |> 
      dplyr::rename_with( ~ paste0("ssGSEA:GLASS:", .x)) |> 
      tibble::rownames_to_column('sid'),
    by=c('sid'='sid'),suffix=c('','')
  )


### more ----



plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p1 = ggplot(plt, aes(x=`ssGSEA:PC1`, y=`ssGSEA:PC2`,col=GBM.transcriptional.subtype.Synapse.2022)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p2 = ggplot(plt, aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`,col=GBM.transcriptional.subtype.Synapse.2022)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p3 = ggplot(plt, aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`,col=GBM.transcriptional.subtype.Synapse.2022)) +
  geom_point() +
  theme_bw()


plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p4 = ggplot(plt, aes(x=`ssGSEA:PC1`, y=`ssGSEA:PC2`,col=`NMF:150:membership`)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p5 = ggplot(plt, aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`,col=`NMF:150:membership`)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p6 = ggplot(plt, aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`,col=`NMF:150:membership`)) +
  geom_point() +
  theme_bw()


(p1 + p2 + p3) /
  (p4 + p5 + p6) 

ggplot(plt , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |>
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021))
ggplot(plt , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`GBM.transcriptional.subtype.Synapse.2021`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`, col=`GBM.transcriptional.subtype.Synapse.2021`)) +
  geom_point() +
  theme_bw()





ggplot(plt , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()

p2 = ggplot(tmp.out , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=ds)) +
  geom_point() +
  theme_bw()

ggplot(plt , aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()


ggplot(plt , aes(x=`NMF:150:1`, y=log(1 - `GBM.transcriptional.subtype.Synapse.Classical.pval` * 0.999), col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:2`, y=`GBM.transcriptional.subtype.Synapse.Proneural.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:3`, y=`GBM.transcriptional.subtype.Synapse.Mesenchymal.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()

ggplot(plt , aes(x=`NMF:7k:1`, y=`GBM.transcriptional.subtype.Synapse.Classical.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:2`, y=`GBM.transcriptional.subtype.Synapse.Classical.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:3`, y=`GBM.transcriptional.subtype.Synapse.Classical.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()



ggplot(plt , aes(x=`NMF:7k:1`, y=log(GBM.transcriptional.subtype.Synapse.Classical.enrichment_score), col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:2`, y=`GBM.transcriptional.subtype.Synapse.Proneural.enrichment_score`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:3`, y=`GBM.transcriptional.subtype.Synapse.Mesenchymal.enrichment_score`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()




ggplot(tmp.out, aes(x=`NMF:7k:PC1`, y=`NMF:150:PC1`,col=ds)) + 
  geom_point() +
  theme_bw()

ggplot(tmp.out, aes(x=`NMF:7k:PC2`, y=`NMF:150:PC2`,col=ds)) + 
  geom_point() +
  theme_bw()


ggplot(tmp.out, aes(x=`NMF:7k:1`, y=`NMF:150:3`,col=ds)) + 
  geom_point() +
  theme_bw()

ggplot(tmp.out, aes(x=`NMF:7k:2`, y=`NMF:150:3`,col=ds)) + 
  geom_point() +
  theme_bw()

ggplot(tmp.out, aes(x=`NMF:7k:3`, y=`NMF:150:1`, col=ds)) + 
  geom_point() +
  theme_bw()


### some plots ----

# 
# plt <- data.frame(PC= paste0("PC",1:3),
#                   var_explained=(tmp.pca$sdev)^2/sum((tmp.pca$sdev)^2)) %>%
#   dplyr::mutate(label = paste0(round(var_explained * 100,1),"%"))
# 
# ggplot(plt, aes(x=PC,y=var_explained, group=1, label=label))+
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   geom_col(fill="gray60", col="black", lwd=0.4) +
#   geom_line()+
#   geom_label(y = 0.35, size=4) + 
#   geom_point(size=4)+
#   labs(title="Scree plot: PCA on NMF (for p=3)", x= NULL, y="Variance in NMF V-matrix explained") + 
#   theme_bw()



p1 = ggplot(tmp.out , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`NMF:7k:membership`)) +
  geom_point() +
  theme_bw()

p2 = ggplot(tmp.out , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=ds)) +
  geom_point() +
  theme_bw()

p1+p2


PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x, varnames='NMF meta-feature 1')
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation) %>%
    dplyr::mutate(varnames = gsub('NMF:123456.','NMF meta-feature ', varnames))
  
  #print(head(datapc$varnames))
  #print(head(datapc))
  
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  
  ggplot(data, aes(x=PC1, y=PC2, col=varnames)) +
    geom_point(size=2.5,  col='black', fill='gray', pch=21) +
    coord_equal() +
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.4,"cm")), lwd=1.3, alpha=0.7) +
    geom_label_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3.5, vjust=1, show.legend=F) +
    youri_gg_theme +
    labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill = "Sub-type (GlioVis)") +
    scale_color_manual(values = subtype_colors_nmfm, labels = c('NMF meta-feature 2'='Classical',
                                                                'NMF meta-feature 1'='Mesenchymal',
                                                                'NMF meta-feature 3'='Proneural')) +
    labs(col = "Associated sub-type")
  #+ theme(legend.position = "none")
}


PCbiplot(p)




## + SVM classification ----
# Aanzienlijk lekkerdere fit dan LDA




if(file.exists("output/tables/gsam_nmf_lda_data.txt")) {
  a = read.table("output/tables/gsam_nmf_lda_data.txt") 
  
  colnames(a)
}  

s150.pca.nmf.subtype.classifier.svm <- readRDS('tmp/s150.pca.nmf.subtype.classifier.svm.Rds')
# s150.pca.nmf.subtype.classifier.svm <- svm(x = plt.single  %>%
#                                              #dplyr::filter(dataset == "GSAM") %>%
#                                              dplyr::select(c('sid' ,'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
#                                              tibble::column_to_rownames('sid') ,
#                                            y = plt.single %>% 
#                                              #dplyr::filter(dataset == "GSAM") %>%
#                                              dplyr::pull(subtype.public), 
#                                            scale = F,
#                                            #type = "C-classification",
# 
#                                            kernel = 'linear',
#                                            tolerance = 0.0001,
#                                            cost = 3
#                                            
#                                            #,probability = T # worse fit, but handy values
#                                            )
#saveRDS(s150.pca.nmf.subtype.classifier.svm, 'tmp/s150.pca.nmf.subtype.classifier.svm.Rds')


# re-fit samples

#attr(train.pred, 'decision.values')
plt.single <- plt.single %>% dplyr::left_join(
  data.frame(
    sid = plt.single$sid,
    
    'NMF:123456.PCA.SVM.class' = as.character(
      predict(object = s150.pca.nmf.subtype.classifier.svm , newdata =  (plt.single  %>%
                                                dplyr::select(c('sid' ,'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
                                                tibble::column_to_rownames('sid')),
              decision.values = T,
              tolerance = 0.0001,
              cost = 3
              )),
    check.names = F
    ) , by = c('sid'='sid') )


plt.single %>% dplyr::filter(dataset == "GSAM") %>% dplyr::mutate(err = subtype.public != `NMF:123456.PCA.SVM.class`) %>% dplyr::pull(err) %>% summary
plt.single %>% dplyr::filter(dataset != "GSAM") %>% dplyr::mutate(err = subtype.public != `NMF:123456.PCA.SVM.class`) %>% dplyr::pull(err) %>% summary
plt.single %>% dplyr::mutate(err = subtype.public != `NMF:123456.PCA.SVM.class`) %>% dplyr::pull(err) %>% summary


### :: redo to estimate probabilities ----
#
# This prediction performs considerable worse
# but probabilities may be useful to other analysis
#

tmp.fit <- svm(
  x = plt.single  %>%
    dplyr::select(c(
      'sid' , 'NMF:123456.PC1', 'NMF:123456.PC2'
    )) %>%
    tibble::column_to_rownames('sid') ,
  y = plt.single %>%
    dplyr::pull(subtype.public),
  scale = F,
  kernel = 'linear',
  tolerance = 0.0001,
  cost = 3,
  probability = T)



tmp <- predict(
  object = tmp.fit ,
  newdata = (
    plt.single %>% dplyr::select(c(
      'sid' , 'NMF:123456.PC1', 'NMF:123456.PC2'
    )) %>% tibble::column_to_rownames('sid')
  ),
  decision.values = T,
  tolerance = 0.0001,
  probability = T,
  cost = 3
) %>% attr('probabilities')  %>%
  as.data.frame () %>%
  `colnames<-`(paste0("NMF:123456.PCA.SVM.", colnames(.), ".p")) %>%
  tibble::rownames_to_column('sid')



plt.single <- plt.single %>%
  dplyr::left_join(tmp, by = c('sid' = 'sid'))



rm(tmp.fit, tmp)


# export without overwriting old SVM
#
# tmp.export <- read.table("output/tables/gsam_nmf_lda_data.txt.old") %>%
#   dplyr::left_join(tmp, by = c('sid' = 'sid'))
# 
# write.table(tmp.export, "output/tables/gsam_nmf_lda_data.txt")


# ~~ Export stats ~~ ----

write.table(
  tmp.out,
  "output/tables/gsam_nmf_lda_data.new.txt")




# Per pair stats ----
## Split table into pairs + eucl dist ----


plt.paired <- plt.single %>%
  #dplyr::filter(resection %in% c("R1", "R2") ) %>%
  dplyr::group_by(pid) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::top_n(1, sid) %>%
  dplyr::select(all_of('pid')) %>%
  as.data.frame() %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R1') %>% `colnames<-`( paste0(colnames(.) , ".R1") ) , by=c('pid'='pid.R1') ) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection %in% c('R2', 'R3', 'R4') ) %>% `colnames<-`( paste0(colnames(.) , ".R2") ) , by=c('pid'='pid.R2') ) %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                        (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2)) %>%
  dplyr::mutate(subtype.public.status = as.factor(ifelse(subtype.public.R1 == subtype.public.R2, "Stable", "Transition"))) %>%
  #dplyr::mutate(`NMF:123456.membership.status` = as.factor(ifelse(`NMF:123456.membership.R1`== `NMF:123456.membership.R2`, "Stable", "Transition"))) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.factor(ifelse(`NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`, "Stable", "Transition")))



plt.paired %>%   dplyr::filter(pid %in% c('G-SM-R056-2', 'GLSS-HF-3081', 'GLSS-HF-2869')  )


plot(density(plt.paired$eucledian.dist))
hist(plt.paired$eucledian.dist,breaks=15)

plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>% 
  nrow()

plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  nrow()


plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  nrow()



dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` != 'Classical') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)




dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` != 'Mesenchymal') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)



dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)



dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)




# Tests ----

## Verschil eucl. distance paren en geshufflede paren ----
## maak een seed voor shufflen voor reproduceerbaarheid of maak alle mogelijke paren
## maak bijbehorende boxplots?

env.test <- env()
env.test$pairs <- plt.paired
env.test$shuffle <- tidyr::crossing(
  plt.paired %>%
    dplyr::select('pid','sid.R1') %>%
    dplyr::rename(pid.R1 = pid) , 
  plt.paired %>%
    dplyr::select('pid','sid.R2') %>%
    dplyr::rename(pid.R2 = pid)
) %>% dplyr::filter(pid.R1 != pid.R2) %>%
  dplyr::mutate(pid.R1 = NULL, pid.R2 = NULL) %>%
  dplyr::left_join(plt.paired %>% dplyr::select(ends_with('.R1')), by = c('sid.R1'='sid.R1')) %>%
  dplyr::left_join(plt.paired %>% dplyr::select(ends_with('.R2')), by = c('sid.R2'='sid.R2')) %>%
  as.data.frame %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                      (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2))



plt.1 <- data.frame(val = env.test$pairs$eucledian.dist , type = 'Actual distances' ) %>% 
  dplyr::mutate(order = (rank(val)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.2 <- data.frame(val = env.test$shuffle$eucledian.dist , type = 'Shuffled distances' ) %>% 
  dplyr::mutate(order = (rank(val)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.5) # band scaling

plt <- rbind(plt.1, plt.2) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type)))

ggplot(plt, aes(x=type, y=val , group = type)) +
  geom_violin() +
  geom_point(aes(x=x),cex=0.8,col="#888888") +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  youri_gg_theme +
  labs(x = NULL, y = "Eucledian distance PNC space")


wilcox.test(env.test$pairs$eucledian.dist, env.test$shuffle$eucledian.dist)


#plot(hist(env.test$shuffle$eucledian.dist,breaks=25))
# library(fitdistrplus)
# plot(fitdist(env.test$shuffle$eucledian.dist , "weibull"))  # meest lijkende
# plot(fitdist(env.test$shuffle$eucledian.dist, "gamma"))
# 
# plot(fitdist(env.test$shuffle$eucledian.dist, "geom"))
# plot(fitdist(env.test$shuffle$eucledian.dist, "cauchy")) # nope!
# #plot(fitdist(env.test$shuffle$eucledian.dist, "t"))
# plot(fitdist(round(env.test$shuffle$eucledian.dist), "pois"))
# #plot(fitdist(round(env.test$shuffle$eucledian.dist), "nbinom"))
# plot( fitdist(env.test$shuffle$eucledian.dist, "lnorm")) # deze zeker niet
# plot( fitdist(env.test$shuffle$eucledian.dist, "exp")) # deze zeker niet
# plot( fitdist(env.test$shuffle$eucledian.dist , "logis")) # deze ook niet



# 
# 
# plt <- data.frame(val = env.test$pairs$eucledian.dist , type = 'real' ) %>% 
#   rbind(data.frame(val = env.test$shuffle %>%
#                      dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`) %>%
#                      dplyr::pull(eucledian.dist)
#                     , type = 'shuffle'))
# 
# 
# ggplot(plt, aes(x=type, y=val )) +
#   geom_violin() +
#   geom_boxplot(width=0.1)
# 

## Verschil distance transities ----

# https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc

t <- seq(0, 180, by = 1) * pi / 180

r <- .5
x <- r * cos(t)
y <- r*4 * sin(t)
y[20:162] <- y[20] # Flattens the arc
y <- y * 0.3
arc.a <- data.frame(x = x, eucledian.dist = y, type="GLASS")


r <- .5
x <- r * cos(t) * 2
y <- r*4 * sin(t)
y[20:162] <- y[20] # Flattens the arc
y <- y * 0.3
arc.b <- data.frame(x = x, eucledian.dist = y, type="GLASS")


r <- .5
x <- r * cos(t) * 1 + 0.5
y <- r*4 * sin(t)
y[20:162] <- y[20] # Flattens the arc
y <- y * 0.3
arc.c <- data.frame(x = x, eucledian.dist = y, type="GLASS")



### Classical ----


plt.c.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Classical\nto\nClassical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.c.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Classical\nto\nProneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.c.p <- wilcox.test(plt.c.c$eucledian.dist , plt.c.p$eucledian.dist )$p.value

plt.c.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Classical\nto\nMesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.c.m <- wilcox.test(plt.c.c$eucledian.dist , plt.c.m$eucledian.dist )$p.value

plt <- rbind(plt.c.c , plt.c.p, plt.c.m) %>%
  dplyr::mutate(type = factor(type, levels = c("Classical\nto\nClassical", 
                                               "Classical\nto\nMesenchymal",
                                               "Classical\nto\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS")))


plt.c <- ggplot(plt, aes(x=type, y=eucledian.dist , group = type , fill = type == "Classical\nto\nClassical")) +
  ylim(0, 4.3) +
  geom_violin(width=1.105) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col = dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset", fill = "Stable") +
  youri_gg_theme + 
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_line(data = arc.a, aes(x = x+1.5, y = eucledian.dist+3.85), lty = 2) +
  geom_line(data = arc.b, aes(x = x+2, y = eucledian.dist+3.85 + 0.2), lty = 2) +
  geom_text(data = data.frame(x = c(2), eucledian.dist = c(4.3), type="GLASS", label = "**" ),
            size=6,
            aes(x=x, y=eucledian.dist, label=label)) +
  geom_text(data = data.frame(x = c(1.5), eucledian.dist = c(4.05), type="GLASS", label = "***" ),
            size=6,
            aes(x=x, y=eucledian.dist, label=label)) +
  scale_fill_manual(values = c('TRUE' = 'gray95', 'FALSE' = 'white'))


rm(plt.c.c , plt.c.p, plt.c.m)




### Mesenchymal ----


plt.m.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Mesenchymal\nto\nMesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.m.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Mesenchymal\nto\nClassical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.m.c <- wilcox.test(plt.m.m$eucledian.dist , plt.m.c$eucledian.dist )$p.value

plt.m.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Mesenchymal\nto\nProneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.m.p <- wilcox.test(plt.m.m$eucledian.dist , plt.m.p$eucledian.dist )$p.value

plt <- rbind(plt.m.c , plt.m.p, plt.m.m) %>%
  dplyr::mutate(type = factor(type, levels = c("Mesenchymal\nto\nClassical",
                                               "Mesenchymal\nto\nMesenchymal",
                                               "Mesenchymal\nto\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS")))


plt.m <- ggplot(plt, aes(x=type, y=eucledian.dist , group = type, fill = type == "Mesenchymal\nto\nMesenchymal")) +
  ylim(0, 4.3) +
  geom_violin(width=1.105) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col = dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset", fill = "Stable") +
  youri_gg_theme + 
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c('TRUE' = 'gray95', 'FALSE' = 'white'))

# no sig diff's



rm(plt.m.c , plt.m.p, plt.m.m)



### Proneural ----



plt.p.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Proneural\nto\nProneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.p.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Proneural\nto\nClassical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.p.c <- wilcox.test(plt.p.p$eucledian.dist , plt.p.c$eucledian.dist )$p.value 

plt.p.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Proneural\nto\nMesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.p.m <- wilcox.test(plt.p.p$eucledian.dist , plt.p.m$eucledian.dist )$p.value


plt <- rbind(plt.p.c , plt.p.p, plt.p.m) %>%
  dplyr::mutate(type = factor(type, levels = c("Proneural\nto\nClassical",
                                               "Proneural\nto\nMesenchymal",
                                               "Proneural\nto\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(type)) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet = "Proneural")


plt.p <- ggplot(plt, aes(x=type, y=eucledian.dist , group = type, fill = type == "Proneural\nto\nProneural")) +
  ylim(0, 4.3) +
  geom_violin(width=1.105) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col=dataset, fill=dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset", fill = "Stable") +
  youri_gg_theme + 
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  #geom_line(data = arc.a, aes(x = x+1.5, y = eucledian.dist+3.85), lty = 2) +
  geom_line(data = arc.c, aes(x = x+2, y = eucledian.dist+3.85 + 0.2), lty = 2) +
  geom_text(data = data.frame(x = c(2.5), eucledian.dist = c(4.3), type="GLASS", label = "*" ),
            size=6,
            aes(x=x, y=eucledian.dist, label=label)) +
  scale_fill_manual(values = c('TRUE' = 'gray95', 'FALSE' = 'white'))

rm(plt.p.c , plt.p.p, plt.p.m)


### combined ----

plt.c / plt.m / plt.p


ggsave('output/figures/paper_subtypes_nmf_eucledian_dists_stable_unstable.pdf',width=8 ,height=6*1.5)


plt.c.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.c.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.c.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.c <- rbind(plt.c.c , plt.c.p, plt.c.m) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) %>%
  dplyr::mutate(type = factor(paste0('to\n', `NMF:123456.PCA.SVM.class.R2`), levels = c("to\nClassical","to\nMesenchymal","to\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet = "Classical")



plt.p.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.p.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.p.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.p <- rbind(plt.p.c , plt.p.p, plt.p.m) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4)  %>%
  dplyr::mutate(type = factor(paste0('to\n', `NMF:123456.PCA.SVM.class.R2`), levels = c("to\nClassical","to\nMesenchymal","to\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet = "Proneural")


plt.m.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.m.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.m.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.m <- rbind(plt.m.c , plt.m.p, plt.m.m) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) %>%
  dplyr::mutate(type = factor(paste0('to\n', `NMF:123456.PCA.SVM.class.R2`), levels = c("to\nClassical","to\nMesenchymal","to\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet ="Mesenchymal")


plt <- rbind(plt.c, plt.m, plt.p) %>%
  dplyr::mutate(stable =   `NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`)


ggplot(plt, aes(x=type, y=eucledian.dist , group = type)) +
  ylim(0, 4.3) +
  geom_violin(width=1.105, aes(fill = stable)) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col = dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset") +
  youri_gg_theme +
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(rows = vars(facet), space="free_x") +
  scale_fill_manual(values = c('TRUE'='gray95', 'FALSE'='white'))


ggsave('output/figures/paper_subtypes_nmf_eucledian_dists_stable_unstable.pdf',width=8 ,height=6*1.5)




# Plots ----


## Determine Contour ----


resolution <- 250 # 1000 x 1000 data points

off_x <- (max(plt.single$`NMF:123456.PC1`) - min(plt.single$`NMF:123456.PC1`)) * 0.025
off_y <- (max(plt.single$`NMF:123456.PC2`) - min(plt.single$`NMF:123456.PC2`)) * 0.025


range_pc1 = seq(from = min(plt.single$`NMF:123456.PC1`) - off_x, to = max(plt.single$`NMF:123456.PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(plt.single$`NMF:123456.PC2`) - off_y, to = max(plt.single$`NMF:123456.PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:123456.PC1' = range_pc1, 'NMF:123456.PC2' = range_pc2)
#nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.lda, newdata = range_df) %>%
nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = range_df) %>% data.frame(class = .) %>%
  #as.data.frame() %>%
  cbind(range_df) %>%
  dplyr::select(c('class', 'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
  dplyr::mutate(type="Contour")

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)


# Figure S1E ----


plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'subtype.public', sid.label, dataset)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = subtype.public),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA))



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col=class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_point(data = subset(plt, type == "Patient Sample" ), aes(fill = class), col="black", 
             pch=21, size=2.5) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features",
       col='Sub-type',fill="Sub-type") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


ggsave("output/figures/Figure_S1_E.png", width=12 * 0.4,height=10 * 0.4)



# Figure S1F ----


plt <- rbind(plt.single %>%
               dplyr::mutate(svm.status = ifelse(subtype.public == `NMF:123456.PCA.SVM.class`, "positive", 'negative')) %>% 
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.SVM.class', sid.label, dataset, svm.status)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.SVM.class`),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA, svm.status = NA))



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col=class, label = sid.label)) + 
  #geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1, na.rm = T) +
  geom_point(data = subset(plt, type == "Patient Sample" ), aes(fill = class), col="black", 
             pch=21, size=2.5) +
  geom_point(data = subset(plt, type == "Patient Sample" & svm.status == 'negative' ),
             col = 'black', pch=19, size=0.7) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(as.factor(class))),
               colour="gray40",
               size=0.25,
               na.rm = T,
               lty=2,breaks=c(1.5,2.5)) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features",
       col='Sub-type',fill="Sub-type") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


ggsave("output/figures/Figure_S1_F.png", width=12 * 0.4,height=10 * 0.4)



# Figure 2 ABC: 3 GITS panels ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.classical = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(primary.proneural = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(complete.pair= pid %in% (plt.single %>% dplyr::filter(resection == "R1") %>% dplyr::pull(pid)) &
                  pid %in% (plt.single %>% dplyr::filter(resection %in% c("R2","R3","R4")) %>% dplyr::pull(pid))                  
  ) %>%
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.character(`NMF:123456.PCA.SVM.status`) ) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = ifelse(is.na(`NMF:123456.PCA.SVM.status`), 'NA', `NMF:123456.PCA.SVM.status`))





ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.classical | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.classical & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_classical.pdf', width=12 * 0.45,height=10 * 0.45)







ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.mesenchymal | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.mesenchymal & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_mesenchymal.pdf', width=12 * 0.45,height=10 * 0.45)





ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.proneural | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.proneural & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_proneural.pdf', width=12 * 0.45,height=10 * 0.45)







## PCA:1+2(NMF:1+2+3) + SVM contours + GlioVis labels ----


plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'subtype.public', sid.label, dataset)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = subtype.public),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA))



plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM")
       , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample" & dataset == "GSAM"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



plt.glass <- ggplot(plt %>% dplyr::filter(dataset != "GSAM")
       , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample" & dataset != "GSAM"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_GlioVis_LDA-countours.png',width=7,height=7.2)




## PCA:1+2(NMF:1+2+3) + SVM contours + SVM labels ----



plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.SVM.class')) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.SVM.class`),
             nmf.pca.lda.countours)

ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, fill = class)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               na.rm = T,
               lty=2,
               breaks=c(1.5,2.5)
               ) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.8, pch=21, lwd=0.5, col="gray20") +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call] / GLASS consortium',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours.png',width=7,height=7.2)
ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours.pdf',width=7,height=7.2)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [classical] ----


plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.classical = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(complete.pair= pid %in% (plt.single %>% dplyr::filter(resection == "R1") %>% dplyr::pull(pid)) &
                               pid %in% (plt.single %>% dplyr::filter(resection %in% c("R2","R3","R4")) %>% dplyr::pull(pid))                  
                ) %>%
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.character(`NMF:123456.PCA.SVM.status`) ) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = ifelse(is.na(`NMF:123456.PCA.SVM.status`), 'NA', `NMF:123456.PCA.SVM.status`))
  # %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4

# 
# plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
#   geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
#   geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
#                colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
#   geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Classical" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1"), size=1.5, alpha=0.8) +
#   geom_path(
#     data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
#       dplyr::filter(primary.classical == T),
#     arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
#   youri_gg_theme +
#   labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
#   scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
#   scale_fill_manual(values = subtype_colors)
# 
# 
# plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
#   geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
#   geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
#                colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
#   geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Classical" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1"), size=1.5, alpha=0.8) +
#   geom_path(
#     data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
#       dplyr::filter(primary.classical == T),
#     arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
#   youri_gg_theme +
#   labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
#   scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
#   scale_fill_manual(values = subtype_colors)
# 
# 
# plt.gsam + plt.glass
# ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_classical.png',width=15,height=6.5)







## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [proneural] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.proneural = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::arrange(pid, resection)
  # %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_proneural.png',width=15,height=6.5)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [mesenchymal] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::arrange(pid, resection)
  # %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_mesenchymal.png',width=15,height=6.5)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [final?] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` ,`NMF:123456.PCA.SVM.class.R1`)) , by=c('pid'='pid')) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::pull(pid))) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(facet = ifelse( (!is.na(`NMF:123456.PCA.SVM.class.R1`) & `NMF:123456.PCA.SVM.class.R1` == "Classical") | (`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") , "Classical" , "?" ) ) %>%
  dplyr::mutate(facet = ifelse( (!is.na(`NMF:123456.PCA.SVM.class.R1`) & `NMF:123456.PCA.SVM.class.R1` == "Mesenchymal") | (`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") , "Mesenchymal" , facet ) ) %>%
  dplyr::mutate(facet = ifelse( (!is.na(`NMF:123456.PCA.SVM.class.R1`) & `NMF:123456.PCA.SVM.class.R1` == "Proneural") | (`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") , "Proneural" , facet ) ) %>%
  dplyr::mutate(facet = ifelse(is.na(`NMF:123456.PCA.SVM.class.R1`) & resection == "R2" , "no R1 present" , facet )) %>%
  dplyr::filter(facet != "no R1 present")



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  #geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  #geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title=NULL) +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors) +
  facet_grid(cols = vars(facet)) +
#  geom_text_repel(data = subset(plt, pid %in% c('CDA','CDD','CBP','CBR','CCZ') )) + 
  theme(strip.background = element_blank(), strip.text = element_blank())






ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_facet.pdf',width=16.7,height=6.5)




# KNN Bootstrapping ~ estimate by chance 

k = 20
df = data.frame()
for(p in 1:nrow(plt.paired) ) {
  #p = 21

  # for pair
  target = plt.paired[p,]
  
  
  # find k closest R1's
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                               (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(k, -ed)
  
  
  #plot(plt.single$`NMF:123456.PC1` , plt.single$`NMF:123456.PC2`, pch=19 , col = as.numeric(plt.single $`NMF:123456.PCA.SVM.class`) + 1 , cex=0.5)
  
  
  #nodes <- data.frame(`NMF:123456.PC1` = target$`NMF:123456.PC1.R1`,
  #                    `NMF:123456.PC2` = target$`NMF:123456.PC2.R1`,
  #                    type="start", check.names=F)
  nodes <- data.frame()
  for(i in 1:k) {
    neighbour <- knn[i,] 
    
    #lines(c(neighbour$`NMF:123456.PC1.R1`,neighbour$`NMF:123456.PC1.R2`), c(neighbour$`NMF:123456.PC2.R1` , neighbour$`NMF:123456.PC2.R2`)  )
    #points(neighbour$`NMF:123456.PC1.R2` , neighbour$`NMF:123456.PC2.R2` , pch=8,cex=0.6  )
    
    delta_PC1 = target$`NMF:123456.PC1.R1` - neighbour$`NMF:123456.PC1.R1`
    delta_PC2 = target$`NMF:123456.PC2.R1` - neighbour$`NMF:123456.PC2.R1`
    
    
    nodes <- rbind(nodes,
              data.frame(`NMF:123456.PC1` = neighbour$`NMF:123456.PC1.R2` + delta_PC1 ,
                         `NMF:123456.PC2` = neighbour$`NMF:123456.PC2.R2` + delta_PC2 ,
                         type = "node", check.names=F))


  }

  nodes <- nodes %>%
    dplyr::mutate(class.svm = predict(s150.pca.nmf.subtype.classifier.svm , newdata = nodes %>% dplyr::mutate(type=NULL)) )
  
  df <- rbind(df, 
              data.frame(n.cl  = nodes %>% dplyr::filter(class.svm == "Classical") %>% nrow(),
                   n.mes = nodes %>% dplyr::filter(class.svm == "Mesenchymal") %>% nrow(),
                   n.pn  = nodes %>% dplyr::filter(class.svm == "Proneural") %>% nrow() ) %>%
    dplyr::mutate(p.cl = n.cl /  rowSums(.),
                  p.mes = n.mes /  rowSums(.),
                  p.pn = n.pn /  rowSums(.) ) %>%
    dplyr::mutate(pid = target$pid)
  )
  
}


df <- df %>%
  #group_by(pid) %>% 
  #dplyr::summarise( p.cl = mean(n.cl),     p.mes = mean(n.mes), p.pn = mean(n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))



df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()



#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()


df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()



#df.2 %>% filter( dataset.R1 == 'GSAM' & `NMF:123456.PCA.SVM.status`  == "Stable") %>% nrow()

## Figure S1E





# Figure 2: STAR / Arrow plots ----



plt <- data.frame()
for(i in 1:nrow(plt.paired)) {
  slice <- plt.paired[i,]
  
  line <- data.frame(x = c(slice$`NMF:123456.PC1.n.R1`,
                           slice$`NMF:123456.PC1.n.R2`),
                    y = c(slice$`NMF:123456.PC2.n.R1`,
                          slice$`NMF:123456.PC2.n.R2`),
                    resection = c('R1','R2'),
                    group = slice$pid,
                    init.subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                    rec.subtype = slice$`NMF:123456.PCA.SVM.class.R2`)
  
  line2 <- line %>%
    dplyr::mutate(x = x - line$x[1] ) %>%
    dplyr::mutate(y = y - line$y[1] )
  
  plt <- plt %>%
    rbind(line2)
  
  rm(slice, line, line2)
}



m <- plt %>%
  dplyr::group_by(init.subtype, resection) %>%
  dplyr::summarise(x = mean(x), y=mean(y)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = "Mean")



ggplot(plt, aes(x=x, y=y, group=group, col=init.subtype)) +
  geom_line(alpha=0.75, lwd=1) +
  geom_point(data = plt %>% dplyr::filter(resection == "R2") , aes(fill=rec.subtype), col="black", pch=21, cex=1.8) +
  geom_line(data=m, col="gray20", lty=1, lwd=0.8, show.legend=F) +
  geom_point(data = m %>% dplyr::filter(resection == "R2") , fill = "white", col="black", pch=21, cex=2.1) +
  labs(x="Centered GITS space traversal", y="Centered GITS space traversal") +
  xlim(-max(c(abs(plt$x),abs(plt$y))),max(c(abs(plt$x),abs(plt$y)))) +
  ylim(-max(c(abs(plt$x),abs(plt$y))),max(c(abs(plt$x),abs(plt$y)))) +
  youri_gg_theme +
  facet_grid(cols = vars(init.subtype)) +
  coord_equal() +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = "Initial subtype", title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_fill_manual(values = subtype_colors, guide = guide_legend(title = "Sub-type recurrence", title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


#ggsave("output/figures/centered_eucledian_distance_traversal_pnc.png",height=10 * 1.1,width=4)
ggsave("output/figures/centered_eucledian_distance_traversal_pnc.pdf",width=10 * 1.1,height=4)


# KNN Bootstrapping ~ estimate by chance ----

# Bij deze bootstrapping zoeken we van de KNN de richtingen/angles
# En gebruiken bootstrappen we de eucledian distances van alle samples


# hoe berekenen we van 2 nodes de angle en verlengen/verkorten we deze
# R1 = ( 1,  1)
# R2 = (-2, -3)

change_length_line <- function(x1, y1, x2, y2, length_new) {
  # From the  line between points (x1, y1) , (x2 ,y2),
  # we want to create a new line with an identical angle
  # but of length `length_new`.
  
  dy <- y2 - y1
  dx <- x2 - x1
  
  #slope <- dy / dx
  angle <- atan2(dy , dx) # in rads
  
  length_x_new <- cos(angle) * length_new
  length_y_new <- sin(angle) * length_new
  
  x2_new <- x1 + length_x_new
  y2_new <- y1 + length_y_new
 
  return (data.frame(x = c(x1, x2_new) ,
             y = c(y1, y2_new) ,
             point = as.factor(c("initial start", "new end"))))
}



nn = round((nrow(plt.paired) - 1) / 10) # 9
bootstrap_n = 2 * nn #1 # 25 # bootstrap iterations per sample-pair

df.A = data.frame()
df.B = data.frame()


## A: test with KNN angles and all distances ----

# plot(c(-4,4), c(-4,4), type="n")
# points(1,1, pch=19)
# points(-2,-3, pch=19)
# lines(c(-2,1), c(-3,1))
# df <- change_length_line(1,1,-2,-3, 1)
# points(df[2,1] , df[2,2], pch=19,col="red")


# 
# plot(plt.single$`NMF:123456.PC1.n` ,
#      plt.single$`NMF:123456.PC2.n`, pch=19 , col = as.numeric(as.factor(plt.single$`NMF:123456.PCA.SVM.class`)) + 1 , cex=0.5)
# 



for(p in 1:nrow(plt.paired) ) { # for pair
  #p = 21 
  #p = 23
  target = plt.paired[p,]
  
  
  # find k closest R1's for angles in close proximity
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                                 (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(nn, -ed)
  

  #lines(c(target$`NMF:123456.PC1.n.R1`,target$`NMF:123456.PC1.n.R2`), c(target$`NMF:123456.PC2.n.R1` , target$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,0,0.15))
  #points(target$`NMF:123456.PC1.n.R1` , target$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,0,0.75)  )
  #points(target$`NMF:123456.PC1.n.R2` , target$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,0,0.75)  )
  
  
  for(i in 1:bootstrap_n) {
    #lines(c(knn$`NMF:123456.PC1.n.R1`,knn$`NMF:123456.PC1.n.R2`), c(knn$`NMF:123456.PC2.n.R1` , knn$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,1,0.15))
    #points(knn$`NMF:123456.PC1.n.R1` , knn$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,1,0.75)  )
    #points(knn$`NMF:123456.PC1.n.R2` , knn$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,1,0.75)  )
    

    # Take the angle from local subsampling
    neighbour <- knn %>% # RANDOM SHUFFLE `MET TERUGLEGGEN` want bootstrappen!!!
      dplyr::slice(sample(1 : dplyr::n() ) [1] )

    # take length from overall subsamples
    random_length <- plt.paired %>% # run once with `knn` and once with `plt.paired`?
      dplyr::slice(sample(1 : dplyr::n() )) %>%
      dplyr::pull(eucledian.dist) %>% 
      purrr::pluck(1)

    # generate new length using samples length and samples angle
    bootstrapped_line <- change_length_line(neighbour$`NMF:123456.PC1.n.R1`,
                                            neighbour$`NMF:123456.PC2.n.R1`,
                                            
                                            neighbour$`NMF:123456.PC1.n.R2`,
                                            neighbour$`NMF:123456.PC2.n.R2`,
                                            
                                            random_length)


    
    # fit to target's R1
    bootstrapped_line$x = bootstrapped_line$x + (target$`NMF:123456.PC1.n.R1` - neighbour$`NMF:123456.PC1.n.R1`)
    bootstrapped_line$y = bootstrapped_line$y + (target$`NMF:123456.PC2.n.R1` - neighbour$`NMF:123456.PC2.n.R1`)
    

    
    #lines(bootstrapped_line$x, bootstrapped_line$y, lwd=2, lty=3, col="red")
    #points(bootstrapped_line$x[1], bootstrapped_line$y[1],  col="red")
    #points(bootstrapped_line$x[2], bootstrapped_line$y[2],  col="red")
    
    
    # @todo UN/SCALE-NORMALISE this before classification
    bootstrapped_line <- bootstrapped_line %>%
      dplyr::mutate(x.orig = (x + attr(`scale.NMF:123456.PC1`,'scaled:center')) * attr(`scale.NMF:123456.PC1`,'scaled:scale') ) %>%
      dplyr::mutate(y.orig = (y + attr(`scale.NMF:123456.PC2`,'scaled:center')) * attr(`scale.NMF:123456.PC2`,'scaled:scale') )

    
    class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = bootstrapped_line %>%
              dplyr::filter(point == "new end") %>%
              dplyr::mutate(point = NULL, x = NULL , y=  NULL) %>%
              dplyr::rename(`NMF:123456.PC1` = x.orig) %>%
              dplyr::rename(`NMF:123456.PC2` = y.orig) )
    
    # still storing co-ordinates in normalised form
    df.A <- rbind(df.A, data.frame(x = bootstrapped_line$x[2],
               y = bootstrapped_line$y[2],
               pid = target$pid,
               neighbour = neighbour$pid,
               class = class ))

  }
}


## @todo check -
## plot original labels (normalised)
## plot new labels (normalised)
## this confirms if un-normalisation went o.k.

# plot(plt.single$`NMF:123456.PC1.n` ,
#      plt.single$`NMF:123456.PC2.n`, pch=19 , col = as.numeric(as.factor(plt.single$`NMF:123456.PCA.SVM.class`)) + 1 , cex=0.5)
# 
# points(df %>% dplyr::filter(class == "Mesenchymal") %>% dplyr::pull(x),
#        df %>% dplyr::filter(class == "Mesenchymal") %>% dplyr::pull(y),
#        pch = 1, col="green", cex=0.2)
# points(df %>% dplyr::filter(class == "Classical") %>% dplyr::pull(x),
#        df %>% dplyr::filter(class == "Classical") %>% dplyr::pull(y),
#        pch = 2, col="red", cex=0.2)
# points(df %>% dplyr::filter(class == "Proneural") %>% dplyr::pull(x),
#        df %>% dplyr::filter(class == "Proneural") %>% dplyr::pull(y),
#        pch = 3, cex=0.2, col="blue")



df.A.integrated <- df.A %>%
  dplyr::group_by(pid) %>% 
  dplyr::summarise(n.cl  = sum(class == "Classical"),
                   n.mes = sum(class == "Mesenchymal"),
                   n.pn  = sum(class == "Proneural")) %>%
  dplyr::mutate(p.cl  = n.cl  / (n.cl + n.mes + n.pn) ,
                p.mes = n.mes / (n.cl + n.mes + n.pn) ,
                p.pn  = n.pn  / (n.cl + n.mes + n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))



## B: test with KNN angles and KNN distances ----


for(p in 1:nrow(plt.paired) ) { # for pair
  #p = 21 
  #p = 23
  target = plt.paired[p,]
  
  
  # find k closest R1's for angles in close proximity
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                                 (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(nn, -ed)
  
  
  #lines(c(target$`NMF:123456.PC1.n.R1`,target$`NMF:123456.PC1.n.R2`), c(target$`NMF:123456.PC2.n.R1` , target$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,0,0.15))
  #points(target$`NMF:123456.PC1.n.R1` , target$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,0,0.75)  )
  #points(target$`NMF:123456.PC1.n.R2` , target$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,0,0.75)  )
  
  
  for(i in 1:bootstrap_n) {
    #lines(c(knn$`NMF:123456.PC1.n.R1`,knn$`NMF:123456.PC1.n.R2`), c(knn$`NMF:123456.PC2.n.R1` , knn$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,1,0.15))
    #points(knn$`NMF:123456.PC1.n.R1` , knn$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,1,0.75)  )
    #points(knn$`NMF:123456.PC1.n.R2` , knn$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,1,0.75)  )
    
    
    # Take the angle from local subsampling
    neighbour <- knn %>% # RANDOM SHUFFLE `MET TERUGLEGGEN` want bootstrappen!!!
      dplyr::slice(sample(1 : dplyr::n() ) [1] )
    
    # take length from overall subsamples
    random_length <- knn %>% # run once with `knn` and once with `plt.paired`?
      dplyr::slice(sample(1 : dplyr::n() )) %>%
      dplyr::pull(eucledian.dist) %>%
      purrr::pluck(1)
    
    # generate new length using samples length and samples angle
    bootstrapped_line <- change_length_line(neighbour$`NMF:123456.PC1.n.R1`,
                                            neighbour$`NMF:123456.PC2.n.R1`,
                                            
                                            neighbour$`NMF:123456.PC1.n.R2`,
                                            neighbour$`NMF:123456.PC2.n.R2`,
                                            
                                            random_length)
    
    
    
    # fit to target's R1
    bootstrapped_line$x = bootstrapped_line$x + (target$`NMF:123456.PC1.n.R1` - neighbour$`NMF:123456.PC1.n.R1`)
    bootstrapped_line$y = bootstrapped_line$y + (target$`NMF:123456.PC2.n.R1` - neighbour$`NMF:123456.PC2.n.R1`)
    
    
    
    #lines(bootstrapped_line$x, bootstrapped_line$y, lwd=2, lty=3, col="red")
    #points(bootstrapped_line$x[1], bootstrapped_line$y[1],  col="red")
    #points(bootstrapped_line$x[2], bootstrapped_line$y[2],  col="red")
    
    
    # @todo UN/SCALE-NORMALISE this before classification
    bootstrapped_line <- bootstrapped_line %>%
      dplyr::mutate(x.orig = (x + attr(`scale.NMF:123456.PC1`,'scaled:center')) * attr(`scale.NMF:123456.PC1`,'scaled:scale') ) %>%
      dplyr::mutate(y.orig = (y + attr(`scale.NMF:123456.PC2`,'scaled:center')) * attr(`scale.NMF:123456.PC2`,'scaled:scale') )
    
    
    class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = bootstrapped_line %>%
                       dplyr::filter(point == "new end") %>%
                       dplyr::mutate(point = NULL, x = NULL , y=  NULL) %>%
                       dplyr::rename(`NMF:123456.PC1` = x.orig) %>%
                       dplyr::rename(`NMF:123456.PC2` = y.orig) )
    
    # still storing co-ordinates in normalised form
    df.B <- rbind(df.B, data.frame(x = bootstrapped_line$x[2],
                                   y = bootstrapped_line$y[2],
                                   pid = target$pid,
                                   neighbour = neighbour$pid,
                                   class = class ))
    
  }
}




df.B.integrated <- df.B %>%
  dplyr::group_by(pid) %>% 
  dplyr::summarise(n.cl  = sum(class == "Classical"),
                   n.mes = sum(class == "Mesenchymal"),
                   n.pn  = sum(class == "Proneural")) %>%
  dplyr::mutate(p.cl  = n.cl  / (n.cl + n.mes + n.pn) ,
                p.mes = n.mes / (n.cl + n.mes + n.pn) ,
                p.pn  = n.pn  / (n.cl + n.mes + n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))


## C: alleen bootstrappen lengtes (Alle) ----

## D: alleen bootstrappen lengtes (binnen zelfde subtype R1) ----



## visualisation ----



help(curvedarrow)

curvedarrow <- function(from, to, lwd=2, lty=1, lcol="black", arr.col=lcol, 
                        arr.pos=0.5, curve=1, dr=0.01, endhead=FALSE, segment = c(0,1), ...)   {
  
  dpos  <- to-from
  angle <- atan(dpos[2]/dpos[1])*180/pi         # angle between both
  if (is.nan(angle)) return()
  mid   <- 0.5*(to+from)                        # midpoint of ellipsoid arrow
  dst   <- dist(rbind(to, from))                # distance from-to
  ry    <- curve*dst                            # small radius of ellepsoid
  aFrom<-0                                      # angle to and from
  aTo<-pi
  if ( from[1] <= to[1]) {
    aFrom <- pi
    aTo <- 2*pi
  }
  
  if (segment [1] != 0)
    From <- segment[1] * aTo + (1-segment[1]) * aFrom
  else
    From <- aFrom
  
  if (segment [2] != 1)
    To <- segment[2] * aTo + (1-segment[2]) * aFrom
  else
    To <- aTo
  
  meanpi <- arr.pos * aTo + (1-arr.pos) * aFrom
  if (endhead) To <- meanpi
  
  
  plotellipse(rx=dst/2,  ry=ry, mid=mid, angle=angle, from = From, to = To,
              lwd=lwd, lty=lty, lcol=lcol)
  ell <- getellipse(rx=dst/2, ry=ry, mid=mid, angle=angle,
                    from=1.001*meanpi, to=0.999*meanpi, dr= 0.002)       #Changed from -0.002
  Arrows(ell[1,1], ell[1,2], ell[nrow(ell),1], ell[nrow(ell),2],
         code=1, lcol=lcol, arr.col=arr.col, ...)
  curvedarrow <- c(ell[nrow(ell),1], ell[nrow(ell),2])
}



O.cc <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Classical" ) %>% nrow()
O.cm <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Mesenchymal" ) %>% nrow()
O.cp <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Proneural" ) %>% nrow()

O.mc <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Classical" ) %>% nrow()
O.mm <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Mesenchymal" ) %>% nrow()
O.mp <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Proneural" ) %>% nrow()

O.pc <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Classical" ) %>% nrow()
O.pm <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Mesenchymal" ) %>% nrow()
O.pp <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Proneural" ) %>% nrow()


A.cc <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
A.cm <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
A.cp <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

A.mm <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
A.mc <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
A.mp <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

A.pp <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)
A.pc <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
A.pm <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)


B.cc <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
B.cm <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
B.cp <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

B.mm <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
B.mc <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
B.mp <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

B.pp <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)
B.pc <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
B.pm <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)





par(mfrow=c(1,3))

### Orig ----
plot(c(-1,1)* 1.25 , c(-1,1) * 1.25 , type='n' , main="Original data")

curvedarrow( c(0 - 0.12, 0.66 ) , c(0 + 0.12, 0.66), curve = -1.5)# CL -> CL
text(0, 1.2, O.cc)
curvedarrow( c(0, 0.66) , c(-0.66, -0.66), curve = -0.001)#        CL -> MES
text(-0.25, -0.1, O.cm)
curvedarrow( c(-0.66, -0.66),  c(0, 0.66)  , curve = -0.2)#      MES -> CL
text(-0.75, 0.1, O.mc)

curvedarrow( c(0.66 - 0.12, -0.66 ) , c(0.66 + 0.12, -0.66), curve = 1.5)# PN -> PN
text(0.66, -1.2, O.pp)
curvedarrow( c(0, 0.66) , c(0.66, -0.66), curve = -0.2)#         CL -> PN
text(0.25, -0.1, O.pc)
curvedarrow( c(0.66, -0.66),  c(0, 0.66)  , curve = -0.001)#       PN -> CL
text(0.75, 0.1, O.cp)

curvedarrow( c(-0.66 - 0.12, -0.66 ) , c(-0.66 + 0.12, -0.66), curve = 1.5)# MES -> MES
text(-0.66, -1.2, O.mm)
curvedarrow( c(-0.66, -0.66) , c(0.66, -0.66), curve = -0.001)#    MES -> PL
text(0, -0.55, O.mp)
curvedarrow( c(0.66, -0.66) , c(-0.66, -0.66), curve = -0.2)#    PL -> MES
text(0, -1.05, O.pm)

# labels
points(0, 0.66, cex = 8, pch=19, col="white")
points(0, 0.66, cex = 8, col="black")
text(0, 0.66, "CL", col="black")

points(-0.66, -0.66, cex = 8, pch=19, col="white")
points(-0.66, -0.66, cex = 8, col="black")
text(-0.66, -0.66, "MES", col="black")

points(0.66, -0.66, cex = 8, pch=19, col="white")
points(0.66, -0.66, cex = 8, col="black")
text(0.66, -0.66, "PN", col="black")


### A ----

plot(c(-1,1)* 1.25 , c(-1,1) * 1.25 , type='n' , main="Bootstrap A (KNN angle, random dist)")

curvedarrow( c(0 - 0.12, 0.66 ) , c(0 + 0.12, 0.66), curve = -1.5)# CL -> CL
text(0, 1.2, A.cc)
curvedarrow( c(0, 0.66) , c(-0.66, -0.66), curve = -0.001)#        CL -> MES
text(-0.25, -0.1, A.cm)
curvedarrow( c(-0.66, -0.66),  c(0, 0.66)  , curve = -0.2)#      MES -> CL
text(-0.75, 0.1, A.mc)

curvedarrow( c(0.66 - 0.12, -0.66 ) , c(0.66 + 0.12, -0.66), curve = 1.5)# PN -> PN
text(0.66, -1.2, A.pp)
curvedarrow( c(0, 0.66) , c(0.66, -0.66), curve = -0.2)#         CL -> PN
text(0.25, -0.1, A.pc)
curvedarrow( c(0.66, -0.66),  c(0, 0.66)  , curve = -0.001)#       PN -> CL
text(0.75, 0.1, A.cp)

curvedarrow( c(-0.66 - 0.12, -0.66 ) , c(-0.66 + 0.12, -0.66), curve = 1.5)# MES -> MES
text(-0.66, -1.2, A.mm)
curvedarrow( c(-0.66, -0.66) , c(0.66, -0.66), curve = -0.001)#    MES -> PL
text(0, -0.55, A.mp)
curvedarrow( c(0.66, -0.66) , c(-0.66, -0.66), curve = -0.2)#    PL -> MES
text(0, -1.05, A.pm)


# labels
points(0, 0.66, cex = 8, pch=19, col="white")
points(0, 0.66, cex = 8, col="black")
text(0, 0.66, "CL", col="black")

points(-0.66, -0.66, cex = 8, pch=19, col="white")
points(-0.66, -0.66, cex = 8, col="black")
text(-0.66, -0.66, "MES", col="black")

points(0.66, -0.66, cex = 8, pch=19, col="white")
points(0.66, -0.66, cex = 8, col="black")
text(0.66, -0.66, "PN", col="black")



### B ----

plot(c(-1,1)* 1.25 , c(-1,1) * 1.25 , type='n' , main="Bootstrap B (KNN angle, KNN dist)")

curvedarrow( c(0 - 0.12, 0.66 ) , c(0 + 0.12, 0.66), curve = -1.5)# CL -> CL
text(0, 1.2, B.cc)
curvedarrow( c(0, 0.66) , c(-0.66, -0.66), curve = -0.001)#        CL -> MES
text(-0.25, -0.1, B.cm)
curvedarrow( c(-0.66, -0.66),  c(0, 0.66)  , curve = -0.2)#      MES -> CL
text(-0.75, 0.1, B.mc)

curvedarrow( c(0.66 - 0.12, -0.66 ) , c(0.66 + 0.12, -0.66), curve = 1.5)# PN -> PN
text(0.66, -1.2, B.pp)
curvedarrow( c(0, 0.66) , c(0.66, -0.66), curve = -0.2)#         CL -> PN
text(0.25, -0.1, B.pc)
curvedarrow( c(0.66, -0.66),  c(0, 0.66)  , curve = -0.001)#       PN -> CL
text(0.75, 0.1, B.cp)

curvedarrow( c(-0.66 - 0.12, -0.66 ) , c(-0.66 + 0.12, -0.66), curve = 1.5)# MES -> MES
text(-0.66, -1.2, B.mm)
curvedarrow( c(-0.66, -0.66) , c(0.66, -0.66), curve = -0.001)#    MES -> PL
text(0, -0.55, B.mp)
curvedarrow( c(0.66, -0.66) , c(-0.66, -0.66), curve = -0.2)#    PL -> MES
text(0, -1.05, B.pm)


# labels
points(0, 0.66, cex = 8, pch=19, col="white")
points(0, 0.66, cex = 8, col="black")
text(0, 0.66, "CL", col="black")

points(-0.66, -0.66, cex = 8, pch=19, col="white")
points(-0.66, -0.66, cex = 8, col="black")
text(-0.66, -0.66, "MES", col="black")

points(0.66, -0.66, cex = 8, pch=19, col="white")
points(0.66, -0.66, cex = 8, col="black")
text(0.66, -0.66, "PN", col="black")


# Figure 2: Proportional distances ----


# 
# plot(plt.single$`NMF:123456.PC1` ,
#      plt.single$`NMF:123456.PC2`, pch=19 , col = as.numeric(as.factor(plt.single$`NMF:123456.PCA.SVM.class`)) + 1 , cex=0.5)



df <- data.frame()
res <- 150 # resolution
for(p in 1:nrow(plt.paired) ) { # for pair
  
  #p <- which(plt.paired$pid == "AZF")
  
  target <- plt.paired[p,]
  d <- sqrt((target$`NMF:123456.PC1.R1` - target$`NMF:123456.PC1.R2`)^2 +
    (target$`NMF:123456.PC2.R1` - target$`NMF:123456.PC2.R2`)^2) # dist
  
  # lines(
  #   c(target$`NMF:123456.PC1.R1`, target$`NMF:123456.PC1.R2`),
  #   c(target$`NMF:123456.PC2.R1`, target$`NMF:123456.PC2.R2`)
  # )
  
  df.classification <- data.frame('pid' = target$pid,
                                  `NMF:123456.PC1` = target$`NMF:123456.PC1.R1`,
                                  `NMF:123456.PC2` = target$`NMF:123456.PC2.R1`,
                                  frac = 0)

  for(r in 1:res) {
    pct <- r / (res+1)

    # generate new length using samples length and samples angle
    bootstrapped_line <- change_length_line(target$`NMF:123456.PC1.R1`,
                                            target$`NMF:123456.PC2.R1`,
                                            
                                            target$`NMF:123456.PC1.R2`,
                                            target$`NMF:123456.PC2.R2`,
                                            
                                            d * pct)
    
    #lines(
    #  c(bootstrapped_line$x[1], bootstrapped_line$x[2]),
    #  c(bootstrapped_line$y[1], bootstrapped_line$y[2])
    #)
    #points(bootstrapped_line$x[2], bootstrapped_line$y[2],col=as.numeric(class) + 1, lwd=2,cex=1.6)
    
    df.classification <- df.classification %>% 
      rbind(data.frame(pid = target$pid,
                       `NMF:123456.PC1` = bootstrapped_line$x[2],
                       `NMF:123456.PC2` = bootstrapped_line$y[2],
                       frac = pct) )
  }

  df.classification <- df.classification %>%
    rbind(data.frame(pid = target$pid,
                     `NMF:123456.PC1` = target$`NMF:123456.PC1.R2`,
                     `NMF:123456.PC2` = target$`NMF:123456.PC2.R2`,
                     frac = 1))
  
  df.classification$class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = df.classification %>% 
                                       dplyr::mutate(frac=NULL, pid=NULL) )
  
  df <- df %>% rbind(df.classification)

}

df <- df %>%
  dplyr::mutate(pid = as.factor(pid)) %>%
  dplyr::mutate(class = as.factor(class))
  
# df %>%
#   dplyr::filter(pid == levels(df$pid)[2]) # == AAC






df.paired <- data.frame()
for(ppid in unique(df$pid)) {
  #ppid = "EBV"
  target <- df %>%
    dplyr::filter(`pid` == ppid) %>%
    dplyr::arrange(frac)
  
  l <- length(unique(target$class))
  
  out <- data.frame(pid = ppid,
                   init = 0,
                   init.class = target$class[1],
                   split = NA,
                   end = 1,
                   end.class = target$class[nrow(target)])

  for(k in 1:nrow(target)) {
    slice = target[k,]
    
    if(is.na(out$split) & slice$class != out$init.class) {
      out$split = (target[k,]$frac + target[k - 1,]$frac) / 2
      #print("!!")
    }
  }
  
  df.paired <- rbind(out, df.paired)
}
df.paired <- df.paired %>%
  dplyr::mutate(init = NULL,
                init.class = NULL,
                end = NULL,
                end.class = NULL) %>%
  dplyr::left_join(plt.paired , by = c('pid' = 'pid')) %>%
  dplyr::left_join(gsam.patient.metadata %>% dplyr::select(studyID, HM), by = c('pid' = 'studyID') ) %>%
  dplyr::mutate(HM = ifelse(is.na(HM),"N/A",HM))





# maak ggplot
plt <- data.frame()
for(i in 1:nrow(df.paired)){
  slice = df.paired[i,]
  

  if(is.na( slice$split )) {
    dist <- slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-0.1,-0.1),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("lpad","lpad")))
    }
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="lpad")
    )
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="lpad")
    )
    
    
    dist <- - slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-2.8,-2.8),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("rpad","rpad")))
    }
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = -slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="rpad")
    )
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="rpad")
    )
    
    
  } else  {
    dist <- slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-0.1,-0.1),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("lpad","lpad")))
    }
    
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = (1 - slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos = (1 - slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos =  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    
    
    
    dist <- (slice$split) * slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-2.8,-2.8),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("rpad","rpad")))
    }
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = - (1 - slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos =  (slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    
  }
}



p1 <- ggplot(plt %>% dplyr::filter(type == "lpad"), aes(x = pos , y = reorder(pid, d), col = subtype, group=group)) +
  geom_point(cex=1) +
  geom_path(lwd = 1.2) + 
  facet_grid(rows = vars(subtype.init), cols=vars(type), scales = "free", space="free_y") + 
  labs(x = "Traversal distance GITS space",
       y=NULL,
       col = "Distance within Subtype space") +
  theme_light() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'bottom',
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    axis.text.y = element_text(size = 6)
  ) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_color_manual(values = c(subtype_colors, 'Hyper-mutant' = 'gray40'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


p2 <- ggplot(plt %>% dplyr::filter(type == "rpad"), aes(x = pos , y = reorder(pid, d), col = subtype, group=group)) +
  geom_point(cex=1.1) +
  geom_path(lwd = 1.2) + 
  facet_grid(rows = vars(subtype.init), cols=vars(type), scales = "free", space="free_y") + 
  labs(x = "Traversal distance GITS space",
       y=NULL,
       col = "Distance within Subtype space") +
  theme_light() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'bottom',
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    axis.text.y = element_text(size = 6)
  ) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_color_manual(values = c(subtype_colors, 'Hyper-mutant' = 'gray40'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))

p1 + p2




#ggsave("output/figures/eucledian_distance_traversal_pnc.png",height=8*1.3,width=8*1.4)
ggsave("output/figures/eucledian_distance_traversal_pnc.pdf",height=8*1.3,width=8*2.0)



## santoesha style plot ----




plt.prop <- data.frame()

tmp <-  df %>% dplyr::left_join(plt.paired %>%
                          dplyr::select(c('pid',`NMF:123456.PCA.SVM.class.R1`, eucledian.dist)),
                          by=c('pid'='pid'))

for(ppid in unique(tmp$pid)) {
  print(ppid)
  
  
  tmp.2 <- tmp %>% dplyr::filter(pid == ppid)

  
  last_class <- tmp.2[1,]$class
  last_frac <- tmp.2[1,]$frac
  
  for(i in 2:nrow(tmp.2))  {
    #print(i)
    if(tmp.2[i - 1,]$class != tmp.2[i,]$class ) {
       split <-  (tmp.2[i -1 ,]$frac + tmp.2[i  ,]$frac) / 2
      
       #print(paste0("(",i,") ", last_class, " -> ", tmp.2[i - 1,]$class))
      
       plt.prop <- plt.prop %>%
         rbind(data.frame(
           pid = tmp.2[i,]$pid,
           
           from_class = last_class,
           to_class = tmp.2[i - 1,]$class,

           from_frac = last_frac  ,
           to_frac = split ,
           delta_frac =  split - last_frac,
           steps = (split - last_frac) / (1/151),
           
           eucl_dist = tmp.2[i,]$eucledian.dist ,
           delta_eucl_dist = tmp.2[i,]$eucledian.dist * (split - last_frac)
         ))
      
      last_class <- tmp.2[i - 1 ,]$class
      last_frac <- split
    }
    
    if ( tmp.2[i,]$frac == 1.0) {
      #print(paste0("(",i,") ", last_class, " -> ", tmp.2[i - 1,]$class))
      
      plt.prop <- plt.prop %>%
        rbind(data.frame(
          pid = tmp.2[i,]$pid,
          
          from_class = last_class,
          to_class = tmp.2[i ,]$class,
          
          from_frac = last_frac  ,
          to_frac = tmp.2[i ,]$frac ,
          delta_frac =  tmp.2[i ,]$frac   - last_frac,
          steps = (tmp.2[i ,]$frac   - last_frac) / (1/151),
          
          eucl_dist = tmp.2[i,]$eucledian.dist ,
          delta_eucl_dist = tmp.2[i,]$eucledian.dist * (tmp.2[i ,]$frac   - last_frac)
        ))
      
    }
  }
  
  #tmp.2
  #plt.prop
  
  test <- plt.prop %>%
    dplyr::filter(pid == ppid)
  
  stopifnot(sum(test$steps ) == 151)
  stopifnot(sum(test$delta_frac) == 1)
}











links <- data.frame() %>%
  rbind(data.frame( from = "CL", to = "CL" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Classical" & to_class == "Classical") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "CL", to = "MES" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Classical" & to_class == "Mesenchymal") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "CL", to = "PN" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Classical" & to_class == "Proneural") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%

  rbind(data.frame( from = "MES", to = "CL" , weight = 
                    plt.prop %>% dplyr::filter(from_class == "Mesenchymal" & to_class == "Classical") %>%
                    dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "MES", to = "MES" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Mesenchymal" & to_class == "Mesenchymal") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "MES", to = "PN" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Mesenchymal" & to_class == "Proneural") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%

  rbind(data.frame( from = "PN", to = "CL" , weight = 
                    plt.prop %>% dplyr::filter(from_class == "Proneural" & to_class == "Classical") %>%
                    dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "PN", to = "MES" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Proneural" & to_class == "Mesenchymal") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "PN", to = "PN" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Proneural" & to_class == "Proneural") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum ))

links[links$from == "CL",]$weight <- links[links$from == "CL",]$weight / (plt.prop %>% dplyr::filter(from_class == "Classical" & from_frac == 0) %>% nrow())
links[links$from == "MES",]$weight <- links[links$from == "MES",]$weight / (plt.prop %>% dplyr::filter(from_class == "Classical" & from_frac == 0) %>% nrow())
links[links$from == "PN",]$weight <- links[links$from == "PN",]$weight / (plt.prop %>% dplyr::filter(from_class == "Classical" & from_frac == 0) %>% nrow())

links <- links %>%
  dplyr::mutate(weight = round(weight,2))

  
    #dplyr::mutate(weight=round(weight/sum(weight),3) * 100)

#Real values
#links <- data.frame(from=c(rep(c("CL"),3),rep(c("MES"),3),rep(c("PN"),3)), to=c(rep(c("CL","MES","PN"),3)))
#links <- links %>% dplyr::group_by(from) %>% dplyr::mutate(thickness=round(weight/sum(weight),2))

nodes <- data.frame(id=c("CL","MES","PN"),
                    subtype=c("Classical","Mesenchymal","Proneural"),
                    size=c(10,8,9))
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
net <- simplify(net, remove.multiple = F, remove.loops = F)
#E(net)$width <- E(net)$thickness*6
#E(net)$size <- E(net)$size*0.75

plot(net,
     edge.arrow.size=.4,
     diag=T,
     edge.label = E(net)$weight,
     arrow.mode=3,
     edge.curved=.2,
     edge.color="gray70",
     vertex.color="coral",
     vertex.label.font=2,
     edge.label.font=1,
     edge.label.cex=1.3,
     edge.label.color="blue",
)


# metadata + some plots ----

tmp.1 <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::select(aliquot_barcode,
                GBM.transcriptional.subtype.Synapse.2021,
                ssGSEA.Synapse.subtype.2022,
                ssGSEA.Synapse.subtype.2022.Proneural.pval,
                ssGSEA.Synapse.subtype.2022.Classical.pval,
                ssGSEA.Synapse.subtype.2022.Mesenchymal.pval,
                ssGSEA.Synapse.subtype.2022.Proneural.enrichment_score,
                ssGSEA.Synapse.subtype.2022.Classical.enrichment_score,
                ssGSEA.Synapse.subtype.2022.Mesenchymal.enrichment_score,
                
                ssGSEA.2022.subtype,
                ssGSEA.2022.Proneural.enrichment_score,
                ssGSEA.2022.Classical.enrichment_score,
                ssGSEA.2022.Mesenchymal.enrichment_score,
                ssGSEA.2022.Proneural_pval,
                ssGSEA.2022.Classical_pval,
                ssGSEA.2022.Mesenchymal_pval
  )  |> 
  dplyr::mutate(ssGSEA.2022.subtype = as.character(ssGSEA.2022.subtype)) |> 
  dplyr::rename(ssGSEA.2022.subtype.from.glass = ssGSEA.2022.subtype,
                ssGSEA.2022.Proneural.enrichment_score.from.glass = ssGSEA.2022.Proneural.enrichment_score,
                ssGSEA.2022.Classical.enrichment_score.from.glass = ssGSEA.2022.Classical.enrichment_score,
                ssGSEA.2022.Mesenchymal.enrichment_score.from.glass = ssGSEA.2022.Mesenchymal.enrichment_score,
                ssGSEA.2022.Proneural_pval.from.glass = ssGSEA.2022.Proneural_pval,
                ssGSEA.2022.Classical_pval.from.glass = ssGSEA.2022.Classical_pval,
                ssGSEA.2022.Mesenchymal_pval.from.glass = ssGSEA.2022.Mesenchymal_pval)

tmp.2 <-    gsam.rna.metadata |>
  dplyr::select(sid,
                gliovis.svm_call, gliovis.knn_call, gliovis.gsea_call, gliovis.equal_call, gliovis.majority_call, 
                svm.Classical, svm.Mesenchymal, svm.Proneural,
                knn.Classical, knn.Mesenchymal, knn.Proneural,
                
                gliovis.ssGSEA.Proneural.score, gliovis.ssGSEA.Classical.score, gliovis.ssGSEA.Mesenchymal.score,
                gliovis.ssGSEA.Proneural, gliovis.ssGSEA.Classical, gliovis.ssGSEA.Mesenchymal,
                
                
                ssGSEA.2022.subtype,
                ssGSEA.2022.Proneural.enrichment_score,
                ssGSEA.2022.Classical.enrichment_score,
                ssGSEA.2022.Mesenchymal.enrichment_score,
                ssGSEA.2022.Proneural_pval,
                ssGSEA.2022.Classical_pval,
                ssGSEA.2022.Mesenchymal_pval
  ) |> 
  dplyr::mutate(ssGSEA.2022.subtype = as.character(ssGSEA.2022.subtype)) |> 
  dplyr::rename(ssGSEA.2022.subtype.from.gsam = ssGSEA.2022.subtype,
                ssGSEA.2022.Proneural.enrichment_score.from.gsam = ssGSEA.2022.Proneural.enrichment_score,
                ssGSEA.2022.Classical.enrichment_score.from.gsam = ssGSEA.2022.Classical.enrichment_score,
                ssGSEA.2022.Mesenchymal.enrichment_score.from.gsam = ssGSEA.2022.Mesenchymal.enrichment_score,
                ssGSEA.2022.Proneural_pval.from.gsam = ssGSEA.2022.Proneural_pval,
                ssGSEA.2022.Classical_pval.from.gsam = ssGSEA.2022.Classical_pval,
                ssGSEA.2022.Mesenchymal_pval.from.gsam = ssGSEA.2022.Mesenchymal_pval)


tmp.out <- tmp.out |> 
  dplyr::left_join(tmp.1, by=c('sid'='aliquot_barcode'),suffix=c('','')) |> 
  dplyr::left_join(tmp.2, by=c('sid'='sid'),suffix=c('','')) |> 
  dplyr::mutate(ssGSEA.2022.subtype = ifelse(is.na(ssGSEA.2022.subtype.from.gsam), ssGSEA.2022.subtype.from.glass, ssGSEA.2022.subtype.from.gsam)) |> 
  dplyr::mutate(ssGSEA.2022.Proneural.enrichment_score = ifelse(is.na(ssGSEA.2022.Proneural.enrichment_score.from.gsam), ssGSEA.2022.Proneural.enrichment_score.from.glass, ssGSEA.2022.Proneural.enrichment_score.from.gsam)) |> 
  dplyr::mutate(ssGSEA.2022.Classical.enrichment_score = ifelse(is.na(ssGSEA.2022.Classical.enrichment_score.from.gsam), ssGSEA.2022.Classical.enrichment_score.from.glass, ssGSEA.2022.Classical.enrichment_score.from.gsam)) |> 
  dplyr::mutate(ssGSEA.2022.Mesenchymal.enrichment_score = ifelse(is.na(ssGSEA.2022.Mesenchymal.enrichment_score.from.gsam), ssGSEA.2022.Mesenchymal.enrichment_score.from.glass, ssGSEA.2022.Mesenchymal.enrichment_score.from.gsam)) |> 
  dplyr::mutate(ssGSEA.2022.Proneural_pval = ifelse(is.na(ssGSEA.2022.Proneural_pval.from.gsam), ssGSEA.2022.Proneural_pval.from.glass, ssGSEA.2022.Proneural_pval.from.gsam)) |> 
  dplyr::mutate(ssGSEA.2022.Classical_pval = ifelse(is.na(ssGSEA.2022.Classical_pval.from.gsam), ssGSEA.2022.Classical_pval.from.glass, ssGSEA.2022.Classical_pval.from.gsam)) |> 
  dplyr::mutate(ssGSEA.2022.Mesenchymal_pval = ifelse(is.na(ssGSEA.2022.Mesenchymal_pval.from.gsam), ssGSEA.2022.Mesenchymal_pval.from.glass, ssGSEA.2022.Mesenchymal_pval.from.gsam)) |> 
  dplyr::select(-contains(".from.g")) |> 
  dplyr::mutate(dataset = ifelse(grepl("GLSS|TCGA", sid),"GLASS","G-SAM"))


rm(tmp.1, tmp.2)





### ssGSEA x 150 ----

# combined 
plt <- tmp.out |>
  dplyr::filter(!is.na(ssGSEA.2022.subtype)) 

p1 = ggplot(plt , aes(x=`NMF:150:1`, y=ssGSEA.2022.Classical.enrichment_score,
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p2 = ggplot(plt , aes(x=`NMF:150:2`, y=`ssGSEA.2022.Proneural.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p3 = ggplot(plt , aes(x=`NMF:150:3`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw()

p1 + p2 + p3


#ggsave("output/figures/")



# # GLASS
# plt <- tmp.out |>
#   dplyr::filter(!is.na(ssGSEA.Synapse.subtype.2022))
# ggplot(plt , aes(x=`NMF:150:1`, y=ssGSEA.Synapse.subtype.2022.Classical.enrichment_score, col=ssGSEA.Synapse.subtype.2022)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:2`, y=ssGSEA.Synapse.subtype.2022.Proneural.enrichment_score, col=ssGSEA.Synapse.subtype.2022)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:3`, y=ssGSEA.Synapse.subtype.2022.Mesenchymal.enrichment_score, col=ssGSEA.Synapse.subtype.2022)) +
#   geom_point() +
#   theme_bw()
# 
# 
# 
# # GSAM
# plt <- tmp.out |>
#   dplyr::filter(!is.na(gliovis.ssGSEA.Mesenchymal.score))
# ggplot(plt , aes(x=`NMF:150:1`, y=gliovis.ssGSEA.Classical.score, col=gliovis.gsea_call)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:2`, y=gliovis.ssGSEA.Proneural.score, col=gliovis.gsea_call)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:150:3`, y=gliovis.ssGSEA.Mesenchymal.score, col=gliovis.gsea_call)) +
#   geom_point() +
#   theme_bw()


### ssGSEA x 7k ----



### ssGSEA x 150 ----

# combined 
plt <- tmp.out |>
  dplyr::filter(!is.na(ssGSEA.2022.subtype)) 

p1 = ggplot(plt , aes(x=`NMF:7k:1`, y=ssGSEA.2022.Classical.enrichment_score,
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p2 = ggplot(plt , aes(x=`NMF:7k:2`, y=`ssGSEA.2022.Proneural.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() 
p3 = ggplot(plt , aes(x=`NMF:7k:4`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, 
                      shape=dataset, col=ssGSEA.2022.subtype, fill=ssGSEA.2022.subtype)) +
  ggpubr::stat_cor(aes(shape=NULL, col=NULL, fill=NULL), method = "pearson") +
  geom_point(cex=2, col="black") +
  scale_fill_manual(values = subtype_colors_ext, guide="none") +
  scale_color_manual(values = subtype_colors_ext) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw()

p1 + p2 + p3



# GLASS - match NMF's w/ subtypes
plt <- tmp.out |>
  dplyr::filter(!is.na(ssGSEA.2022.subtype))

# GLASS - match NMF's w/ subtypes
ggplot(plt , aes(x=`NMF:7k:1`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
  geom_point() +
  theme_bw()
# ggplot(plt , aes(x=`NMF:7k:2`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:3`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:4`, y=ssGSEA.2022.Classical.enrichment_score, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()

ggplot(plt , aes(x=`NMF:7k:2`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
  geom_point() +
  theme_bw()
# ggplot(plt , aes(x=`NMF:7k:1`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:3`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:4`, y=`ssGSEA.2022.Proneural.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()

ggplot(plt , aes(x=`NMF:7k:4`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
  geom_point() +
  theme_bw()
# ggplot(plt , aes(x=`NMF:7k:1`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:2`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()
# ggplot(plt , aes(x=`NMF:7k:3`, y=`ssGSEA.2022.Mesenchymal.enrichment_score`, col=ssGSEA.2022.subtype)) +
#   geom_point() +
#   theme_bw()


### ssGSEA [glass] ----

### PCA ----


tmp.pca <- tmp.out |>
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022)) |> 
  dplyr::select(sid, 
                GBM.transcriptional.subtype.Synapse.Proneural.enrichment_score,
                GBM.transcriptional.subtype.Synapse.Classical.enrichment_score,
                GBM.transcriptional.subtype.Synapse.Mesenchymal.enrichment_score
  ) |> 
  tibble::column_to_rownames('sid') |> 
  prcomp()


tmp.out <- tmp.out |> 
  dplyr::left_join(
    tmp.pca |> 
      purrr::pluck('x') |> 
      as.data.frame() |> 
      dplyr::rename_with( ~ paste0("ssGSEA:GLASS:", .x)) |> 
      tibble::rownames_to_column('sid'),
    by=c('sid'='sid'),suffix=c('','')
  )


### more ----



plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p1 = ggplot(plt, aes(x=`ssGSEA:PC1`, y=`ssGSEA:PC2`,col=GBM.transcriptional.subtype.Synapse.2022)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p2 = ggplot(plt, aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`,col=GBM.transcriptional.subtype.Synapse.2022)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p3 = ggplot(plt, aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`,col=GBM.transcriptional.subtype.Synapse.2022)) +
  geom_point() +
  theme_bw()


plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p4 = ggplot(plt, aes(x=`ssGSEA:PC1`, y=`ssGSEA:PC2`,col=`NMF:150:membership`)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p5 = ggplot(plt, aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`,col=`NMF:150:membership`)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |> dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2022))
p6 = ggplot(plt, aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`,col=`NMF:150:membership`)) +
  geom_point() +
  theme_bw()


(p1 + p2 + p3) /
  (p4 + p5 + p6) 

ggplot(plt , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()

plt <- tmp.out |>
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021))
ggplot(plt , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`GBM.transcriptional.subtype.Synapse.2021`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`, col=`GBM.transcriptional.subtype.Synapse.2021`)) +
  geom_point() +
  theme_bw()





ggplot(plt , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()

p2 = ggplot(tmp.out , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=ds)) +
  geom_point() +
  theme_bw()

ggplot(plt , aes(x=`NMF:150:PC1`, y=`NMF:150:PC2`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()


ggplot(plt , aes(x=`NMF:150:1`, y=log(1 - `GBM.transcriptional.subtype.Synapse.Classical.pval` * 0.999), col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:2`, y=`GBM.transcriptional.subtype.Synapse.Proneural.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:150:3`, y=`GBM.transcriptional.subtype.Synapse.Mesenchymal.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()

ggplot(plt , aes(x=`NMF:7k:1`, y=`GBM.transcriptional.subtype.Synapse.Classical.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:2`, y=`GBM.transcriptional.subtype.Synapse.Classical.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:3`, y=`GBM.transcriptional.subtype.Synapse.Classical.pval`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()



ggplot(plt , aes(x=`NMF:7k:1`, y=log(GBM.transcriptional.subtype.Synapse.Classical.enrichment_score), col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:2`, y=`GBM.transcriptional.subtype.Synapse.Proneural.enrichment_score`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()
ggplot(plt , aes(x=`NMF:7k:3`, y=`GBM.transcriptional.subtype.Synapse.Mesenchymal.enrichment_score`, col=`GBM.transcriptional.subtype.Synapse.2022`)) +
  geom_point() +
  theme_bw()




ggplot(tmp.out, aes(x=`NMF:7k:PC1`, y=`NMF:150:PC1`,col=ds)) + 
  geom_point() +
  theme_bw()

ggplot(tmp.out, aes(x=`NMF:7k:PC2`, y=`NMF:150:PC2`,col=ds)) + 
  geom_point() +
  theme_bw()


ggplot(tmp.out, aes(x=`NMF:7k:1`, y=`NMF:150:3`,col=ds)) + 
  geom_point() +
  theme_bw()

ggplot(tmp.out, aes(x=`NMF:7k:2`, y=`NMF:150:3`,col=ds)) + 
  geom_point() +
  theme_bw()

ggplot(tmp.out, aes(x=`NMF:7k:3`, y=`NMF:150:1`, col=ds)) + 
  geom_point() +
  theme_bw()


### some plots ----

# 
# plt <- data.frame(PC= paste0("PC",1:3),
#                   var_explained=(tmp.pca$sdev)^2/sum((tmp.pca$sdev)^2)) %>%
#   dplyr::mutate(label = paste0(round(var_explained * 100,1),"%"))
# 
# ggplot(plt, aes(x=PC,y=var_explained, group=1, label=label))+
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   geom_col(fill="gray60", col="black", lwd=0.4) +
#   geom_line()+
#   geom_label(y = 0.35, size=4) + 
#   geom_point(size=4)+
#   labs(title="Scree plot: PCA on NMF (for p=3)", x= NULL, y="Variance in NMF V-matrix explained") + 
#   theme_bw()



p1 = ggplot(tmp.out , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=`NMF:7k:membership`)) +
  geom_point() +
  theme_bw()

p2 = ggplot(tmp.out , aes(x=`NMF:7k:PC1`, y=`NMF:7k:PC2`, col=ds)) +
  geom_point() +
  theme_bw()

p1+p2


PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x, varnames='NMF meta-feature 1')
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation) %>%
    dplyr::mutate(varnames = gsub('NMF:123456.','NMF meta-feature ', varnames))
  
  #print(head(datapc$varnames))
  #print(head(datapc))
  
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  
  ggplot(data, aes(x=PC1, y=PC2, col=varnames)) +
    geom_point(size=2.5,  col='black', fill='gray', pch=21) +
    coord_equal() +
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.4,"cm")), lwd=1.3, alpha=0.7) +
    geom_label_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3.5, vjust=1, show.legend=F) +
    youri_gg_theme +
    labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", fill = "Sub-type (GlioVis)") +
    scale_color_manual(values = subtype_colors_nmfm, labels = c('NMF meta-feature 2'='Classical',
                                                                'NMF meta-feature 1'='Mesenchymal',
                                                                'NMF meta-feature 3'='Proneural')) +
    labs(col = "Associated sub-type")
  #+ theme(legend.position = "none")
}


PCbiplot(p)




## + SVM classification ----
# Aanzienlijk lekkerdere fit dan LDA




if(file.exists("output/tables/gsam_nmf_lda_data.txt")) {
  a = read.table("output/tables/gsam_nmf_lda_data.txt") 
  
  colnames(a)
}  

s150.pca.nmf.subtype.classifier.svm <- readRDS('tmp/s150.pca.nmf.subtype.classifier.svm.Rds')
# s150.pca.nmf.subtype.classifier.svm <- svm(x = plt.single  %>%
#                                              #dplyr::filter(dataset == "GSAM") %>%
#                                              dplyr::select(c('sid' ,'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
#                                              tibble::column_to_rownames('sid') ,
#                                            y = plt.single %>% 
#                                              #dplyr::filter(dataset == "GSAM") %>%
#                                              dplyr::pull(subtype.public), 
#                                            scale = F,
#                                            #type = "C-classification",
# 
#                                            kernel = 'linear',
#                                            tolerance = 0.0001,
#                                            cost = 3
#                                            
#                                            #,probability = T # worse fit, but handy values
#                                            )
#saveRDS(s150.pca.nmf.subtype.classifier.svm, 'tmp/s150.pca.nmf.subtype.classifier.svm.Rds')


# re-fit samples

#attr(train.pred, 'decision.values')
plt.single <- plt.single %>% dplyr::left_join(
  data.frame(
    sid = plt.single$sid,
    
    'NMF:123456.PCA.SVM.class' = as.character(
      predict(object = s150.pca.nmf.subtype.classifier.svm , newdata =  (plt.single  %>%
                                                                           dplyr::select(c('sid' ,'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
                                                                           tibble::column_to_rownames('sid')),
              decision.values = T,
              tolerance = 0.0001,
              cost = 3
      )),
    check.names = F
  ) , by = c('sid'='sid') )


plt.single %>% dplyr::filter(dataset == "GSAM") %>% dplyr::mutate(err = subtype.public != `NMF:123456.PCA.SVM.class`) %>% dplyr::pull(err) %>% summary
plt.single %>% dplyr::filter(dataset != "GSAM") %>% dplyr::mutate(err = subtype.public != `NMF:123456.PCA.SVM.class`) %>% dplyr::pull(err) %>% summary
plt.single %>% dplyr::mutate(err = subtype.public != `NMF:123456.PCA.SVM.class`) %>% dplyr::pull(err) %>% summary


### :: redo to estimate probabilities ----
#
# This prediction performs considerable worse
# but probabilities may be useful to other analysis
#

tmp.fit <- svm(
  x = plt.single  %>%
    dplyr::select(c(
      'sid' , 'NMF:123456.PC1', 'NMF:123456.PC2'
    )) %>%
    tibble::column_to_rownames('sid') ,
  y = plt.single %>%
    dplyr::pull(subtype.public),
  scale = F,
  kernel = 'linear',
  tolerance = 0.0001,
  cost = 3,
  probability = T)



tmp <- predict(
  object = tmp.fit ,
  newdata = (
    plt.single %>% dplyr::select(c(
      'sid' , 'NMF:123456.PC1', 'NMF:123456.PC2'
    )) %>% tibble::column_to_rownames('sid')
  ),
  decision.values = T,
  tolerance = 0.0001,
  probability = T,
  cost = 3
) %>% attr('probabilities')  %>%
  as.data.frame () %>%
  `colnames<-`(paste0("NMF:123456.PCA.SVM.", colnames(.), ".p")) %>%
  tibble::rownames_to_column('sid')



plt.single <- plt.single %>%
  dplyr::left_join(tmp, by = c('sid' = 'sid'))



rm(tmp.fit, tmp)


# export without overwriting old SVM
#
# tmp.export <- read.table("output/tables/gsam_nmf_lda_data.txt.old") %>%
#   dplyr::left_join(tmp, by = c('sid' = 'sid'))
# 
# write.table(tmp.export, "output/tables/gsam_nmf_lda_data.txt")


# ~~ Export stats ~~ ----

write.table(
  tmp.out,
  "output/tables/gsam_nmf_lda_data.new.txt")




# Per pair stats ----
## Split table into pairs + eucl dist ----


plt.paired <- plt.single %>%
  #dplyr::filter(resection %in% c("R1", "R2") ) %>%
  dplyr::group_by(pid) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::top_n(1, sid) %>%
  dplyr::select(all_of('pid')) %>%
  as.data.frame() %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection == 'R1') %>% `colnames<-`( paste0(colnames(.) , ".R1") ) , by=c('pid'='pid.R1') ) %>%
  dplyr::left_join(plt.single %>% dplyr::filter(resection %in% c('R2', 'R3', 'R4') ) %>% `colnames<-`( paste0(colnames(.) , ".R2") ) , by=c('pid'='pid.R2') ) %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                        (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2)) %>%
  dplyr::mutate(subtype.public.status = as.factor(ifelse(subtype.public.R1 == subtype.public.R2, "Stable", "Transition"))) %>%
  #dplyr::mutate(`NMF:123456.membership.status` = as.factor(ifelse(`NMF:123456.membership.R1`== `NMF:123456.membership.R2`, "Stable", "Transition"))) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.factor(ifelse(`NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`, "Stable", "Transition")))



plt.paired %>%   dplyr::filter(pid %in% c('G-SM-R056-2', 'GLSS-HF-3081', 'GLSS-HF-2869')  )


plot(density(plt.paired$eucledian.dist))
hist(plt.paired$eucledian.dist,breaks=15)

plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>% 
  nrow()

plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  nrow()


plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  nrow()



dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` != 'Classical') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)




dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` != 'Mesenchymal') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)



dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)



dist.stable <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::pull(eucledian.dist)


dist.switch <- plt.paired %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural')  %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::pull(eucledian.dist)


wilcox.test(dist.stable, dist.switch)




# Tests ----

## Verschil eucl. distance paren en geshufflede paren ----
## maak een seed voor shufflen voor reproduceerbaarheid of maak alle mogelijke paren
## maak bijbehorende boxplots?

env.test <- env()
env.test$pairs <- plt.paired
env.test$shuffle <- tidyr::crossing(
  plt.paired %>%
    dplyr::select('pid','sid.R1') %>%
    dplyr::rename(pid.R1 = pid) , 
  plt.paired %>%
    dplyr::select('pid','sid.R2') %>%
    dplyr::rename(pid.R2 = pid)
) %>% dplyr::filter(pid.R1 != pid.R2) %>%
  dplyr::mutate(pid.R1 = NULL, pid.R2 = NULL) %>%
  dplyr::left_join(plt.paired %>% dplyr::select(ends_with('.R1')), by = c('sid.R1'='sid.R1')) %>%
  dplyr::left_join(plt.paired %>% dplyr::select(ends_with('.R2')), by = c('sid.R2'='sid.R2')) %>%
  as.data.frame %>%
  dplyr::mutate(eucledian.dist = sqrt((`NMF:123456.PC1.n.R1` - `NMF:123456.PC1.n.R2`)^2 +
                                        (`NMF:123456.PC2.n.R1` - `NMF:123456.PC2.n.R2`)^2))



plt.1 <- data.frame(val = env.test$pairs$eucledian.dist , type = 'Actual distances' ) %>% 
  dplyr::mutate(order = (rank(val)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.2 <- data.frame(val = env.test$shuffle$eucledian.dist , type = 'Shuffled distances' ) %>% 
  dplyr::mutate(order = (rank(val)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.5) # band scaling

plt <- rbind(plt.1, plt.2) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type)))

ggplot(plt, aes(x=type, y=val , group = type)) +
  geom_violin() +
  geom_point(aes(x=x),cex=0.8,col="#888888") +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  youri_gg_theme +
  labs(x = NULL, y = "Eucledian distance PNC space")


wilcox.test(env.test$pairs$eucledian.dist, env.test$shuffle$eucledian.dist)


#plot(hist(env.test$shuffle$eucledian.dist,breaks=25))
# library(fitdistrplus)
# plot(fitdist(env.test$shuffle$eucledian.dist , "weibull"))  # meest lijkende
# plot(fitdist(env.test$shuffle$eucledian.dist, "gamma"))
# 
# plot(fitdist(env.test$shuffle$eucledian.dist, "geom"))
# plot(fitdist(env.test$shuffle$eucledian.dist, "cauchy")) # nope!
# #plot(fitdist(env.test$shuffle$eucledian.dist, "t"))
# plot(fitdist(round(env.test$shuffle$eucledian.dist), "pois"))
# #plot(fitdist(round(env.test$shuffle$eucledian.dist), "nbinom"))
# plot( fitdist(env.test$shuffle$eucledian.dist, "lnorm")) # deze zeker niet
# plot( fitdist(env.test$shuffle$eucledian.dist, "exp")) # deze zeker niet
# plot( fitdist(env.test$shuffle$eucledian.dist , "logis")) # deze ook niet



# 
# 
# plt <- data.frame(val = env.test$pairs$eucledian.dist , type = 'real' ) %>% 
#   rbind(data.frame(val = env.test$shuffle %>%
#                      dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`) %>%
#                      dplyr::pull(eucledian.dist)
#                     , type = 'shuffle'))
# 
# 
# ggplot(plt, aes(x=type, y=val )) +
#   geom_violin() +
#   geom_boxplot(width=0.1)
# 

## Verschil distance transities ----

# https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc

t <- seq(0, 180, by = 1) * pi / 180

r <- .5
x <- r * cos(t)
y <- r*4 * sin(t)
y[20:162] <- y[20] # Flattens the arc
y <- y * 0.3
arc.a <- data.frame(x = x, eucledian.dist = y, type="GLASS")


r <- .5
x <- r * cos(t) * 2
y <- r*4 * sin(t)
y[20:162] <- y[20] # Flattens the arc
y <- y * 0.3
arc.b <- data.frame(x = x, eucledian.dist = y, type="GLASS")


r <- .5
x <- r * cos(t) * 1 + 0.5
y <- r*4 * sin(t)
y[20:162] <- y[20] # Flattens the arc
y <- y * 0.3
arc.c <- data.frame(x = x, eucledian.dist = y, type="GLASS")



### Classical ----


plt.c.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Classical\nto\nClassical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.c.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Classical\nto\nProneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.c.p <- wilcox.test(plt.c.c$eucledian.dist , plt.c.p$eucledian.dist )$p.value

plt.c.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Classical\nto\nMesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.c.m <- wilcox.test(plt.c.c$eucledian.dist , plt.c.m$eucledian.dist )$p.value

plt <- rbind(plt.c.c , plt.c.p, plt.c.m) %>%
  dplyr::mutate(type = factor(type, levels = c("Classical\nto\nClassical", 
                                               "Classical\nto\nMesenchymal",
                                               "Classical\nto\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS")))


plt.c <- ggplot(plt, aes(x=type, y=eucledian.dist , group = type , fill = type == "Classical\nto\nClassical")) +
  ylim(0, 4.3) +
  geom_violin(width=1.105) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col = dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset", fill = "Stable") +
  youri_gg_theme + 
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_line(data = arc.a, aes(x = x+1.5, y = eucledian.dist+3.85), lty = 2) +
  geom_line(data = arc.b, aes(x = x+2, y = eucledian.dist+3.85 + 0.2), lty = 2) +
  geom_text(data = data.frame(x = c(2), eucledian.dist = c(4.3), type="GLASS", label = "**" ),
            size=6,
            aes(x=x, y=eucledian.dist, label=label)) +
  geom_text(data = data.frame(x = c(1.5), eucledian.dist = c(4.05), type="GLASS", label = "***" ),
            size=6,
            aes(x=x, y=eucledian.dist, label=label)) +
  scale_fill_manual(values = c('TRUE' = 'gray95', 'FALSE' = 'white'))


rm(plt.c.c , plt.c.p, plt.c.m)




### Mesenchymal ----


plt.m.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Mesenchymal\nto\nMesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.m.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Mesenchymal\nto\nClassical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.m.c <- wilcox.test(plt.m.m$eucledian.dist , plt.m.c$eucledian.dist )$p.value

plt.m.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Mesenchymal\nto\nProneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.m.p <- wilcox.test(plt.m.m$eucledian.dist , plt.m.p$eucledian.dist )$p.value

plt <- rbind(plt.m.c , plt.m.p, plt.m.m) %>%
  dplyr::mutate(type = factor(type, levels = c("Mesenchymal\nto\nClassical",
                                               "Mesenchymal\nto\nMesenchymal",
                                               "Mesenchymal\nto\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS")))


plt.m <- ggplot(plt, aes(x=type, y=eucledian.dist , group = type, fill = type == "Mesenchymal\nto\nMesenchymal")) +
  ylim(0, 4.3) +
  geom_violin(width=1.105) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col = dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset", fill = "Stable") +
  youri_gg_theme + 
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c('TRUE' = 'gray95', 'FALSE' = 'white'))

# no sig diff's



rm(plt.m.c , plt.m.p, plt.m.m)



### Proneural ----



plt.p.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Proneural\nto\nProneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

plt.p.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Proneural\nto\nClassical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.p.c <- wilcox.test(plt.p.p$eucledian.dist , plt.p.c$eucledian.dist )$p.value 

plt.p.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::select(c(dataset.R1, eucledian.dist)) %>%
  dplyr::mutate(type = 'Proneural\nto\nMesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1)) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) # band scaling

sig.p.m <- wilcox.test(plt.p.p$eucledian.dist , plt.p.m$eucledian.dist )$p.value


plt <- rbind(plt.p.c , plt.p.p, plt.p.m) %>%
  dplyr::mutate(type = factor(type, levels = c("Proneural\nto\nClassical",
                                               "Proneural\nto\nMesenchymal",
                                               "Proneural\nto\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(type)) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet = "Proneural")


plt.p <- ggplot(plt, aes(x=type, y=eucledian.dist , group = type, fill = type == "Proneural\nto\nProneural")) +
  ylim(0, 4.3) +
  geom_violin(width=1.105) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col=dataset, fill=dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset", fill = "Stable") +
  youri_gg_theme + 
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  #geom_line(data = arc.a, aes(x = x+1.5, y = eucledian.dist+3.85), lty = 2) +
  geom_line(data = arc.c, aes(x = x+2, y = eucledian.dist+3.85 + 0.2), lty = 2) +
  geom_text(data = data.frame(x = c(2.5), eucledian.dist = c(4.3), type="GLASS", label = "*" ),
            size=6,
            aes(x=x, y=eucledian.dist, label=label)) +
  scale_fill_manual(values = c('TRUE' = 'gray95', 'FALSE' = 'white'))

rm(plt.p.c , plt.p.p, plt.p.m)


### combined ----

plt.c / plt.m / plt.p


ggsave('output/figures/paper_subtypes_nmf_eucledian_dists_stable_unstable.pdf',width=8 ,height=6*1.5)


plt.c.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.c.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.c.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Classical') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.c <- rbind(plt.c.c , plt.c.p, plt.c.m) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) %>%
  dplyr::mutate(type = factor(paste0('to\n', `NMF:123456.PCA.SVM.class.R2`), levels = c("to\nClassical","to\nMesenchymal","to\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet = "Classical")



plt.p.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.p.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.p.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Proneural') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.p <- rbind(plt.p.c , plt.p.p, plt.p.m) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4)  %>%
  dplyr::mutate(type = factor(paste0('to\n', `NMF:123456.PCA.SVM.class.R2`), levels = c("to\nClassical","to\nMesenchymal","to\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet = "Proneural")


plt.m.m <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Mesenchymal') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.m.c <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Classical') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.m.p <- env.test$pairs %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == 'Mesenchymal') %>%
  dplyr::filter(`NMF:123456.PCA.SVM.class.R2` == 'Proneural') %>%
  dplyr::mutate(order = (rank(eucledian.dist)) - 1) %>%
  dplyr::mutate(x = order / (nrow(.) - 1))

plt.m <- rbind(plt.m.c , plt.m.p, plt.m.m) %>%
  dplyr::mutate(x = x - 0.5) %>%
  dplyr::mutate(x = x *  0.4) %>%
  dplyr::mutate(type = factor(paste0('to\n', `NMF:123456.PCA.SVM.class.R2`), levels = c("to\nClassical","to\nMesenchymal","to\nProneural"))) %>%
  dplyr::mutate(x = x + as.numeric(as.factor(type))) %>%
  dplyr::mutate(dataset = as.factor(ifelse(dataset.R1 == "GSAM", "G-SAM" , "GLASS"))) %>%
  dplyr::mutate(facet ="Mesenchymal")


plt <- rbind(plt.c, plt.m, plt.p) %>%
  dplyr::mutate(stable =   `NMF:123456.PCA.SVM.class.R1` == `NMF:123456.PCA.SVM.class.R2`)


ggplot(plt, aes(x=type, y=eucledian.dist , group = type)) +
  ylim(0, 4.3) +
  geom_violin(width=1.105, aes(fill = stable)) +
  geom_boxplot(width=0.1,outlier.shape = NA,alpha=0.5)  +
  geom_point(aes(x=x, col = dataset),cex=2) +
  labs(x = NULL, y = "Eucledian distance in PNC space", col = "Dataset") +
  youri_gg_theme +
  theme(  axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(rows = vars(facet), space="free_x") +
  scale_fill_manual(values = c('TRUE'='gray95', 'FALSE'='white'))


ggsave('output/figures/paper_subtypes_nmf_eucledian_dists_stable_unstable.pdf',width=8 ,height=6*1.5)




# Plots ----


## Determine Contour ----


resolution <- 250 # 1000 x 1000 data points

off_x <- (max(plt.single$`NMF:123456.PC1`) - min(plt.single$`NMF:123456.PC1`)) * 0.025
off_y <- (max(plt.single$`NMF:123456.PC2`) - min(plt.single$`NMF:123456.PC2`)) * 0.025


range_pc1 = seq(from = min(plt.single$`NMF:123456.PC1`) - off_x, to = max(plt.single$`NMF:123456.PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(plt.single$`NMF:123456.PC2`) - off_y, to = max(plt.single$`NMF:123456.PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:123456.PC1' = range_pc1, 'NMF:123456.PC2' = range_pc2)
#nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.lda, newdata = range_df) %>%
nmf.pca.lda.countours <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = range_df) %>% data.frame(class = .) %>%
  #as.data.frame() %>%
  cbind(range_df) %>%
  dplyr::select(c('class', 'NMF:123456.PC1', 'NMF:123456.PC2')) %>%
  dplyr::mutate(type="Contour")

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)


# Figure S1E ----


plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'subtype.public', sid.label, dataset)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = subtype.public),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA))



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col=class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_point(data = subset(plt, type == "Patient Sample" ), aes(fill = class), col="black", 
             pch=21, size=2.5) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features",
       col='Sub-type',fill="Sub-type") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


ggsave("output/figures/Figure_S1_E.png", width=12 * 0.4,height=10 * 0.4)



# Figure S1F ----


plt <- rbind(plt.single %>%
               dplyr::mutate(svm.status = ifelse(subtype.public == `NMF:123456.PCA.SVM.class`, "positive", 'negative')) %>% 
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.SVM.class', sid.label, dataset, svm.status)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.SVM.class`),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA, svm.status = NA))



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col=class, label = sid.label)) + 
  #geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1, na.rm = T) +
  geom_point(data = subset(plt, type == "Patient Sample" ), aes(fill = class), col="black", 
             pch=21, size=2.5) +
  geom_point(data = subset(plt, type == "Patient Sample" & svm.status == 'negative' ),
             col = 'black', pch=19, size=0.7) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(as.factor(class))),
               colour="gray40",
               size=0.25,
               na.rm = T,
               lty=2,breaks=c(1.5,2.5)) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features",
       col='Sub-type',fill="Sub-type") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


ggsave("output/figures/Figure_S1_F.png", width=12 * 0.4,height=10 * 0.4)



# Figure 2 ABC: 3 GITS panels ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.classical = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(primary.proneural = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(complete.pair= pid %in% (plt.single %>% dplyr::filter(resection == "R1") %>% dplyr::pull(pid)) &
                  pid %in% (plt.single %>% dplyr::filter(resection %in% c("R2","R3","R4")) %>% dplyr::pull(pid))                  
  ) %>%
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.character(`NMF:123456.PCA.SVM.status`) ) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = ifelse(is.na(`NMF:123456.PCA.SVM.status`), 'NA', `NMF:123456.PCA.SVM.status`))





ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.classical | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.classical & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_classical.pdf', width=12 * 0.45,height=10 * 0.45)







ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.mesenchymal | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.mesenchymal & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_mesenchymal.pdf', width=12 * 0.45,height=10 * 0.45)





ggplot(plt , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, group=pid,  col = `NMF:123456.PCA.SVM.status`, label=pid)) +
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) ,  aes(fill = factor(class)), alpha=0.05) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)), colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(!primary.proneural | !complete.pair), size=2.5 * 0.65, alpha=0.15, col="black",pch=21, fill='gray80') +
  #geom_path( data = plt %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>% dplyr::filter(primary.classical == T), arrow = arrow(ends = "last",  type = "closed", angle=15, length = unit(0.135, "inches")), col="red", alpha = 0.8 ) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & `NMF:123456.PCA.SVM.class.R2` == "Mesenchymal") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#eab509", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & `NMF:123456.PCA.SVM.class.R2` == "Proneural") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#ff5f68", col=rgb(0,0,0,0.6),
               lwd = 0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_segment(data=plt.paired %>%
                 dplyr::filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & `NMF:123456.PCA.SVM.class.R2` == "Classical") %>% 
                 dplyr::mutate(`NMF:123456.PC1`=`NMF:123456.PC1.R1`, `NMF:123456.PC2`=`NMF:123456.PC2.R1`),
               aes(xend = `NMF:123456.PC1.R2`, yend=`NMF:123456.PC2.R2`, fill = `NMF:123456.PCA.SVM.class.R2`), arrow.fill="#6ba6e5", col=rgb(0,0,0,0.6), 
               lwd=0.35,
               arrow = arrow(ends = "last", type = "closed",  angle=15, length = unit(0.135 * 0.65, "inches"))) +
  geom_point(data = plt %>% dplyr::filter(primary.proneural & complete.pair & resection == "R1"), aes(fill = `NMF:123456.PCA.SVM.class`), size=2.5 * 0.65, alpha=0.8, col="black",pch=21) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="Sub-type") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'gray30','NA'='#cc00cc')) +
  scale_fill_manual(values = subtype_colors)


ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_proneural.pdf', width=12 * 0.45,height=10 * 0.45)







## PCA:1+2(NMF:1+2+3) + SVM contours + GlioVis labels ----


plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'subtype.public', sid.label, dataset)) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = subtype.public),
             nmf.pca.lda.countours %>% dplyr::mutate(sid.label = NA, dataset=NA))



plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM")
                   , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample" & dataset == "GSAM"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



plt.glass <- ggplot(plt %>% dplyr::filter(dataset != "GSAM")
                    , aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, col = class, label = sid.label)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = subset(plt, type == "Patient Sample" & dataset != "GSAM"), size=1.5) +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call]',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_GlioVis_LDA-countours.png',width=7,height=7.2)




## PCA:1+2(NMF:1+2+3) + SVM contours + SVM labels ----



plt <- rbind(plt.single %>%
               dplyr::select(c(`NMF:123456.PC1`, `NMF:123456.PC2`, 'NMF:123456.PCA.SVM.class')) %>%
               dplyr::mutate(type = "Patient Sample") %>%
               dplyr::rename(class = `NMF:123456.PCA.SVM.class`),
             nmf.pca.lda.countours)

ggplot(plt, aes(x = `NMF:123456.PC1`, y = `NMF:123456.PC2`, fill = class)) + 
  geom_raster(data = subset(plt, type == "Contour"), aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data= subset(plt, type == "Contour"), aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               na.rm = T,
               lty=2,
               breaks=c(1.5,2.5)
  ) +
  geom_point(data = subset(plt, type == "Patient Sample"), size=1.8, pch=21, lwd=0.5, col="gray20") +
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col='Subtype by GlioVis [majority call] / GLASS consortium',fill="LDA countour") +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) + 
  scale_fill_manual(values = subtype_colors)



ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours.png',width=7,height=7.2)
ggsave('output/figures/paper_GLASS_subtypes_nmf_S150G_PC1_PC2_LDA-reclassification_LDA-countours.pdf',width=7,height=7.2)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [classical] ----


plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.classical = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(complete.pair= pid %in% (plt.single %>% dplyr::filter(resection == "R1") %>% dplyr::pull(pid)) &
                  pid %in% (plt.single %>% dplyr::filter(resection %in% c("R2","R3","R4")) %>% dplyr::pull(pid))                  
  ) %>%
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = as.character(`NMF:123456.PCA.SVM.status`) ) %>%
  dplyr::mutate(`NMF:123456.PCA.SVM.status` = ifelse(is.na(`NMF:123456.PCA.SVM.status`), 'NA', `NMF:123456.PCA.SVM.status`))
# %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4

# 
# plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
#   geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
#   geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
#                colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
#   geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Classical" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1"), size=1.5, alpha=0.8) +
#   geom_path(
#     data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
#       dplyr::filter(primary.classical == T),
#     arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
#   youri_gg_theme +
#   labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
#   scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
#   scale_fill_manual(values = subtype_colors)
# 
# 
# plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
#   geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
#   geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
#                colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
#   geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Classical" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
#   geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1"), size=1.5, alpha=0.8) +
#   geom_path(
#     data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
#       dplyr::filter(primary.classical == T),
#     arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
#   youri_gg_theme +
#   labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
#   scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
#   scale_fill_manual(values = subtype_colors)
# 
# 
# plt.gsam + plt.glass
# ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_classical.png',width=15,height=6.5)







## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [proneural] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::mutate(primary.proneural = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::arrange(pid, resection)
# %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Proneural" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.proneural == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_proneural.png',width=15,height=6.5)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [mesenchymal] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` )) , by=c('pid'='pid')) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") %>% dplyr::pull(pid))) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::arrange(pid, resection)
# %>%  dplyr::filter(pid %in% c('GLSS-SM-R056', 'GLSS-HF-3081', 'GLSS-HF-2869')  ) << R3 R4


plt.gsam <- ggplot(plt %>% dplyr::filter(dataset == "GSAM") , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="G-SAM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.glass <- ggplot(plt  %>% dplyr::filter(dataset != "GSAM"), aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset != "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset != "GSAM")%>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title="GLASS GBM") +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors)


plt.gsam + plt.glass



ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_mesenchymal.png',width=15,height=6.5)



## PCA:1+2(NMF:1+2+3) + LDA countors + LDA labels [final?] ----



plt <- plt.single %>%
  dplyr::left_join(plt.paired %>% dplyr::select(c("pid", `NMF:123456.PCA.SVM.status` ,`NMF:123456.PCA.SVM.class.R1`)) , by=c('pid'='pid')) %>%
  dplyr::mutate(primary.mesenchymal = pid %in% (plt.single %>% dplyr::pull(pid))) %>%
  dplyr::mutate(resection = gsub("R3","R2",resection)) %>% # of these 2 patients [G-SM-R056-2 & GLSS-HF-3081] the 2nd was not available in the seq data []
  dplyr::mutate(resection = gsub("R4","R2",resection)) %>% # of this patient [GLSS-HF-2869] the 2nd and 3rd were not available in the seq data
  dplyr::arrange(pid, resection) %>%
  dplyr::mutate(facet = ifelse( (!is.na(`NMF:123456.PCA.SVM.class.R1`) & `NMF:123456.PCA.SVM.class.R1` == "Classical") | (`NMF:123456.PCA.SVM.class` == "Classical" & resection == "R1") , "Classical" , "?" ) ) %>%
  dplyr::mutate(facet = ifelse( (!is.na(`NMF:123456.PCA.SVM.class.R1`) & `NMF:123456.PCA.SVM.class.R1` == "Mesenchymal") | (`NMF:123456.PCA.SVM.class` == "Mesenchymal" & resection == "R1") , "Mesenchymal" , facet ) ) %>%
  dplyr::mutate(facet = ifelse( (!is.na(`NMF:123456.PCA.SVM.class.R1`) & `NMF:123456.PCA.SVM.class.R1` == "Proneural") | (`NMF:123456.PCA.SVM.class` == "Proneural" & resection == "R1") , "Proneural" , facet ) ) %>%
  dplyr::mutate(facet = ifelse(is.na(`NMF:123456.PCA.SVM.class.R1`) & resection == "R2" , "no R1 present" , facet )) %>%
  dplyr::filter(facet != "no R1 present")



ggplot(plt  , aes(x = `NMF:123456.PC1`, y = -`NMF:123456.PC2`, group=pid, col = `NMF:123456.PCA.SVM.status`, label=pid)) + 
  geom_raster(data = nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(fill = factor(class)), alpha=0.1) +
  geom_contour(data=  nmf.pca.lda.countours %>% dplyr::mutate('NMF:123456.PCA.SVM.status' = NA, pid=NA, GBM.transcriptional.subtype.Synapse.status=NA) , aes(z=as.numeric(class)),
               colour="gray40", size=0.25, lty=2,breaks=c(1.5,2.5)) +
  #geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(`NMF:123456.PCA.SVM.class` != "Mesenchymal" ), size=1.0, col="gray80") +
  #geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(resection == "R2" & `NMF:123456.PCA.SVM.status` != "stable" ), size=1.0, col="gray80") +
  geom_point(data = plt %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(resection == "R1"), size=1.5, alpha=0.8) +
  geom_path(
    data = plt  %>% dplyr::filter(dataset == "GSAM") %>% dplyr::filter(pid %in% .$pid[duplicated(.$pid)]) %>%
      dplyr::filter(primary.mesenchymal == T),
    arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 ) + 
  youri_gg_theme +
  labs(x="PC1 on NMF meta-features", y="PC2 on NMF meta-features", col=NULL,fill="SVM class", title=NULL) +
  scale_color_manual(values=c('Stable'= rgb(0,0.75,0), 'Transition' = 'black','NA'='purple')) +
  scale_fill_manual(values = subtype_colors) +
  facet_grid(cols = vars(facet)) +
  #  geom_text_repel(data = subset(plt, pid %in% c('CDA','CDD','CBP','CBR','CCZ') )) + 
  theme(strip.background = element_blank(), strip.text = element_blank())






ggsave('output/figures/paper_subtypes_nmf_S150G_PC1_PC2_SVM-reclassification_SVM-countours_facet.pdf',width=16.7,height=6.5)




# KNN Bootstrapping ~ estimate by chance 

k = 20
df = data.frame()
for(p in 1:nrow(plt.paired) ) {
  #p = 21
  
  # for pair
  target = plt.paired[p,]
  
  
  # find k closest R1's
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                                 (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(k, -ed)
  
  
  #plot(plt.single$`NMF:123456.PC1` , plt.single$`NMF:123456.PC2`, pch=19 , col = as.numeric(plt.single $`NMF:123456.PCA.SVM.class`) + 1 , cex=0.5)
  
  
  #nodes <- data.frame(`NMF:123456.PC1` = target$`NMF:123456.PC1.R1`,
  #                    `NMF:123456.PC2` = target$`NMF:123456.PC2.R1`,
  #                    type="start", check.names=F)
  nodes <- data.frame()
  for(i in 1:k) {
    neighbour <- knn[i,] 
    
    #lines(c(neighbour$`NMF:123456.PC1.R1`,neighbour$`NMF:123456.PC1.R2`), c(neighbour$`NMF:123456.PC2.R1` , neighbour$`NMF:123456.PC2.R2`)  )
    #points(neighbour$`NMF:123456.PC1.R2` , neighbour$`NMF:123456.PC2.R2` , pch=8,cex=0.6  )
    
    delta_PC1 = target$`NMF:123456.PC1.R1` - neighbour$`NMF:123456.PC1.R1`
    delta_PC2 = target$`NMF:123456.PC2.R1` - neighbour$`NMF:123456.PC2.R1`
    
    
    nodes <- rbind(nodes,
                   data.frame(`NMF:123456.PC1` = neighbour$`NMF:123456.PC1.R2` + delta_PC1 ,
                              `NMF:123456.PC2` = neighbour$`NMF:123456.PC2.R2` + delta_PC2 ,
                              type = "node", check.names=F))
    
    
  }
  
  nodes <- nodes %>%
    dplyr::mutate(class.svm = predict(s150.pca.nmf.subtype.classifier.svm , newdata = nodes %>% dplyr::mutate(type=NULL)) )
  
  df <- rbind(df, 
              data.frame(n.cl  = nodes %>% dplyr::filter(class.svm == "Classical") %>% nrow(),
                         n.mes = nodes %>% dplyr::filter(class.svm == "Mesenchymal") %>% nrow(),
                         n.pn  = nodes %>% dplyr::filter(class.svm == "Proneural") %>% nrow() ) %>%
                dplyr::mutate(p.cl = n.cl /  rowSums(.),
                              p.mes = n.mes /  rowSums(.),
                              p.pn = n.pn /  rowSums(.) ) %>%
                dplyr::mutate(pid = target$pid)
  )
  
}


df <- df %>%
  #group_by(pid) %>% 
  #dplyr::summarise( p.cl = mean(n.cl),     p.mes = mean(n.mes), p.pn = mean(n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))



df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()



#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
#df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()


df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.mes) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.cl) %>% sum()
df.2 %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" & dataset.R1 == 'GSAM') %>% dplyr::pull(p.pn) %>% sum()



#df.2 %>% filter( dataset.R1 == 'GSAM' & `NMF:123456.PCA.SVM.status`  == "Stable") %>% nrow()

## Figure S1E





# Figure 2: STAR / Arrow plots ----



plt <- data.frame()
for(i in 1:nrow(plt.paired)) {
  slice <- plt.paired[i,]
  
  line <- data.frame(x = c(slice$`NMF:123456.PC1.n.R1`,
                           slice$`NMF:123456.PC1.n.R2`),
                     y = c(slice$`NMF:123456.PC2.n.R1`,
                           slice$`NMF:123456.PC2.n.R2`),
                     resection = c('R1','R2'),
                     group = slice$pid,
                     init.subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                     rec.subtype = slice$`NMF:123456.PCA.SVM.class.R2`)
  
  line2 <- line %>%
    dplyr::mutate(x = x - line$x[1] ) %>%
    dplyr::mutate(y = y - line$y[1] )
  
  plt <- plt %>%
    rbind(line2)
  
  rm(slice, line, line2)
}



m <- plt %>%
  dplyr::group_by(init.subtype, resection) %>%
  dplyr::summarise(x = mean(x), y=mean(y)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = "Mean")



ggplot(plt, aes(x=x, y=y, group=group, col=init.subtype)) +
  geom_line(alpha=0.75, lwd=1) +
  geom_point(data = plt %>% dplyr::filter(resection == "R2") , aes(fill=rec.subtype), col="black", pch=21, cex=1.8) +
  geom_line(data=m, col="gray20", lty=1, lwd=0.8, show.legend=F) +
  geom_point(data = m %>% dplyr::filter(resection == "R2") , fill = "white", col="black", pch=21, cex=2.1) +
  labs(x="Centered GITS space traversal", y="Centered GITS space traversal") +
  xlim(-max(c(abs(plt$x),abs(plt$y))),max(c(abs(plt$x),abs(plt$y)))) +
  ylim(-max(c(abs(plt$x),abs(plt$y))),max(c(abs(plt$x),abs(plt$y)))) +
  youri_gg_theme +
  facet_grid(cols = vars(init.subtype)) +
  coord_equal() +
  scale_color_manual(values = subtype_colors, guide = guide_legend(title = "Initial subtype", title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) +
  scale_fill_manual(values = subtype_colors, guide = guide_legend(title = "Sub-type recurrence", title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


#ggsave("output/figures/centered_eucledian_distance_traversal_pnc.png",height=10 * 1.1,width=4)
ggsave("output/figures/centered_eucledian_distance_traversal_pnc.pdf",width=10 * 1.1,height=4)


# KNN Bootstrapping ~ estimate by chance ----

# Bij deze bootstrapping zoeken we van de KNN de richtingen/angles
# En gebruiken bootstrappen we de eucledian distances van alle samples


# hoe berekenen we van 2 nodes de angle en verlengen/verkorten we deze
# R1 = ( 1,  1)
# R2 = (-2, -3)

change_length_line <- function(x1, y1, x2, y2, length_new) {
  # From the  line between points (x1, y1) , (x2 ,y2),
  # we want to create a new line with an identical angle
  # but of length `length_new`.
  
  dy <- y2 - y1
  dx <- x2 - x1
  
  #slope <- dy / dx
  angle <- atan2(dy , dx) # in rads
  
  length_x_new <- cos(angle) * length_new
  length_y_new <- sin(angle) * length_new
  
  x2_new <- x1 + length_x_new
  y2_new <- y1 + length_y_new
  
  return (data.frame(x = c(x1, x2_new) ,
                     y = c(y1, y2_new) ,
                     point = as.factor(c("initial start", "new end"))))
}



nn = round((nrow(plt.paired) - 1) / 10) # 9
bootstrap_n = 2 * nn #1 # 25 # bootstrap iterations per sample-pair

df.A = data.frame()
df.B = data.frame()


## A: test with KNN angles and all distances ----

# plot(c(-4,4), c(-4,4), type="n")
# points(1,1, pch=19)
# points(-2,-3, pch=19)
# lines(c(-2,1), c(-3,1))
# df <- change_length_line(1,1,-2,-3, 1)
# points(df[2,1] , df[2,2], pch=19,col="red")


# 
# plot(plt.single$`NMF:123456.PC1.n` ,
#      plt.single$`NMF:123456.PC2.n`, pch=19 , col = as.numeric(as.factor(plt.single$`NMF:123456.PCA.SVM.class`)) + 1 , cex=0.5)
# 



for(p in 1:nrow(plt.paired) ) { # for pair
  #p = 21 
  #p = 23
  target = plt.paired[p,]
  
  
  # find k closest R1's for angles in close proximity
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                                 (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(nn, -ed)
  
  
  #lines(c(target$`NMF:123456.PC1.n.R1`,target$`NMF:123456.PC1.n.R2`), c(target$`NMF:123456.PC2.n.R1` , target$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,0,0.15))
  #points(target$`NMF:123456.PC1.n.R1` , target$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,0,0.75)  )
  #points(target$`NMF:123456.PC1.n.R2` , target$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,0,0.75)  )
  
  
  for(i in 1:bootstrap_n) {
    #lines(c(knn$`NMF:123456.PC1.n.R1`,knn$`NMF:123456.PC1.n.R2`), c(knn$`NMF:123456.PC2.n.R1` , knn$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,1,0.15))
    #points(knn$`NMF:123456.PC1.n.R1` , knn$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,1,0.75)  )
    #points(knn$`NMF:123456.PC1.n.R2` , knn$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,1,0.75)  )
    
    
    # Take the angle from local subsampling
    neighbour <- knn %>% # RANDOM SHUFFLE `MET TERUGLEGGEN` want bootstrappen!!!
      dplyr::slice(sample(1 : dplyr::n() ) [1] )
    
    # take length from overall subsamples
    random_length <- plt.paired %>% # run once with `knn` and once with `plt.paired`?
      dplyr::slice(sample(1 : dplyr::n() )) %>%
      dplyr::pull(eucledian.dist) %>% 
      purrr::pluck(1)
    
    # generate new length using samples length and samples angle
    bootstrapped_line <- change_length_line(neighbour$`NMF:123456.PC1.n.R1`,
                                            neighbour$`NMF:123456.PC2.n.R1`,
                                            
                                            neighbour$`NMF:123456.PC1.n.R2`,
                                            neighbour$`NMF:123456.PC2.n.R2`,
                                            
                                            random_length)
    
    
    
    # fit to target's R1
    bootstrapped_line$x = bootstrapped_line$x + (target$`NMF:123456.PC1.n.R1` - neighbour$`NMF:123456.PC1.n.R1`)
    bootstrapped_line$y = bootstrapped_line$y + (target$`NMF:123456.PC2.n.R1` - neighbour$`NMF:123456.PC2.n.R1`)
    
    
    
    #lines(bootstrapped_line$x, bootstrapped_line$y, lwd=2, lty=3, col="red")
    #points(bootstrapped_line$x[1], bootstrapped_line$y[1],  col="red")
    #points(bootstrapped_line$x[2], bootstrapped_line$y[2],  col="red")
    
    
    # @todo UN/SCALE-NORMALISE this before classification
    bootstrapped_line <- bootstrapped_line %>%
      dplyr::mutate(x.orig = (x + attr(`scale.NMF:123456.PC1`,'scaled:center')) * attr(`scale.NMF:123456.PC1`,'scaled:scale') ) %>%
      dplyr::mutate(y.orig = (y + attr(`scale.NMF:123456.PC2`,'scaled:center')) * attr(`scale.NMF:123456.PC2`,'scaled:scale') )
    
    
    class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = bootstrapped_line %>%
                       dplyr::filter(point == "new end") %>%
                       dplyr::mutate(point = NULL, x = NULL , y=  NULL) %>%
                       dplyr::rename(`NMF:123456.PC1` = x.orig) %>%
                       dplyr::rename(`NMF:123456.PC2` = y.orig) )
    
    # still storing co-ordinates in normalised form
    df.A <- rbind(df.A, data.frame(x = bootstrapped_line$x[2],
                                   y = bootstrapped_line$y[2],
                                   pid = target$pid,
                                   neighbour = neighbour$pid,
                                   class = class ))
    
  }
}


## @todo check -
## plot original labels (normalised)
## plot new labels (normalised)
## this confirms if un-normalisation went o.k.

# plot(plt.single$`NMF:123456.PC1.n` ,
#      plt.single$`NMF:123456.PC2.n`, pch=19 , col = as.numeric(as.factor(plt.single$`NMF:123456.PCA.SVM.class`)) + 1 , cex=0.5)
# 
# points(df %>% dplyr::filter(class == "Mesenchymal") %>% dplyr::pull(x),
#        df %>% dplyr::filter(class == "Mesenchymal") %>% dplyr::pull(y),
#        pch = 1, col="green", cex=0.2)
# points(df %>% dplyr::filter(class == "Classical") %>% dplyr::pull(x),
#        df %>% dplyr::filter(class == "Classical") %>% dplyr::pull(y),
#        pch = 2, col="red", cex=0.2)
# points(df %>% dplyr::filter(class == "Proneural") %>% dplyr::pull(x),
#        df %>% dplyr::filter(class == "Proneural") %>% dplyr::pull(y),
#        pch = 3, cex=0.2, col="blue")



df.A.integrated <- df.A %>%
  dplyr::group_by(pid) %>% 
  dplyr::summarise(n.cl  = sum(class == "Classical"),
                   n.mes = sum(class == "Mesenchymal"),
                   n.pn  = sum(class == "Proneural")) %>%
  dplyr::mutate(p.cl  = n.cl  / (n.cl + n.mes + n.pn) ,
                p.mes = n.mes / (n.cl + n.mes + n.pn) ,
                p.pn  = n.pn  / (n.cl + n.mes + n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))



## B: test with KNN angles and KNN distances ----


for(p in 1:nrow(plt.paired) ) { # for pair
  #p = 21 
  #p = 23
  target = plt.paired[p,]
  
  
  # find k closest R1's for angles in close proximity
  knn <- plt.paired[- p,] %>%
    dplyr::mutate(ed =  sqrt(  (`NMF:123456.PC1.n.R1` - target$`NMF:123456.PC1.n.R1`)^2 +
                                 (`NMF:123456.PC2.n.R1` - target$`NMF:123456.PC2.n.R1`)^2 ) ) %>%
    dplyr::arrange(ed) %>%
    top_n(nn, -ed)
  
  
  #lines(c(target$`NMF:123456.PC1.n.R1`,target$`NMF:123456.PC1.n.R2`), c(target$`NMF:123456.PC2.n.R1` , target$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,0,0.15))
  #points(target$`NMF:123456.PC1.n.R1` , target$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,0,0.75)  )
  #points(target$`NMF:123456.PC1.n.R2` , target$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,0,0.75)  )
  
  
  for(i in 1:bootstrap_n) {
    #lines(c(knn$`NMF:123456.PC1.n.R1`,knn$`NMF:123456.PC1.n.R2`), c(knn$`NMF:123456.PC2.n.R1` , knn$`NMF:123456.PC2.n.R2`) , lwd=1.5 , col = rgb(0,0,1,0.15))
    #points(knn$`NMF:123456.PC1.n.R1` , knn$`NMF:123456.PC2.n.R1` , pch=1,cex=0.6 , col = rgb(0,0,1,0.75)  )
    #points(knn$`NMF:123456.PC1.n.R2` , knn$`NMF:123456.PC2.n.R2` , pch=4,cex=0.6 ,  col = rgb(0,0,1,0.75)  )
    
    
    # Take the angle from local subsampling
    neighbour <- knn %>% # RANDOM SHUFFLE `MET TERUGLEGGEN` want bootstrappen!!!
      dplyr::slice(sample(1 : dplyr::n() ) [1] )
    
    # take length from overall subsamples
    random_length <- knn %>% # run once with `knn` and once with `plt.paired`?
      dplyr::slice(sample(1 : dplyr::n() )) %>%
      dplyr::pull(eucledian.dist) %>%
      purrr::pluck(1)
    
    # generate new length using samples length and samples angle
    bootstrapped_line <- change_length_line(neighbour$`NMF:123456.PC1.n.R1`,
                                            neighbour$`NMF:123456.PC2.n.R1`,
                                            
                                            neighbour$`NMF:123456.PC1.n.R2`,
                                            neighbour$`NMF:123456.PC2.n.R2`,
                                            
                                            random_length)
    
    
    
    # fit to target's R1
    bootstrapped_line$x = bootstrapped_line$x + (target$`NMF:123456.PC1.n.R1` - neighbour$`NMF:123456.PC1.n.R1`)
    bootstrapped_line$y = bootstrapped_line$y + (target$`NMF:123456.PC2.n.R1` - neighbour$`NMF:123456.PC2.n.R1`)
    
    
    
    #lines(bootstrapped_line$x, bootstrapped_line$y, lwd=2, lty=3, col="red")
    #points(bootstrapped_line$x[1], bootstrapped_line$y[1],  col="red")
    #points(bootstrapped_line$x[2], bootstrapped_line$y[2],  col="red")
    
    
    # @todo UN/SCALE-NORMALISE this before classification
    bootstrapped_line <- bootstrapped_line %>%
      dplyr::mutate(x.orig = (x + attr(`scale.NMF:123456.PC1`,'scaled:center')) * attr(`scale.NMF:123456.PC1`,'scaled:scale') ) %>%
      dplyr::mutate(y.orig = (y + attr(`scale.NMF:123456.PC2`,'scaled:center')) * attr(`scale.NMF:123456.PC2`,'scaled:scale') )
    
    
    class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = bootstrapped_line %>%
                       dplyr::filter(point == "new end") %>%
                       dplyr::mutate(point = NULL, x = NULL , y=  NULL) %>%
                       dplyr::rename(`NMF:123456.PC1` = x.orig) %>%
                       dplyr::rename(`NMF:123456.PC2` = y.orig) )
    
    # still storing co-ordinates in normalised form
    df.B <- rbind(df.B, data.frame(x = bootstrapped_line$x[2],
                                   y = bootstrapped_line$y[2],
                                   pid = target$pid,
                                   neighbour = neighbour$pid,
                                   class = class ))
    
  }
}




df.B.integrated <- df.B %>%
  dplyr::group_by(pid) %>% 
  dplyr::summarise(n.cl  = sum(class == "Classical"),
                   n.mes = sum(class == "Mesenchymal"),
                   n.pn  = sum(class == "Proneural")) %>%
  dplyr::mutate(p.cl  = n.cl  / (n.cl + n.mes + n.pn) ,
                p.mes = n.mes / (n.cl + n.mes + n.pn) ,
                p.pn  = n.pn  / (n.cl + n.mes + n.pn) ) %>%
  dplyr::left_join(plt.paired %>%
                     dplyr::select(c('pid','dataset.R1','NMF:123456.PCA.SVM.status','NMF:123456.PCA.SVM.class.R1','NMF:123456.PCA.SVM.class.R2') )
                   , by=c('pid'='pid'))


## C: alleen bootstrappen lengtes (Alle) ----

## D: alleen bootstrappen lengtes (binnen zelfde subtype R1) ----



## visualisation ----



help(curvedarrow)

curvedarrow <- function(from, to, lwd=2, lty=1, lcol="black", arr.col=lcol, 
                        arr.pos=0.5, curve=1, dr=0.01, endhead=FALSE, segment = c(0,1), ...)   {
  
  dpos  <- to-from
  angle <- atan(dpos[2]/dpos[1])*180/pi         # angle between both
  if (is.nan(angle)) return()
  mid   <- 0.5*(to+from)                        # midpoint of ellipsoid arrow
  dst   <- dist(rbind(to, from))                # distance from-to
  ry    <- curve*dst                            # small radius of ellepsoid
  aFrom<-0                                      # angle to and from
  aTo<-pi
  if ( from[1] <= to[1]) {
    aFrom <- pi
    aTo <- 2*pi
  }
  
  if (segment [1] != 0)
    From <- segment[1] * aTo + (1-segment[1]) * aFrom
  else
    From <- aFrom
  
  if (segment [2] != 1)
    To <- segment[2] * aTo + (1-segment[2]) * aFrom
  else
    To <- aTo
  
  meanpi <- arr.pos * aTo + (1-arr.pos) * aFrom
  if (endhead) To <- meanpi
  
  
  plotellipse(rx=dst/2,  ry=ry, mid=mid, angle=angle, from = From, to = To,
              lwd=lwd, lty=lty, lcol=lcol)
  ell <- getellipse(rx=dst/2, ry=ry, mid=mid, angle=angle,
                    from=1.001*meanpi, to=0.999*meanpi, dr= 0.002)       #Changed from -0.002
  Arrows(ell[1,1], ell[1,2], ell[nrow(ell),1], ell[nrow(ell),2],
         code=1, lcol=lcol, arr.col=arr.col, ...)
  curvedarrow <- c(ell[nrow(ell),1], ell[nrow(ell),2])
}



O.cc <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Classical" ) %>% nrow()
O.cm <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Mesenchymal" ) %>% nrow()
O.cp <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Proneural" ) %>% nrow()

O.mc <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Classical" ) %>% nrow()
O.mm <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Mesenchymal" ) %>% nrow()
O.mp <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Proneural" ) %>% nrow()

O.pc <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Classical" ) %>% nrow()
O.pm <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Mesenchymal" ) %>% nrow()
O.pp <- plt.paired %>%  filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>%  filter(`NMF:123456.PCA.SVM.class.R2` == "Proneural" ) %>% nrow()


A.cc <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
A.cm <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
A.cp <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

A.mm <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
A.mc <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
A.mp <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

A.pp <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)
A.pc <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
A.pm <- df.A.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)


B.cc <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
B.cm <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
B.cp <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Classical" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

B.mm <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)
B.mc <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
B.mp <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Mesenchymal" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)

B.pp <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.pn) %>% sum() %>% round(1)
B.pc <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.cl) %>% sum() %>% round(1)
B.pm <- df.B.integrated %>% filter(`NMF:123456.PCA.SVM.class.R1` == "Proneural" ) %>% dplyr::pull(p.mes) %>% sum() %>% round(1)





par(mfrow=c(1,3))

### Orig ----
plot(c(-1,1)* 1.25 , c(-1,1) * 1.25 , type='n' , main="Original data")

curvedarrow( c(0 - 0.12, 0.66 ) , c(0 + 0.12, 0.66), curve = -1.5)# CL -> CL
text(0, 1.2, O.cc)
curvedarrow( c(0, 0.66) , c(-0.66, -0.66), curve = -0.001)#        CL -> MES
text(-0.25, -0.1, O.cm)
curvedarrow( c(-0.66, -0.66),  c(0, 0.66)  , curve = -0.2)#      MES -> CL
text(-0.75, 0.1, O.mc)

curvedarrow( c(0.66 - 0.12, -0.66 ) , c(0.66 + 0.12, -0.66), curve = 1.5)# PN -> PN
text(0.66, -1.2, O.pp)
curvedarrow( c(0, 0.66) , c(0.66, -0.66), curve = -0.2)#         CL -> PN
text(0.25, -0.1, O.pc)
curvedarrow( c(0.66, -0.66),  c(0, 0.66)  , curve = -0.001)#       PN -> CL
text(0.75, 0.1, O.cp)

curvedarrow( c(-0.66 - 0.12, -0.66 ) , c(-0.66 + 0.12, -0.66), curve = 1.5)# MES -> MES
text(-0.66, -1.2, O.mm)
curvedarrow( c(-0.66, -0.66) , c(0.66, -0.66), curve = -0.001)#    MES -> PL
text(0, -0.55, O.mp)
curvedarrow( c(0.66, -0.66) , c(-0.66, -0.66), curve = -0.2)#    PL -> MES
text(0, -1.05, O.pm)

# labels
points(0, 0.66, cex = 8, pch=19, col="white")
points(0, 0.66, cex = 8, col="black")
text(0, 0.66, "CL", col="black")

points(-0.66, -0.66, cex = 8, pch=19, col="white")
points(-0.66, -0.66, cex = 8, col="black")
text(-0.66, -0.66, "MES", col="black")

points(0.66, -0.66, cex = 8, pch=19, col="white")
points(0.66, -0.66, cex = 8, col="black")
text(0.66, -0.66, "PN", col="black")


### A ----

plot(c(-1,1)* 1.25 , c(-1,1) * 1.25 , type='n' , main="Bootstrap A (KNN angle, random dist)")

curvedarrow( c(0 - 0.12, 0.66 ) , c(0 + 0.12, 0.66), curve = -1.5)# CL -> CL
text(0, 1.2, A.cc)
curvedarrow( c(0, 0.66) , c(-0.66, -0.66), curve = -0.001)#        CL -> MES
text(-0.25, -0.1, A.cm)
curvedarrow( c(-0.66, -0.66),  c(0, 0.66)  , curve = -0.2)#      MES -> CL
text(-0.75, 0.1, A.mc)

curvedarrow( c(0.66 - 0.12, -0.66 ) , c(0.66 + 0.12, -0.66), curve = 1.5)# PN -> PN
text(0.66, -1.2, A.pp)
curvedarrow( c(0, 0.66) , c(0.66, -0.66), curve = -0.2)#         CL -> PN
text(0.25, -0.1, A.pc)
curvedarrow( c(0.66, -0.66),  c(0, 0.66)  , curve = -0.001)#       PN -> CL
text(0.75, 0.1, A.cp)

curvedarrow( c(-0.66 - 0.12, -0.66 ) , c(-0.66 + 0.12, -0.66), curve = 1.5)# MES -> MES
text(-0.66, -1.2, A.mm)
curvedarrow( c(-0.66, -0.66) , c(0.66, -0.66), curve = -0.001)#    MES -> PL
text(0, -0.55, A.mp)
curvedarrow( c(0.66, -0.66) , c(-0.66, -0.66), curve = -0.2)#    PL -> MES
text(0, -1.05, A.pm)


# labels
points(0, 0.66, cex = 8, pch=19, col="white")
points(0, 0.66, cex = 8, col="black")
text(0, 0.66, "CL", col="black")

points(-0.66, -0.66, cex = 8, pch=19, col="white")
points(-0.66, -0.66, cex = 8, col="black")
text(-0.66, -0.66, "MES", col="black")

points(0.66, -0.66, cex = 8, pch=19, col="white")
points(0.66, -0.66, cex = 8, col="black")
text(0.66, -0.66, "PN", col="black")



### B ----

plot(c(-1,1)* 1.25 , c(-1,1) * 1.25 , type='n' , main="Bootstrap B (KNN angle, KNN dist)")

curvedarrow( c(0 - 0.12, 0.66 ) , c(0 + 0.12, 0.66), curve = -1.5)# CL -> CL
text(0, 1.2, B.cc)
curvedarrow( c(0, 0.66) , c(-0.66, -0.66), curve = -0.001)#        CL -> MES
text(-0.25, -0.1, B.cm)
curvedarrow( c(-0.66, -0.66),  c(0, 0.66)  , curve = -0.2)#      MES -> CL
text(-0.75, 0.1, B.mc)

curvedarrow( c(0.66 - 0.12, -0.66 ) , c(0.66 + 0.12, -0.66), curve = 1.5)# PN -> PN
text(0.66, -1.2, B.pp)
curvedarrow( c(0, 0.66) , c(0.66, -0.66), curve = -0.2)#         CL -> PN
text(0.25, -0.1, B.pc)
curvedarrow( c(0.66, -0.66),  c(0, 0.66)  , curve = -0.001)#       PN -> CL
text(0.75, 0.1, B.cp)

curvedarrow( c(-0.66 - 0.12, -0.66 ) , c(-0.66 + 0.12, -0.66), curve = 1.5)# MES -> MES
text(-0.66, -1.2, B.mm)
curvedarrow( c(-0.66, -0.66) , c(0.66, -0.66), curve = -0.001)#    MES -> PL
text(0, -0.55, B.mp)
curvedarrow( c(0.66, -0.66) , c(-0.66, -0.66), curve = -0.2)#    PL -> MES
text(0, -1.05, B.pm)


# labels
points(0, 0.66, cex = 8, pch=19, col="white")
points(0, 0.66, cex = 8, col="black")
text(0, 0.66, "CL", col="black")

points(-0.66, -0.66, cex = 8, pch=19, col="white")
points(-0.66, -0.66, cex = 8, col="black")
text(-0.66, -0.66, "MES", col="black")

points(0.66, -0.66, cex = 8, pch=19, col="white")
points(0.66, -0.66, cex = 8, col="black")
text(0.66, -0.66, "PN", col="black")


# Figure 2: Proportional distances ----


# 
# plot(plt.single$`NMF:123456.PC1` ,
#      plt.single$`NMF:123456.PC2`, pch=19 , col = as.numeric(as.factor(plt.single$`NMF:123456.PCA.SVM.class`)) + 1 , cex=0.5)



df <- data.frame()
res <- 150 # resolution
for(p in 1:nrow(plt.paired) ) { # for pair
  
  #p <- which(plt.paired$pid == "AZF")
  
  target <- plt.paired[p,]
  d <- sqrt((target$`NMF:123456.PC1.R1` - target$`NMF:123456.PC1.R2`)^2 +
              (target$`NMF:123456.PC2.R1` - target$`NMF:123456.PC2.R2`)^2) # dist
  
  # lines(
  #   c(target$`NMF:123456.PC1.R1`, target$`NMF:123456.PC1.R2`),
  #   c(target$`NMF:123456.PC2.R1`, target$`NMF:123456.PC2.R2`)
  # )
  
  df.classification <- data.frame('pid' = target$pid,
                                  `NMF:123456.PC1` = target$`NMF:123456.PC1.R1`,
                                  `NMF:123456.PC2` = target$`NMF:123456.PC2.R1`,
                                  frac = 0)
  
  for(r in 1:res) {
    pct <- r / (res+1)
    
    # generate new length using samples length and samples angle
    bootstrapped_line <- change_length_line(target$`NMF:123456.PC1.R1`,
                                            target$`NMF:123456.PC2.R1`,
                                            
                                            target$`NMF:123456.PC1.R2`,
                                            target$`NMF:123456.PC2.R2`,
                                            
                                            d * pct)
    
    #lines(
    #  c(bootstrapped_line$x[1], bootstrapped_line$x[2]),
    #  c(bootstrapped_line$y[1], bootstrapped_line$y[2])
    #)
    #points(bootstrapped_line$x[2], bootstrapped_line$y[2],col=as.numeric(class) + 1, lwd=2,cex=1.6)
    
    df.classification <- df.classification %>% 
      rbind(data.frame(pid = target$pid,
                       `NMF:123456.PC1` = bootstrapped_line$x[2],
                       `NMF:123456.PC2` = bootstrapped_line$y[2],
                       frac = pct) )
  }
  
  df.classification <- df.classification %>%
    rbind(data.frame(pid = target$pid,
                     `NMF:123456.PC1` = target$`NMF:123456.PC1.R2`,
                     `NMF:123456.PC2` = target$`NMF:123456.PC2.R2`,
                     frac = 1))
  
  df.classification$class <- predict(s150.pca.nmf.subtype.classifier.svm , newdata = df.classification %>% 
                                       dplyr::mutate(frac=NULL, pid=NULL) )
  
  df <- df %>% rbind(df.classification)
  
}

df <- df %>%
  dplyr::mutate(pid = as.factor(pid)) %>%
  dplyr::mutate(class = as.factor(class))

# df %>%
#   dplyr::filter(pid == levels(df$pid)[2]) # == AAC






df.paired <- data.frame()
for(ppid in unique(df$pid)) {
  #ppid = "EBV"
  target <- df %>%
    dplyr::filter(`pid` == ppid) %>%
    dplyr::arrange(frac)
  
  l <- length(unique(target$class))
  
  out <- data.frame(pid = ppid,
                    init = 0,
                    init.class = target$class[1],
                    split = NA,
                    end = 1,
                    end.class = target$class[nrow(target)])
  
  for(k in 1:nrow(target)) {
    slice = target[k,]
    
    if(is.na(out$split) & slice$class != out$init.class) {
      out$split = (target[k,]$frac + target[k - 1,]$frac) / 2
      #print("!!")
    }
  }
  
  df.paired <- rbind(out, df.paired)
}
df.paired <- df.paired %>%
  dplyr::mutate(init = NULL,
                init.class = NULL,
                end = NULL,
                end.class = NULL) %>%
  dplyr::left_join(plt.paired , by = c('pid' = 'pid')) %>%
  dplyr::left_join(gsam.patient.metadata %>% dplyr::select(studyID, HM), by = c('pid' = 'studyID') ) %>%
  dplyr::mutate(HM = ifelse(is.na(HM),"N/A",HM))





# maak ggplot
plt <- data.frame()
for(i in 1:nrow(df.paired)){
  slice = df.paired[i,]
  
  
  if(is.na( slice$split )) {
    dist <- slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-0.1,-0.1),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("lpad","lpad")))
    }
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="lpad")
    )
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="lpad")
    )
    
    
    dist <- - slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-2.8,-2.8),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("rpad","rpad")))
    }
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = -slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="rpad")
    )
    
    plt <- rbind(plt, 
                 data.frame(
                   pid = slice$pid,
                   group = slice$pid,
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1` , 
                   d = dist,
                   type="rpad")
    )
    
    
  } else  {
    dist <- slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-0.1,-0.1),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("lpad","lpad")))
    }
    
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = (1 - slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos = (1 - slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos =  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="lpad"
                 ))
    
    
    
    dist <- (slice$split) * slice$eucledian.dist
    
    # Add HM status
    if(slice$HM[1] == "Yes") {
      plt <- rbind(plt,
                   data.frame(
                     pid = c(slice$pid,slice$pid),
                     group = paste0(c(slice$pid,slice$pid),".HM"),
                     pos = c(-2.8,-2.8),
                     subtype = c('Hyper-mutant','Hyper-mutant'),
                     subtype.init = c(slice$`NMF:123456.PCA.SVM.class.R1`,slice$`NMF:123456.PCA.SVM.class.R1`),
                     d = c(dist,dist),
                     type=c("rpad","rpad")))
    }
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = - (1 - slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R1"),
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R1`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos = 0,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    plt <- rbind(plt,
                 data.frame(
                   pid = slice$pid,
                   group = paste0(slice$pid , "-R2"),
                   pos =  (slice$split) *  slice$eucledian.dist,
                   subtype = slice$`NMF:123456.PCA.SVM.class.R2`,
                   subtype.init = slice$`NMF:123456.PCA.SVM.class.R1`,
                   d = dist,
                   type="rpad"
                 ))
    
  }
}



p1 <- ggplot(plt %>% dplyr::filter(type == "lpad"), aes(x = pos , y = reorder(pid, d), col = subtype, group=group)) +
  geom_point(cex=1) +
  geom_path(lwd = 1.2) + 
  facet_grid(rows = vars(subtype.init), cols=vars(type), scales = "free", space="free_y") + 
  labs(x = "Traversal distance GITS space",
       y=NULL,
       col = "Distance within Subtype space") +
  theme_light() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'bottom',
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    axis.text.y = element_text(size = 6)
  ) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_color_manual(values = c(subtype_colors, 'Hyper-mutant' = 'gray40'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))


p2 <- ggplot(plt %>% dplyr::filter(type == "rpad"), aes(x = pos , y = reorder(pid, d), col = subtype, group=group)) +
  geom_point(cex=1.1) +
  geom_path(lwd = 1.2) + 
  facet_grid(rows = vars(subtype.init), cols=vars(type), scales = "free", space="free_y") + 
  labs(x = "Traversal distance GITS space",
       y=NULL,
       col = "Distance within Subtype space") +
  theme_light() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'bottom',
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    axis.text.y = element_text(size = 6)
  ) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_color_manual(values = c(subtype_colors, 'Hyper-mutant' = 'gray40'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75))

p1 + p2




#ggsave("output/figures/eucledian_distance_traversal_pnc.png",height=8*1.3,width=8*1.4)
ggsave("output/figures/eucledian_distance_traversal_pnc.pdf",height=8*1.3,width=8*2.0)



## santoesha style plot ----




plt.prop <- data.frame()

tmp <-  df %>% dplyr::left_join(plt.paired %>%
                                  dplyr::select(c('pid',`NMF:123456.PCA.SVM.class.R1`, eucledian.dist)),
                                by=c('pid'='pid'))

for(ppid in unique(tmp$pid)) {
  print(ppid)
  
  
  tmp.2 <- tmp %>% dplyr::filter(pid == ppid)
  
  
  last_class <- tmp.2[1,]$class
  last_frac <- tmp.2[1,]$frac
  
  for(i in 2:nrow(tmp.2))  {
    #print(i)
    if(tmp.2[i - 1,]$class != tmp.2[i,]$class ) {
      split <-  (tmp.2[i -1 ,]$frac + tmp.2[i  ,]$frac) / 2
      
      #print(paste0("(",i,") ", last_class, " -> ", tmp.2[i - 1,]$class))
      
      plt.prop <- plt.prop %>%
        rbind(data.frame(
          pid = tmp.2[i,]$pid,
          
          from_class = last_class,
          to_class = tmp.2[i - 1,]$class,
          
          from_frac = last_frac  ,
          to_frac = split ,
          delta_frac =  split - last_frac,
          steps = (split - last_frac) / (1/151),
          
          eucl_dist = tmp.2[i,]$eucledian.dist ,
          delta_eucl_dist = tmp.2[i,]$eucledian.dist * (split - last_frac)
        ))
      
      last_class <- tmp.2[i - 1 ,]$class
      last_frac <- split
    }
    
    if ( tmp.2[i,]$frac == 1.0) {
      #print(paste0("(",i,") ", last_class, " -> ", tmp.2[i - 1,]$class))
      
      plt.prop <- plt.prop %>%
        rbind(data.frame(
          pid = tmp.2[i,]$pid,
          
          from_class = last_class,
          to_class = tmp.2[i ,]$class,
          
          from_frac = last_frac  ,
          to_frac = tmp.2[i ,]$frac ,
          delta_frac =  tmp.2[i ,]$frac   - last_frac,
          steps = (tmp.2[i ,]$frac   - last_frac) / (1/151),
          
          eucl_dist = tmp.2[i,]$eucledian.dist ,
          delta_eucl_dist = tmp.2[i,]$eucledian.dist * (tmp.2[i ,]$frac   - last_frac)
        ))
      
    }
  }
  
  #tmp.2
  #plt.prop
  
  test <- plt.prop %>%
    dplyr::filter(pid == ppid)
  
  stopifnot(sum(test$steps ) == 151)
  stopifnot(sum(test$delta_frac) == 1)
}











links <- data.frame() %>%
  rbind(data.frame( from = "CL", to = "CL" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Classical" & to_class == "Classical") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "CL", to = "MES" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Classical" & to_class == "Mesenchymal") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "CL", to = "PN" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Classical" & to_class == "Proneural") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  
  rbind(data.frame( from = "MES", to = "CL" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Mesenchymal" & to_class == "Classical") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "MES", to = "MES" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Mesenchymal" & to_class == "Mesenchymal") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "MES", to = "PN" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Mesenchymal" & to_class == "Proneural") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  
  rbind(data.frame( from = "PN", to = "CL" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Proneural" & to_class == "Classical") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "PN", to = "MES" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Proneural" & to_class == "Mesenchymal") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum )) %>%
  rbind(data.frame( from = "PN", to = "PN" , weight = 
                      plt.prop %>% dplyr::filter(from_class == "Proneural" & to_class == "Proneural") %>%
                      dplyr::pull(delta_eucl_dist) %>% sum ))

links[links$from == "CL",]$weight <- links[links$from == "CL",]$weight / (plt.prop %>% dplyr::filter(from_class == "Classical" & from_frac == 0) %>% nrow())
links[links$from == "MES",]$weight <- links[links$from == "MES",]$weight / (plt.prop %>% dplyr::filter(from_class == "Classical" & from_frac == 0) %>% nrow())
links[links$from == "PN",]$weight <- links[links$from == "PN",]$weight / (plt.prop %>% dplyr::filter(from_class == "Classical" & from_frac == 0) %>% nrow())

links <- links %>%
  dplyr::mutate(weight = round(weight,2))


#dplyr::mutate(weight=round(weight/sum(weight),3) * 100)

#Real values
#links <- data.frame(from=c(rep(c("CL"),3),rep(c("MES"),3),rep(c("PN"),3)), to=c(rep(c("CL","MES","PN"),3)))
#links <- links %>% dplyr::group_by(from) %>% dplyr::mutate(thickness=round(weight/sum(weight),2))

nodes <- data.frame(id=c("CL","MES","PN"),
                    subtype=c("Classical","Mesenchymal","Proneural"),
                    size=c(10,8,9))
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
net <- simplify(net, remove.multiple = F, remove.loops = F)
#E(net)$width <- E(net)$thickness*6
#E(net)$size <- E(net)$size*0.75

plot(net,
     edge.arrow.size=.4,
     diag=T,
     edge.label = E(net)$weight,
     arrow.mode=3,
     edge.curved=.2,
     edge.color="gray70",
     vertex.color="coral",
     vertex.label.font=2,
     edge.label.font=1,
     edge.label.cex=1.3,
     edge.label.color="blue",
)





# 〰 © Dr. Youri Hoogstrate 〰 ----
