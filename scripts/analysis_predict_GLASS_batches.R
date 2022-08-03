#!/usr/bin/env R


if(!exists('glass.gbm.rnaseq.expression.all.samples.vst')) {
  source('scripts/load_glass_expression_data.R')
}


# batches ~ PCA ~ plot ----


plt <- glass.gbm.rnaseq.expression.all.samples.vst |> 
  t() |> 
  prcomp() |> 
  purrr::pluck('x') |> 
  as.data.frame(stringsAsFactors=F) |> 
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8) |> 
  tibble::rownames_to_column('aliquot_barcode') |> 
  dplyr::left_join(glass.gbm.rnaseq.metadata.all.samples, by=c('aliquot_barcode'='aliquot_barcode'), suffix=c('','')) |> 
  dplyr::mutate(project = gsub("^([^\\-]+)-[^\\-]+-.+$","\\1",aliquot_barcode)) |>  # project
  dplyr::mutate(tissue.source = gsub("^[^\\-]+-([^\\-]+)-.+$","\\1",aliquot_barcode)) |>  # tissue.source
  dplyr::mutate(is.sm.r099 = grepl("SM-R099",aliquot_barcode)) |> 
  dplyr::left_join(
    data.frame(colSums = colSums(glass.gbm.rnaseq.expression.all.samples)) |> 
      tibble::rownames_to_column('aliquot_barcode'),
    by=c('aliquot_barcode'='aliquot_barcode'), suffix=c('','')
  ) |> 
  dplyr::mutate(depth = case_when(colSums < 18000000 ~ "low",
                                  colSums > 52500000 ~ "high",
                                  T ~ "average")) |> 
  dplyr::mutate(deduced.batch = case_when(
    grepl("CU-R004-P", sid.label) ~ "outlier",
    grepl("G-SM-R099-1", sid.label) ~ "outlier",
    grepl("G-SM-R111-1", sid.label) ~ "outlier",
    T ~ "t.b.d"
  ))


# check how well the provided batches are ----
## combined ----
ggplot(plt, aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  #facet_grid(cols = vars(aliquot_batch_synapse)) + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-19-RNA ----
#'@details Batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-19-RNA"| grepl('-19-',aliquot_barcode,fixed=T)),
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-CU-RNA ----
#'@details Two distinct batches - of which some are presumably GLSS-PD-RNA
#'@details "GLSS-CU-R006-TP-01R-RNA-VRCJGW","GLSS-CU-R018-R1-01R-RNA-0OVRM9",
#'@details "GLSS-CU-R019-TP-01R-RNA-S5OA7D","GLSS-CU-R019-R1-01R-RNA-BGF7KM",
#'@details "GLSS-CU-R006-R1-01R-RNA-4EGGFS","GLSS-CU-R017-R1-01R-RNA-S3WCU4" moved to GLSS-PM-RNA and fixed it
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-CU-RNA" | grepl('-CU-',aliquot_barcode,fixed=T))
       , aes(x=PC1, y=PC2, label=aliquot_barcode, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)



## GLSS-H2-RNA ----
#'@details Batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-H2-RNA"  | grepl('-HF-',aliquot_barcode,fixed=T)), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-HF-RNA ----
#'@details Batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-HF-RNA" | grepl('-HF-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-HK-RNA ----
#'@details Batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-HK-RNA" | grepl('-HK-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-JG-RNA ----
#'@details Batch = presumably part of -MD-
#'@details "GLSS-MD-LP04-TP-01R-RNA-LTGS2F" "GLSS-MD-LP04-R1-01R-RNA-TPSKGU" moved to GLASS-MD-RNA
# batch is now gone
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-JG-RNA" | grepl('-MD-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)




## GLSS-LU-RNA ----
#'@details Batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-LU-RNA" | grepl('-LU-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-LU-RNA ----
#'@details two of JG presumably should be included
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-MD-RNA" | grepl('-MD-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)

## GLSS-PD-RNA ----
#'@details six of CU presumably should be included
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-PD-RNA" | grepl('-CU-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)


## GLSS-SM-RNA ----
#'@details should be split in two batches, with 2 outliers

plt |>
  dplyr::filter(aliquot_batch_synapse == "GLSS-SM-RNA") |> 
  dplyr::filter(PC1 > 150) |> 
  dplyr::pull(aliquot_barcode)

plt |>
  dplyr::filter(aliquot_batch_synapse == "GLSS-SM-RNA") |> 
  dplyr::filter(PC1 <= 150) |> 
  dplyr::pull(aliquot_barcode)


ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-SM-RNA" | grepl('-SM-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)

## GLSS-SN-RNA ----
#'@details batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "GLSS-SN-RNA" | grepl('-xx-',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)

## TCGA-GB-RNA ----
#'@details batch = good
ggplot(plt |> dplyr::filter(aliquot_batch_synapse == "TCGA-GB-RNA" | grepl('TCGA',aliquot_barcode,fixed=T) ), 
       aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse,group=case_barcode)) +
  geom_line(alpha=0.15,col="black") +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  xlim(-120,270) +
  ylim(-100,200)






# further analysis? ----
# extreme batch effects
ggplot(plt, aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  facet_grid(cols = vars(tissue.source))

ggplot(plt, aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw() + 
  facet_grid(cols = vars(aliquot_batch_synapse))


# extreme batch effects
ggplot(plt, aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()



# find sub batches in tumor source 06;
ggplot(plt, aes(x=PC1, y=PC2, label=sid.label, col=tissue.source == "06")) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source 19;
ggplot(plt |> dplyr::filter(tissue.source == "19"), aes(x=PC1, y=PC2, label=sid.label, col=batch)) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source CU;
ggplot(plt |> dplyr::filter(tissue.source == "CU"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source HF;
ggplot(plt |> dplyr::filter(tissue.source == "HF"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source HK;
ggplot(plt |> dplyr::filter(tissue.source == "HK"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source LX;
ggplot(plt |> dplyr::filter(tissue.source == "LX"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source SM;
ggplot(plt |> dplyr::filter(tissue.source == "SM"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source SN;
ggplot(plt |> dplyr::filter(tissue.source == "SN"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source MD
ggplot(plt |> dplyr::filter(tissue.source == "MD"), aes(x=PC1, y=PC2, label=sid.label, col=deduced.batch, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()

# find sub batches in tumor source MD
ggplot(plt |> dplyr::filter(tissue.source == "MD"), aes(x=PC1, y=PC2, label=sid.label, col=aliquot_batch_synapse, group=case_barcode)) +
  geom_line(col="gray60", lwd=0.5) +
  geom_point() +
  ggrepel::geom_text_repel(cex=3) +
  theme_bw()




