#!/usr/bin/env R


if(!exists('glass.gbm.rnaseq.expression.all.samples.vst')) {
  source('scripts/load_glass_expression_data.R')
}


## batches ~ PCA ~ plot ----


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
    tissue.source == "06" ~ "TCGA",
    tissue.source == "14" ~ "TCGA",
    tissue.source == "19" & project == "TCGA" & grepl("19-4065|19-1389", aliquot_barcode) ~ "TCGA",
    tissue.source == "19" & project == "TCGA" & grepl("19-0957", aliquot_barcode) ~ "GLASS:19",
    tissue.source == "19" & project == "GLSS" ~ "GLASS:19",
    tissue.source == "CU" & project == "GLSS" & PC1 < -15 ~ "GLASS:CU [A]",
    tissue.source == "CU" & project == "GLSS" & PC1 >= -15 ~ "GLASS:CU [B]",
    tissue.source == "HF" & project == "GLSS" & PC2 < -10 ~ "GLASS:HF [numeric]",
    tissue.source == "HF" & project == "GLSS" & PC2 >= -10 ~ "GLASS:HF [alpha-numeric]",
    tissue.source == "HK" & project == "GLSS" ~ "GLASS:HK",
    tissue.source == "LU" & project == "GLSS" ~ "GLASS:LU",
    tissue.source == "LX" & project == "GLSS" ~ "GLASS:LX",
    tissue.source == "MD" & project == "GLSS" ~ "GLASS:MD",
    tissue.source == "SM" & project == "GLSS" & PC1 > 150 ~ "GLASS:SM [A]",
    tissue.source == "SM" & project == "GLSS" & PC1 <= 150 ~ "GLASS:SM [B]",
    tissue.source == "SN" & project == "GLSS" ~ "GLASS:SN",
    T ~ "t.b.d"
  ))


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



# export ----



saveRDS(plt |>  dplyr::select(aliquot_barcode, deduced.batch) |> dplyr::rename(predicted.GLASS.batch = deduced.batch),"cache/analysis_predict_GLASS_batches.R")



