#!/usr/bin/env R

# settings ----


options(warnPartialMatchDollar = TRUE) # https://stackoverflow.com/questions/32854683/data-frames-in-r-name-autocompletion


# load libs ----


# load data ----


source("scripts/load_G-SAM_metadata.R")
source("scripts/load_G-SAM_expression_data.R")

source("scripts/load_GLASS_data.R") # glass & tcga validation set + metedata

source('scripts/load_results.out.R')



# prepare data ----

## take relevant subsets ----




tmp.sel.gsam <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(old.batch == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
  dplyr::filter(tumour.percentage.dna >= 15) |> 
  dplyr::select(sid) |> 
  dplyr::mutate(dataset = "G-SAM")
stopifnot(nrow(tmp.sel.gsam) == 287)


tmp.sel.glass <- glass.gbm.rnaseq.metadata.all.samples |>
  dplyr::filter(`tumour.percentage.2022` >= 15) |>
  dplyr::rename(sid = aliquot_barcode) |>
  dplyr::select(sid) |>
  dplyr::mutate(dataset = "GLASS")
stopifnot(nrow(tmp.sel.glass) == 216)


#tmp.gsam.gene.expression.all <- gsam.rnaseq.expression %>%
#  dplyr::select(tmp.gsam.metadata.all$sid)
#stopifnot(colnames(tmp.gsam.gene.expression.all) == tmp.gsam.metadata.all$sid)



## Combined ----


tmp.combined.metadata <- rbind(tmp.sel.gsam, tmp.sel.glass) |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |>
      dplyr::select(aliquot_barcode, aliquot_batch_synapse),
    by=c('sid'='aliquot_barcode')
  ) |> 
  dplyr::mutate(batch = ifelse(dataset == "G-SAM","G-SAM",aliquot_batch_synapse)) |> 
  dplyr::mutate(aliquot_batch_synapse = NULL ) |> 
  dplyr::mutate(batch = gsub("[",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = gsub("]",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = gsub(" ",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = gsub("-",".",batch,fixed=T)) |> 
  dplyr::mutate(batch = as.factor(batch))


stopifnot(nrow(tmp.combined.metadata) == 287+216)




tmp.combined.gene.expression <- dplyr::inner_join(
  gsam.rnaseq.expression %>%
    dplyr::select(
      tmp.combined.metadata |>
        dplyr::filter(dataset == "G-SAM") |>
        dplyr::pull(sid)
    ) %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::mutate(gid = gsub('^(ENSG[0-9]+).+$','\\1', gid))
  ,
  glass.gbm.rnaseq.expression.all.samples %>%
    dplyr::select(
      tmp.combined.metadata |>
        dplyr::filter(dataset == "GLASS") |>
        dplyr::pull(sid)
    ) %>%
    tibble::rownames_to_column('gid'),
  by=c('gid'='gid') ) %>%
  tibble::column_to_rownames('gid') %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  limma::removeBatchEffect(tmp.combined.metadata$batch) %>% # remove batch effect :)
  as.data.frame()
stopifnot(tmp.combined.metadata$sid == colnames(tmp.combined.gene.expression))


rm(tmp.sel.glass, tmp.sel.gsam)



# C0: Fuzzy ----


tmp.c0.fuzzy <- tmp.combined.gene.expression |> 
  tibble::rownames_to_column('ens') |> 
  dplyr::filter(ens %in% (results.out |> dplyr::filter(C0.2022) |> dplyr::pull(ensembl_id)) ) |> 
  tibble::column_to_rownames('ens') |> 
  t() |> 
  prcomp()


plot(tmp.c0.fuzzy)
ggbiplot::ggbiplot(tmp.c0.fuzzy)
# tmp.c0.fuzzy$rotation[1,] << see if most are + or -


signature.c0.fuzzy <- tmp.c0.fuzzy$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::select(c('sid','PC1')) %>%
  dplyr::mutate(inverse = sum(tmp.c0.fuzzy$rotation[,1] < 0) > sum(tmp.c0.fuzzy$rotation[,1] > 0)) |> 
  dplyr::mutate(component = ifelse(inverse , -PC1, PC1) ) |> 
  dplyr::mutate(inverse = NULL, PC1 = NULL) |> 
  dplyr::rename(rna.signature.C0.fuzzy.2022 = component)

rm(tmp.c0.fuzzy)



## plt ----


plt <- gsam.rna.metadata |> 
  dplyr::mutate(rna.signature.C0.fuzzy.2022 = NULL) |> 
  dplyr::left_join(signature.c0.fuzzy, by=c('sid'='sid'), suffix=c('','')) |> 
  dplyr::select(sid, rna.signature.C0.fuzzy.2022, C4.signature.2021, C5.signature.2021) |> 
  tidyr::pivot_longer(cols=c(C4.signature.2021, C5.signature.2021)) |>
  dplyr::mutate(name = gsub(".signature.2021","",name)) |> 
  dplyr::rename(signature.2021 = value)

n.gsam <- plt |>
  dplyr::filter(!is.na(`rna.signature.C0.fuzzy.2022`) & !is.na(`signature.2021`)) |> 
  dplyr::pull(sid) |>
  unique() |> 
  length()


ggplot(plt, aes(x = rna.signature.C0.fuzzy.2022, y = signature.2021)) +
  facet_grid(cols = vars(name), scales = "free") +
  geom_point() +
  ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..)) +
  theme_bw() +
  labs(
    x = "C0 signature (current revision)",
    y = "<cluster> signature (initial revision)",
    caption = paste0("G-SAM: n=", n.gsam, " samples")
  )
ggsave("output/figures/2022_reviewers_comment__comparison_rna_signature_C0.pdf", width = 10, height = 5)


rm(plt, n.gsam)



# C1: Collagen ----


tmp.c1.collagen <- tmp.combined.gene.expression |> 
  tibble::rownames_to_column('ens') |> 
  dplyr::filter(ens %in% (results.out |> dplyr::filter(C1.2022) |> dplyr::pull(ensembl_id)) ) |> 
  tibble::column_to_rownames('ens') |> 
  t() |> 
  prcomp()


plot(tmp.c1.collagen)
ggbiplot::ggbiplot(tmp.c1.collagen)
# tmp.c1.collagen$rotation[1,] << see if most are + or -


signature.c1.collagen <- tmp.c1.collagen$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::select(c('sid','PC1')) %>%
  dplyr::mutate(inverse = sum(tmp.c1.collagen$rotation[,1] < 0) > sum(tmp.c1.collagen$rotation[,1] > 0)) |> 
  dplyr::mutate(component = ifelse(inverse , -PC1, PC1) ) |> 
  dplyr::mutate(inverse = NULL, PC1 = NULL) |> 
  dplyr::rename(rna.signature.C1.collagen.2022 = component)

rm(tmp.c1.collagen)


## plt ----


plt <- gsam.rna.metadata |> 
  dplyr::mutate(rna.signature.C1.collagen.2022 = NULL) |> 
  dplyr::left_join(signature.c1.collagen, by=c('sid'='sid'), suffix=c('',''))

n.gsam <- plt |>
  dplyr::filter(!is.na(`rna.signature.C1.collagen.2022`) & !is.na(`extracellular.matrix.component.2021`)) |> 
  dplyr::pull(sid) |>
  unique() |> 
  length()



ggplot(plt, aes(x = rna.signature.C1.collagen.2022, y = extracellular.matrix.component.2021)) + 
  geom_point() +
  ggpubr::stat_cor(method = "pearson",  aes(label = ..r.label..)) +
  theme_bw() +
  labs(x = "C1 signature (current revision)", 
       y = "C6 [Collagen] signature (initial revision)",
       caption = paste0("G-SAM: n=",n.gsam, " samples"))
ggsave("output/figures/2022_reviewers_comment__comparison_rna_signature_C1.pdf", width=5,height=5)


rm(plt, n.gsam)



# C2: Endothelial ----


tmp.c2.endothelial <- tmp.combined.gene.expression |> 
  tibble::rownames_to_column('ens') |> 
  dplyr::filter(ens %in% (results.out |> dplyr::filter(C2.2022) |> dplyr::pull(ensembl_id)) ) |> 
  tibble::column_to_rownames('ens') |> 
  t() |> 
  prcomp()


plot(tmp.c2.endothelial)
ggbiplot::ggbiplot(tmp.c2.endothelial)
# tmp.c2.endothelial$rotation[1,] << see if most are + or -


signature.c2.endothelial <- tmp.c2.endothelial$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::select(c('sid','PC1')) %>%
  dplyr::mutate(inverse = sum(tmp.c2.endothelial$rotation[,1] < 0) > sum(tmp.c2.endothelial$rotation[,1] > 0)) |> 
  dplyr::mutate(component = ifelse(inverse , -PC1, PC1) ) |> 
  dplyr::mutate(inverse = NULL, PC1 = NULL) |> 
  dplyr::rename(rna.signature.C2.endothelial.2022 = component)

rm(tmp.c2.endothelial)



## plt ----

plt <- gsam.rna.metadata |> 
  dplyr::left_join(signature.c2.endothelial, by=c('sid'='sid'), suffix=c('',''))


n.gsam <- plt |>
  dplyr::filter(!is.na(`rna.signature.C2.endothelial.2022`) & !is.na(`endothelial.component.2021`)) |> 
  dplyr::pull(sid) |>
  unique() |> 
  length()



ggplot(plt, aes(x = rna.signature.C2.endothelial.2022, y = - endothelial.component.2021)) + 
  geom_point() +
  ggpubr::stat_cor(method = "pearson",  aes(label = ..r.label..)) +
  theme_bw() +
  labs(x = "C2 signature (current revision)",
       y = "- C3 [Endothelial] signature (initial revision)",
       caption = paste0("G-SAM: n=", n.gsam, " samples"))
ggsave("output/figures/2022_reviewers_comment__comparison_rna_signature_C2.pdf", width=5,height=5)


rm(plt, n.gsam)





# C3: Oligodendrocyte ----


tmp.c3.oligodendrocyte <- tmp.combined.gene.expression |> 
  tibble::rownames_to_column('ens') |> 
  dplyr::filter(ens %in% (results.out |> dplyr::filter(C3.2022) |> dplyr::pull(ensembl_id)) ) |> 
  tibble::column_to_rownames('ens') |> 
  t() |> 
  prcomp()


plot(tmp.c3.oligodendrocyte)
ggbiplot::ggbiplot(tmp.c3.oligodendrocyte)
# tmp.c3.oligodendrocyte$rotation[1,] << see if most are + or -


signature.c3.oligodendrocyte <- tmp.c3.oligodendrocyte$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::select(c('sid','PC1')) %>%
  dplyr::mutate(inverse = sum(tmp.c3.oligodendrocyte$rotation[,1] < 0) > sum(tmp.c3.oligodendrocyte$rotation[,1] > 0)) |> 
  dplyr::mutate(component = ifelse(inverse , -PC1, PC1) ) |> 
  dplyr::mutate(inverse = NULL, PC1 = NULL) |> 
  dplyr::rename(rna.signature.C3.oligodendrocyte.2022 = component)

rm(tmp.c3.oligodendrocyte)



## plt ----


plt <- gsam.rna.metadata |> 
  dplyr::left_join(signature.c3.oligodendrocyte, by=c('sid'='sid'), suffix=c('',''))

n.gsam <- plt |>
  dplyr::filter(!is.na(`rna.signature.C3.oligodendrocyte.2022`) & !is.na(`oligodendrocyte.component.2021`)) |> 
  dplyr::pull(sid) |>
  unique() |> 
  length()




ggplot(plt, aes(x = rna.signature.C3.oligodendrocyte.2022, y = oligodendrocyte.component.2021)) + 
  geom_point() +
  ggpubr::stat_cor(method = "pearson",  aes(label = ..r.label..)) +
  theme_bw() +
  labs(x = "C3 signature (current revision)", 
       y = "C2 [Oligodendrocyte] signature (initial revision)",
       caption = paste0("G-SAM: n=", n.gsam, " samples"))
ggsave("output/figures/2022_reviewers_comment__comparison_rna_signature_C3.pdf", width=5,height=5)



rm(plt, n.gsam)



# C4: Neuron ----


tmp.c4.neuron <- tmp.combined.gene.expression |> 
  tibble::rownames_to_column('ens') |> 
  dplyr::filter(ens %in% (results.out |> dplyr::filter(C4.2022) |> dplyr::pull(ensembl_id)) ) |> 
  tibble::column_to_rownames('ens') |> 
  t() |> 
  prcomp()




plot(tmp.c4.neuron)
ggbiplot::ggbiplot(tmp.c4.neuron)
# tmp.c4.neuron$rotation[1,] << see if most are + or -


signature.c4.neuron <- tmp.c4.neuron$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::select(c('sid','PC1')) %>%
  dplyr::mutate(inverse = sum(tmp.c4.neuron$rotation[,1] < 0) > sum(tmp.c4.neuron$rotation[,1] > 0)) |> 
  dplyr::mutate(component = ifelse(inverse , -PC1, PC1) ) |> 
  dplyr::mutate(inverse = NULL, PC1 = NULL) |> 
  dplyr::rename(rna.signature.C4.neuron.2022 = component)

rm(tmp.c4.neuron)



## plt ----


plt <- gsam.rna.metadata |>
  dplyr::left_join(signature.c4.neuron, by = c("sid" = "sid"), suffix = c("", ""))

n.gsam <- plt |>
  dplyr::filter(!is.na(`rna.signature.C4.neuron.2022`) & !is.na(`neuron.component.2021`)) |> 
  dplyr::pull(sid) |>
  unique() |> 
  length()



ggplot(plt, aes(x = rna.signature.C4.neuron.2022, y = -neuron.component.2021)) +
  geom_point() +
  ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..)) +
  theme_bw() +
  labs(
    x = "C4 signature (current revision)",
    y = "- C1 [Neuron] signature (initial revision)",
    caption = paste0("G-SAM: n=", n.gsam, " samples")  )
ggsave("output/figures/2022_reviewers_comment__comparison_rna_signature_C4.pdf", width = 5, height = 5)




rm(plt, n.gsam)



# combine & export ----


tmp.export <- signature.c0.fuzzy |> 
  dplyr::left_join(signature.c1.collagen, by=c('sid'='sid'), suffix=c('','')) |> 
  dplyr::left_join(signature.c2.endothelial, by=c('sid'='sid'), suffix=c('','')) |> 
  dplyr::left_join(signature.c3.oligodendrocyte, by=c('sid'='sid'), suffix=c('','')) |> 
  dplyr::left_join(signature.c4.neuron, by=c('sid'='sid'), suffix=c('','')) 


write.table(tmp.export , "output/tables/principal_DE_cluster_components_2022.txt")


rm(tmp.export, n.gsam)



# clean-up ----


rm(signature.c0.fuzzy, 
   signature.c1.collagen, 
   signature.c2.endothelial, 
   signature.c3.oligodendrocyte, 
   signature.c4.neuron)

rm(tmp.combined.gene.expression, tmp.combined.metadata)



## indep val ----
# 
# tmp.1 <- read.table("output/tables/principal_DE_cluster_components_2022.txt") |> 
#   dplyr::rename_with( ~ paste0(.x,".new"))
# tmp.2 <- read.table("output/tables/principal_DE_cluster_components_2022.txt.bak") |> 
#   dplyr::rename_with( ~ paste0(.x,".old"))
# 
# 
# plt <- tmp.1 |> 
#   dplyr::left_join(tmp.2, by=c('sid.new'='sid.old'),suffix=c('',''))
# 
# 
# ggplot(plt, aes(x=`rna.signature.C1.collagen.2022.new`,`rna.signature.C1.collagen.2022.old`)) +
#   geom_point()
# 


