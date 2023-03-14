#!/usr/bin/env R


# Select and merge samples ----


tmp.combined.metadata <- rbind(
  gsam.rna.metadata %>%
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |>
    dplyr::select(sid) |>
    dplyr::mutate(dataset = "G-SAM"),
  glass.gbm.rnaseq.metadata.all.samples |>
    dplyr::filter(tumour.percentage.2022 >= 15) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::select(sid) |>
    dplyr::mutate(dataset = "GLASS")
) |> 
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


# 7k run ----

## Select raw expression data ----


tmp.gsam.gene.expression.all <- gsam.rnaseq.expression |> 
  dplyr::select(tmp.combined.metadata |> dplyr::filter(dataset == "G-SAM") |> dplyr::pull(sid)) |> 
  dplyr::filter(rowSums(dplyr::across()) > ncol(dplyr::across()) * 3)
stopifnot(colnames(tmp.gsam.gene.expression.all) == tmp.combined.metadata |> dplyr::filter(dataset == "G-SAM") |> dplyr::pull(sid))


tmp.glass.gene.expression.all <- glass.gbm.rnaseq.expression.all.samples |> 
  dplyr::select(tmp.combined.metadata |> dplyr::filter(dataset == "GLASS") |> dplyr::pull(sid)) |> 
  dplyr::filter(rowSums(dplyr::across()) > ncol(dplyr::across()) * 3)
stopifnot(colnames(tmp.glass.gene.expression.all) == tmp.combined.metadata |> dplyr::filter(dataset == "GLASS") |> dplyr::pull(sid) )




### combine expression ----


tmp.combined.gene.expression <- dplyr::inner_join(
  tmp.gsam.gene.expression.all |> tibble::rownames_to_column('gid') |> dplyr::mutate(gid = gsub('^(ENSG[0-9]+).+$','\\1', gid)), # change to ENS id's
  tmp.glass.gene.expression.all |> tibble::rownames_to_column('gid') |> dplyr::mutate(gid = ifelse(gid == "ENSG00000165659", "ENSG00000276644", gid)), # DACH1 equivalent
  by=c('gid'='gid') ) |> 
  tibble::column_to_rownames('gid') %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  limma::removeBatchEffect(tmp.combined.metadata$batch) %>%  # remove batch effects
  as.data.frame()
stopifnot(tmp.combined.metadata$sid == colnames(tmp.combined.gene.expression))


rm(tmp.gsam.gene.expression.all, tmp.glass.gene.expression.all)




### collect metadata ----





tmp.plt <- rbind(
  gsam.rna.metadata |>
    
    dplyr::filter(blacklist.pca == F) %>%
    dplyr::filter(pat.with.IDH == F) %>%
    dplyr::filter(
      sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F
    ) %>%
    dplyr::filter(tumour.percentage.dna >= 15) |> # avoid NA values
    
    dplyr::mutate(is.primary = resection == "r1") |> 
    dplyr::select(
      sid,
      pid,
      is.primary,
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
      GITS.150.svm.2022.subtype
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
) |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype)) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup()


tmp.plt <- tmp.plt |> 
  dplyr::left_join(
    tmp.plt |>
      dplyr::filter(is.primary) |> 
      dplyr::select(pid, `GITS.150.svm.2022.subtype`) |> 
      dplyr::rename(subtype.r1 = GITS.150.svm.2022.subtype),
    by=c('pid'='pid'), suffix=c('','')
  ) |> 
  dplyr::left_join(
    tmp.plt |>
      dplyr::filter(!is.primary) |> 
      dplyr::select(pid, `GITS.150.svm.2022.subtype`) |> 
      dplyr::rename(subtype.r2 = GITS.150.svm.2022.subtype),
    by=c('pid'='pid'), suffix=c('','')
  )





tmp.expr <- tmp.combined.gene.expression |> 
  as.data.frame() |> 
  dplyr::select(tmp.plt$sid) |> 
  tibble::rownames_to_column('gid') |> 
  dplyr::filter(grepl("ENSG00000116641|ENSG00000135905|ENSG00000143603", gid)) |> # DOCK7 | DOCK10
  tibble::column_to_rownames('gid') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('sid')



plt <- tmp.plt |> 
  dplyr::left_join(tmp.expr, by=c('sid'='sid'), suffix=c('','')) |> 
  dplyr::mutate(mes.transition = case_when(
    subtype.r1 == "Proneural" & subtype.r2 == "Mesenchymal" ~ "PN -> MES",
    subtype.r1 == "Classical" & subtype.r2 == "Mesenchymal" ~ "CL -> MES",
    #subtype.r1 == "Mesenchymal" & subtype.r2 == "Mesenchymal" ~ "MES -> MES",
    T ~ "other"
  )) |> 
  dplyr::mutate(mes.transition = factor(mes.transition, levels=c("other"    ,"CL -> MES",  "PN -> MES"))) |> 
  dplyr::mutate(resection = dplyr::recode(as.character(is.primary),'TRUE'='Primary','FALSE'='Recurrence'))



# DOCK7
ggplot(plt, aes(x=resection, y=`ENSG00000116641`, group=pid)) +
  facet_grid(cols = vars(mes.transition), scales = "free",space="free_x") +
  geom_line(alpha=0.5) +
  geom_point(size=3) + 
  labs(y="DOCK7 - top DE gene MES subtype in snRNA")


t.test(plt |>  dplyr::filter(mes.transition == "other" & resection == "Primary") |> dplyr::pull(ENSG00000135905),
       plt |>  dplyr::filter(mes.transition == "other" & resection == "Recurrence") |> dplyr::pull(ENSG00000135905),
       paired = TRUE, alternative = "two.sided")
t.test(plt |>  dplyr::filter(mes.transition == "CL -> MES" & resection == "Primary") |> dplyr::pull(ENSG00000135905),
       plt |>  dplyr::filter(mes.transition == "CL -> MES" & resection == "Recurrence") |> dplyr::pull(ENSG00000135905),
       paired = TRUE, alternative = "two.sided")
t.test(plt |>  dplyr::filter(mes.transition == "PN -> MES" & resection == "Primary") |> dplyr::pull(ENSG00000135905),
       plt |>  dplyr::filter(mes.transition == "PN -> MES" & resection == "Recurrence") |> dplyr::pull(ENSG00000135905),
       paired = TRUE, alternative = "two.sided")



# DOCK10 - marker of non MES cells but environment
ggplot(plt, aes(x=resection, y=`ENSG00000135905`, group=pid)) +
  facet_grid(cols = vars(mes.transition), scales = "free",space="free_x") +
  geom_line() +
  geom_point(size=3)

t.test(plt |>  dplyr::filter(mes.transition == "other" & resection == "Primary") |> dplyr::pull(ENSG00000116641),
       plt |>  dplyr::filter(mes.transition == "other" & resection == "Recurrence") |> dplyr::pull(ENSG00000116641),
       paired = TRUE, alternative = "two.sided")
t.test(plt |>  dplyr::filter(mes.transition == "CL -> MES" & resection == "Primary") |> dplyr::pull(ENSG00000116641),
       plt |>  dplyr::filter(mes.transition == "CL -> MES" & resection == "Recurrence") |> dplyr::pull(ENSG00000116641),
       paired = TRUE, alternative = "two.sided")
t.test(plt |>  dplyr::filter(mes.transition == "PN -> MES" & resection == "Primary") |> dplyr::pull(ENSG00000116641),
       plt |>  dplyr::filter(mes.transition == "PN -> MES" & resection == "Recurrence") |> dplyr::pull(ENSG00000116641),
       paired = TRUE, alternative = "two.sided")




# # KCNN3
# ggplot(plt, aes(x=resection, y=`ENSG00000143603`, group=pid)) +
#   facet_grid(cols = vars(mes.transition), scales = "free",space="free_x") +
#   geom_line() +
#   geom_point(size=3)
# 
# 
# 

