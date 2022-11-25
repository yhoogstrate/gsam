#!/usr/bin/env R

# load libs ----


library(e1071)

# #install.packages('NMF')
# library(scales)
# 
# library(DESeq2)
# library(MASS)
# library(limma)
# library(fitdistrplus)
# library(patchwork)
# library(circlize)
# library(ggrepel)
# 
# library(splancs) # https://stackoverflow.com/questions/17571602/r-filter-coordinates
# library(combinat)
# 
# library(rlang) # https://adv-r.hadley.nz/environments.html
# 
# library(limma)
# library(igraph)
# 
# library(diagram)


# load data ----


source('scripts/load_G-SAM_expression_data.R')

source('scripts/load_GLASS_data.R')


source('scripts/load_results.out.R')


source('data/wang/msig.library.12.R')  # no license w/ code provided, can't include it in source tree


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



### combine ----


tmp.combined.gene.expression <- dplyr::inner_join(
  tmp.gsam.gene.expression.all |> tibble::rownames_to_column('gid') |> dplyr::mutate(gid = gsub('^(ENSG[0-9]+).+$','\\1', gid)), # change to ENS id's
  tmp.glass.gene.expression.all |> tibble::rownames_to_column('gid') |> dplyr::mutate(gid = ifelse(gid == "ENSG00000165659", "ENSG00000276644", gid)), # DACH1 equivalent
  by=c('gid'='gid') ) |> 
  tibble::column_to_rownames('gid') %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  limma::removeBatchEffect(tmp.combined.metadata$batch) %>%  # remove batch effects
  as.data.frame() |> 
  tibble::rownames_to_column('gid') |>  
  dplyr::filter(gid %in% (results.out |> dplyr::filter(wang.glioma.intrinsic.gene) |> dplyr::pull(ensembl_id))) %>% # take the 150 subtype genes
  tibble::column_to_rownames('gid')
stopifnot(tmp.combined.metadata$sid == colnames(tmp.combined.gene.expression))


rm(tmp.gsam.gene.expression.all, tmp.glass.gene.expression.all)




## NMF:1-4, using Wang-code ----




#if(!file.exists("tmp/analysis_NMF_7k_2022.Rds" ) & F ) {
# avoid re-calculating
# NMF is not a pure mathematical solution
# it (randomly) converges into more or less the same answers
# use the cached data when

#analysis_NMF_7k_2022 <- NMF(as.matrix(tmp.combined.gene.expression ), 4, seed = 123456) # seed 123456 is the one used by their paper
#saveRDS(analysis_NMF_7k_2022, "tmp/analysis_NMF_7k_2022.Rds")
analysis_NMF_7k_2022 <- readRDS("tmp/analysis_NMF_7k_2022.Rds")


tmp.out <- data.frame(
  analysis_NMF_7k_2022  %>%
    purrr::pluck('H') %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`(c('NMF:7k:1','NMF:7k:2','NMF:7k:3','NMF:7k:4')) %>%
    tibble::rownames_to_column('sid'),
  analysis_NMF_7k_2022 |> 
    purrr::pluck('error.v.per.M') |> 
    as.data.frame() |> 
    `colnames<-`('NMF:7k:error.v.per.M'),
  `NMF:7k:membership` = analysis_NMF_7k_2022 |> 
    purrr::pluck('membership') 
) |> 
  dplyr::mutate(NMF.7k.membership = paste0('NMF:7k:', NMF.7k.membership)) |> 
  dplyr::rename("NMF:7k:1" = "NMF.7k.1",
                "NMF:7k:2" = "NMF.7k.2",
                "NMF:7k:3" = "NMF.7k.3",
                "NMF:7k:4" = "NMF.7k.4",
                "NMF:7k:error.v.per.M" = "NMF.7k.error.v.per.M",
                "NMF:7k:membership" = "NMF.7k.membership")
stopifnot(rownames(tmp.out) == tmp.out$sid)

tmp.out <- tmp.out |> 
  tibble::remove_rownames()



## PCA:1-2 over NMF:1,2,4 ----

tmp.pca <- tmp.out %>%
  dplyr::select(c("sid", "NMF:7k:1", "NMF:7k:2", "NMF:7k:4")) %>% # 1 = CL, 2 = PN, 4 = MES
  tibble::column_to_rownames('sid') |> 
  prcomp()


tmp.out <- tmp.out |> 
  dplyr::left_join(
    tmp.pca |> 
      purrr::pluck('x') |> 
      as.data.frame() |> 
      dplyr::rename_with( ~ paste0("NMF:7k:", .x)) |> 
      tibble::rownames_to_column('sid'),
    by=c('sid'='sid'),suffix=c('','')
  )




# 150 run ----

## Select raw expression data ----


tmp.gsam.gene.expression.all <- gsam.rnaseq.expression %>%
  dplyr::select(
    tmp.combined.metadata %>%
      dplyr::filter(dataset == "G-SAM") |>
      dplyr::pull(sid)
  ) %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3)
stopifnot(colnames(tmp.gsam.gene.expression.all) == 
            tmp.combined.metadata %>%
            dplyr::filter(dataset == "G-SAM") |>
            dplyr::pull(sid)
)


tmp.glass.gene.expression.all <- glass.gbm.rnaseq.expression.all.samples %>%
  dplyr::select(
    tmp.combined.metadata %>%
      dplyr::filter(dataset == "GLASS") |>
      dplyr::pull(sid)
  ) %>%
  dplyr::filter(rowSums(.) > ncol(.) * 3)
stopifnot(colnames(tmp.glass.gene.expression.all) ==     tmp.combined.metadata %>%
            dplyr::filter(dataset == "GLASS") |>
            dplyr::pull(sid)
)



### combine ----

#"ENSG00000165659" %in% (results.out |> dplyr::filter(!is.na(TCGA.subtype.marker ))  |> dplyr::pull(ensembl_id))
#"ENSG00000276644" %in% (results.out |> dplyr::filter(!is.na(TCGA.subtype.marker ))  |> dplyr::pull(ensembl_id))

tmp.combined.gene.expression <- dplyr::inner_join(
  tmp.gsam.gene.expression.all %>%
    tibble::rownames_to_column('gid') %>%
    dplyr::mutate(gid = gsub('^(ENSG[0-9]+).+$','\\1', gid)) # change to ENS id's
  ,
  tmp.glass.gene.expression.all %>%
    tibble::rownames_to_column('gid') |> 
    dplyr::mutate(gid = ifelse(gid == "ENSG00000165659", "ENSG00000276644", gid)) # DACH1 equivalent
  ,
  by=c('gid'='gid') ) %>%
  tibble::column_to_rownames('gid') %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::filter(gid %in% (results.out |> dplyr::filter(!is.na(TCGA.subtype.marker)) |> dplyr::pull(ensembl_id))) |> 
  tibble::column_to_rownames('gid') %>%
  limma::removeBatchEffect(tmp.combined.metadata$batch) # remove batch effects
stopifnot(tmp.combined.metadata$sid == colnames(tmp.combined.gene.expression))
stopifnot(nrow(tmp.combined.gene.expression) == 150)

rm(tmp.gsam.gene.expression.all, tmp.glass.gene.expression.all)



#'@todo compare error w/ purirty?


## NMF:1-3 ----



#analysis_NMF_150_2022 <- NMF(as.matrix(tmp.combined.gene.expression), 3, seed = 123456) # seed 123456 is the one used by their paper
#saveRDS(analysis_NMF_150_2022, "tmp/analysis_NMF_150_2022.Rds")
analysis_NMF_150_2022 <- readRDS("tmp/analysis_NMF_150_2022.Rds")

tmp.out <- tmp.out |> 
  dplyr::left_join(
    data.frame(
      analysis_NMF_150_2022  %>%
        purrr::pluck('H') %>%
        t() %>%
        as.data.frame() %>%
        `colnames<-`(c('NMF:150:1','NMF:150:2','NMF:150:3')) %>%
        tibble::rownames_to_column('sid'),
      analysis_NMF_150_2022 |> 
        purrr::pluck('error.v.per.M') |> 
        as.data.frame() |> 
        `colnames<-`('NMF:150:error.v.per.M'),
      `NMF:150:membership` = analysis_NMF_150_2022 |> 
        purrr::pluck('membership') 
    ) |> 
      dplyr::mutate(NMF.150.membership = paste0('NMF:150:', NMF.150.membership)) |> 
      dplyr::rename("NMF:150:1" = "NMF.150.1",
                    "NMF:150:2" = "NMF.150.2",
                    "NMF:150:3" = "NMF.150.3",
                    "NMF:150:error.v.per.M" = "NMF.150.error.v.per.M",
                    "NMF:150:membership" = "NMF.150.membership") |> 
      tibble::remove_rownames(),
    by=c('sid'='sid'), suffix=c('','')
  ) 



#'@todo compare error w/ purirty?


## PCA:1-2 over NMF:1-3 ----



tmp.pca <- tmp.out %>%
  dplyr::select(c("sid", "NMF:150:1", "NMF:150:2", "NMF:150:3")) %>%
  tibble::column_to_rownames('sid') |> 
  prcomp()



# split them as PCA object needs to be stored as well
tmp <- tmp.pca |> 
  purrr::pluck('x') |> 
  as.data.frame() |> 
  dplyr::rename_with( ~ paste0("NMF:150:", .x)) |> 
  tibble::rownames_to_column('sid') |> 
  dplyr::mutate(`NMF:150:PC1.n` = as.numeric(scale(`NMF:150:PC1`))) |> 
  dplyr::mutate(`NMF:150:PC2.n` = as.numeric(scale(`NMF:150:PC2`)))


tmp.out <- tmp.out |> 
  dplyr::left_join( tmp, by=c('sid'='sid'),suffix=c('','') )


# saveRDS(tmp.pca, file="tmp/analysis_GITS_space_GITS_PCA_150.Rds")
rm(tmp.pca, tmp)



## Re-classification within GITs space ----

# train
tmp.labels.1 <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::select(aliquot_barcode, ssGSEA.2022.subtype)  |> 
  dplyr::mutate(ssGSEA.2022.subtype = as.character(ssGSEA.2022.subtype)) |> 
  dplyr::rename(ssGSEA.2022.subtype.from.glass = ssGSEA.2022.subtype)

tmp.labels.2 <-    gsam.rna.metadata |>
  dplyr::select(sid, ssGSEA.2022.subtype) |> 
  dplyr::mutate(ssGSEA.2022.subtype = as.character(ssGSEA.2022.subtype)) |> 
  dplyr::rename(ssGSEA.2022.subtype.from.gsam = ssGSEA.2022.subtype)


tmp.train <- tmp.out |> 
  dplyr::select(c('sid' ,'NMF:150:PC1', 'NMF:150:PC2')) |> 
  dplyr::left_join(tmp.labels.1, by=c('sid'='aliquot_barcode'), suffix=c('','')) |> 
  dplyr::left_join(tmp.labels.2, by=c('sid'='sid'), suffix=c('','')) |> 
  dplyr::mutate(ssGSEA.2022.subtype = ifelse(is.na(ssGSEA.2022.subtype.from.glass), ssGSEA.2022.subtype.from.gsam, ssGSEA.2022.subtype.from.glass)) |> 
  dplyr::select(-contains('.from.g')) |> 
  tibble::column_to_rownames('sid') |> 
  dplyr::filter(grepl("|",ssGSEA.2022.subtype,fixed=T) == F) |> 
  dplyr::mutate(ssGSEA.2022.subtype = as.factor(ssGSEA.2022.subtype))

rm(tmp.labels.1, tmp.labels.2)


# Classifying over training data is not an issue, it's not about performance
# but about finding a nice straight fitting line
set.seed(123456)
# tmp.model <- e1071::svm(x = tmp.train |> dplyr::select(c('NMF:150:PC1', 'NMF:150:PC2')) ,
#                         y = tmp.train |> dplyr::pull('ssGSEA.2022.subtype')
#                      ,
#                      scale = F,
#                      #type = "C-classification",
#                      kernel = 'linear',
#                      tolerance = 0.0001,
#                      cost = 3
#                      #,probability = T # worse fit, but handy values
#                      )
# saveRDS(tmp.model, 'tmp/analysis_GITS_space.svm.model.Rds')
tmp.model <- readRDS('tmp/analysis_GITS_space.svm.model.Rds')
rm(tmp.train)


# re-fit
tmp.fit <- predict(object = tmp.model,
                   newdata = tmp.out |>
                     dplyr::select(c('sid', 'NMF:150:PC1', 'NMF:150:PC2')) |> 
                     tibble::column_to_rownames('sid'),
                   decision.values = T,
                   tolerance = 0.0001,
                   cost = 3)

tmp.fit.df <- data.frame(
  sid = attr(tmp.fit,"names"),
  GITS.150.svm.2022.subtype = as.character(tmp.fit)) |> 
  dplyr::left_join(
    attr(tmp.fit,"decision.values") |> 
      as.data.frame() |> 
      dplyr::rename_with( ~ paste0("GITS.150.svm.2022.", .x)) |> 
      tibble::rownames_to_column('sid'),
    by=c('sid'='sid'),suffix=c('','')
  )


tmp.out <- tmp.out |> 
  dplyr::left_join(tmp.fit.df, by=c('sid'='sid'),suffix=c('',''))


rm(tmp.fit, tmp.fit.df)



## define contours for SVM ----

resolution <- 250 # 1000 x 1000 data points

off_x <- (max(-tmp.out$`NMF:150:PC1`) - min(-tmp.out$`NMF:150:PC1`)) * 0.025
off_y <- (max(-tmp.out$`NMF:150:PC2`) - min(-tmp.out$`NMF:150:PC2`)) * 0.025


range_pc1 = seq(from = min(-tmp.out$`NMF:150:PC1`) - off_x, to = max(-tmp.out$`NMF:150:PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(-tmp.out$`NMF:150:PC2`) - off_y, to = max(-tmp.out$`NMF:150:PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:150:PC1' = range_pc1, 'NMF:150:PC2' = range_pc2)
gits.contours <- data.frame(class = predict(tmp.model , newdata = range_df)) |> 
  cbind(range_df) |> 
  dplyr::select(c('class', 'NMF:150:PC1', 'NMF:150:PC2')) |> 
  dplyr::mutate(type="Contour")

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)


# saveRDS(gits.contours, file="cache/analysis_GITS_space_GITS_contours.Rds")

## eucl dist for pairs ----


tmp.out <- tmp.out |> 
  dplyr::mutate(pid = case_when(
    grepl("GLSS|TCGA", sid) ~ gsub("^(............).+$","\\1",sid),
    T ~ gsub("^(...).+$","\\1",sid)
  ))
  


tmp.paired <- tmp.out |> 
  dplyr::select(sid, pid, `NMF:150:PC1.n`, `NMF:150:PC2.n`, GITS.150.svm.2022.subtype) |> 
  dplyr::mutate(resection = case_when(
    grepl("GLSS|TCGA", sid) ~ ifelse(gsub("^.............(..).+$","\\1",sid) == "TP", "primary", "recurrence"),
    T ~ ifelse(gsub("^...(.).*?$","R\\1",sid) == "R1", "primary", "recurrence")
  )) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() > 1) |> 
  dplyr::ungroup() |> 
  tidyr::pivot_wider(id_cols=  pid,
                     names_from = resection, 
                     values_from = c(sid, `NMF:150:PC1.n`, `NMF:150:PC2.n`, GITS.150.svm.2022.subtype )) |> 
  dplyr::mutate(`NMF:150:PCA:eucledian.dist` =
                  sqrt((`NMF:150:PC1.n_primary` - `NMF:150:PC1.n_recurrence`)^2 +
                       (`NMF:150:PC2.n_primary` - `NMF:150:PC2.n_recurrence`)^2)
                  ) |> 
  #dplyr::mutate(GITS.150.svm.2022.subtype.status = as.factor(ifelse(subtype.public.R1 == subtype.public.R2, "Stable", "Transition"))) %>% # should go into subtype loaders
  ###dplyr::mutate(`NMF:123456.membership.status` = as.factor(ifelse(`NMF:123456.membership.R1`== `NMF:123456.membership.R2`, "Stable", "Transition"))) %>% # not needed
  dplyr::mutate(`GITS.150.svm.2022.subtype.status` = as.factor(ifelse(`GITS.150.svm.2022.subtype_primary` == `GITS.150.svm.2022.subtype_recurrence`, "Stable", "Transition"))) |> 
  dplyr::mutate(
    sid_primary = NULL, 
    sid_recurrence  = NULL,
    `NMF:150:PC1.n_primary` = NULL,
    `NMF:150:PC1.n_recurrence` = NULL,
    `NMF:150:PC2.n_primary` = NULL,
    `NMF:150:PC2.n_recurrence` = NULL,
    GITS.150.svm.2022.subtype_primary = NULL,
    GITS.150.svm.2022.subtype_recurrence = NULL
  )



tmp.out <- tmp.out |> 
  dplyr::left_join(tmp.paired, by=c('pid'='pid'),suffix=c('','')) |> 
  dplyr::mutate(pid = NULL)


rm(tmp.paired)


## eucl dist for shuffled pairs ----


tmp.out <- tmp.out |> 
  dplyr::mutate(pid = case_when(
    grepl("GLSS|TCGA", sid) ~ gsub("^(............).+$","\\1",sid),
    T ~ gsub("^(...).+$","\\1",sid)
  ))


rep <- tmp.out |>
  dplyr::select(sid, pid, `NMF:150:PC1.n`, `NMF:150:PC2.n`) |>
  dplyr::mutate(resection = case_when(
    grepl("GLSS|TCGA", sid) ~ ifelse(gsub("^.............(..).+$", "\\1", sid) == "TP", "primary", "recurrence"),
    T ~ ifelse(gsub("^...(.).*?$", "R\\1", sid) == "R1", "primary", "recurrence")
  )) |>
  dplyr::group_by(pid) |>
  dplyr::filter(n() > 1) |>
  dplyr::ungroup()



shuffle <- tidyr::crossing(
  rep |>
    dplyr::filter(.data$resection == "primary") |> 
    dplyr::select(.data$sid, .data$pid) |>
    dplyr::rename(`sid.R1` = .data$sid) |>
    dplyr::rename(`pid.R1` = .data$pid)
  ,
  rep |>
    dplyr::filter(.data$resection == "recurrence") |> 
    dplyr::select(.data$sid, .data$pid) |>
    dplyr::rename(`sid.R2` = .data$sid) |>
    dplyr::rename(`pid.R2` = .data$pid)
) |> dplyr::filter(.data$pid.R1 != .data$pid.R2) |> 
  dplyr::mutate(pid.R1 = NULL) |> 
  dplyr::mutate(pid.R2 = NULL) |> 
  dplyr::left_join(
    rep |>
      dplyr::filter(.data$resection == "primary") |> 
      dplyr::select(.data$sid, .data$`NMF:150:PC1.n`, .data$`NMF:150:PC2.n`) |> 
      dplyr::rename_with( ~ paste0(.x, ".R1"))
    ,
    by=c('sid.R1'='sid.R1'), suffix=c('','')
  ) |> 
  dplyr::left_join(
    rep |>
      dplyr::filter(.data$resection == "recurrence") |> 
      dplyr::select(.data$sid, .data$`NMF:150:PC1.n`, .data$`NMF:150:PC2.n`) |> 
      dplyr::rename_with( ~ paste0(.x, ".R2"))
    ,
    by=c('sid.R2'='sid.R2'), suffix=c('','')
  )

shuffle <- shuffle |> 
  dplyr::mutate(`NMF:150:PCA:eucledian.dist` =
                  sqrt((`NMF:150:PC1.n.R1` - `NMF:150:PC1.n.R2`)^2 +
                         (`NMF:150:PC2.n.R1` - `NMF:150:PC2.n.R2`)^2))


#saveRDS(shuffle, file = "tmp/PCA.eucledian.distances.shuffled.Rds")


rm(rep, shuffle)


## calc proprtional segments ----

# 1. vind de lijn tussen R1 en R2
# 2. reken de hoek uit
# 3. maak k=150 lijnen met deze hoek, startend vanuit R1, en maak deze structureel langer
#  - predict de klasse hiervan


# x-check GLSS−HF−3050
tmp.out <- tmp.out |> 
  dplyr::mutate(pid = case_when(
    grepl("GLSS|TCGA", sid) ~ gsub("^(............).+$","\\1",sid),
    T ~ gsub("^(...).+$","\\1",sid)
  ))


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

  return(c("x"=x2_new, "y"=y2_new))
    
  # return (data.frame(x = c(x1, x2_new) ,
  #                    y = c(y1, y2_new) ,
  #                    point = as.factor(c("initial start", "new end"))))
}



calc_transitions <- function(pid, x1, y1, x2, y2) {
  # for testing:
  #x1 = 1
  #y1 = 1
  #x2 = 4
  #y2 = 3.5

  d = sqrt((x1 - x2)^2 + (y1 - y2)^2) # dist
  
    
  res <- 150
  df <- data.frame(pct = (0:(res-1) / (res-1))) |>
    dplyr::mutate(line = lapply(pct * d, change_length_line, x1=x1,y1=y1,x2=x2,y2=y2)) |> 
    dplyr::mutate(x2 = unlist(lapply(line, function(x){return(x['x'])} ))) |> 
    dplyr::mutate(y2 = unlist(lapply(line, function(x){return(x['y'])} ))) |> 
    dplyr::mutate(line = NULL) |> 
    tibble::remove_rownames()

  df$pred <- as.character(predict(tmp.model, newdata = df |>
                    dplyr::mutate(pct =  NULL) |>
                    dplyr::rename(`NMF:150:PC1` = x2) |>
                    dplyr::rename(`NMF:150:PC2` = y2)
                  ))

  df <- rbind(
    df[1,] |> dplyr::mutate(pred = "*init"), 
    df,
    df[res,] |> dplyr::mutate(pred = "end*")
  ) |>
    tibble::remove_rownames()
  

  df.transitions <- df |> 
    dplyr::mutate(pct_from = lag(pct),
                  x2_from = lag(x2),
                  y2_from = lag(y2),
                  pred_from = lag(pred)) |> 
    dplyr::rename(pct_to = pct,
                  x2_to = x2,
                  y2_to = y2,
                  pred_to = pred) |> 
    dplyr::mutate(transition.status = ifelse(pred_from != pred_to, T ,F)) |> 
    dplyr::mutate(x2_transition = (x2_from + x2_to) / 2) |> 
    dplyr::mutate(y2_transition = (y2_from + y2_to) / 2) |> 
    dplyr::mutate(pct_transition = (pct_from + pct_to) / 2) |> 
    dplyr::filter(transition.status)
  
  df.transitions <- df.transitions |> 
    tidyr::pivot_longer(cols = c(`pct_to`,`x2_to`,`y2_to`,`pred_to`,`pct_from`,`x2_from`,`y2_from`,`pred_from`),
                        names_to = c('.value','status'),
                        names_pattern="(.+)_(.+)") |> 
    dplyr::filter(pred %in% c('*init', 'end*') == F) |> 
    dplyr::arrange(pct) |>
    dplyr::mutate(segment= paste0( "segment.",pid,".",((1:n())+1) %/% 2)) |> 
    dplyr::mutate(pid = pid) |> 
    dplyr::mutate(
      transition.status = NULL,
      status = NULL,
      pct = NULL,  
      x2 = NULL,
      y2 = NULL
    )
  
  
  return(df.transitions)
}
#calc_transitions("id", 1,1,4,3.5)
calc_transitions("id", -1.237046,   0.5114374,
                       -0.2791471, -0.8142251) # GLSS-HF-3050



tmp.paired <- tmp.out |> 
  dplyr::select(sid, pid,
                `NMF:150:PC1`, `NMF:150:PC2`,
                `NMF:150:PCA:eucledian.dist`,
                GITS.150.svm.2022.subtype) |> 
  dplyr::mutate(resection = case_when(
    grepl("GLSS|TCGA", sid) ~ ifelse(gsub("^.............(..).+$","\\1",sid) == "TP", "primary", "recurrence"),
    T ~ ifelse(gsub("^...(.).*?$","R\\1",sid) == "R1", "primary", "recurrence")
  )) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() > 1) |> 
  dplyr::ungroup() |> 
  tidyr::pivot_wider(id_cols=  c(pid),
                     names_from = resection, 
                     values_from = c(sid, `NMF:150:PC1`, `NMF:150:PC2`, GITS.150.svm.2022.subtype, `NMF:150:PCA:eucledian.dist` )) |> 
  dplyr::rename(`NMF:150:PCA:eucledian.dist` = `NMF:150:PCA:eucledian.dist_primary`) |>  # manual has no clear way of keeping columns in pivot_wider?
  dplyr::mutate(`NMF:150:PCA:eucledian.dist_recurrence` = NULL) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(transition.segments = list(calc_transitions(pid, `NMF:150:PC1_primary`, `NMF:150:PC2_primary`,
                            `NMF:150:PC1_recurrence`,`NMF:150:PC2_recurrence`))) |> 
  dplyr::ungroup()


segments <- data.table::rbindlist(tmp.paired$transition.segments) |> 
  as.data.frame()



# check if first matches metadata
tmp <- segments |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(pct_transition == min(pct_transition)) |> 
  dplyr::select(pid, pred) |> 
  dplyr::left_join(tmp.paired |> dplyr::select(pid, GITS.150.svm.2022.subtype_primary), by=c('pid'='pid'))
stopifnot(tmp$pred == tmp$GITS.150.svm.2022.subtype_primary)
rm(tmp)

tmp <- segments |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(pct_transition == max(pct_transition)) |> 
  dplyr::select(pid, pred) |> 
  dplyr::left_join(tmp.paired |> dplyr::select(pid, GITS.150.svm.2022.subtype_recurrence), by=c('pid'='pid'))
stopifnot(tmp$pred == tmp$GITS.150.svm.2022.subtype_recurrence)
rm(tmp)


saveRDS(segments, 'tmp/analysis_GITS_space.transition.segments.Rds')
rm(segments, tmp.paired)




# export ----


saveRDS(
  tmp.out |> dplyr::select(
    `sid`,
    `NMF:150:1`,
    `NMF:150:2`,
    `NMF:150:3`,
    `NMF:150:error.v.per.M`,
    `NMF:150:membership`,
    
    `NMF:150:PC1`,
    `NMF:150:PC2`,
    `NMF:150:PC3`,
    
    `NMF:150:PC1.n`,
    `NMF:150:PC2.n`,
    
    `NMF:7k:1`,
    `NMF:7k:2`,
    `NMF:7k:3`,
    `NMF:7k:4`,
    `NMF:7k:error.v.per.M`,
    `NMF:7k:membership`,
    
    `NMF:7k:PC1`,
    `NMF:7k:PC2`,
    `NMF:7k:PC3`,
    
    `GITS.150.svm.2022.subtype`,
    `GITS.150.svm.2022.Classical/Mesenchymal`,
    `GITS.150.svm.2022.Classical/Proneural`,
    `GITS.150.svm.2022.Mesenchymal/Proneural`,
    
    `NMF:150:PCA:eucledian.dist`,
    `GITS.150.svm.2022.subtype.status`
  ),
  file = "cache/analysis_GITS_space.Rds"
)




# contour diff revisions ----

# 1. subtypes as in initial manuscript
# - G-SAM: gliovis.majority_call
# - GLASS: GBM.transcriptional.subtype.Synapse.2021

stopifnot(nrow(tmp.combined.metadata) == (287+216))


fv <- Vectorize(function(x,y){
  x <- unique(strsplit(x, "\\|")[[1]])
  x[x != "Mesenchymal"] <- "Non-Mesenchymal"
  
  y <- unique(strsplit(y, "\\|")[[1]])
  y[y != "Mesenchymal"] <- "Non-Mesenchymal"
  
  isct <- intersect(x,y)
  
  x <- setdiff(x, isct)
  y <- setdiff(y, isct)
  
  if(length(x) > 0) {
    x <- paste0(x, " (first submission)")
  }
  if(length(y) > 0) {
    y <- paste0(y, " (revision)")
  }
  
  z <- c(isct, sort(unique(c(x,y))))
  
  return(as.character(paste0(z, collapse="|")))
})

fv("ASD","DEF")
fv("ASD","ASD")
fv("ASD|DEF","DEF")
fv("ASD","DEF|DEF")
fv("ASD|DEF","ASD|DEF")


## metadata table ----

tmp.train <- tmp.combined.metadata |>
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |>
      dplyr::select(
        aliquot_barcode,
        case_barcode,
        resection,

        # subtype classes:
        ssGSEA.2022.subtype,
        GBM.transcriptional.subtype.Synapse.2021,

        # coordinates:
        "NMF:150:PC1", "NMF:150:PC2"
      ) |>
      dplyr::rename(`NMF:150:PC1_glass` = `NMF:150:PC1`) |>
      dplyr::rename(`NMF:150:PC2_glass` = `NMF:150:PC2`) |> 
      dplyr::rename(`subtype.first.draft_glass` = GBM.transcriptional.subtype.Synapse.2021) |> 
      dplyr::rename(`subtype.revision_glass` = ssGSEA.2022.subtype)  |> 
      dplyr::rename(pid_glass = case_barcode) |> 
      dplyr::mutate(resection_glass = ifelse(resection == "TP", "primary","recurrence"), resection=NULL)
    ,
    by = c("sid" = "aliquot_barcode"), suffix = c("", "")
  ) |>
  dplyr::mutate(dataset = NULL) |>
  dplyr::mutate(batch = NULL) |>
  dplyr::mutate(ssGSEA.2022.subtype = NULL) |>
  dplyr::left_join(
    gsam.rna.metadata |>
      dplyr::select(
        sid,
        pid,
        resection, 
        
        # subtype classes:
        ssGSEA.2022.subtype,
        gliovis.majority_call,

        # coordinates:
        "NMF:150:PC1", "NMF:150:PC2"
      ) |>
      dplyr::rename(pid_gsam = pid) |> 
      dplyr::rename(`NMF:150:PC1_gsam` = `NMF:150:PC1`) |>
      dplyr::rename(`NMF:150:PC2_gsam` = `NMF:150:PC2`) |> 
      dplyr::rename(`subtype.first.draft_gsam` = gliovis.majority_call) |> 
      dplyr::rename(`subtype.revision_gsam` = ssGSEA.2022.subtype) |> 
      dplyr::mutate(resection_gsam = ifelse(resection == "r1", "primary","recurrence"), resection=NULL)
    ,
    by = c("sid" = "sid"), suffix = c("", "")
  ) |>
  dplyr::mutate(subtype.first.draft_gsam = as.character(subtype.first.draft_gsam)) |> 
  dplyr::mutate(subtype.revision_gsam = as.character(subtype.revision_gsam)) |> 
  dplyr::mutate(subtype.first.draft_glass = as.character(subtype.first.draft_glass)) |> 
  dplyr::mutate(subtype.revision_glass = as.character(subtype.revision_glass)) |> 
  
  dplyr::mutate(subtype.first.draft = ifelse(is.na(subtype.first.draft_gsam), subtype.first.draft_glass, subtype.first.draft_gsam)) |>
  dplyr::mutate(subtype.revision = ifelse(is.na(subtype.revision_gsam), subtype.revision_glass, subtype.revision_gsam)) |>
  dplyr::mutate(resection = ifelse(is.na(resection_gsam), resection_glass, resection_gsam)) |>
  dplyr::mutate(pid = ifelse(is.na(pid_gsam), pid_glass, pid_gsam)) |>
  
  dplyr::mutate(`NMF:150:PC1` = ifelse(is.na(`NMF:150:PC1_gsam`), `NMF:150:PC1_glass`, `NMF:150:PC1_gsam`)) |>
  dplyr::mutate(`NMF:150:PC2` = ifelse(is.na(`NMF:150:PC2_gsam`), `NMF:150:PC2_glass`, `NMF:150:PC2_gsam`)) |> 
  dplyr::select(!contains("_gsam") & !contains("_glass")) |> 
  #dplyr::mutate(subtype_label = fv(ssGSEA.2022.subtype, subtypeme.gliovis.majority.call.2022)) |> 
  
  dplyr::mutate(subtype = fv(subtype.first.draft, subtype.revision)) |> 

  dplyr::mutate(subtype.first.draft = as.factor(subtype.first.draft)) |> 
  dplyr::mutate(subtype.revision = as.factor(subtype.revision)) |> 
  dplyr::mutate(revisions = ifelse(is.na(subtype.first.draft),"only revised manuscript","both")) |> 

  tibble::column_to_rownames('sid')


## contours first revision ----

tmp.model.first.draft <- e1071::svm(x = tmp.train |> dplyr::select(c('NMF:150:PC1', 'NMF:150:PC2')) ,
                                    y = tmp.train |> dplyr::pull('subtype.first.draft')
                                    ,
                                    scale = F,
                                    #type = "C-classification",
                                    kernel = 'linear',
                                    tolerance = 0.0001,
                                    cost = 3
                                    )


resolution <- 250 # 1000 x 1000 data points

off_x <- (max(-tmp.train$`NMF:150:PC1`) - min(-tmp.train$`NMF:150:PC1`)) * 0.025
off_y <- (max(-tmp.train$`NMF:150:PC2`) - min(-tmp.train$`NMF:150:PC2`)) * 0.025


range_pc1 = seq(from = min(-tmp.train$`NMF:150:PC1`) - off_x, to = max(-tmp.train$`NMF:150:PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(-tmp.train$`NMF:150:PC2`) - off_y, to = max(-tmp.train$`NMF:150:PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:150:PC1' = range_pc1, 'NMF:150:PC2' = range_pc2)
gits.contours.first.draft <- data.frame(class = predict(tmp.model.first.draft , newdata = range_df)) |> 
  cbind(range_df) |> 
  dplyr::select(c('class', 'NMF:150:PC1', 'NMF:150:PC2')) |> 
  dplyr::mutate(type="Contour") |> 
  dplyr::mutate(class = as.character(class))

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)


## contours revision ----

tmp.model.revision <- e1071::svm(x = tmp.train |> dplyr::select(c('NMF:150:PC1', 'NMF:150:PC2')) ,
                                 y = tmp.train |> dplyr::pull('subtype.revision')
                                 ,
                                 scale = F,
                                 #type = "C-classification",
                                 kernel = 'linear',
                                 tolerance = 0.0001,
                                 cost = 3)

resolution <- 250 # 1000 x 1000 data points

off_x <- (max(-tmp.train$`NMF:150:PC1`) - min(-tmp.train$`NMF:150:PC1`)) * 0.025
off_y <- (max(-tmp.train$`NMF:150:PC2`) - min(-tmp.train$`NMF:150:PC2`)) * 0.025


range_pc1 = seq(from = min(-tmp.train$`NMF:150:PC1`) - off_x, to = max(-tmp.train$`NMF:150:PC1`) + off_x, length.out = resolution)
range_pc2 = seq(from = min(-tmp.train$`NMF:150:PC2`) - off_y, to = max(-tmp.train$`NMF:150:PC2`) + off_y, length.out = resolution)

range_df = expand.grid('NMF:150:PC1' = range_pc1, 'NMF:150:PC2' = range_pc2)
gits.contours.revision <- data.frame(class = predict(tmp.model.revision , newdata = range_df)) |> 
  cbind(range_df) |> 
  dplyr::select(c('class', 'NMF:150:PC1', 'NMF:150:PC2')) |> 
  dplyr::mutate(type="Contour") |> 
  dplyr::mutate(class = as.character(class))

rm(resolution, off_x, off_y, range_pc1, range_pc2, range_df)




## integrate / diff contours ----

stopifnot(gits.contours.first.draft$`NMF:150:PC1` == gits.contours.revision$`NMF:150:PC1`)
stopifnot(gits.contours.first.draft$`NMF:150:PC2` == gits.contours.revision$`NMF:150:PC2`)
stopifnot(sum(gits.contours.first.draft$class != gits.contours.revision$class) > 0)


gits.contours <- data.frame(
  `NMF:150:PC1` = gits.contours.first.draft$`NMF:150:PC1`,
  `NMF:150:PC2` = gits.contours.first.draft$`NMF:150:PC2`,
  class.first.draft = gits.contours.first.draft$class,
  class.revision = gits.contours.revision$class,
  check.names=F
) |> 
  dplyr::mutate(class.first.draft = as.character(class.first.draft)) |> 
  dplyr::mutate(class.revision = as.character (class.revision)) |> 
  #dplyr::filter(class.first.draft != class.revision) |> 
  dplyr::mutate(class = fv(class.first.draft, class.revision))



## plt ----



primary.mes.first.draft <- tmp.train |> 
  dplyr::filter(resection == "primary" & subtype.first.draft == "Mesenchymal") |> 
  dplyr::pull(pid)

primary.mes.revision <- tmp.train |> 
  dplyr::filter(resection == "primary" & subtype.revision == "Mesenchymal") |> 
  dplyr::pull(pid)


plt <- tmp.train |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(primary.mes.first.draft = pid %in% primary.mes.first.draft) |> 
  dplyr::mutate(primary.mes.revision = pid %in% primary.mes.revision)




plt.contours <- gits.contours |> 
  dplyr::mutate(`subtype` = class, `pid` = NA)



p1 <- ggplot(plt |> 
               dplyr::filter(resection == "primary" &
                               revisions == "both"
               ), aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`)) +
  theme_bw() +
  geom_raster(aes(fill=subtype), data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)),
               colour="black",
               size=1.25,
               lty=2,
               breaks=c(1.5,2.5)*0.00008) +
  # geom_point(size=3, 
  #            fill = "white",
  #            col='black',
  #            pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)"
  ) +
  scale_fill_manual(values = c(
    subtype_colors_ext['Mesenchymal'],
    'Non-Mesenchymal'= as.character("white"),
    "Mesenchymal (revision)|Non-Mesenchymal (first submission)" = 'red',
    'TRUE'='blue',
    'FALSE'='gray90'),
    labels=c('TRUE'='primary tumor MES', 'FALSE'='primary tumor non-MES')
  ) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  labs(subtitle = "GITS space at revision, classifiers by GlioVis + Synapse (first submission) or ssGSEA (revision)")




p2 <- ggplot(plt |> 
               dplyr::filter(resection == "recurrence" &
                               revisions == "both"
               ), aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`)) +
  theme_bw() +
  geom_raster(aes(fill=subtype), data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)),
               colour="black",
               size=1.25,
               lty=2,
               breaks=c(1.5,2.5)*0.00008) +
  geom_point(size=3, 
             aes(fill = primary.mes.first.draft), 
             #fill = "gray",
             col='black',
             pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)"
       #caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
  ) +
  scale_fill_manual(values = c(
    subtype_colors_ext['Mesenchymal'],
    'Non-Mesenchymal'= as.character("white"),
    "Mesenchymal (revision)|Non-Mesenchymal (first submission)" = 'red',
    'TRUE'='blue',
    'FALSE'='gray90'),
    labels=c('TRUE'='primary tumor MES', 'FALSE'='primary tumor non-MES')
  ) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  labs(subtitle = "Recurrent tumors, both manuscript versions, MES status @ first submission")





p3 <- ggplot(plt |> 
               dplyr::filter(resection == "recurrence" &
                               revisions == "both"
               ), aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=subtype)) +
  
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,
               breaks=c(1.5,2.5)) +
  geom_point(size=3, 
             aes(fill = primary.mes.revision), 
             #fill = "gray",
             col='black',
             pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)"
  ) +
  scale_fill_manual(values = c(
    subtype_colors_ext['Mesenchymal'],
    'Non-Mesenchymal'= as.character("white"),
    "Mesenchymal (revision)|Non-Mesenchymal (first submission)" = 'red',
    'TRUE'='blue',
    'FALSE'='gray90'),
    labels=c('TRUE'='primary tumor MES', 'FALSE'='primary tumor non-MES')
  ) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  labs(subtitle = "Recurrent tumors, both manuscript versions, MES status @ revision")



p4 <- ggplot(plt |> 
               dplyr::filter(resection == "primary" &
                               revisions == "only revised manuscript"
               ), aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=subtype)) +
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,
               breaks=c(1.5,2.5)) +
  geom_point(size=3, 
             #aes(fill = primary.mes.first.draft), 
             fill = "white",
             col='black',
             pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)"
       #caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
  ) +
  scale_fill_manual(values = c(
    subtype_colors_ext['Mesenchymal'],
    'Non-Mesenchymal'= as.character("white"),
    "Mesenchymal (revision)|Non-Mesenchymal (first submission)" = 'red',
    'TRUE'='blue',
    'FALSE'='gray90'),
    labels=c('TRUE'='primary tumor MES', 'FALSE'='primary tumor non-MES')
  ) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  labs(subtitle = "Primary tumors, appended in revision")



p5 <- ggplot(plt |> 
               dplyr::filter(resection == "recurrence" &
                               revisions == "only revised manuscript"
               ), aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=subtype)) +
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,
               breaks=c(1.5,2.5)) +
  geom_point(size=3, 
             aes(fill = primary.mes.first.draft), 
             #fill = "gray",
             col='black',
             pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)"
       #caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
  ) +
  scale_fill_manual(values = c(
    subtype_colors_ext['Mesenchymal'],
    'Non-Mesenchymal'= as.character("white"),
    "Mesenchymal (revision)|Non-Mesenchymal (first submission)" = 'red',
    'TRUE'='blue',
    'FALSE'='gray90'),
    labels=c('TRUE'='primary tumor MES', 'FALSE'='primary tumor non-MES')
  ) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  labs(subtitle = "Primary tumors, appended in revision")



p6 <- ggplot(plt |> 
               dplyr::filter(resection == "recurrence" & revisions == "only revised manuscript"),
             
             
             aes(x=-`NMF:150:PC1`, y=-`NMF:150:PC2`, fill=subtype)) +
  theme_bw() +
  geom_raster(data = plt.contours, alpha=0.15) +
  geom_contour(data=  plt.contours,
               aes(z=as.numeric(class)),
               colour="gray40",
               size=0.25,
               lty=2,
               breaks=c(1.5,2.5)) +
  geom_point(size=3, 
             aes(fill = primary.mes.revision), 
             #fill = "gray",
             col='black',
             pch=21, alpha=0.7) +
  coord_equal() +
  labs(x="PC1 on NMF meta-features", 
       y="PC2 on NMF meta-features",
       fill = "Subtype (ssGSEA)"
       #caption = paste0("G-SAM: n=",n.gsam, "  -  GLASS: n=",n.glass," samples")
  ) +
  scale_fill_manual(values = c(
    subtype_colors_ext['Mesenchymal'],
    'Non-Mesenchymal'= as.character("white"),
    "Mesenchymal (revision)|Non-Mesenchymal (first submission)" = 'red',
    'TRUE'='blue',
    'FALSE'='gray90'),
    labels=c('TRUE'='primary tumor MES', 'FALSE'='primary tumor non-MES')
  ) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  guides(fill=guide_legend(ncol=2)) +
  labs(subtitle = "Recurrent tumors, appended in revision")



(p1 + p2 + p3) / (p4 + p5 + p6)




# vis 2e ----


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


gits.contours.dual <- readRDS("tmp/analysis_GITS_space_GITS_contours_dual.Rds") |> 
  dplyr::mutate(class = ifelse(
    
  ))
plt.contours <- gits.contours.dual |> 
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
  #scale_fill_manual(values = mixcol(subtype_colors_ext,rep("black",length(subtype_colors_ext)),0.1),
  #                  label=c('Mesenchymal'='MES','Proneural'='PN','Classical'='CL','Proneural|Classical'='PN|CL')) +
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



### show the transitions


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
    dplyr::mutate(dataset = "G-SAM") |> 
    dplyr::mutate(GBM.transcriptional.subtype.Synapse.2021 = NA)
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
      GITS.150.svm.2022.subtype,
      GBM.transcriptional.subtype.Synapse.2021
    ) |>
    dplyr::rename(sid = aliquot_barcode) |>
    dplyr::rename(pid = case_barcode) |> 
    dplyr::mutate(dataset = "GLASS")
)  |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup() |> 
  dplyr::select(-`NMF:150:PC1`, -`NMF:150:PC2`, -`ssGSEA.2022.subtype`)

plt.paired <- plt |> 
  tidyr::pivot_wider(id_cols = pid, names_from = is.primary, values_from = c(GITS.150.svm.2022.subtype)) |> 
  dplyr::left_join(
    plt |> dplyr::select(pid, dataset, GBM.transcriptional.subtype.Synapse.2021) |> 
      dplyr::distinct(),
    by=c('pid'='pid'),suffix=c('','')
  ) |> 
  dplyr::rename(subtype_primary = `TRUE`)|> 
  dplyr::rename(subtype_recurrence = `FALSE`) |> 
  dplyr::filter(grepl("Mes",subtype_primary)) |> 
  dplyr::mutate(facet = case_when(
    dataset == "G-SAM" ~ "G-SAM",
    dataset != "G-SAM" & !is.na(GBM.transcriptional.subtype.Synapse.2021) ~ "GLASS (old batch)",
    dataset != "G-SAM" & is.na(GBM.transcriptional.subtype.Synapse.2021) ~ "GLASS (new samples)"
  ))


ggplot(plt.paired, aes(x=subtype_recurrence)) +
  facet_grid(cols = vars(facet)) +
  geom_bar(position = "stack") +
  labs(x = "GITS subtype at recurrence")




