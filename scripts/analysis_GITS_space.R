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


#source('scripts/load_G-SAM_metadata.R')
source('scripts/load_G-SAM_expression_data.R')

source('scripts/load_GLASS_data.R')


source('scripts/load_results.out.R')


source('data/wang/msig.library.12.R')  # no license w/ code provided, can't include it in source tree

#source('scripts/R/subtype_genes.R')
#source('scripts/R/wang_glioma_intrinsic_genes.R')


#'@todo move to vis_GITS_space.R
# source('scripts/R/job_gg_theme.R')
# source('scripts/R/youri_gg_theme.R')
# source('scripts/R/palette.R')


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


#saveRDS(gits.contours, file="cache/analysis_GITS_space_GITS_contours.Rds")

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


## calc proprtional segments ----

# 1. vind de lijn tussen R1 en R2
# 2. reken de hoek uit
# 3. maak k=150 lijnen met deze hoek, startend vanuit R1, en maak deze structureel langer
#  - predict de klasse hiervan



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



calc_transitions <- function(x1, y1, x2, y2) {
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
    df[1,] |> 
      dplyr::mutate(pred = "*init")
    ,
    df,
    df[res,]|> 
      dplyr::mutate(pred = "end*")
  ) |>
    tibble::remove_rownames() 
  

  df.transitions <- df |> 
    dplyr::mutate(transition = ifelse(pred != lag(pred), T ,F))|> 
    dplyr::mutate(transition.previous_pred = lag(pred)) |> 
    dplyr::mutate(transition.boundary_x2 = (lag(x2) + x2) / 2) |> 
    dplyr::mutate(transition.boundary_y2 = (lag(y2) + y2) / 2) |> 
    dplyr::mutate(transition.pct = (lag(pct) + pct) / 2) |> 
    dplyr::mutate(transition.type = paste0(transition.previous_pred, "->",pred)) |> 
    dplyr::select(transition, transition.boundary_x2, transition.boundary_y2, transition.pct, transition.type) |> 
    dplyr::filter(transition == T) |> 
    dplyr::mutate(transition = NULL)

  return(df.transitions)
}
#calc_transitions(1,1,4,3.5)



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
  dplyr::mutate(transition.segments = list(calc_transitions(`NMF:150:PC1_primary`, `NMF:150:PC2_primary`,
                            `NMF:150:PC1_recurrence`,`NMF:150:PC2_recurrence`))) |> 
  dplyr::ungroup()




# tmp.out <- tmp.out |> 
#   dplyr::left_join(tmp.paired, by=c('pid'='pid'),suffix=c('','')) |> 
#   dplyr::mutate(pid = NULL)
# 
# 
# rm(tmp.paired)
# 
# 


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


# 〰 © Dr. Youri Hoogstrate 〰 ----

