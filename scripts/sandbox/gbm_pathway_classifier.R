#!/usr/bin/env

# lib

library(tidyverse)

# export ----

tmp <- read.csv('/tmp/gsam.rnaseq.expression-322.csv',check.names = F)
colnames(tmp)[1] <- 'gid'
tmp <- tmp |> 
  dplyr::mutate(gid.clean = gsub("^.+\\|(.+)\\|.+$","\\1",gid))


excl <- tmp |>
  dplyr::filter(duplicated(gid.clean)) |> 
  dplyr::pull(gid.clean) |> 
  unique()


tmp <- tmp |> 
  dplyr::filter(gid.clean %in% excl == F) |> 
  tibble::column_to_rownames('gid.clean') |>
  dplyr::mutate(gid = NULL ) |> 
  dplyr::filter(rowSums(dplyr::across()) > ncol(dplyr::across()) * 3)


tmp.vst <- tmp %>%
  DESeq2::DESeqDataSetFromMatrix(data.frame(cond = as.factor(paste0("r",round(runif( ncol(.) , 1, 2))))) , ~cond) %>%
  DESeq2::vst(blind=T) %>%
  SummarizedExperiment::assay() %>%
  as.data.frame()


dim(tmp) == dim(tmp.vst)

write.csv(tmp.vst, "output/tables/G-SAM.expression.322.vst.csv")

View(tmp.vst)



slice_size <- 20
for(i in 1:ceiling(322/slice_size)) {
  slice <- tmp.vst[,] |> 
    tibble::rownames_to_column('GeneName') |> 
    dplyr::arrange(GeneName) |> 
    tibble::column_to_rownames('GeneName')
  
  s_pos <- 1 + ((i - 1) * slice_size)
  e_pos <- ((i) * slice_size)
  e_pos <- min(e_pos,  ncol(tmp.vst) )
  
  #print(c(s_pos, e_pos))
  
  slice <- slice[,s_pos:e_pos] |> 
    tibble::rownames_to_column('GeneName') 
  
  print(dim(slice))
  
  write.table(slice, file=paste0("/tmp/GSAM-322-vst_",s_pos,"-",e_pos,".txt"), sep="\t",quote=F, row.names=F)
}


# import results ----
# https://lucgar88.shinyapps.io/GBMclassifier/
# set to FFPE

source('scripts/load_G-SAM_metadata.R')

included <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::filter(tumour.percentage.dna >= 15)  |> 
  dplyr::pull(sid)

gbm.class <- do.call(rbind, lapply(sort(Sys.glob("output/tables/gbm-pathway-based-classifier/Table*.txt")), read.table, header=T))



plt <- gbm.class |> 
  dplyr::mutate(sid = gsub(".","-",TumorID, fixed=T), TumorID = NULL) |> 
  dplyr::filter(sid %in% included) |> 
  dplyr::mutate(pat = gsub("^(...).+$","\\1",sid)) |> 
  dplyr::mutate(res = gsub("^...(.).*$","R\\1",sid)) |> 
  dplyr::arrange(pat,sid)

stopifnot(nrow(plt) == 287)
plt |> dplyr::filter(pat=="AAS")


plt <- plt |> 
  dplyr::group_by(pat) |> 
  dplyr::filter(n() == 2)


plt.wide <- plt |> 
  #head(n=12) |> 
  tidyr::pivot_wider(id_cols = pat, 
                     names_from = res, 
                     values_from = c(GPM_score, MTC_score, NEU_score, PPR_score, Simplicity_score, Subtype, MTC_logOddsRatio, sid)) |> 
  as.data.frame()


# plt scores per resection type ----


ggplot(plt |> dplyr::mutate(col=Subtype == "GPM"), aes(x=res, y=GPM_score, group=pat, col=col)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="Classified as GPM")


ggplot(plt |> dplyr::mutate(col=Subtype == "NEU"), aes(x=res, y=NEU_score, group=pat, col=col)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="Classified as NEU")


ggplot(plt |> dplyr::mutate(col=Subtype == "PPR"), aes(x=res, y=PPR_score, group=pat, col=col)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="Classified as PPR")


ggplot(plt |> dplyr::mutate(col=Subtype == "MTC"), aes(x=res, y=MTC_score, group=pat, col=col)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="Classified as MTC")



plt.GPM.NEU <- rbind(
  plt |>
    dplyr::filter(res == "R1") |> 
    dplyr::mutate(score = GPM_score) |> 
    dplyr::mutate(score_name = "GPM")
  ,
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = NEU_score) |> 
    dplyr::mutate(score_name = "NEU")
)


ggplot(plt.GPM.NEU, aes(x=res, y=score, group=pat, col=score_name)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="score")




plt.PPR.NEU <- rbind(
  plt |>
    dplyr::filter(res == "R1") |> 
    dplyr::mutate(score = PPR_score) |> 
    dplyr::mutate(score_name = "PPR")
  ,
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = NEU_score) |> 
    dplyr::mutate(score_name = "NEU")
  
)

ggplot(plt.PPR.NEU, aes(x=res, y=score, group=pat, col=score_name)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="score")




plt.MTC.NEU <- rbind(
  plt |>
    dplyr::filter(res == "R1") |> 
    dplyr::mutate(score = MTC_score) |> 
    dplyr::mutate(score_name = "MTC")
  ,
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = NEU_score) |> 
    dplyr::mutate(score_name = "NEU")
)


ggplot(plt.GPM.NEU, aes(x=res, y=score, group=pat, col=score_name)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="score")








plt.GPM.NEU <- rbind(
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = GPM_score) |> 
    dplyr::mutate(score_name = "GPM")
  ,
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = NEU_score) |> 
    dplyr::mutate(score_name = "NEU")
  
)

ggplot(plt.GPM.NEU, aes(x=score_name, y=score, group=pat, col=score_name)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="score")




plt.PPR.NEU <- rbind(
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = PPR_score) |> 
    dplyr::mutate(score_name = "PPR")
  ,
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = NEU_score) |> 
    dplyr::mutate(score_name = "NEU")
  
)

ggplot(plt.PPR.NEU, aes(x=score_name, y=score, group=pat, col=score_name)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="score")




plt.MTC.NEU <- rbind(
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = MTC_score) |> 
    dplyr::mutate(score_name = "MTC")
  ,
  plt |>
    dplyr::filter(res == "R2") |> 
    dplyr::mutate(score = NEU_score) |> 
    dplyr::mutate(score_name = "NEU")
)


ggplot(plt.GPM.NEU, aes(x=score_name, y=score, group=pat, col=score_name)) + 
  geom_violin(aes(group=res)) +
  theme_bw() +
  geom_line(aes(group=pat),col="gray80",alpha=0.6) +
  ggbeeswarm::geom_quasirandom() +
  labs(col="score")


