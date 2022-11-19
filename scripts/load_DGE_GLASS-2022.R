#!/usr/bin/env R 


stopifnot(exists("results.out"))
stopifnot(file.exists("cache/analysis_DGE_GLASS-2020.Rds"))



# load ----
analysis.dge.glass.2022 <- readRDS("cache/analysis.dge.glass.2022.Rds")


# merge ----

results.out <- results.out |> 
  dplyr::left_join(glass.gene.res.res.2022, by = c('ensembl_id' = 'ensembl_id'), suffix=c('',''))



# cleanup ----


rm(analysis.dge.glass.2022, metadata)


# smome test plots ----

tmp <- results.out |> 
  dplyr::filter(!is.na(`stat.glass.res`) & !is.na(`stat.glass-2022.res`))
plot(tmp$`stat.glass.res` ,
     tmp$`stat.glass-2022.res`)


cor(x=tmp$`stat.glass.res` ,
    y=tmp$`stat.glass-2022.res`)

cor(x=c(1,2,3,4,5),
    y=c(2,4,5,NA,10),
    rm.na=TRUE)


