#!/usr/bin/env

# load data ----


if(!exists('gsam.rna.metadata')) {
  source('scripts/load_G-SAM_metadata.R')
}

if(!exists('glass.gbm.rnaseq.metadata.all.samples')) {
  source('scripts/load_GLASS_data.R')
}



# input labels ----


## G-SAM ----

### GlioVis MJ -> ssGSEA 2022 ----

tmp <- gsam.rna.metadata  |> 
  dplyr::filter(!is.na(`gliovis.majority_call`) & !is.na(`ssGSEA.2022.subtype`)) |> 
  dplyr::select(`gliovis.majority_call` , `ssGSEA.2022.subtype`) |> 
  dplyr::mutate(textual = paste0(`gliovis.majority_call` ," => ", `ssGSEA.2022.subtype`))

stopifnot(nrow(tmp) == 287)


table(tmp$gliovis.majority_call)
table(tmp$ssGSEA.2022.subtype)


table(tmp$textual)


sum(tmp$gliovis.majority_call == "Classical")
sum(tmp$ssGSEA.2022.subtype == "Classical")

sum(tmp$gliovis.majority_call == "Mesenchymal")
sum(tmp$ssGSEA.2022.subtype == "Mesenchymal")

sum(tmp$gliovis.majority_call == "Proneural")
sum(tmp$ssGSEA.2022.subtype == "Proneural")



rm(tmp)


## GLASS ----
### Synapse 2022 x ssGSEA 2022 ----

tmp <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |> 
  dplyr::filter(!is.na(ssGSEA.Synapse.subtype.2022) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::select(ssGSEA.Synapse.subtype.2022, ssGSEA.2022.subtype) |> 
  dplyr::mutate(textual = paste0(`ssGSEA.Synapse.subtype.2022` ," => ", `ssGSEA.2022.subtype`))

stopifnot(nrow(tmp) == 216)

table(tmp$textual)


### Synapse 2021 x ssGSEA 2022  ----

tmp <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |> 
  dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::select(GBM.transcriptional.subtype.Synapse.2021, ssGSEA.2022.subtype) |> 
  dplyr::mutate(textual = paste0(`GBM.transcriptional.subtype.Synapse.2021` ," => ", `ssGSEA.2022.subtype`))

stopifnot(nrow(tmp) == 52)

table(tmp$textual)

sum(tmp$GBM.transcriptional.subtype.Synapse.2021 == "Classical")
sum(tmp$ssGSEA.2022.subtype == "Classical")

sum(tmp$GBM.transcriptional.subtype.Synapse.2021 == "Mesenchymal")
sum(tmp$ssGSEA.2022.subtype == "Mesenchymal")

sum(tmp$GBM.transcriptional.subtype.Synapse.2021 == "Proneural")
sum(tmp$ssGSEA.2022.subtype == "Proneural")




### Synapse 2021 x Synapse 2022 ----

# tmp <- glass.gbm.rnaseq.metadata.all.samples |> 
#   dplyr::filter(tumour.percentage.2022 >= 15) |> 
#   dplyr::filter(!is.na(GBM.transcriptional.subtype.Synapse.2021) & !is.na(ssGSEA.Synapse.subtype.2022)) |> 
#   dplyr::select(GBM.transcriptional.subtype.Synapse.2021, ssGSEA.Synapse.subtype.2022) |> 
#   dplyr::mutate(textual = paste0(`GBM.transcriptional.subtype.Synapse.2021` ," => ", `ssGSEA.Synapse.subtype.2022`))
# 
# stopifnot(nrow(tmp) == 52)
# 
# table(tmp$textual)





# output labels ----
## G-SAM ----
# stats over verandering in finale GITS subtype labels

tmp <- gsam.rna.metadata  |> 
  dplyr::filter(!is.na(`NMF.123456.PCA.SVM.class_2021`) & !is.na(`GITS.150.svm.2022.subtype`)) |> 
  dplyr::select(`NMF.123456.PCA.SVM.class_2021` , `GITS.150.svm.2022.subtype`) |> 
  dplyr::mutate(textual = paste0(`NMF.123456.PCA.SVM.class_2021` ," => ", `GITS.150.svm.2022.subtype`))

stopifnot(nrow(tmp) == 287)

table(tmp$textual)


rm(tmp)


## GITS x ssGSEA ----


tmp.1 <- gsam.rna.metadata |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::mutate(textual = paste0(`ssGSEA.2022.subtype` ," => ", `GITS.150.svm.2022.subtype`))
table(tmp.1$textual)

(2+1+4+2+3) / (125+2+1+4+90+2+3+56)


tmp.2 <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(!is.na(GITS.150.svm.2022.subtype) & !is.na(ssGSEA.2022.subtype)) |> 
  dplyr::mutate(textual = paste0(`ssGSEA.2022.subtype` ," => ", `GITS.150.svm.2022.subtype`))
table(tmp.2$textual)


# stats ----


tmp <- gsam.rna.metadata |> 
  dplyr::filter(!is.na(ssGSEA.2022.subtype)) |> 
  
  dplyr::filter(!is.na(gliovis.majority_call)) |> 
  dplyr::filter(!is.na(NMF.123456.PCA.SVM.class_2021)) |> 
  dplyr::filter(!is.na(gliovis.gsea_call)) |> 
  dplyr::filter(!is.na(gliovis.knn_call)) |>  
  dplyr::filter(!is.na(gliovis.svm_call))
  
  


sum(tmp$NMF.123456.PCA.SVM.class_2021 != tmp$gliovis.majority_call) # 24

sum(tmp$gliovis.gsea_call != tmp$gliovis.majority_call)
sum(tmp$gliovis.knn_call != tmp$gliovis.majority_call)
sum(tmp$gliovis.svm_call != tmp$gliovis.majority_call)




sum(tmp$gliovis.majority_call != tmp$gliovis.gsea_call)
sum(tmp$gliovis.knn_call != tmp$gliovis.gsea_call)
sum(tmp$gliovis.svm_call != tmp$gliovis.gsea_call)


