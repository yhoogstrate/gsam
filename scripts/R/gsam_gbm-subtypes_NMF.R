#!/usr/bin/env R

# load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/gsam_rna-seq_expression.R')

source('data/wang/msig.library.12.R') # license is incomplete regarding sharing, can't include it in source tree

source('scripts/R/subtype_genes.R')
source('scripts/R/wang_glioma_intrinsic_genes.R')


# NMF ----

seeds <- c(123456, 12345, 1337, 909)

gsam_nmf_150 <- {{}}
gsam_nmf_15k <- {{}}
gsam_nmf_7k <- {{}}

## Normalize for NMF ----

### normalize the VST values ----
tmp <- gsam.rnaseq.expression.vst

if(min(tmp) < 0) { # not the case
  tmp <- tmp - min(tmp) + .Machine$double.eps
}


### make 15k set ----
### make 7k set ----
### make 150 set ----


gsam.rnaseq.expression.vst.150 <- tmp %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::mutate(gid = gsub("\\..+$","",gid)) %>%
  dplyr::filter(gid %in%
                  (wang.glioma.intrinsic.genes %>%
                    dplyr::filter(Subtyping_Signature_Gene. != "") %>%
                    dplyr::pull(ENSG.short))
                  ) %>%
  dplyr::mutate(gid = NULL)

metadata <- gsam.rna.metadata %>%
  dplyr::filter(sid %in% colnames(gsam.rnaseq.expression.vst.150))

stopifnot(ncol(gsam.rnaseq.expression.vst.150) == nrow(metadata))





### exclude IDH muts [after normalisation, so they can be mapped back into the coordinates later] ----

## Perform NMF by Wang-code [seeds: ] ----


for(seed in seeds) {
  gsam_nmf_150[[as.character(seed)]] <- NMF(as.matrix(gsam.rnaseq.expression.vst.150), 3, seed = seed)
}


gsam_nmf_150$`909`$membership
gsam_nmf_150$`12345`$membership
gsam_nmf_150$`123456`$membership
gsam_nmf_150$`1337`$membership


plot(
  gsam_nmf_150$`909`$W[,1],
  gsam_nmf_150$`909`$W[,2],
  col = as.factor(gsam_nmf_150$`909`$membership)
)


plot(
  gsam_nmf_150$`1337`$W[,1],
  gsam_nmf_150$`1337`$W[,2],
)



plot(
  gsam_nmf_150$`909`$W[,1],
  gsam_nmf_150$`909`$W[,2],
)



plot(
  gsam_nmf_150$`909`$W[,1],
  gsam_nmf_150$`909`$W[,2],
)




## Perform PCA on NMF coordinates ? ----

## Export ----



