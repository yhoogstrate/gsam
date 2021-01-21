#!/usr/bin/env R

# read counts from GLASS samples ----

# possibly some match w/ Wang dataset?

# Kallisto counts somehow may harbour floating point digits..
glass.gbm.rnaseq.expression <- 'data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt' %>%
  read.delim(stringsAsFactors = F) %>%
  dplyr::mutate(length = NULL) %>% # effective length from Kalliosto used to quantify
  tibble::column_to_rownames('target_id')


# metadata ----


glass.gbm.rnaseq.metadata <- data.frame(sid = colnames(glass.gbm.rnaseq.expression),
                                        stringsAsFactors = F) %>%
  dplyr::mutate(pid = gsub("^(............).+$","\\1",sid)) %>%
  dplyr::arrange(pid) %>%
  dplyr::mutate(resection = as.factor(gsub("^.............(..).+$","\\1",sid))) %>% # TP is primary tumour? https://github.com/fpbarthel/GLASS
  dplyr::mutate(dataset =  as.factor(gsub("^(....).+$","\\1",sid)) )


