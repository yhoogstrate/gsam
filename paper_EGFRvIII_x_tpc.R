#!/usr/bin/env R

# load libs ----


# load data ----


source('scripts/R/gsam_metadata.R')
source('scripts/R/gsam_rnaseq_egfrviii_expression.R')



# plt ----


plt <- gsam.patient.metadata %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::filter(batch != "old" & resection == 'r1') %>%
      dplyr::filter(blacklist.pca == F) %>%
      dplyr::filter(pat.with.IDH == F) %>%
      dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
      dplyr::select(c('pid','tumour.percentage.dna')) %>%
      dplyr::rename(tumour.percentage.dna.R1 = tumour.percentage.dna)
    , by = c('studyID'='pid') ) %>%
  dplyr::left_join(
    gsam.rna.metadata %>% dplyr::filter(batch != "old" & resection == 'r2') %>%
      dplyr::filter(blacklist.pca == F) %>%
      dplyr::filter(pat.with.IDH == F) %>%
      dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
      dplyr::select(c('pid','tumour.percentage.dna')) %>%
      dplyr::rename(tumour.percentage.dna.R2 = tumour.percentage.dna) ,
    by = c('studyID'='pid') ) %>%
  dplyr::filter(!is.na(tumour.percentage.dna.R1) & !is.na(tumour.percentage.dna.R2) ) %>%
  dplyr::mutate(tpc.lfc = log(tumour.percentage.dna.R2 / tumour.percentage.dna.R1))  %>%
  dplyr::mutate(extentOfResectionFirstSurgery = ifelse(grepl("Complete resection",extentOfResectionFirstSurgery) , "Complete resection" , extentOfResectionFirstSurgery)  ) %>%
  dplyr::mutate(extentOfResectionSecondSurgery = ifelse(grepl("Complete resection",extentOfResectionSecondSurgery) , "Complete resection" , extentOfResectionSecondSurgery)  ) %>%
  dplyr::mutate(extentOfResectionState = paste0('R1:',extentOfResectionFirstSurgery, '->R2:' , extentOfResectionSecondSurgery) ) %>%
  dplyr::mutate(change.extend.of.resection = extentOfResectionFirstSurgery != extentOfResectionSecondSurgery )




plt <-  gsam.rna.metadata %>% 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) %>%
  dplyr::filter(blacklist.pca == F) %>%
  dplyr::filter(pat.with.IDH == F) %>%
  dplyr::select(c('sid','pid','resection','wt.reads.v3','vIII.reads.v3','vIII.percentage','tumour.percentage.dna')) %>%
  dplyr::filter(!is.na(vIII.percentage)) %>% dplyr::arrange(pid, resection)

tmp <- plt %>% filter(duplicated(pid)) %>% pull(pid) %>% as.character() # ids with pairs

plt <- plt %>%
  dplyr::filter(pid %in% tmp)

rm(tmp)
  

stopifnot  ( (plt %>% dplyr::filter(resection == 'r1') %>% pull(pid) %>% duplicated() %>% sum()) == 0  )
stopifnot  ( (plt %>% dplyr::filter(resection == 'r2') %>% pull(pid) %>% duplicated() %>% sum()) == 0  )


ggplot(plt, aes(x = vIII.percentage , y = tumour.percentage.dna, group=pid, col=resection)) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8 )  +
  geom_line(col="gray") +
  geom_point(data = subset(plt, resection == 'r1'))



# aggregate ----


plt.r1 <- plt %>% dplyr::filter(resection == 'r1')
plt.r2 <- plt %>% dplyr::filter(resection == 'r2')

stopifnot(plt.r1$pid  == plt.r2$pid)


plt.ag <- dplyr::full_join(
  plt.r1 %>%  `colnames<-`(paste0(colnames(.), '.R1')) %>% dplyr::rename(pid = pid.R1) ,
  plt.r2 %>%  `colnames<-`(paste0(colnames(.), '.R2')) %>% dplyr::rename(pid = pid.R2),
  by = c('pid'='pid') ) %>%
  dplyr::mutate(lfc.egfrviii.pct =  log(vIII.percentage.R2 / vIII.percentage.R1 )) %>%
  dplyr::mutate(lfc.tpc = log(tumour.percentage.dna.R2 / tumour.percentage.dna.R1) ) %>%
  dplyr::filter(vIII.percentage.R2 > 0.01 & vIII.percentage.R1 > 0.01)


ggplot(plt.ag, aes(x = lfc.egfrviii.pct , y = lfc.tpc)) + 
  geom_point()




