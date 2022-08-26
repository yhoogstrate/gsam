#!/usr/bin/env R 

# load libs ----


#library(EPIC)
#library(patchwork)



# load data ----


source("scripts/load_G-SAM_metadata.R")
source("scripts/R/palette.R")


# fig 2a ----


tmp.paired <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |>  # replicates
  dplyr::filter(!is.na(`tumour.percentage.dna`)) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |>  # only pairs - with low purity included for statistical analysis regarding purity'
  dplyr::mutate(names_from = as.character(dplyr::recode(resection,`r1`="primary",`r2`="recurrence"))) |> 
  tidyr::pivot_wider(names_from = names_from,
                     id_cols = pid,
                     values_from = c(`tumour.percentage.dna`,
                                     `EPIC: CD4 T-cells`,
                                     `EPIC: Macrophages`,
                                     `NMF:150:3`)
                     ) |> 
  dplyr::mutate(has.low.purity.sample = tumour.percentage.dna_recurrence < 15 | tumour.percentage.dna_primary < 15) |> 
  dplyr::mutate(`tumour.percentage.dna.log.odds` = log((tumour.percentage.dna_recurrence + 1) / (tumour.percentage.dna_primary + 1) )) |> 
  dplyr::mutate(`tumour.percentage.dna.delta` = tumour.percentage.dna_recurrence - tumour.percentage.dna_primary) |> 
  
  dplyr::mutate(`EPIC: CD4 T-cells log-odds` = log((`EPIC: CD4 T-cells_recurrence` + 1 ) /  ( `EPIC: CD4 T-cells_primary` + 1 ))) |> 
  dplyr::mutate(`EPIC: CD4 T-cells delta` = `EPIC: CD4 T-cells_recurrence` - `EPIC: CD4 T-cells_primary` ) |> 
  dplyr::mutate(`EPIC: Macrophages log-odds` = log((`EPIC: Macrophages_recurrence` + 1 ) /  ( `EPIC: Macrophages_primary` + 1 ))) |> 
  dplyr::mutate(`EPIC: Macrophages delta` = `EPIC: Macrophages_recurrence` - `EPIC: Macrophages_primary`) 



cor.logodds.purtiy.logodds.epic.tcell <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`EPIC: CD4 T-cells log-odds`)
cor.logodds.purtiy.logodds.epic.macrophages <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`EPIC: Macrophages log-odds`)

cor.delta.purtiy.delta.epic.tcell <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`EPIC: CD4 T-cells delta`)
cor.delta.purtiy.delta.epic.macrophages <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`EPIC: Macrophages delta`)


# seems odd to go from single -> paired -> single, but this way we can easily preserve the log odds stats to single
tmp.single <- tmp.paired |>
  tidyr::pivot_longer(
    cols = c(
      "tumour.percentage.dna_primary", "tumour.percentage.dna_recurrence" ,
      "EPIC: CD4 T-cells_primary", "EPIC: CD4 T-cells_recurrence",
      "EPIC: Macrophages_primary", "EPIC: Macrophages_recurrence"  ,
      "NMF:150:3_primary", "NMF:150:3_recurrence"
    ),  names_to = c(".value", "resection"), names_pattern = "(.+)_(.+)"
  )





## panel top ----


plt <- tmp.single |> 
  tidyr::pivot_longer(cols=c('tumour.percentage.dna',
                             'EPIC: CD4 T-cells',
                             'EPIC: Macrophages',
                             'NMF:150:3'
                             ), values_to = 'y', names_to = 'type') |>  # separate the panels for CD4, Macrrophage and NMF W3 (MES)
  dplyr::mutate(panel = ifelse(has.low.purity.sample,"Patients with purity < 15% sample(s)", "Patients with purity >= 15% for both samples")) |> 
  dplyr::filter((has.low.purity.sample & type == "NMF:150:3") == F) |> 
  dplyr::group_by(pid, type) |> 
  dplyr::mutate(`sign` = ifelse((y == max(y) & resection == "recurrence") | (y == min(y) & resection == "primary") , "increase", "decrease")) |> 
  dplyr::ungroup()
  


plt |>
  dplyr::filter(type == 'tumour.percentage.dna') |> 
  dplyr::filter(pid %in% c('FAG','HAE')) |> 
  dplyr::group_by(pid, type) |> 
  dplyr::mutate(`sign` = ifelse((y == max(y) & resection == "recurrence") | (y == min(y) & resection == "primary") , "increase", "decrease")) |> 
  dplyr::ungroup()



#p1 <- 
ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=y, col=sign)) +
  geom_point(data = subset(plt, resection == "primary"), pch=19, cex=1.2, alpha=0.8) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8, lwd=1.05 )  +
  scale_color_manual(values = c('increase'='#bb5f6c', 'decrease'='#79b1b1')) +
  facet_grid(rows=vars(type),cols=vars(panel), scales="free", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  labs(y=NULL,x="G-SAM patient",col="Longitudinal direction")








plt <- out.sample %>%
  dplyr::filter(pid %in% out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds,
                                                   tumour.percentage.dna.R1,
                                                   tumour.percentage.dna.R2,
                                                   EPIC.Macrophages.log.odds,EPIC.CD4_Tcells.log.odds
                                                   ,`NMF:123456.1.log.odds`
                                                   #,`NMF:123456.2.log.odds`
                                                   #,`NMF:123456.3.log.odds`
  )), by=c('pid'='pid')) %>%
  dplyr::mutate(tumour.percentage.dna.delta = tumour.percentage.dna.R2 - tumour.percentage.dna.R1 ) %>%
  dplyr::mutate(panel = ifelse(
    pid %in% (out.sample %>% dplyr::filter(!is.na(tumour.percentage.dna) & tumour.percentage.dna < 15) %>% dplyr::pull(pid)),
    "< 15%", ">= 15%"
  )) %>% 
  dplyr::mutate(panel = factor(panel, levels=c("< 15%", ">= 15%"))) # fix order



plt <- rbind(plt %>% dplyr::mutate(y = tumour.percentage.dna, type="Tumour purity",
                                   sign = ifelse(tumour.percentage.dna.log.odds > 0, "increase", "decrease")),
             plt %>% dplyr::mutate(y = EPIC.Macrophages, type="EPIC Macrophages score",
                                   sign = ifelse(EPIC.Macrophages.log.odds > 0, "increase", "decrease")),
             plt %>% dplyr::mutate(y = EPIC.CD4_Tcells, type="EPIC CD4 T-cells score",
                                   sign = ifelse(EPIC.CD4_Tcells.log.odds > 0, "increase", "decrease"))
             ,plt %>% dplyr::mutate(y = `NMF:123456.1`, type="NMF W1 (MES)",
                                    sign = ifelse(`NMF:123456.1.log.odds` > 0, "increase", "decrease")) 
             # ,plt %>% dplyr::mutate(y = `NMF:123456.2`, type="NMF W2 (MES)",
             #                        sign = ifelse(`NMF:123456.2.log.odds` > 0, "increase", "decrease"))
             # ,plt %>% dplyr::mutate(y = `NMF:123456.3`, type="NMF W3 (MES)",
             #                        sign = ifelse(`NMF:123456.3.log.odds` > 0, "increase", "decrease"))
)


p1 <- ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=y, col=sign)) +
  #p1 <- ggplot(plt, aes(x = reorder(pid, `NMF:123456.1.log.odds`), y=y, col=sign)) +
  geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2, alpha=0.8) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8, lwd=1.05 )  +
  scale_color_manual(values = c('increase'='#bb5f6c', 'decrease'='#79b1b1')) +
  facet_grid(rows=vars(type),cols=vars(panel), scales="free", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  labs(y=NULL,x="G-SAM patient",col="Longitudinal direction")

p1




plt <- out.sample %>%
  dplyr::filter(pid %in% out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds)), by=c('pid'='pid')) %>%
  dplyr::mutate(panel = ifelse(
    pid %in% (out.sample %>% dplyr::filter(!is.na(tumour.percentage.dna) & tumour.percentage.dna < 15) %>% dplyr::pull(pid)),
    "< 15%", ">= 15%"
  )) %>% 
  dplyr::select(pid, panel, tumour.percentage.dna.log.odds) %>%
  dplyr::distinct() %>%
  dplyr::left_join(out.patient %>% dplyr::select(pid, NMF.123456.PCA.SVM.class.R1, NMF.123456.PCA.SVM.class.R2), by = c('pid' = 'pid')) %>%
  dplyr::rename(`Sub-type R1` = NMF.123456.PCA.SVM.class.R1) %>%
  dplyr::rename(`Sub-type R2` = NMF.123456.PCA.SVM.class.R2) %>%
  dplyr::mutate(panel = factor(panel, levels=c("< 15%", ">= 15%"))) %>%
  dplyr::left_join(gsam.patient.metadata %>%
                     dplyr::select(studyID, 
                                   initialMGMT,HM,
                                   AXIN2,APC,JAK2,
                                   RB1,MSH2,BRCA1,
                                   BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,
                                   NF1,ERBB3,
                                   EGFR,TP53BP1,TP53,PIK3R1,PIK3CA,TSC2,
                                   # PI3K, TP53Signaling, Wnt, Telomere, RTK, RAS, DNADamage, primary7, primary10, recurrent7, recurrent10, cnStatusEGFRs,
                                   
                                   cnStatusRB1s,
                                   cnStatusNF1s,
                                   cnStatusCDK4s,
                                   cnStatusMDM2s,
                                   
                                   SETD2,PDGFRA
                     ),by = c('pid' = 'studyID')
  ) %>%
  reshape2::melt(id = c('pid','panel', 'tumour.percentage.dna.log.odds')) %>% 
  dplyr::mutate(value = gsub('Methylated','Yes',value)) %>% 
  dplyr::mutate(value = gsub('Unmethylated','No',value)) %>%
  dplyr::mutate(value = gsub('Normal','No',value)) %>%
  dplyr::mutate(value = gsub('Loss','Yes',value)) %>%
  dplyr::mutate(value = gsub('^Gain$','Yes',value)) %>%
  dplyr::mutate(value = gsub('^Gained$','Gained/increased',value)) %>%
  dplyr::mutate(value = gsub('^Lost$','Lost/decreased',value)) %>%
  dplyr::mutate(ypanel = case_when(
    grepl("NMF|Sub-type",variable) ~ "A",
    variable %in% c("EGFR",  "cnStatusEGFRs", "NF1") ~ "B",
    T ~ "C"
  )) %>% 
  dplyr::mutate(ypanel = factor(ypanel, levels=c('A','B','C')))



p2 <- ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=variable, fill = value)) +
  facet_grid(cols=vars(panel), rows=vars(ypanel), scales="free", space="free") + 
  geom_tile(colour = "black", size = 0.3) +
  theme_bw() +
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  scale_fill_manual(values = c( subtype_colors ,
                                "Wildtype"="white",
                                "No"="white",
                                
                                "Gained/increased"="#bb5f6c", # rood #bb5f6c
                                "Lost/decreased"="#79b1b1", # lichtlauw #79b1b1
                                "Stable"="#2e415e", # donker blauw #2e415e
                                "Yes" = "#2e415e", # zelfde als stable #2e415e
                                
                                "NA"="grey"  )) + 
  labs(y=NULL,x=NULL)




p1 / p2 +  plot_layout(heights = c(2, 1))





ggsave("output/figures/epic_tumor_percentage_traversal_vertical.pdf", width = 16 * 1.05, height=11.5 * 1.05)







rm(plt)



# Additional changes  ----






out.patient <- out.patient %>%  
  dplyr::filter(pid != 'ECD') %>% # excessive outlier, R1 or R2 RNA odd?
  dplyr::left_join(
    gsam.rna.metadata %>%
      dplyr::filter(resection == 'r1') %>%
      dplyr::select(c('pid','NMF.123456.PCA.SVM.class', 'NMF.123456.PCA.SVM.Classical.p','NMF.123456.PCA.SVM.Proneural.p','NMF.123456.PCA.SVM.Mesenchymal.p','NMF:123456.1')) %>%
      dplyr::filter(!is.na(NMF.123456.PCA.SVM.class) ) %>%
      dplyr::rename(NMF.123456.PCA.SVM.class.R1 = NMF.123456.PCA.SVM.class) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Classical.p.R1 = NMF.123456.PCA.SVM.Classical.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Proneural.p.R1 = NMF.123456.PCA.SVM.Proneural.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Mesenchymal.p.R1 = NMF.123456.PCA.SVM.Mesenchymal.p) %>% 
      dplyr::rename('NMF:123456.1.R1' = 'NMF:123456.1')
      
    , by=c('pid'='pid')) %>%
  dplyr::left_join(
    gsam.rna.metadata %>%
      dplyr::filter(resection == 'r2') %>%
      dplyr::select(c('pid','NMF.123456.PCA.SVM.class', 'NMF.123456.PCA.SVM.Classical.p','NMF.123456.PCA.SVM.Proneural.p','NMF.123456.PCA.SVM.Mesenchymal.p','NMF:123456.1')) %>%
      dplyr::filter(!is.na(NMF.123456.PCA.SVM.class) ) %>%
      dplyr::rename(NMF.123456.PCA.SVM.class.R2 = NMF.123456.PCA.SVM.class) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Classical.p.R2 = NMF.123456.PCA.SVM.Classical.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Proneural.p.R2 = NMF.123456.PCA.SVM.Proneural.p) %>%
      dplyr::rename(NMF.123456.PCA.SVM.Mesenchymal.p.R2 = NMF.123456.PCA.SVM.Mesenchymal.p) %>% 
      dplyr::rename('NMF:123456.1.R2' = 'NMF:123456.1')
    , by=c('pid'='pid')) %>%
  dplyr::mutate(NMF.123456.PCA.SVM.Classical.p.log.odds = log( (NMF.123456.PCA.SVM.Classical.p.R2+1) / (NMF.123456.PCA.SVM.Classical.p.R1+1) )) %>%
  dplyr::mutate(NMF.123456.PCA.SVM.Mesenchymal.p.log.odds = log( (NMF.123456.PCA.SVM.Mesenchymal.p.R2+1) / (NMF.123456.PCA.SVM.Mesenchymal.p.R1+1) )) %>%
  dplyr::mutate(NMF.123456.PCA.SVM.Proneural.p.log.odds = log( (NMF.123456.PCA.SVM.Proneural.p.R2+1) / (NMF.123456.PCA.SVM.Proneural.p.R1+1 ))) %>% 
  dplyr::mutate(`NMF:123456.1.log.odds` = log( (`NMF:123456.1.R1` + 1) / (`NMF:123456.1.R2` +1)))
  
  
  #dplyr::mutate(NMF.123456.PCA.SVM.Classical.p.log.odds2 = NMF.123456.PCA.SVM.Classical.p.R2 - NMF.123456.PCA.SVM.Classical.p.R1 ) %>%
  #dplyr::mutate(NMF.123456.PCA.SVM.Mesenchymal.p.log.odds2 = NMF.123456.PCA.SVM.Mesenchymal.p.R2 - NMF.123456.PCA.SVM.Mesenchymal.p.R1 ) %>%
  #dplyr::mutate(NMF.123456.PCA.SVM.Proneural.p.log.odds2 = NMF.123456.PCA.SVM.Proneural.p.R2 - NMF.123456.PCA.SVM.Proneural.p.R1 )







# plots ----

## Figure 2fg ----

plt <- out.sample %>%
  dplyr::filter(pid %in%  out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds,
                                                   EPIC.Macrophages.log.odds,
                                                   EPIC.CD4_Tcells.log.odds
                                                   )), by=c('pid'='pid')) %>%
  dplyr::mutate(tumour.percentage.status = ifelse(tumour.percentage.dna.log.odds > 0, "increase", "decrease"))


plt <- plt %>%
  dplyr::select(pid, resection, tumour.percentage.dna, tumour.percentage.status, EPIC.Macrophages, EPIC.CD4_Tcells) %>% 
  reshape2::melt(id = c('pid','resection', 'tumour.percentage.dna', 'tumour.percentage.status'))


ggplot(plt, aes(x=tumour.percentage.dna, y=value, group=pid, col= tumour.percentage.status)) + 
  facet_grid(cols = vars(variable), scales = "free",space="free") +
  geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
  geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
  labs(x = "Tumor cell percentage (DNA-seq estimate)",y='EPIC deconvolution score') +
  scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1" )) +
  xlim(0, 100) +
  ylim(0, 0.5) +
  youri_gg_theme


# p1 <- ggplot(plt, aes(x=tumour.percentage.dna, y=EPIC.Macrophages, group=pid, col= tumour.percentage.status)) + 
#   geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
#   geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
#   labs(x = "Tumor cell percentage (DNA-seq estimate)",y="EPIC Macrophages score") +
#   scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1" )) +
#   xlim(0, 100) +
#   ylim(0, 0.5) +
#   youri_gg_theme
# 
# p2 <- ggplot(plt, aes(x=tumour.percentage.dna, y=EPIC.CD4_Tcells, group=pid, col= tumour.percentage.status)) + 
#   geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
#   geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
#   labs(x = "Tumor cell percentage (DNA-seq estimate)",y="EPIC CD4 T-cells score") +
#   scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1" )) +
#   xlim(0, 100) +
#   ylim(0, 0.5) +
#   youri_gg_theme
# 
# 
# p1 + p2
# 


# ggsave("output/figures/epic_tumor_percentage_traversal.pdf", width = 12, height=4.8)
# ggsave("output/figures/epic_tumor_percentage_traversal_4.2.pdf", width = 12, height=4.2)



## Figure 2a ----



plt <- out.sample %>%
  dplyr::filter(pid %in% out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds,
                                                   tumour.percentage.dna.R1,
                                                   tumour.percentage.dna.R2,
                                                   EPIC.Macrophages.log.odds,EPIC.CD4_Tcells.log.odds
                                                   ,`NMF:123456.1.log.odds`
                                                   #,`NMF:123456.2.log.odds`
                                                   #,`NMF:123456.3.log.odds`
                                                   )), by=c('pid'='pid')) %>%
  dplyr::mutate(tumour.percentage.dna.delta = tumour.percentage.dna.R2 - tumour.percentage.dna.R1 ) %>%
  dplyr::mutate(panel = ifelse(
    pid %in% (out.sample %>% dplyr::filter(!is.na(tumour.percentage.dna) & tumour.percentage.dna < 15) %>% dplyr::pull(pid)),
    "< 15%", ">= 15%"
  )) %>% 
  dplyr::mutate(panel = factor(panel, levels=c("< 15%", ">= 15%"))) # fix order



plt <- rbind(plt %>% dplyr::mutate(y = tumour.percentage.dna, type="Tumour purity",
                                   sign = ifelse(tumour.percentage.dna.log.odds > 0, "increase", "decrease")),
             plt %>% dplyr::mutate(y = EPIC.Macrophages, type="EPIC Macrophages score",
                                   sign = ifelse(EPIC.Macrophages.log.odds > 0, "increase", "decrease")),
             plt %>% dplyr::mutate(y = EPIC.CD4_Tcells, type="EPIC CD4 T-cells score",
                                   sign = ifelse(EPIC.CD4_Tcells.log.odds > 0, "increase", "decrease"))
             ,plt %>% dplyr::mutate(y = `NMF:123456.1`, type="NMF W1 (MES)",
                                    sign = ifelse(`NMF:123456.1.log.odds` > 0, "increase", "decrease")) 
             # ,plt %>% dplyr::mutate(y = `NMF:123456.2`, type="NMF W2 (MES)",
             #                        sign = ifelse(`NMF:123456.2.log.odds` > 0, "increase", "decrease"))
             # ,plt %>% dplyr::mutate(y = `NMF:123456.3`, type="NMF W3 (MES)",
             #                        sign = ifelse(`NMF:123456.3.log.odds` > 0, "increase", "decrease"))
             )


p1 <- ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=y, col=sign)) +
#p1 <- ggplot(plt, aes(x = reorder(pid, `NMF:123456.1.log.odds`), y=y, col=sign)) +
  geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2, alpha=0.8) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.8, lwd=1.05 )  +
  scale_color_manual(values = c('increase'='#bb5f6c', 'decrease'='#79b1b1')) +
  facet_grid(rows=vars(type),cols=vars(panel), scales="free", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  labs(y=NULL,x="G-SAM patient",col="Longitudinal direction")

p1




plt <- out.sample %>%
  dplyr::filter(pid %in% out.patient$pid) %>%
  dplyr::arrange(pid, resection, sid) %>%
  dplyr::left_join(out.patient %>% dplyr::select(c(pid, tumour.percentage.dna.log.odds)), by=c('pid'='pid')) %>%
  dplyr::mutate(panel = ifelse(
    pid %in% (out.sample %>% dplyr::filter(!is.na(tumour.percentage.dna) & tumour.percentage.dna < 15) %>% dplyr::pull(pid)),
    "< 15%", ">= 15%"
  )) %>% 
  dplyr::select(pid, panel, tumour.percentage.dna.log.odds) %>%
  dplyr::distinct() %>%
  dplyr::left_join(out.patient %>% dplyr::select(pid, NMF.123456.PCA.SVM.class.R1, NMF.123456.PCA.SVM.class.R2), by = c('pid' = 'pid')) %>%
  dplyr::rename(`Sub-type R1` = NMF.123456.PCA.SVM.class.R1) %>%
  dplyr::rename(`Sub-type R2` = NMF.123456.PCA.SVM.class.R2) %>%
  dplyr::mutate(panel = factor(panel, levels=c("< 15%", ">= 15%"))) %>%
  dplyr::left_join(gsam.patient.metadata %>%
      dplyr::select(studyID, 
                    initialMGMT,HM,
                    AXIN2,APC,JAK2,
                    RB1,MSH2,BRCA1,
                    BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,
                    NF1,ERBB3,
                    EGFR,TP53BP1,TP53,PIK3R1,PIK3CA,TSC2,
                    # PI3K, TP53Signaling, Wnt, Telomere, RTK, RAS, DNADamage, primary7, primary10, recurrent7, recurrent10, cnStatusEGFRs,
                    
                    cnStatusRB1s,
                    cnStatusNF1s,
                    cnStatusCDK4s,
                    cnStatusMDM2s,
                    
                    SETD2,PDGFRA
                    ),by = c('pid' = 'studyID')
  ) %>%
  reshape2::melt(id = c('pid','panel', 'tumour.percentage.dna.log.odds')) %>% 
  dplyr::mutate(value = gsub('Methylated','Yes',value)) %>% 
  dplyr::mutate(value = gsub('Unmethylated','No',value)) %>%
  dplyr::mutate(value = gsub('Normal','No',value)) %>%
  dplyr::mutate(value = gsub('Loss','Yes',value)) %>%
  dplyr::mutate(value = gsub('^Gain$','Yes',value)) %>%
  dplyr::mutate(value = gsub('^Gained$','Gained/increased',value)) %>%
  dplyr::mutate(value = gsub('^Lost$','Lost/decreased',value)) %>%
  dplyr::mutate(ypanel = case_when(
    grepl("NMF|Sub-type",variable) ~ "A",
    variable %in% c("EGFR",  "cnStatusEGFRs", "NF1") ~ "B",
    T ~ "C"
    )) %>% 
  dplyr::mutate(ypanel = factor(ypanel, levels=c('A','B','C')))



p2 <- ggplot(plt, aes(x = reorder(pid, tumour.percentage.dna.log.odds), y=variable, fill = value)) +
  facet_grid(cols=vars(panel), rows=vars(ypanel), scales="free", space="free") + 
  geom_tile(colour = "black", size = 0.3) +
  theme_bw() +
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  scale_fill_manual(values = c( subtype_colors ,
"Wildtype"="white",
"No"="white",

"Gained/increased"="#bb5f6c", # rood #bb5f6c
"Lost/decreased"="#79b1b1", # lichtlauw #79b1b1
"Stable"="#2e415e", # donker blauw #2e415e
"Yes" = "#2e415e", # zelfde als stable #2e415e
  
"NA"="grey"  )) + 
  labs(y=NULL,x=NULL)




p1 / p2 +  plot_layout(heights = c(2, 1))





ggsave("output/figures/epic_tumor_percentage_traversal_vertical.pdf", width = 16 * 1.05, height=11.5 * 1.05)



