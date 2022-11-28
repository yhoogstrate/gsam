#!/usr/bin/env R 

# load libs ----

# load data ----


source("scripts/load_G-SAM_metadata.R")
source("scripts/R/palette.R")


# F] Figure M2A ----


tmp.paired <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |>  # replicates
  dplyr::filter(!is.na(`tumour.percentage.dna`)) |> 
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |>  # only pairs - with low purity included for statistical analysis regarding purity'
  dplyr::ungroup() |> 
  dplyr::mutate(names_from = as.character(dplyr::recode(resection,`r1`="primary",`r2`="recurrence"))) |> 
  tidyr::pivot_wider(names_from = names_from,
                     id_cols = pid,
                     values_from = c(`tumour.percentage.dna`,
                                     
                                     
                                     
                                     `EPIC: CD4 T-cells`,
                                     `EPIC: Macrophages`,
                                     
                                     `NMF:150:1`,
                                     `NMF:150:2`,
                                     `NMF:150:3`,
                                     
                                     `GITS.150.svm.2022.subtype`,
                                     extent,
                                     
                                     rna.signature.C0.fuzzy.2022,
                                     rna.signature.C1.collagen.2022,
                                     rna.signature.C2.endothelial.2022,
                                     rna.signature.C3.oligodendrocyte.2022,
                                     rna.signature.C4.neuron.2022
                                     
                                     )
                     ) |> 
  dplyr::mutate(`has.low.purity.sample` = tumour.percentage.dna_recurrence < 15 | tumour.percentage.dna_primary < 15) |> 
  dplyr::mutate(
    panel = ifelse(
      has.low.purity.sample,
      "Patients with purity < 15% sample(s)",
      "Patients with purity \u2265 15% for both samples"
    )
  ) |>
  
  dplyr::mutate(`tumour.percentage.dna.log.odds` = log((tumour.percentage.dna_recurrence + 1) / (tumour.percentage.dna_primary + 1) )) |> 
  dplyr::mutate(`tumour.percentage.dna.delta` = tumour.percentage.dna_recurrence - tumour.percentage.dna_primary) |> 
  
  dplyr::mutate(`NMF:150:1 log-odds` = log((`NMF:150:1_recurrence` + 1 ) /  ( `NMF:150:1_primary` + 1 ))) |> 
  dplyr::mutate(`NMF:150:1 delta` = `NMF:150:1_recurrence` - `NMF:150:1_primary`) |> 
  dplyr::mutate(`NMF:150:2 log-odds` = log((`NMF:150:2_recurrence` + 1 ) /  ( `NMF:150:2_primary` + 1 ))) |> 
  dplyr::mutate(`NMF:150:2 delta` = `NMF:150:2_recurrence` - `NMF:150:2_primary`) |> 
  dplyr::mutate(`NMF:150:3 log-odds` = log((`NMF:150:3_recurrence` + 1 ) /  ( `NMF:150:3_primary` + 1 ))) |> 
  dplyr::mutate(`NMF:150:3 delta` = `NMF:150:3_recurrence` - `NMF:150:3_primary`) |> 
  
  dplyr::mutate(`rna.signature.C0.fuzzy.2022 delta` = `rna.signature.C0.fuzzy.2022_recurrence` - `rna.signature.C0.fuzzy.2022_primary`) |> 
  dplyr::mutate(`rna.signature.C1.collagen.2022 delta` = `rna.signature.C1.collagen.2022_recurrence` - `rna.signature.C1.collagen.2022_primary`) |> 
  dplyr::mutate(`rna.signature.C2.endothelial.2022 delta` = `rna.signature.C2.endothelial.2022_recurrence` - `rna.signature.C2.endothelial.2022_primary`) |> 
  dplyr::mutate(`rna.signature.C3.oligodendrocyte.2022 delta` = `rna.signature.C3.oligodendrocyte.2022_recurrence` - `rna.signature.C3.oligodendrocyte.2022_primary`) |> 
  dplyr::mutate(`rna.signature.C4.neuron.2022 delta` = `rna.signature.C4.neuron.2022_recurrence` - `rna.signature.C4.neuron.2022_primary`) |> 
  
  dplyr::mutate(`EPIC: CD4 T-cells log-odds` = log((`EPIC: CD4 T-cells_recurrence` + 1 ) /  ( `EPIC: CD4 T-cells_primary` + 1 ))) |> 
  dplyr::mutate(`EPIC: CD4 T-cells delta` = `EPIC: CD4 T-cells_recurrence` - `EPIC: CD4 T-cells_primary` ) |> 
  dplyr::mutate(`EPIC: Macrophages log-odds` = log((`EPIC: Macrophages_recurrence` + 1 ) /  ( `EPIC: Macrophages_primary` + 1 ))) |> 
  dplyr::mutate(`EPIC: Macrophages delta` = `EPIC: Macrophages_recurrence` - `EPIC: Macrophages_primary`) |> 
  
  dplyr::left_join(
    gsam.patient.metadata |> 
      dplyr::select(
        studyID,
        
        # svvl
        survivalDays,
        survivalFromSecondSurgeryDays,
        progressionFreeDays,
        status,
        
        # asl
        age,
        gender,
        performanceAtSecondSurgery,
        tumorLocation,
        
        # treat
        treatedWithTMZ,
        treatedWithRT,
        bevacizumab.before.recurrence,
        PTK787.before.recurrence,
        
        # general genetics
        mgmtStability, HM,
        
        # gene muts
        AXIN2, APC, JAK2,
        RB1, MSH2, BRCA1,
        BRCA2, ATM, SETD2, ARID2, KMT2C, KMT2D,
        NF1, ERBB3,
        EGFR, TP53BP1, TP53, PIK3R1, PIK3CA, TSC2,
        SETD2, PDGFRA,
        
        # cn
        cnStatusCDKN2ABs,
        cnStatusEGFRs,
        cnStatusRB1s,
        cnStatusNF1s,
        cnStatusCDK4s,
        cnStatusMDM2s
      ), by=c('pid'='studyID'), suffix=c('','')
  ) |> 
  dplyr::mutate(time.to.progression = survivalDays -  survivalFromSecondSurgeryDays) |> 
  dplyr::rename(event = status) |> 
  dplyr::mutate(rank = order(order(tumour.percentage.dna.log.odds, pid))) |>  # order(order()) equals rank() but allows more than 1 factor 
  dplyr::mutate(event = ifelse(.data$event == "Deceased", 1, 0)) |>
  dplyr::mutate(Deceased = dplyr::recode(event, "1" = "Yes", "0" = "No")) |> 
  dplyr::mutate(`Age above 50` = ifelse(age > 50, "Yes", "No")) |>
  dplyr::mutate(gender = as.character(gender)) |>
  dplyr::rename(Sex = gender) |>
  
  dplyr::mutate(`MGMT meth` = dplyr::recode(mgmtStability,
                                            "Stable methylated" = "Stable",
                                            "Stable unmethylated" = "Wildtype"
  ), mgmtStability = NULL) |>
  dplyr::mutate(`Treatment: Beva` = case_when(
    is.na(bevacizumab.before.recurrence) ~ "NA",
    bevacizumab.before.recurrence == "Trial participant" ~ "Randomized trial",
    bevacizumab.before.recurrence == "Yes" ~ "Yes",
    T ~ "No"
  )) |>
  
  dplyr::mutate(`Treatment: Beva` = factor(`Treatment: Beva`, levels = c("No", "Yes", "Randomized trial"))) |>
  dplyr::mutate(bevacizumab.before.recurrence = NULL) |>
  dplyr::rename(`Treatment: TMZ` = treatedWithTMZ) |>
  dplyr::rename(`Treatment: RT` = treatedWithRT) |>
  dplyr::rename(`Treatment: PTK787` = PTK787.before.recurrence) |>
  dplyr::rename(`Resection/Biopsy R1` = extent_primary) |>
  dplyr::rename(`Resection/Biopsy R2` = extent_recurrence) |>
  dplyr::mutate(`KPS 70 or above` = factor(ifelse(is.na(performanceAtSecondSurgery) | performanceAtSecondSurgery >= 70, "Yes", "No"), levels = c("Yes", "No"))) 



cor.logodds.purtiy.logodds.epic.tcell <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`EPIC: CD4 T-cells log-odds`)
cor.logodds.purtiy.logodds.epic.macrophages <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`EPIC: Macrophages log-odds`)
cor.logodds.purtiy.logodds.epic.nmf1 <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`NMF:150:1 log-odds`)
cor.logodds.purtiy.logodds.epic.nmf2 <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`NMF:150:2 log-odds`)
cor.logodds.purtiy.logodds.epic.nmf3 <- cor.test(tmp.paired$`tumour.percentage.dna.log.odds`, tmp.paired$`NMF:150:3 log-odds`)

cor.delta.purtiy.delta.epic.tcell <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`EPIC: CD4 T-cells delta`)
cor.delta.purtiy.delta.epic.macrophages <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`EPIC: Macrophages delta`)
cor.delta.purtiy.delta.nmf1 <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`NMF:150:1 delta`)
cor.delta.purtiy.delta.nmf2 <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`NMF:150:2 delta`)
cor.delta.purtiy.delta.nmf3 <- cor.test(tmp.paired$`tumour.percentage.dna.delta`, tmp.paired$`NMF:150:3 delta`)



# seems odd to go from single -> paired -> single, but this way we can easily preserve the log odds stats to single
tmp.single <- tmp.paired |>
  tidyr::pivot_longer(
    cols = c(
      "tumour.percentage.dna_primary", "tumour.percentage.dna_recurrence",
      "EPIC: CD4 T-cells_primary", "EPIC: CD4 T-cells_recurrence",
      "EPIC: Macrophages_primary", "EPIC: Macrophages_recurrence",
      "GITS.150.svm.2022.subtype_primary", "GITS.150.svm.2022.subtype_recurrence", 
      "NMF:150:3_primary", "NMF:150:3_recurrence"
    ),  names_to = c(".value", "resection"), names_pattern = "(.+)_(.+)"
  )




## panel top ----


plt <- tmp.single |>
  tidyr::pivot_longer(
    cols = c(
      'tumour.percentage.dna',
      'EPIC: CD4 T-cells',
      'EPIC: Macrophages',
      'NMF:150:3'
    ),
    values_to = 'y',
    names_to = 'type'
  ) |>  # separate the panels for CD4, Macrrophage and NMF W3 (MES)
  dplyr::filter((has.low.purity.sample &
                   type == "NMF:150:3") == F) |>
  dplyr::group_by(pid, type) |>
  dplyr::mutate(`sign` = ifelse((y == max(y) & resection == "recurrence") |
                                (y == min(y) & resection == "primary") ,
                                "increase", "decrease"
  )) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    type = dplyr::recode(
      type,
      `tumour.percentage.dna` = 'Tumor purtiy',
      `NMF:150:3` = 'NMF meta 3 (MES)'
    )
  ) |>
  dplyr::mutate(type = factor(
    type,
    levels = c(
      "Tumor purtiy"          ,
      "EPIC: Macrophages"  ,
      "EPIC: CD4 T-cells"   ,
      "NMF meta 3 (MES)"
    ))) |> 
  dplyr::mutate(label = NA)

# add cor labels
tmp.tcell <- plt[1,] |> 
  dplyr::mutate(type = "EPIC: CD4 T-cells") |> 
  dplyr::mutate(panel = "Patients with purity \u2265 15% for both samples") |> 
  dplyr::mutate(y = plt |> dplyr::filter(type ==  "EPIC: CD4 T-cells") |> dplyr::pull(y) |> max() ) |> 
  dplyr::mutate(label = paste0("R = ",round(cor.delta.purtiy.delta.epic.tcell$estimate,2)))
tmp.mac <- plt[1,] |> 
  dplyr::mutate(type = "EPIC: Macrophages") |> 
  dplyr::mutate(panel = "Patients with purity \u2265 15% for both samples") |> 
  dplyr::mutate(y = plt |> dplyr::filter(type ==  "EPIC: Macrophages") |> dplyr::pull(y) |> max()  ) |>
  dplyr::mutate(label = paste0("R = ",round(cor.delta.purtiy.delta.epic.macrophages$estimate,2)))
tmp.nmf3 <- plt[1,] |> 
  dplyr::mutate(type = "NMF meta 3 (MES)") |> 
  dplyr::mutate(panel = "Patients with purity \u2265 15% for both samples") |> 
  dplyr::mutate(y = plt |> dplyr::filter(type ==  "NMF meta 3 (MES)") |> dplyr::pull(y) |> max()  ) |>
  dplyr::mutate(label = paste0("R = ",round(cor.delta.purtiy.delta.nmf3$estimate,2)))

tmp.labels <- rbind(tmp.mac, tmp.tcell, tmp.nmf3)
rm(tmp.mac, tmp.tcell, tmp.nmf3)


tmp.n.pairs.below.15 <- tmp.paired |>
  dplyr::filter(has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()
tmp.n.pairs.leq.15 <-  tmp.paired |>
  dplyr::filter(!has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()


p1 <- 
ggplot(plt, aes(x = reorder(pid, rank), y=y, col=sign, label=label)) +
  ggplot2::geom_point(data = subset(plt, resection == "primary"), pch=19, cex=1.2 * 2) +
  ggplot2::geom_path(
    lineend = "butt",
    linejoin = "mitre",
    lwd=1.25,
   arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125*0.65 , "inches")))  +
  ggplot2::scale_color_manual(values = c('increase'='#bb5f6c', 'decrease'='#79b1b1'),guide="none") +
  ggplot2::facet_grid(rows=vars(type), cols=vars(panel), scales="free", space="free_x") +
  ggplot2::geom_text(x=6, data=tmp.labels, col="black",fontface = "italic",hjust = 0, vjust=1) +
  ggplot2::labs(
    y = NULL,
    x = NULL, #"G-SAM patient",
    col = NULL, #$"Longitudinal direction",
    caption = paste0("Pairs = ",(tmp.n.pairs.below.15 + tmp.n.pairs.leq.15), " (",tmp.n.pairs.leq.15," >= 15%, ",tmp.n.pairs.below.15," < 15%)")
  ) +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) 



rm(plt, tmp.n.pairs.below.15, tmp.n.pairs.leq.15, tmp.labels)



## panel bottom ----



plt <- tmp.paired |> 
  dplyr::select(
    pid,
    panel,
    Deceased,
    tumour.percentage.dna.log.odds,
    `KPS 70 or above`,
    `Age above 50`,
    Sex,
    `Resection/Biopsy R1`,
    `Resection/Biopsy R2`,
    
    `Treatment: Beva`,
    `Treatment: TMZ`,
    `Treatment: RT`,
    
    GITS.150.svm.2022.subtype_primary, GITS.150.svm.2022.subtype_recurrence,
    
    `MGMT meth`,#mgmtStability,
    HM,

    # gene muts
    AXIN2, APC, JAK2,
    RB1, MSH2, BRCA1,
    BRCA2, ATM, SETD2, ARID2, KMT2C, KMT2D,
    NF1, ERBB3,
    EGFR, TP53BP1, TP53, PIK3R1, PIK3CA, TSC2,
    SETD2, PDGFRA,
    
    # cn
    cnStatusCDKN2ABs,
    cnStatusEGFRs,
    cnStatusRB1s,
    cnStatusNF1s,
    cnStatusCDK4s,
    cnStatusMDM2s,
    
    rank
  ) |> 
  dplyr::rename(`GITS subtype R1` = GITS.150.svm.2022.subtype_primary) |> 
  dplyr::rename(`GITS subtype R2` = GITS.150.svm.2022.subtype_recurrence) |> 
  #dplyr::rename(`MGMT meth` = mgmtStability) |> 
  dplyr::rename(`Hyper mut` = HM) |> 
  reshape2::melt(id = c('pid','panel', 'tumour.percentage.dna.log.odds', 'rank')) |> 

  dplyr::mutate(value = gsub('Unmethylated|Stable unmethylated|Normal|No|Wildtype|Biopsy','No / wildtype / biopsy',value)) |> 
  dplyr::mutate(value = gsub('^Methylated|Stable methylated|Loss|Yes|Stable|Resection$','Yes / stable / resection',value)) |> 
  dplyr::mutate(value = gsub('^Gained|Female$','Gained / increased / female',value)) |>
  dplyr::mutate(value = gsub('^Lost|Male$','Lost / decreased / male',value)) |> 
  
  dplyr::mutate(value = ifelse(is.na(value) | value == "Randomized trial", "NA / Beva: random trial" ,value)) |> 
  
  dplyr::mutate(ypanel = case_when(
    grepl("KPS|Age|eceased|Sex|Biopsy",variable) ~ "A",
    grepl("reatment",variable) ~ "B",
    grepl("ubtype",variable) ~ "C",
    grepl("MGMT|Hyper m",variable) ~ "E",
    T ~ "D"
  )) %>% 
  dplyr::mutate(ypanel = factor(ypanel, levels=c('A','B','C','D','E'))) |> 
  dplyr::filter(
    (
      (panel == "Patients with purity < 15% sample(s)" & (variable) %in% c(
        "Hyper mut", "MGMT meth", "GITS subtype R2", "GITS subtype R1",

        "AXIN2", "APC", "JAK2",
        "RB1", "MSH2", "BRCA1",
        "BRCA2", "ATM", "SETD2", "ARID2", "KMT2C", "KMT2D",
        "NF1", "ERBB3",
        "EGFR", "TP53BP1", "TP53", "PIK3R1", "PIK3CA", "TSC2",
        "SETD2", "PDGFRA",

        "cnStatusCDKN2ABs",
        "cnStatusEGFRs",
        "cnStatusRB1s",
        "cnStatusNF1s",
        "cnStatusCDK4s",
        "cnStatusMDM2s"
      )) == F
    ))

  

tmp.n.pairs.below.15 <- tmp.paired |>
  dplyr::filter(has.low.purity.sample) |>
  dplyr::pull(pid) |>
  unique() |>
  length()
tmp.n.pairs.leq.15 <- tmp.paired |>
  dplyr::filter(!has.low.purity.sample) |>
  dplyr::pull(pid) |>
  unique() |>
  length()


p2 <- ggplot(plt, aes(x = reorder(pid, rank), y=variable, fill = value)) +
  facet_grid(cols=vars(panel), rows=vars(ypanel), scales="free", space="free") + 
  geom_tile(colour = "black", size = 0.3) +
  theme_bw() +
  #coord_equal() +
  scale_fill_manual(values = c( subtype_colors ,
                                "No / wildtype / biopsy"="white",
                                
                                "Gained / increased / female"="#bb5f6c", # rood #bb5f6c
                                "Lost / decreased / male"="#79b1b1", # lichtlauw #79b1b1
                                "Yes / stable / resection"="#2e415e", # donker blauw #2e415e
                                
                                "NA / Beva: random trial"="grey"  )) + 
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  labs(y=NULL,x=NULL,fill = NULL,
       caption = paste0("Pairs = ",(tmp.n.pairs.below.15 + tmp.n.pairs.leq.15), " (",tmp.n.pairs.leq.15," >= 15% ",tmp.n.pairs.below.15," < 15%)")
       ) +
  guides(fill=guide_legend(ncol=10))

table(plt$value)


p1 / p2 + patchwork::plot_layout(heights = c(1.8, 1.6))

ggsave("output/figures/2022_Figure_M2A.pdf", width=8.3 ,height=8.3/1.5 , scale=2)
rm(plt, tmp.n.pairs.below.15, tmp.n.pairs.leq.15, p1, p2)




# F] Figure M2G, M2H - EPIC ~ purity ----


plt <- tmp.single |> 
  dplyr::group_by(pid) |> 
  dplyr::mutate(`tumour.percentage.status` = ifelse((tumour.percentage.dna == max(tumour.percentage.dna) & resection == "recurrence") |
                                  (tumour.percentage.dna == min(tumour.percentage.dna) & resection == "primary") ,
                                "increase", "decrease"
  )) |>
  dplyr::select(pid, resection, tumour.percentage.dna, tumour.percentage.status, `EPIC: Macrophages`, `EPIC: CD4 T-cells`) %>% 
  reshape2::melt(id = c('pid','resection', 'tumour.percentage.dna', 'tumour.percentage.status')) |> 
  dplyr::mutate(label = NA) |> 
  dplyr::mutate(variable = as.character(variable )) |> 
  dplyr::mutate(variable = dplyr::recode(variable, 'EPIC: Macrophages'='Macrophages','EPIC: CD4 T-cells'='CD4 T-cells'))


cor.delta.purity.delta.epic.tcells <- cor(
  tmp.single |>
    dplyr::filter(resection == "primary") |>  # either select primary or recurrence, but just one per patient
    dplyr::pull(tumour.percentage.dna.delta),
  tmp.single |>
    dplyr::filter(resection == "primary") |>  # either select primary or recurrence, but just one per patient
    dplyr::pull(`EPIC: CD4 T-cells delta`)
)
cor.delta.purity.delta.epic.macrophages <- cor(
  tmp.single |>
    dplyr::filter(resection == "primary") |>  # either select primary or recurrence, but just one per patient
    dplyr::pull(tumour.percentage.dna.delta),
  tmp.single |>
    dplyr::filter(resection == "primary") |>  # either select primary or recurrence, but just one per patient
    dplyr::pull(`EPIC: Macrophages delta`)
)



tmp.labels <- rbind(
  plt[1,] |> 
    dplyr::mutate(variable = 'Macrophages') |> 
    dplyr::mutate(label = paste0("R = ",round(cor.delta.purity.delta.epic.macrophages,2), "")) |> 
    dplyr::mutate(tumour.percentage.dna = 75) |> 
    dplyr::mutate(value = 0.45)
  ,
  plt[1,] |> 
    dplyr::mutate(variable = 'CD4 T-cells') |> 
    dplyr::mutate(label = paste0("R = ",round(cor.delta.purity.delta.epic.tcells,2), "")) |> 
    dplyr::mutate(tumour.percentage.dna = 75) |> 
    dplyr::mutate(value = 0.25)
  )


tmp.n.pairs.below.15 <- tmp.paired |>
  dplyr::filter(has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()
tmp.n.pairs.leq.15 <-  tmp.paired |>
  dplyr::filter(!has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()


ggplot(plt, aes(x=tumour.percentage.dna, y=value, group=pid, col= tumour.percentage.status, label=label)) + 
  facet_wrap(~variable, scales = "free") +
  geom_point(data = subset(plt, resection == "r1"), pch=19, cex=1.2 * 1.5, alpha=0.5) +
  geom_text(data = tmp.labels, col="black",fontface = "italic",hjust = 0, vjust=0.5) + 
  geom_path(  arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 ) + 
  scale_color_manual(values = c("increase"="#bb5f6c", "decrease"="#79b1b1"),name="") +
  xlim(0, 100) +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  )  +
  labs(x = "Tumor purity (%)", 
       y = 'EPIC deconvolution score',
       caption = paste0( "Pairs = ", (tmp.n.pairs.below.15 + tmp.n.pairs.leq.15), " (", tmp.n.pairs.leq.15, " >= 15%, ", tmp.n.pairs.below.15, " < 15%)"))
ggsave("output/figures/2022_Figure_M2G_M2H.pdf", width=8.3 * 2/3,height=8.3/3.5 , scale=2)





# figure 2ij ----
## figure 2i [t-cell] ----



tmp.n.pairs.below.15 <- tmp.paired |>
  dplyr::filter(has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()
tmp.n.pairs.leq.15 <-  tmp.paired |>
  dplyr::filter(!has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()



rho <- cor(tmp.paired$`EPIC: CD4 T-cells delta` , log(tmp.paired$time.to.progression),method="pearson")
#cor(tmp.paired$`EPIC: CD4 T-cells delta` , log(tmp.paired$time.to.progression),method="spearman")

df <- data.frame('EPIC: CD4 T-cells delta' = 0.15, time.to.progression=1250,label=paste0("R = ",round(rho,2)),check.names=F)



p1 <- ggplot(tmp.paired, aes(x=`EPIC: CD4 T-cells delta` , y=time.to.progression)) +
  geom_vline(xintercept=0,  color = "gray",lty='dotted') +
  stat_smooth(method="lm", se=FALSE,lwd=0.25,col="red",lty=1) +
  geom_point() +
  ggpubr::stat_cor(method = "pearson") +
  #geom_text(data=df,aes(label=label)) +
  scale_y_continuous(breaks=NULL,trans='log2') +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  labs(y = "Time between resections (log)",
       x="Increase EPIC CD4 T-cell score between resections",
       caption = paste0( "G-SAM: pairs = ", (tmp.n.pairs.below.15 + tmp.n.pairs.leq.15), " (", tmp.n.pairs.leq.15, " >= 15%, ", tmp.n.pairs.below.15, " < 15%)")
)








## figure 2j [macrophage] ----





#cor(tmp.paired$`EPIC: Macrophages delta` , tmp.paired$time.to.progression,method="pearson")
#cor(tmp.paired$`EPIC: Macrophages delta` , tmp.paired$time.to.progression,method="spearman")
rho <- cor(tmp.paired$`EPIC: Macrophages delta` , log(tmp.paired$time.to.progression),method="pearson")
#cor(tmp.paired$`EPIC: Macrophages delta` , log(tmp.paired$time.to.progression),method="spearman")


df <- data.frame('EPIC: Macrophages delta' = 0.20, time.to.progression=1250,label=paste0("R = ",round(rho,2)),check.names=F)



tmp.n.pairs.below.15 <- tmp.paired |>
  dplyr::filter(has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()
tmp.n.pairs.leq.15 <-  tmp.paired |>
  dplyr::filter(!has.low.purity.sample) |> 
  dplyr::pull(pid) |> 
  unique() |> 
  length()



p2 <- ggplot(tmp.paired , aes(x=`EPIC: Macrophages delta` , y=time.to.progression)) +
  geom_vline(xintercept=0,  color = "gray",lty='dotted') +
  stat_smooth(method="lm", se=FALSE,lwd=0.25,col="red",lty=1) +
  geom_point() +
  #geom_text(data=df,aes(label=label)) +
  ggpubr::stat_cor(method = "pearson") +
  scale_y_continuous(breaks=NULL, trans='log2') +
  ggplot2::theme_bw()  +
  ggplot2::theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    legend.position = 'bottom',
    
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.1)
  ) +
  labs(y = "Time between resections (log)",
       x="Increase EPIC Macrophage score between resections",
       caption = paste0( "G-SAM: pairs = ", (tmp.n.pairs.below.15 + tmp.n.pairs.leq.15), " (", tmp.n.pairs.leq.15, " >= 15%, ", tmp.n.pairs.below.15, " < 15%)")
       )



## export ----

p1 + p2


ggsave("output/figures/2022_figure_S4ij.pdf", width=8.3 / 2,height=8.3/4, scale=2)



#' # survival style C0-C4 ----
#' # too much redundancy with regular SVVL analysis
#' 
#' plt <- tmp.paired |> 
#'   dplyr::filter(!is.na(`rna.signature.C0.fuzzy.2022 delta`))
#' 
#' plt$ `rna.signature.C0.fuzzy.2022 delta` = scale(plt$`rna.signature.C0.fuzzy.2022 delta`)[,1]
#' plt$ `rna.signature.C1.collagen.2022 delta` = scale(plt$`rna.signature.C1.collagen.2022 delta`)[,1]
#' plt$ `rna.signature.C2.endothelial.2022 delta` = scale(plt$`rna.signature.C2.endothelial.2022 delta`)[,1]
#' plt$ `rna.signature.C3.oligodendrocyte.2022 delta` = scale(plt$`rna.signature.C3.oligodendrocyte.2022 delta`)[,1]
#' plt$ `rna.signature.C4.neuron.2022 delta` = scale(plt$`rna.signature.C4.neuron.2022 delta`)[,1]
#' 
#' plt <- plt |> 
#'   dplyr::rename( `delta C0 (fuzzy) signature` = `rna.signature.C0.fuzzy.2022 delta` ) |> 
#'   dplyr::rename( `delta C1 (collagen) signature` = `rna.signature.C1.collagen.2022 delta` ) |> 
#'   dplyr::rename( `delta C2 (endothelial) signature` = `rna.signature.C2.endothelial.2022 delta` ) |> 
#'   dplyr::rename( `delta C3 (oligodendrocyte) signature` = `rna.signature.C3.oligodendrocyte.2022 delta` ) |> 
#'   dplyr::rename( `delta C4 (neuron) signature` = `rna.signature.C4.neuron.2022 delta` )
#' 
#' 
#' surv_object <- survival::Surv(time = plt$time.to.progression)
#' fit.cox <- survival::coxph(surv_object ~
#'                            `delta C0 (fuzzy) signature` +
#'                            `delta C1 (collagen) signature` +
#'                            `delta C2 (endothelial) signature` +
#'                            `delta C3 (oligodendrocyte) signature` +
#'                            `delta C4 (neuron) signature` 
#'                            ,
#'                            data = plt)
#' survminer::ggforest(fit.cox, data = plt)
#' 
#' 
#' ggsave("output/figures/2022_figure_SXXG.pdf", width=8.3 / 1.5,height=8.3/6, scale=2)



# F] Figure S5 - Svvl ~ EPIC ----


plt <- tmp.paired
plt$ `EPIC: CD4 T-cells delta` <- scale(plt$ `EPIC: CD4 T-cells delta`)[, 1]
plt$ `EPIC: Macrophages delta` <- scale(plt$ `EPIC: Macrophages delta`)[, 1]
plt$`delta.purity` <- scale(tmp.paired$ `tumour.percentage.dna.delta`)[, 1]

plt <- plt |>
  dplyr::rename(`delta 'EPIC: CD4 T-cells' score` = `EPIC: CD4 T-cells delta`) |>
  dplyr::rename(`delta 'EPIC: Macrophages' score` = `EPIC: Macrophages delta`)


surv_object <- survival::Surv(time = plt$time.to.progression)
fit.cox <- survival::coxph(surv_object ~
  `delta 'EPIC: CD4 T-cells' score` +
  `delta 'EPIC: Macrophages' score`,
data = plt
)
survminer::ggforest(fit.cox, data = plt)
ggsave("output/figures/2022_Figure_S5.pdf", width = 8.3 / 1.5, height = 8.3 / 8, scale = 2)

