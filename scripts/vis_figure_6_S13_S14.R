#!/usr/bin/env R 

# load libs ----


library(patchwork) # avoid calling infix functions (+ and /) pls


#  load data ----


source('scripts/R/palette.R')

source('scripts/load_G-SAM_metadata.R')
source('scripts/load_G-SAM_expression_data.R')



delta <- function(x) {
  n <- length(x)
  
  return ( x[2:n] - x[1:(n-1)] )
}


## ggsurvplot export bugfix - https://github.com/kassambara/survminer/issues/152#issuecomment-938941051
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}


## prep data table ----


tmp.metadata <- gsam.rna.metadata |>
  dplyr::filter(.data$blacklist.pca == F) |>
  dplyr::filter(.data$pat.with.IDH == F) |>
  dplyr::filter(.data$sid %in% c("BAI2", "CAO1-replicate", "FAB2", "GAS2-replicate") == F) |>
  dplyr::filter(.data$batch != "old") |>
  dplyr::filter(.data$tumour.percentage.dna >= 15) |>
  dplyr::select(dplyr::contains("rna.signature") | `sid` | `pid` | `GITS.150.svm.2022.subtype` | `resection` | `extent` | `MGMT`) |>
  dplyr::rename(`Resection or Biopsy` = extent) |>
  dplyr::mutate(resection = ifelse(resection == "r1", "primary", "recurrence")) |>
  dplyr::filter(!is.na(.data$rna.signature.C1.collagen.2022))


tmp.metadata.paired <- tmp.metadata |>
  dplyr::mutate(MGMT = NULL) |>
  tidyr::pivot_wider(
    id_cols = pid,
    names_from = resection,
    values_from = -c(pid, resection)
  ) |>
  as.data.frame() |>
  dplyr::filter(!is.na(sid_primary) & !is.na(sid_recurrence)) |> # only complete pairs for these stats

  dplyr::mutate(delta.rna.signature.C0.fuzzy.2022 = rna.signature.C0.fuzzy.2022_recurrence - rna.signature.C0.fuzzy.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C1.collagen.2022 = rna.signature.C1.collagen.2022_recurrence - rna.signature.C1.collagen.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C2.endothelial.2022 = rna.signature.C2.endothelial.2022_recurrence - rna.signature.C2.endothelial.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C3.oligodendrocyte.2022 = rna.signature.C3.oligodendrocyte.2022_recurrence - rna.signature.C3.oligodendrocyte.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C4.neuron.2022 = rna.signature.C4.neuron.2022_recurrence - rna.signature.C4.neuron.2022_primary) |>
  dplyr::left_join(gsam.patient.metadata |>
    dplyr::select(
      studyID,

      # svvl
      survivalDays,
      survivalFromSecondSurgeryDays,
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
    ) |>
    dplyr::rename(event = status),
  by = c("pid" = "studyID"), suffix = c("", "")
  ) |>
  dplyr::mutate(event = ifelse(.data$event == "Deceased", 1, 0)) |>
  dplyr::mutate(Deceased = dplyr::recode(event, "1" = "Yes", "0" = "No")) |>
  dplyr::mutate(rank = order(order(delta.rna.signature.C1.collagen.2022, delta.rna.signature.C1.collagen.2022, pid))) |>
  dplyr::mutate(em.pc.status = ifelse(.data$`rna.signature.C1.collagen.2022_recurrence` > .data$`rna.signature.C1.collagen.2022_primary`, "increase", "decrease")) |>
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
  dplyr::rename(`GITS subtype R1` = .data$GITS.150.svm.2022.subtype_primary) |>
  dplyr::rename(`GITS subtype R2` = .data$GITS.150.svm.2022.subtype_recurrence) |>
  dplyr::mutate(`Age above 50` = ifelse(age > 50, "Yes", "No")) |>
  dplyr::mutate(gender = as.character(gender)) |>
  dplyr::rename(Sex = gender) |>
  dplyr::mutate(`KPS 70 or above` = factor(ifelse(is.na(performanceAtSecondSurgery) | performanceAtSecondSurgery >= 70, "Yes", "No"), levels = c("Yes", "No"))) |>
  dplyr::mutate(performanceAtSecondSurgery = as.factor(as.character(performanceAtSecondSurgery))) |>
  dplyr::mutate(daysToProgression = (survivalDays - survivalFromSecondSurgeryDays)) |>
  dplyr::mutate(progression.event = 1) |>
  dplyr::rename(`Resection/Biopsy R1` = `Resection or Biopsy_primary`) |>
  dplyr::rename(`Resection/Biopsy R2` = `Resection or Biopsy_recurrence`) |>
  dplyr::mutate(`tumorLocation` = factor(`tumorLocation`, levels = c("Temporal", "Frontal", "Parietal", "Fossa posterior", "Occipital"))) |>
  dplyr::rename(`Tumor location` = `tumorLocation`)





tmp.metadata.paired.breed <- tmp.metadata |>
  dplyr::mutate(MGMT = NULL) |>
  tidyr::pivot_wider(
    id_cols = pid,
    names_from = resection,
    values_from = -c(pid, resection)
  ) |>
  as.data.frame() |>
  #dplyr::filter(!is.na(sid_primary) & !is.na(sid_recurrence)) |> # only complete pairs for these stats
  
  dplyr::mutate(delta.rna.signature.C0.fuzzy.2022 = rna.signature.C0.fuzzy.2022_recurrence - rna.signature.C0.fuzzy.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C1.collagen.2022 = rna.signature.C1.collagen.2022_recurrence - rna.signature.C1.collagen.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C2.endothelial.2022 = rna.signature.C2.endothelial.2022_recurrence - rna.signature.C2.endothelial.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C3.oligodendrocyte.2022 = rna.signature.C3.oligodendrocyte.2022_recurrence - rna.signature.C3.oligodendrocyte.2022_primary) |>
  dplyr::mutate(delta.rna.signature.C4.neuron.2022 = rna.signature.C4.neuron.2022_recurrence - rna.signature.C4.neuron.2022_primary) |>
  dplyr::left_join(gsam.patient.metadata |>
                     dplyr::select(
                       studyID,
                       
                       # svvl
                       survivalDays,
                       survivalFromSecondSurgeryDays,
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
                     ) |>
                     dplyr::rename(event = status),
                   by = c("pid" = "studyID"), suffix = c("", "")
  ) |>
  dplyr::mutate(event = ifelse(.data$event == "Deceased", 1, 0)) |>
  dplyr::mutate(Deceased = dplyr::recode(event, "1" = "Yes", "0" = "No")) |>
  dplyr::mutate(rank = order(order(delta.rna.signature.C1.collagen.2022, delta.rna.signature.C1.collagen.2022, pid))) |>
  dplyr::mutate(em.pc.status = ifelse(.data$`rna.signature.C1.collagen.2022_recurrence` > .data$`rna.signature.C1.collagen.2022_primary`, "increase", "decrease")) |>
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
  dplyr::rename(`GITS subtype R1` = .data$GITS.150.svm.2022.subtype_primary) |>
  dplyr::rename(`GITS subtype R2` = .data$GITS.150.svm.2022.subtype_recurrence) |>
  dplyr::mutate(`Age above 50` = ifelse(age > 50, "Yes", "No")) |>
  dplyr::mutate(gender = as.character(gender)) |>
  dplyr::rename(Sex = gender) |>
  dplyr::mutate(`KPS 70 or above` = factor(ifelse(is.na(performanceAtSecondSurgery) | performanceAtSecondSurgery >= 70, "Yes", "No"), levels = c("Yes", "No"))) |>
  dplyr::mutate(performanceAtSecondSurgery = as.factor(as.character(performanceAtSecondSurgery))) |>
  dplyr::mutate(daysToProgression = (survivalDays - survivalFromSecondSurgeryDays)) |>
  dplyr::mutate(progression.event = 1) |>
  dplyr::rename(`Resection/Biopsy R1` = `Resection or Biopsy_primary`) |>
  dplyr::rename(`Resection/Biopsy R2` = `Resection or Biopsy_recurrence`) |>
  dplyr::mutate(`tumorLocation` = factor(`tumorLocation`, levels = c("Temporal", "Frontal", "Parietal", "Fossa posterior", "Occipital"))) |>
  dplyr::rename(`Tumor location` = `tumorLocation`)







tmp.metadata <- tmp.metadata |>
  dplyr::left_join(
    tmp.metadata.paired |> dplyr::select(pid, rank, em.pc.status,
    ),# stats only accessible through pairing
    by=c('pid'='pid'), suffix=c('','')
  ) |> 

dplyr::left_join(
  gsam.patient.metadata |>
    dplyr::select(
      studyID,

      # svvl
      survivalDays,
      survivalFromSecondSurgeryDays,
      status,

      # treat
      treatedWithTMZ,
      treatedWithRT,
      bevacizumab.before.recurrence,
      PTK787.before.recurrence, # also angio

      # loc
      tumorLocation,
      
      # mgmt
      mgmtStability,
      
      # asl
      age,
      gender,
      performanceAtSecondSurgery
    ) |>
    dplyr::rename(event = status), by = c("pid" = "studyID"), suffix = c("", "")) |>
  dplyr::mutate(event = ifelse(.data$event == "Deceased", 1, 0)) |>
  dplyr::mutate(`MGMT meth` = dplyr::recode(mgmtStability,
    "Stable methylated" = "Stable",
    "Stable unmethylated" = "Wildtype"
  ), mgmtStability = NULL) |>
  dplyr::mutate(`Treatment: Beva` = case_when(
    is.na(bevacizumab.before.recurrence) ~ "NA",
    bevacizumab.before.recurrence == "Trial participant" ~ "NA",
    bevacizumab.before.recurrence == "Yes" ~ "Yes",
    T ~ "No"
  )) |>
  dplyr::mutate(`Treatment: Beva (randomized)` = case_when(
    is.na(bevacizumab.before.recurrence) ~ "NA",
    bevacizumab.before.recurrence == "Trial participant" ~ "Yes",
    T ~ "No"
  )) |>
  dplyr::mutate(bevacizumab.before.recurrence = NULL) |>
  dplyr::rename(`Treatment: PT787` = PTK787.before.recurrence) |>
  dplyr::rename(`Treatment: TMZ` = treatedWithTMZ) |>
  dplyr::rename(`Treatment: RT` = treatedWithRT) |> 
  dplyr::mutate(daysToProgression = (survivalDays - survivalFromSecondSurgeryDays)) |>
  dplyr::mutate(progression.event = 1) |> 
  dplyr::mutate(`Age above 50` = ifelse(age > 50, "Yes", "No")) |> 
  dplyr::mutate(gender = as.character(gender)) |>
  dplyr::rename(Sex = gender) |> 
  dplyr::mutate(`KPS 70 or above` = factor(ifelse(is.na(performanceAtSecondSurgery) | performanceAtSecondSurgery >= 70, "Yes", "No"), levels = c("Yes", "No"))) |> 
  dplyr::mutate(`tumorLocation` = factor(`tumorLocation`, levels = c("Temporal", "Frontal", "Parietal", "Fossa posterior", "Occipital"))) |>
  dplyr::rename(`Tumor location` = `tumorLocation`)


## figure S13a: determine cut-off C0 / fuz ----


d.prim <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C0.fuzzy.2022_primary) |>
  dplyr::pull(rna.signature.C0.fuzzy.2022_primary)
d.rec <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C0.fuzzy.2022_recurrence) |>
  dplyr::pull(rna.signature.C0.fuzzy.2022_recurrence)


#c0.fuz.cutoff.p <- d[(k-1):k] |> sum() |> (\(x){x/2})()
c0.fuz.cutoff.p <- d.rec |> median()


plt.1 <- rbind(
  data.frame(y = d.rec) |> 
    dplyr::mutate(type = 'recurrence') |> 
    dplyr::mutate(x = order(order(1:n()))),
  data.frame(y = d.prim) |> 
    dplyr::mutate(type = 'primary') |> 
    dplyr::mutate(x = order(order(1:n())))
  )

plt.2 <- data.frame(y = delta(d.rec)) |> 
  dplyr::rename(`delta C0/fuzzy at Rec.` = y ) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=y, col=type, x=x)) +
  geom_hline(yintercept=c0.fuz.cutoff.p, lwd=1.5, color = "red") +
  annotate(geom="text", x=25, y=9, label="\nCutoff: median\nat recurrence",size=3) +
  geom_line(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_line(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  xlim(1,max(plt.1$x)) +
  labs(y='C0/fuzzy') + 
  scale_color_manual(name = "ordered at", values = c('primary'='gray60','recurrence'='black')) +
  #scale_color_manual(name = "ordered at", values = resection_colors[c('primary','recurrence')]) +
  #scale_alpha_manual(name=NULL, values=c('primary'=0.4,'recurrence'=1),guide="none") +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL)

p2 <- ggplot(plt.2, aes(y=`delta C0/fuzzy at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_line() +
  geom_point() +
  xlim(1,max(plt.1$x)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL, y="1st derivative")


p1 / p2



#ggsave("output/figures/2022_figure_S13a.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
#ggsave("output/figures/2022_figure_S13a.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


rm(p1, p2, plt.1, plt.2, d.prim, d.rec)



## figure S13b: determine cut-off C1 / col ----


d.prim <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C1.collagen.2022_primary) |>
  dplyr::pull(rna.signature.C1.collagen.2022_primary)
d.rec <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C1.collagen.2022_recurrence) |>
  dplyr::pull(rna.signature.C1.collagen.2022_recurrence)


k <- 82+1
c1.em.cutoff.p <- d.rec[(k-1):k] |>
  sum() |>
  (\(x){x/2})() #equivalent to: (function(x){return(x/2)})()


plt.1 <- rbind(
  data.frame(y = d.rec) |> 
    dplyr::mutate(type = 'recurrence') |> 
    dplyr::mutate(x = order(order(1:n()))),
  data.frame(y = d.prim) |> 
    dplyr::mutate(type = 'primary') |> 
    dplyr::mutate(x = order(order(1:n())))
)

plt.2 <- data.frame(y = delta(d.rec)) |> 
  dplyr::rename(`delta C1/col at Rec.` = y) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=y, col=type, x=x)) +
  geom_hline(yintercept=c1.em.cutoff.p, lwd=1.5, color = "red") +
  annotate(geom="text", x=40, y=14, label="\nCutoff: change in rate\nat recurrence", size=3) +
  geom_line(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_line(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  xlim(1,max(plt.1$x)) +
  labs(y='C1/collagen') + 
  scale_color_manual(name = "ordered at", values = c('primary'='gray60','recurrence'='black')) +
  #scale_color_manual(name = "ordered at", values = resection_colors[c('primary','recurrence')]) +
  #scale_alpha_manual(name=NULL, values=c('primary'=0.5,'recurrence'=1),guide="none") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL)

p2 <- ggplot(plt.2, aes(y=`delta C1/col at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,max(plt.1$x)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )  +
  labs(x=NULL, y="1st derivative")


p1 / p2



#ggsave("output/figures/2022_figure_S13b.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
#ggsave("output/figures/2022_figure_S13b.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


rm(p1, p2, plt.1, plt.2, k, d.prim, d.rec)




## figure S13c: determine cut-off C2 / endo ----


d.prim <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C2.endothelial.2022_primary) |>
  dplyr::pull(rna.signature.C2.endothelial.2022_primary)
d.rec <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C2.endothelial.2022_recurrence) |>
  dplyr::pull(rna.signature.C2.endothelial.2022_recurrence)


#k <- 52
#c2.endo.cutoff.p <- d[(k-1):k] |> sum() |> (\(x){x/2})()
c2.endo.cutoff.p <- d.rec |>
  median()


plt.1 <- rbind(
  data.frame(y = d.rec) |> 
    dplyr::mutate(type = 'recurrence') |> 
    dplyr::mutate(x = order(order(1:n()))),
  data.frame(y = d.prim) |> 
    dplyr::mutate(type = 'primary') |> 
    dplyr::mutate(x = order(order(1:n())))
)

plt.2 <- data.frame(y = delta(d.rec)) |> 
  dplyr::rename(`delta C2/endo at Rec.` = y) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=y, col=type,x=x)) +
  geom_hline(yintercept=c2.endo.cutoff.p, lwd=1.5, color = "red") +
  annotate(geom="text", x=25, y=6, label="\nCutoff: median\nat recurrence", size=3) +
  geom_line(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_line(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  xlim(1,max(plt.1$x)) +
  labs(y='C2/endo') + 
  scale_color_manual(name = "ordered at", values = c('primary'='gray60','recurrence'='black')) +
  #scale_color_manual(name = "ordered at", values = resection_colors[c('primary','recurrence')]) +
  #scale_alpha_manual(name=NULL, values=c('primary'=0.5,'recurrence'=1),guide="none") +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL)

p2 <- ggplot(plt.2, aes(y=`delta C2/endo at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  #geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,max(plt.1$x)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )  +
  labs(x=NULL, y="1st derivative")


p1 / p2



#ggsave("output/figures/2022_figure_S13c.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
#ggsave("output/figures/2022_figure_S13c.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


rm(p1, p2, plt.1, plt.2, d.prim, d.rec)




## figure S13d: determine cut-off C3 / olig ----


d.prim <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C3.oligodendrocyte.2022_primary) |>
  dplyr::pull(rna.signature.C3.oligodendrocyte.2022_primary)
d.rec <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C3.oligodendrocyte.2022_recurrence) |>
  dplyr::pull(rna.signature.C3.oligodendrocyte.2022_recurrence)



#k <- 76
#c3.olig.cutoff.p <- d[(k-1):k] |> sum() |> (\(x){x/2})()
rm(k)
c3.olig.cutoff.p <- d.rec |>
  median()


plt.1 <- rbind(
  data.frame(y = d.rec) |> 
    dplyr::mutate(type = 'recurrence') |> 
    dplyr::mutate(x = order(order(1:n()))),
  data.frame(y = d.prim) |> 
    dplyr::mutate(type = 'primary') |> 
    dplyr::mutate(x = order(order(1:n())))
)

plt.2 <- data.frame(y = delta(d.rec)) |> 
  dplyr::rename(`delta C3/olig at Rec.` = y) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)



p1 <- ggplot(plt.1, aes(y=y, col=type, x=x)) +
  geom_hline(yintercept=c3.olig.cutoff.p, lwd=1.5, color = "red") +
  annotate(geom="text", x=25, y=25, label="\nCutoff: median\nat recurrence", size=3) +
  geom_line(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_line(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  xlim(1,max(plt.1$x)) +
  labs(y='C3/olig') + 
  scale_color_manual(name = "ordered at", values = c('primary'='gray60','recurrence'='black')) +
  #scale_color_manual(name = "ordered at", values = resection_colors[c('primary','recurrence')]) +
  #scale_alpha_manual(name=NULL, values=c('primary'=0.5,'recurrence'=1),guide="none") +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL)

p2 <- ggplot(plt.2, aes(y=`delta C3/olig at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  #geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,max(plt.1$x)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL, y="1st derivative")


p1 / p2




#ggsave("output/figures/2022_figure_S13d.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
#ggsave("output/figures/2022_figure_S13d.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


rm(p1, p2, plt.1, plt.2, d.prim, d.rec)





## figure S13e: determine cut-off C4 / neur ----


d.prim <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C4.neuron.2022_primary) |>
  dplyr::pull(rna.signature.C4.neuron.2022_primary)
d.rec <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C4.neuron.2022_recurrence) |>
  dplyr::pull(rna.signature.C4.neuron.2022_recurrence)


k <- 83
c4.neur.cutoff.p <- d.rec[(k-1):k] |>
  sum() |>
  (\(x){x/2})()


plt.1 <- rbind(
  data.frame(y = d.rec) |> 
    dplyr::mutate(type = 'recurrence') |> 
    dplyr::mutate(x = order(order(1:n()))),
  data.frame(y = d.prim) |> 
    dplyr::mutate(type = 'primary') |> 
    dplyr::mutate(x = order(order(1:n())))
)

plt.2 <- data.frame(y = delta(d.rec)) |> 
  dplyr::rename(`delta C4/neu at Rec.` = y) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)




p1 <- ggplot(plt.1, aes(y=y, col=type,x=x)) +
  geom_hline(yintercept=c4.neur.cutoff.p, lwd=1.5, color = "red") +
  annotate(geom="text", x=40, y=38, label="\nCutoff: change in rate\nat recurrence", size=3) +
  geom_line(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_line(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  xlim(1,max(plt.1$x)) +
  labs(y='C4/Neuron') + 
  scale_color_manual(name = "ordered at", values = c('primary'='gray60','recurrence'='black')) +
  #scale_color_manual(name = "ordered at", values = resection_colors[c('primary','recurrence')]) +
  #scale_alpha_manual(name=NULL, values=c('primary'=0.5,'recurrence'=1),guide="none") +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL)

p2 <- ggplot(plt.2, aes(y=`delta C4/neu at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,max(plt.1$x)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL, y="1st derivative") 


p1 / p2




#ggsave("output/figures/2022_figure_S13e.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
#ggsave("output/figures/2022_figure_S13e.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


rm(p1, p2, plt.1, plt.2, k, d.prim, d.rec)



## figure S13f: determine cut-off C1 GLASS ----


d.prim <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |>
  dplyr::filter(resection == "TP") |> 
  dplyr::arrange(rna.signature.C1.collagen.2022) |>
  dplyr::pull(rna.signature.C1.collagen.2022)
d.rec <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |>
  dplyr::filter(resection != "TP") |> 
  dplyr::arrange(rna.signature.C1.collagen.2022) |>
  dplyr::pull(rna.signature.C1.collagen.2022)



k <- 97
c1.em.cutoff.p.glass <- d.rec[(k-1):k] |>
  sum() |>
  (\(x){x/2})() #equivalent to: (function(x){return(x/2)})()

plt.1 <- rbind(
  data.frame(y = d.rec) |> 
    dplyr::mutate(type = 'recurrence') |> 
    dplyr::mutate(x = order(order(1:n()))),
  data.frame(y = d.prim) |> 
    dplyr::mutate(type = 'primary') |> 
    dplyr::mutate(x = order(order(1:n())))  |> 
    dplyr::mutate(x = x * length(d.rec) / length(d.prim))
)

plt.2 <- data.frame(y = delta(d.rec)) |> 
  dplyr::rename(`delta C1/col at Rec.` = y) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=y, col=type, x=x)) +
  geom_hline(yintercept=c1.em.cutoff.p.glass, lwd=1.5, color = "red") +
  annotate(geom="text", x=40, y=14, label="\nCutoff: change in rate\nat recurrence", size=3) +
  geom_line(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "primary")) +
  geom_line(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  geom_point(data = plt.1 |> dplyr::filter(type == "recurrence")) +
  xlim(1,max(plt.1$x)) +
  labs(y='C1/collagen') + 
  scale_color_manual(name = "ordered at", values = c('primary'='gray60','recurrence'='black')) +
  #scale_color_manual(name = "ordered at", values = resection_colors[c('primary','recurrence')]) +
  #scale_alpha_manual(name=NULL, values=c('primary'=0.5,'recurrence'=1),guide="none") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  labs(x=NULL)

p2 <- ggplot(plt.2, aes(y=`delta C1/col at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,max(plt.1$x)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_blank(),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )  +
  labs(x=NULL, y="1st derivative") + 
  coord_cartesian(ylim=c(0, 2.0))


p1 / p2



#ggsave("output/figures/2022_figure_S13f.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)




## plot(s) ----

### figure 6a: svvl ----


plt <- tmp.metadata.paired |> 
  dplyr::select(pid, rank, event, daysToProgression, survivalFromSecondSurgeryDays) |> 
  dplyr::mutate(daysToProgression =  -1 * daysToProgression) |> 
  tidyr::pivot_longer(cols=c(daysToProgression,survivalFromSecondSurgeryDays))


p1 <- ggplot(plt, aes(x=reorder(pid, rank),y=value, group=pid)) +
  geom_line(lwd=2) +
  geom_hline(yintercept=0,col="white") +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.45, size=5.5)  ) +
  labs(x = NULL, y="Survival")


#ggsave("output/figures/2022_figure_6a.pdf", width=8.3 / 2,height=8.3/4 * 0.5, scale=2)





### figure 6b: signature arrows ----


plt <- tmp.metadata |>  # only resections from complete pairs
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup()


p2 <- ggplot(plt, aes(x = reorder(pid, rank), y=`rna.signature.C1.collagen.2022`, col=em.pc.status, group=pid))  +
  geom_hline(yintercept=2.5, lty=1, color = "#FFFFFF44",lwd=3) +
  geom_point(data = subset(plt, resection == "primary"), pch=19, cex=1.1, alpha=0.65) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125*0.65, "inches")) , alpha = 0.75 )  + 
  labs(x=NULL, col=NULL, y="C0 signature"
       #caption="G-SAM: n=122 pairs"
       ) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c(
    "increase"="#bb5f6c", # rood #bb5f6c
    "decrease"="#79b1b1" # lichtlauw #79b1b
  )
    , guide="none" # zit in c
  ) +
  geom_hline(yintercept=c1.em.cutoff.p, lty=3, color = "red",lwd=0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.45, size=5.5)  )


rm(plt)



### figure 6c: squares ----

#'@todo x-check @ fig2 fgh -- mgmt


plt <- tmp.metadata.paired |>
  dplyr::select(-em.pc.status, -sid_primary, -sid_recurrence,
                -survivalDays, -survivalFromSecondSurgeryDays, -daysToProgression, -event, -progression.event,
                -age, -performanceAtSecondSurgery) |> 
  dplyr::select(-contains(".signature")) |> 
  reshape2::melt(id.vars=c("pid", "rank")) |> 
  dplyr::mutate(value = gsub('^Gained$','Gained / increased / female',value)) %>% 
  dplyr::mutate(value = gsub('^Lost$','Lost / decreased / male',value)) %>% 
  
  dplyr::mutate(value = gsub('^Female$','Gained / increased / female',value)) %>% 
  dplyr::mutate(value = gsub('^Male$','Lost / decreased / male',value)) %>% 
  
  dplyr::mutate(value = gsub('^Yes|Stable|Resection$','Yes / stable / resection',value)) %>% 
  dplyr::mutate(value = gsub('^No|Wildtype|Biopsy$','No / wildtype / biopsy',value)) %>% 
  
  dplyr::mutate(panel = case_when(
    grepl("Deceased|Age|Sex|KPS|Biops", variable) ~ "A",
    grepl("Treatment", variable) ~ "B",
    grepl("subtype", variable) ~ "C",
    grepl("cnStatus", variable) ~ "D",
    variable %in% c("HM",  "MGMT meth") ~ "F",
    T ~ "E"
  ))


p3 <- ggplot(plt, aes(x = reorder(pid, rank), y = variable, fill = value)) +
  facet_grid(rows=vars(panel), scales="free", space="free") + 
  geom_tile(colour = "black", size = 0.3) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
        axis.title.y=element_blank(),
        axis.ticks.y =element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size=5),
        legend.key.size = unit(0.3, 'cm')) +
  scale_fill_manual(values = c( subtype_colors ,
                                #"Wildtype"="white",
                                #"No"="white",
                                "No / wildtype / biopsy"="white",
                                
                                "Gained / increased / female"="#bb5f6c", # rood #bb5f6c
                                "Lost / decreased / male"="#79b1b1", # lichtlauw #79b1b1
                                "Yes / stable / resection"="#2e415e", # donker blauw #2e415e
                                #"Yes" = "#2e415e", # zelfde als stable #2e415e
                                
                                "NA"="grey")) + 
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  labs(caption="G-SAM: n=122 pairs")

#p3



#### export ----

p1 / p2 / p3 + patchwork::plot_layout(heights = c(0.7,0.7,1.95))


ggsave("output/figures/2022_figure_6abc.pdf", width=8.3 / 2,height=8.3/2 * 0.75 * 0.93, scale=2.2)




### stats ? ----





# svvl analysis ----


tmp.metadata.paired <- tmp.metadata.paired |>
  
  dplyr::mutate(`C0/fuzzy signature at Prim.` = factor(ifelse(rna.signature.C0.fuzzy.2022_primary >           c0.fuz.cutoff.p,  "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C1/col signature at Prim.`   = factor(ifelse(rna.signature.C1.collagen.2022_primary >        c1.em.cutoff.p,   "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C2/endo signature at Prim.`  = factor(ifelse(rna.signature.C2.endothelial.2022_primary >     c2.endo.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C3/olig signature at Prim.`  = factor(ifelse(rna.signature.C3.oligodendrocyte.2022_primary > c3.olig.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C4/neu signature at Prim.`   = factor(ifelse(rna.signature.C4.neuron.2022_primary >          c4.neur.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  
  dplyr::mutate(`C0.fuzzy.signature.prim` = `C0/fuzzy signature at Prim.`) |> 
  dplyr::mutate(`C1.col.signature.prim`   = `C1/col signature at Prim.`  ) |> 
  dplyr::mutate(`C2.endo.signature.prim`  = `C2/endo signature at Prim.` ) |> 
  dplyr::mutate(`C3.olig.signature.prim`  = `C3/olig signature at Prim.` ) |> 
  dplyr::mutate(`C4.neu.signature.prim`   = `C4/neu signature at Prim.`  ) |> 
  
  
  dplyr::mutate(`C0/fuzzy signature at Rec.` = factor(ifelse(rna.signature.C0.fuzzy.2022_recurrence >           c0.fuz.cutoff.p,  "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C1/col signature at Rec.`   = factor(ifelse(rna.signature.C1.collagen.2022_recurrence >        c1.em.cutoff.p,   "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C2/endo signature at Rec.`  = factor(ifelse(rna.signature.C2.endothelial.2022_recurrence >     c2.endo.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C3/olig signature at Rec.`  = factor(ifelse(rna.signature.C3.oligodendrocyte.2022_recurrence > c3.olig.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C4/neu signature at Rec.`   = factor(ifelse(rna.signature.C4.neuron.2022_recurrence >          c4.neur.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  
  dplyr::mutate(`C0.fuzzy.signature.rec` = `C0/fuzzy signature at Rec.`) |> 
  dplyr::mutate(`C1.col.signature.rec`   = `C1/col signature at Rec.`  ) |> 
  dplyr::mutate(`C2.endo.signature.rec`  = `C2/endo signature at Rec.` ) |> 
  dplyr::mutate(`C3.olig.signature.rec`  = `C3/olig signature at Rec.` ) |> 
  dplyr::mutate(`C4.neu.signature.rec`   = `C4/neu signature at Rec.`  ) |> 
  
  dplyr::mutate(`MGMT at Rec.` = factor(
    dplyr::recode(`MGMT meth`,
                  'Gained'='Methylated',
                  'Stable'='Methylated',
                  'Wildtype'='Unmethylated',
                  'Lost'='Unmethylated'
    ) ,levels=c('Unmethylated','Methylated'))) |> 
  dplyr::mutate(`MGMT at Prim.` = factor(
  dplyr::recode(`MGMT meth`,
                'Gained'='Unmethylated',
                'Stable'='Methylated',
                'Wildtype'='Unmethylated',
                'Lost'='Methylated'
  ), levels=c('Unmethylated','Methylated')))




tmp.metadata.paired.breed <- tmp.metadata.paired.breed |>
  
  dplyr::mutate(`C0/fuzzy signature at Prim.` = factor(ifelse(rna.signature.C0.fuzzy.2022_primary >           c0.fuz.cutoff.p,  "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C1/col signature at Prim.`   = factor(ifelse(rna.signature.C1.collagen.2022_primary >        c1.em.cutoff.p,   "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C2/endo signature at Prim.`  = factor(ifelse(rna.signature.C2.endothelial.2022_primary >     c2.endo.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C3/olig signature at Prim.`  = factor(ifelse(rna.signature.C3.oligodendrocyte.2022_primary > c3.olig.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C4/neu signature at Prim.`   = factor(ifelse(rna.signature.C4.neuron.2022_primary >          c4.neur.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  
  dplyr::mutate(`C0.fuzzy.signature.prim` = `C0/fuzzy signature at Prim.`) |> 
  dplyr::mutate(`C1.col.signature.prim`   = `C1/col signature at Prim.`  ) |> 
  dplyr::mutate(`C2.endo.signature.prim`  = `C2/endo signature at Prim.` ) |> 
  dplyr::mutate(`C3.olig.signature.prim`  = `C3/olig signature at Prim.` ) |> 
  dplyr::mutate(`C4.neu.signature.prim`   = `C4/neu signature at Prim.`  ) |> 
  
  
  dplyr::mutate(`C0/fuzzy signature at Rec.` = factor(ifelse(rna.signature.C0.fuzzy.2022_recurrence >           c0.fuz.cutoff.p,  "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C1/col signature at Rec.`   = factor(ifelse(rna.signature.C1.collagen.2022_recurrence >        c1.em.cutoff.p,   "high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C2/endo signature at Rec.`  = factor(ifelse(rna.signature.C2.endothelial.2022_recurrence >     c2.endo.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C3/olig signature at Rec.`  = factor(ifelse(rna.signature.C3.oligodendrocyte.2022_recurrence > c3.olig.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  dplyr::mutate(`C4/neu signature at Rec.`   = factor(ifelse(rna.signature.C4.neuron.2022_recurrence >          c4.neur.cutoff.p ,"high", "low"), levels=c('low','high'))) |> 
  
  dplyr::mutate(`C0.fuzzy.signature.rec` = `C0/fuzzy signature at Rec.`) |> 
  dplyr::mutate(`C1.col.signature.rec`   = `C1/col signature at Rec.`  ) |> 
  dplyr::mutate(`C2.endo.signature.rec`  = `C2/endo signature at Rec.` ) |> 
  dplyr::mutate(`C3.olig.signature.rec`  = `C3/olig signature at Rec.` ) |> 
  dplyr::mutate(`C4.neu.signature.rec`   = `C4/neu signature at Rec.`  ) |> 
  
  dplyr::mutate(`MGMT at Rec.` = factor(
    dplyr::recode(`MGMT meth`,
                  'Gained'='Methylated',
                  'Stable'='Methylated',
                  'Wildtype'='Unmethylated',
                  'Lost'='Unmethylated'
    ) ,levels=c('Unmethylated','Methylated'))) |> 
  dplyr::mutate(`MGMT at Prim.` = factor(
    dplyr::recode(`MGMT meth`,
                  'Gained'='Unmethylated',
                  'Stable'='Methylated',
                  'Wildtype'='Unmethylated',
                  'Lost'='Methylated'
    ), levels=c('Unmethylated','Methylated')))


## stats prevalence MES @ C1 high ----


tmp.stats <- tmp.metadata |> 
  dplyr::filter(resection == "recurrence") |> 
  #dplyr::mutate(`C1/col signature` = factor(ifelse(rna.signature.C1.collagen.2022 > c1.em.cutoff.p, "high", "low"), levels=c('low','high'))) |> 
  
  dplyr::mutate(GITS.150.svm.2022.subtype = ifelse(GITS.150.svm.2022.subtype != "Mesenchymal", "Non-MES", "MES")) |> 
  dplyr::select(GITS.150.svm.2022.subtype, rna.signature.C1.collagen.2022 )

wilcox.test(
  tmp.stats |>
    dplyr::filter(GITS.150.svm.2022.subtype == "MES") |> 
    dplyr::pull(rna.signature.C1.collagen.2022),
  tmp.stats |>
    dplyr::filter(GITS.150.svm.2022.subtype == "Non-MES") |> 
    dplyr::pull(rna.signature.C1.collagen.2022)
  )



# figure 6d: signature x svvl ----


plt <- tmp.metadata.paired |> 
  dplyr::filter(event == 1) |> 
  dplyr::select(pid, `C1/col signature at Rec.`,
                rna.signature.C1.collagen.2022_primary, rna.signature.C1.collagen.2022_recurrence,
                survivalDays, survivalFromSecondSurgeryDays) |> 
  dplyr::rename(survivalDays_primary = survivalDays) |> 
  dplyr::rename(survivalDays_recurrence = survivalFromSecondSurgeryDays) |> 
  #dplyr::mutate(survival.from.first.surgery = ifelse(survivalDays > 365* 2.5, "> 2.5yr", "<= 2.5yr") ) %>%
  #dplyr::mutate(survival.from.second.surgery = factor(survival.from.second.surgery, levels=c("> 1yr", "<= 1yr")))
  tidyr::pivot_longer(cols=-c(pid, `C1/col signature at Rec.`),  names_to = c(".value", "resection"), names_pattern = "(.+)_(.+)")


#ggplot(plt, aes(y = extracellular.matrix.component , x=tts, group=pid, fill=em.status.res2 )) 

ggplot(plt, aes(y = rna.signature.C1.collagen.2022 , x=survivalDays, group=pid, fill=`C1/col signature at Rec.`)) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.3, col="gray50" )  + 
  geom_point(data = subset(plt, resection == "primary"), pch=21, cex=1.2, alpha=0.7, fill="gray50", col="gray50") +
  geom_point(data = subset(plt, resection == "recurrence"), pch=19, cex=4.5, alpha=0.25, col="white") +
  geom_point(data = subset(plt, resection == "recurrence"), pch=21, cex=2.0, alpha=0.7) +
  labs(y = "C0 / col signature", x = "Survival") +
  geom_hline(yintercept = c1.em.cutoff.p, lty=1, lwd=3.5, col="#FFFFFFBB") +
  geom_hline(yintercept = c1.em.cutoff.p, lty=2, col="red", lwd=0.25) +
  geom_vline(xintercept = 0, lty=1, col="black", lwd=0.35) +
  geom_vline(xintercept = -6, lty=1, col="black", lwd=0.35) +
  scale_x_reverse(breaks = (0:15) * 365) +
  scale_fill_manual(values = c('high'='#009E74', 'low'='#CB75A4') ) +
  theme_bw()  +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    #strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    #axis.text.x = element_blank(),
    legend.position = 'bottom',
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )


ggsave("output/figures/2022_figure_6d.pdf", width=8.3 / 2 * 0.65,height=8.3/5 *0.85, scale=2.2)





## R2 col high ----
### figure 6e: KM R2 -> death (ppt) ----
# unpaired 



surv_object <- survival::Surv(time = tmp.metadata.paired.breed$survivalFromSecondSurgeryDays, event=tmp.metadata.paired.breed$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired.breed)
pval.ppt.r2 <- survminer::surv_pvalue(fit1)$pval # post progression time
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired.breed, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#CB75A4',0.7),
                              'C1/col signature: low'=alpha('#009E74',0.7)
                            ),
                            legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                          'C1.col.signature=low'='C1/col signature: low'),
                            xlab="Survival time from recurrence (days)")
p1


ggsave("output/figures/2022_figure_6e.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)





# paired
# surv_object <- survival::Surv(time = tmp.metadata.paired$survivalFromSecondSurgeryDays, event=tmp.metadata.paired$event)
# fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired)
# p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
#                             palette = c(
#                               'C1/col signature: high'=alpha('#CB75A4',0.7),
#                               'C1/col signature: low'=alpha('#009E74',0.7)
#                             ),
#                             legend.labs=c('C1.col.signature=high'='C1/col signature: high',
#                                           'C1.col.signature=low'='C1/col signature: low'),
#                             xlab="Survival time from recurrence")
# p1



### figure 6f: KM R1 -> R2 (ttp) ----
# unpaired 

surv_object <- survival::Surv(time = tmp.metadata.paired.breed$daysToProgression, 
                              event= tmp.metadata.paired.breed$progression.event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired.breed)
pval.ttp.r2 <- survminer::surv_pvalue(fit1)$pval
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired.breed, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#CB75A4',0.7),
                              'C1/col signature: low'=alpha('#009E74',0.7)
                            ),
                            legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                          'C1.col.signature=low'='C1/col signature: low'),
                            xlab="Time to progression (days)")
p1


ggsave("output/figures/2022_figure_6f.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)



# ## paired
# 
# surv_object <- survival::Surv(time = tmp.metadata.paired$daysToProgression, event=tmp.metadata.paired$progression.event)
# fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired)
# p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
#                             palette = c(
#                               'C1/col signature: high'=alpha('#CB75A4',0.7),
#                               'C1/col signature: low'=alpha('#009E74',0.7)
#                             ),
#                             legend.labs=c('C1.col.signature=high'='C1/col signature: high',
#                                           'C1.col.signature=low'='C1/col signature: low'),
#                             xlab="Time to progression")
# p1
# 




### figure 6g: KM R1 -> R2 (os) ----
# unpaired 

surv_object <- survival::Surv(time = tmp.metadata.paired.breed$survivalDays, event=tmp.metadata.paired.breed$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired.breed)
pval.os.r2 <- survminer::surv_pvalue(fit1)$pval
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired.breed, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#CB75A4',0.7),
                              'C1/col signature: low'=alpha('#009E74',0.7)
                            ),
                            legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                          'C1.col.signature=low'='C1/col signature: low'),
                            xlab="Overall survival (days)")
p1

ggsave("output/figures/2022_figure_6g.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)



# paired ---
# surv_object <- survival::Surv(time = tmp.metadata.paired$survivalDays, event=tmp.metadata.paired$event)
# fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired)
# p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
#                             palette = c(
#                               'C1/col signature: high'=alpha('#CB75A4',0.7),
#                               'C1/col signature: low'=alpha('#009E74',0.7)
#                             ),
#                             legend.labs=c('C1.col.signature=high'='C1/col signature: high',
#                                           'C1.col.signature=low'='C1/col signature: low'),
#                             xlab="Overall survival")
# p1




### figure 6efg+S14abc padj stats ----

df = data.frame(pval = c(pval.ppt.r2,
                         pval.ttp.r2,
                         pval.os.r2,
                         
                         pval.ppt.r1,
                         pval.ttp.r1,
                         pval.os.r1
                         )) |>
  dplyr::mutate(padj = p.adjust(pval,method="fdr"))
df


rm(df)


### figure S14d OS GLASS ----



svvl <- glass.gbm.rnaseq.metadata.all.samples |> 
  dplyr::filter(tumour.percentage.2022 >= 15) |>
  dplyr::filter(resection != "TP") |> 
  dplyr::mutate(`col.sig` = factor(ifelse(rna.signature.C1.collagen.2022 > c1.em.cutoff.p.glass , "high", "low"), levels=c('low','high')))
  #dplyr::arrange(rna.signature.C1.collagen.2022) |>
  #dplyr::pull(rna.signature.C1.collagen.2022)

plot(sort(svvl$rna.signature.C1.collagen.2022))
#abline(h=5.6)
abline(h=4.3)



surv_object <- survival::Surv(time = svvl$case_overall_survival_mo, event=svvl$os.event)
fit1 <- survival::survfit(surv_object ~  col.sig , data = svvl)
print(survminer::surv_pvalue(fit1)$pval)
p1 <- survminer::ggsurvplot(fit1, data = svvl, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#f346b3',0.7),
                              'C1/col signature: low'=alpha('#4384fd',0.7)
                            ),
                            legend.labs=c('col.sig=high'='C1/col signature: high',
                                          'col.sig=low'='C1/col signature: low'),
                            xlab="Overall survival (months)")
p1

ggsave("output/figures/2022_figure_S14d.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)







### figure 6h: ggforest R2 -> death ----
# unpaired


svvl <- tmp.metadata.paired.breed |> 
  dplyr::filter(!is.na(`C1/col signature at Rec.`)) |> 
  dplyr::filter(!is.na(`Tumor location`))



surv_object <- survival::Surv(time = svvl$survivalFromSecondSurgeryDays,
                              event= svvl$event)
fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             `Sex` +
                             `KPS 70 or above` +
                             `Treatment: Beva` +
                             `Treatment: TMZ` +
                             `Tumor location` +

                             #`MGMT meth` + # incomplete data

                             `C0/fuzzy signature at Rec.` +
                             `C1/col signature at Rec.` +
                             `C2/endo signature at Rec.` +
                             `C3/olig signature at Rec.` +
                             `C4/neu signature at Rec.`
                           ,
                           data = svvl)
survminer::ggforest(fit.cox, data = svvl)



data.frame(pval = summary(fit.cox)$coefficients[,5]) |>
  dplyr::mutate(padj = p.adjust(pval, method="fdr")) |>
  dplyr::mutate(padj.f = format.pval(padj, digits=2))



ggsave("output/figures/2022_figure_6h.pdf", width=8.3 / 2,height=8.3/3.4, scale=2)


sum(is.na(tmp.metadata.paired.breed$`Age above 50` ))
sum(is.na(tmp.metadata.paired.breed$`Sex` ))
sum(is.na(tmp.metadata.paired.breed$`KPS 70 or above` ))
sum(is.na(tmp.metadata.paired.breed$`Treatment: Beva` ))
sum(is.na(tmp.metadata.paired.breed$`Treatment: TMZ` ))
sum(is.na(tmp.metadata.paired.breed$`Tumor location` ))

sum(is.na(tmp.metadata.paired.breed$`C1/col signature at Rec.` ))
sum(is.na(tmp.metadata.paired.breed$`C2/endo signature at Rec.` ))
sum(is.na(tmp.metadata.paired.breed$`C3/olig signature at Rec.` ))
sum(is.na(tmp.metadata.paired.breed$`C4/neu signature at Rec.`))



# paired samples only
# 
# surv_object <- survival::Surv(time = tmp.metadata.paired$survivalFromSecondSurgeryDays, event=tmp.metadata.paired$event)
# fit.cox <- survival::coxph(surv_object ~
#                              `Age above 50` +
#                              `Sex` +
#                              `KPS 70 or above` +
#                              `Treatment: Beva` +
#                              `Treatment: TMZ` +
#                              `Tumor location` +
#                              
#                              #`MGMT meth` + # incomplete data
#                              
#                              #`C0/fuzzy signature at Rec.` +
#                              `C1/col signature at Rec.` +
#                              `C2/endo signature at Rec.` +
#                              `C3/olig signature at Rec.` +
#                              `C4/neu signature at Rec.`
#                            ,
#                            data = tmp.metadata.paired)
# survminer::ggforest(fit.cox, data = tmp.metadata.paired)
# 
# 
# 
# data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
#   dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
#   dplyr::mutate(padj.f = format.pval(padj, digits=2))








### figure S14e: ggforest R2 -> death [MGMT] ----



tmp.metadata.mgmt <- tmp.metadata.paired.breed |> 
  dplyr::filter(!is.na(`C1/col signature at Rec.`)) |> 
  dplyr::filter(!is.na(`MGMT meth`) & !is.na(`C1/col signature at Rec.`))
#dplyr::mutate(`MGMT at Rec.` = as.character(`MGMT at Rec.`)) |> 
#dplyr::mutate(`MGMT at Rec.` = ifelse(is.na(`MGMT at Rec.`), "unknown", `MGMT at Rec.`))
#dplyr::mutate(`MGMT meth` = ifelse(is.na(`MGMT meth`), "unknown", `MGMT meth`))


surv_object <- survival::Surv(time = tmp.metadata.mgmt$survivalFromSecondSurgeryDays, 
                              event = tmp.metadata.mgmt$event)
fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             `Sex` +
                             `KPS 70 or above` +
                             `Treatment: Beva` +
                             `Treatment: TMZ` +
                             `Tumor location` +
                             
                             `MGMT meth` +
                             
                             `C1/col signature at Rec.`
                           
                           ,
                           data = tmp.metadata.mgmt)
survminer::ggforest(fit.cox, data = tmp.metadata.mgmt)


data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
  dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
  dplyr::mutate(padj.f = format.pval(padj, digits=3))



ggsave("output/figures/2022_figure_S14e.pdf", width=8.3 / 2,height=8.3/3.4, scale=2)




# # paired
# 
# tmp.metadata.paired.mgmt <- tmp.metadata.paired |> 
#   dplyr::filter(!is.na(`MGMT meth`)) |> 
#   #dplyr::mutate(`MGMT at Rec.` = as.character(`MGMT at Rec.`)) |> 
#   dplyr::mutate(`MGMT at Rec.` = ifelse(is.na(`MGMT at Rec.`), "unknown", `MGMT at Rec.`))
# #dplyr::mutate(`MGMT meth` = ifelse(is.na(`MGMT meth`), "unknown", `MGMT meth`))
# 
# 
# surv_object <- survival::Surv(time = tmp.metadata.paired.mgmt$survivalFromSecondSurgeryDays, 
#                               event = tmp.metadata.paired.mgmt$event)
# fit.cox <- survival::coxph(surv_object ~
#                              `Age above 50` +
#                              `Sex` +
#                              `KPS 70 or above` +
#                              `Treatment: Beva` +
#                              `Treatment: TMZ` +
#                              `Tumor location` +
#                              
#                              `MGMT meth` +
#                              
#                              `C1/col signature at Rec.`
#                            
#                            ,
#                            data = tmp.metadata.paired.mgmt)
# survminer::ggforest(fit.cox, data = tmp.metadata.paired.mgmt)
# 
# 
# data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
#   dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
#   dplyr::mutate(padj.f = format.pval(padj, digits=1))



### figure 6i: ggforest R1 -> R2 (ttp) ----
# unpaired


svvl <- tmp.metadata.paired.breed |> 
  dplyr::filter(!is.na(`C1/col signature at Rec.`)) |> 
  dplyr::filter(!is.na(`Tumor location`))



surv_object <- survival::Surv(time = svvl$daysToProgression, 
                              event= svvl$progression.event)
fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             `Sex` +
                             `KPS 70 or above` +
                             `Treatment: Beva` +
                             `Treatment: TMZ` +
                             `Tumor location` +
                             
                             
                             #`MGMT meth` + # incomplete data
                             
                             `C0/fuzzy signature at Rec.` +
                             `C1/col signature at Rec.` +
                             `C2/endo signature at Rec.` +
                             `C3/olig signature at Rec.` +
                             `C4/neu signature at Rec.`
                           ,
                           data = svvl)
survminer::ggforest(fit.cox, data = svvl)


data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
  dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
  dplyr::mutate(padj.f = format.pval(padj, digits=1))



ggsave("output/figures/2022_figure_6i.pdf", width=8.3 / 2,height=8.3/3.4, scale=2)


# paired
# surv_object <- survival::Surv(time = tmp.metadata.paired$daysToProgression, event=tmp.metadata.paired$progression.event)
# fit.cox <- survival::coxph(surv_object ~
#                              `Age above 50` +
#                              `Sex` +
#                              `KPS 70 or above` +
#                              `Treatment: Beva` +
#                              `Treatment: TMZ` +
#                              `Tumor location` +
#                              
# 
#                              #`MGMT meth` + # incomplete data
#                              
#                              `C0/fuzzy signature at Rec.` +
#                              `C1/col signature at Rec.` +
#                              `C2/endo signature at Rec.` +
#                              `C3/olig signature at Rec.` +
#                              `C4/neu signature at Rec.`
#                            ,
#                            data = tmp.metadata.paired)
# survminer::ggforest(fit.cox, data = tmp.metadata.paired)
# 
# 
# 
# data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
#   dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
#   dplyr::mutate(padj.f = format.pval(padj, digits=1))
# 







### figure S14f: ggforest R1 -> R2 (ttp) [MGMT] ----

# unpaired
tmp.metadata.mgmt <- tmp.metadata.paired.breed |> 
  dplyr::filter(!is.na(`MGMT meth`)) |> 
  dplyr::filter(!is.na(`C1/col signature at Rec.`)) |> 
  #dplyr::mutate(`MGMT at Rec.` = as.character(`MGMT at Rec.`)) |> 
  dplyr::mutate(`MGMT at Rec.` = ifelse(is.na(`MGMT at Rec.`), "unknown", `MGMT at Rec.`))
#dplyr::mutate(`MGMT meth` = ifelse(is.na(`MGMT meth`), "unknown", `MGMT meth`))


surv_object <- survival::Surv(time = tmp.metadata.mgmt$daysToProgression, event=tmp.metadata.mgmt$progression.event)
fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             `Sex` +
                             `KPS 70 or above` +
                             `Treatment: Beva` +
                             `Treatment: TMZ` +
                             `Tumor location` +
                             
                             
                             `MGMT meth` +
                             
                             `C1/col signature at Rec.`
                           
                           ,
                           data = tmp.metadata.mgmt)
survminer::ggforest(fit.cox, data = tmp.metadata.mgmt)

data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
  dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
  dplyr::mutate(padj.f = format.pval(padj, digits=2))


sum(is.na( tmp.metadata.mgmt$`Age above 50`  ))
sum(is.na( tmp.metadata.mgmt$`Sex`  ))
sum(is.na( tmp.metadata.mgmt$`KPS 70 or above`  ))
sum(is.na( tmp.metadata.mgmt$`Treatment: Beva`  ))
sum(is.na( tmp.metadata.mgmt$`Treatment: TMZ`  ))
sum(is.na( tmp.metadata.mgmt$`Tumor location`  ))
sum(is.na( tmp.metadata.mgmt$`MGMT meth`  ))
sum(is.na( tmp.metadata.mgmt$`C1/col signature at Rec.` ))


ggsave("output/figures/2022_figure_S14f.pdf", width=8.3 / 2,height=8.3/3.4, scale=2)







# paired
# tmp.metadata.paired.mgmt <- tmp.metadata.paired |> 
#   dplyr::filter(!is.na(`MGMT meth`)) |> 
#   #dplyr::mutate(`MGMT at Rec.` = as.character(`MGMT at Rec.`)) |> 
#   dplyr::mutate(`MGMT at Rec.` = ifelse(is.na(`MGMT at Rec.`), "unknown", `MGMT at Rec.`))
# #dplyr::mutate(`MGMT meth` = ifelse(is.na(`MGMT meth`), "unknown", `MGMT meth`))
# 
# 
# surv_object <- survival::Surv(time = tmp.metadata.paired.mgmt$daysToProgression, event=tmp.metadata.paired.mgmt$progression.event)
# fit.cox <- survival::coxph(surv_object ~
#                              `Age above 50` +
#                              `Sex` +
#                              `KPS 70 or above` +
#                              `Treatment: Beva` +
#                              `Treatment: TMZ` +
#                              `Tumor location` +
#                              
# 
#                              `MGMT meth` +
#                              
#                              `C1/col signature at Rec.`
#                            
#                            ,
#                            data = tmp.metadata.paired.mgmt)
# survminer::ggforest(fit.cox, data = tmp.metadata.paired.mgmt)
# 
# data.frame(pval = summary(fit.cox)$coefficients[,5]) |> 
#   dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
#   dplyr::mutate(padj.f = format.pval(padj, digits=1))
# 
# 
# 
# 



## R1 col high ----
### figure 14a: KM R2 -> death (ppt) ----


# unpaired
surv_object <- survival::Surv(time = tmp.metadata.paired.breed$survivalFromSecondSurgeryDays, 
                              event= tmp.metadata.paired.breed$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.prim` , data = tmp.metadata.paired.breed)
pval.ppt.r1 <- survminer::surv_pvalue(fit1)$pval # post progression time
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired.breed, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#CB75A4',0.7),
                              'C1/col signature: low'=alpha('#009E74',0.7)
                            ),
                            legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                          'C1.col.signature=low'='C1/col signature: low'),
                            xlab="Survival time from recurrence (days)")
p1


ggsave("output/figures/2022_figure_14a.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)




# paired
# surv_object <- survival::Surv(time = tmp.metadata.paired$survivalFromSecondSurgeryDays, event=tmp.metadata.paired$event)
# fit1 <- survival::survfit(surv_object ~  `C1.col.signature.prim` , data = tmp.metadata.paired)
# p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
#                             palette = c(
#                               'C1/col signature: high'=alpha('#CB75A4',0.7),
#                               'C1/col signature: low'=alpha('#009E74',0.7)
#                             ),
#                             legend.labs=c('C1.col.signature=high'='C1/col signature: high',
#                                           'C1.col.signature=low'='C1/col signature: low'),
#                             xlab="Survival time from recurrence")
# p1



### figure 14b: KM R1 -> R2 (ttp) ----


surv_object <- survival::Surv(time = tmp.metadata.paired.breed$daysToProgression,
                              event= tmp.metadata.paired.breed$progression.event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.prim` , data = tmp.metadata.paired.breed)
pval.ttp.r1 <- survminer::surv_pvalue(fit1)$pval
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired.breed, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#CB75A4',0.7),
                              'C1/col signature: low'=alpha('#009E74',0.7)
                            ),
                            legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                          'C1.col.signature=low'='C1/col signature: low'),
                            xlab="Time to progression (days)")
p1


ggsave("output/figures/2022_figure_14b.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)



### figure 14c: KM R1 -> R2 (os) ----


surv_object <- survival::Surv(time = tmp.metadata.paired.breed$survivalDays,
                              event=tmp.metadata.paired.breed$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.prim` , data = tmp.metadata.paired.breed)
pval.os.r1 <- survminer::surv_pvalue(fit1)$pval
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired.breed, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                            palette = c(
                              'C1/col signature: high'=alpha('#CB75A4',0.7),
                              'C1/col signature: low'=alpha('#009E74',0.7)
                            ),
                            legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                          'C1.col.signature=low'='C1/col signature: low'),
                            xlab="Overall survival (days)")
p1



ggsave("output/figures/2022_figure_14c.pdf", width=8.3 / 2 * 0.8,height=8.3/3.4, scale=2, plot=p1)








# figure S13g: C1 x TMZ ----
stopifnot(nrow(tmp.metadata) == 287)
unique(length(tmp.metadata$pid))

plt_r1 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "primary") |>
  dplyr::select(sid, pid, `Treatment: TMZ`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "primary")
unique(length(plt_r1$pid))

plt_r2 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "recurrence") |>
  dplyr::select(sid, pid, `Treatment: TMZ`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "recurrence")
unique(length(plt_r2$pid))


plt_both <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::select(sid, pid, `Treatment: TMZ`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "combined")
unique(length(plt_both$pid))


plt <- rbind(plt_r1, plt_r2, plt_both) |>
  #dplyr::filter(pid %in% tmp.metadata.paired$pid) |> 
  dplyr::mutate(type = factor(type, levels = c("primary", "recurrence", "combined")))



ggplot(plt, aes(x = `Treatment: TMZ`, y = rna.signature.C1.collagen.2022)) +
  facet_grid(cols = vars(type)) +
  ggbeeswarm::geom_quasirandom() +
  ggsignif::geom_signif(
    comparisons = list(c("Yes", "No")),
    test = "wilcox.test",
    col = "black",
    tip_length = 0
  ) +
  theme_bw() +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(caption = paste0("G-SAM: n=", nrow(plt_both), " samples"), y = "C1/col signature") +
  scale_y_continuous(expand = expansion(mult = .075))


ggsave("output/figures/2022_figure_S13g.pdf", width=8.3 / 2,height=8.3/4.5, scale=2)



# figure S13h: C1 x Beva ----


tmp.metadata <- tmp.metadata |> 
  dplyr::mutate(`Treatment: Beva: factor` = case_when(
    `Treatment: Beva (randomized)` == "Yes" ~ "Rand. trial",
    `Treatment: Beva` == "Yes" ~ "Yes",
    `Treatment: Beva` == "No" ~ "No"
  )) |> 
  dplyr::mutate(`Treatment: Beva: factor`= factor(`Treatment: Beva: factor`, levels=c("No", "Yes", "Rand. trial")))

plt_r1 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "primary") |>
  dplyr::select(sid, pid, `Treatment: Beva: factor`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "primary")


plt_r2 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "recurrence") |>
  dplyr::select(sid, pid, `Treatment: Beva: factor`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "recurrence")


plt_both <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::select(sid, pid, `Treatment: Beva: factor`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "combined")


plt <- rbind(plt_r1, plt_r2, plt_both) |>
  #dplyr::filter(pid %in% tmp.metadata.paired$pid) |> 
  dplyr::mutate(type = factor(type, levels = c("primary", "recurrence", "combined")))



ggplot(plt, aes(x = `Treatment: Beva: factor`, y = rna.signature.C1.collagen.2022, col=`Treatment: Beva: factor`)) +
  facet_grid(cols = vars(type)) +
  ggbeeswarm::geom_quasirandom() +
  ggsignif::geom_signif(
    comparisons = list(c("Yes", "No")),
    test = "wilcox.test",
    col = "black",
    tip_length = 0
  ) +
  theme_bw() +
  scale_color_manual(values=c('Yes'='black', 'No'='black','Rand. trial'='gray60'), guide="none") +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(caption = paste0("G-SAM: n=", nrow(plt_both), " samples"), y = "C1/col signature",x="Treatment: Beva") +
  scale_y_continuous(expand = expansion(mult = .075))



ggsave("output/figures/2022_figure_S13h.pdf", width=8.3 / 2,height=8.3/4.5, scale=2)


# figure S13i: C1 x Biopsy~Resection ----



plt_r1 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "primary") |>
  dplyr::select(sid, pid, `Resection or Biopsy`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "primary")


plt_r2 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "recurrence") |>
  dplyr::select(sid, pid, `Resection or Biopsy`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "recurrence")


plt_both <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::select(sid, pid, `Resection or Biopsy`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "combined")


plt <- rbind(plt_r1, plt_r2, plt_both) |>
  dplyr::filter(pid %in% tmp.metadata.paired$pid) |> 
  dplyr::mutate(`Resection or Biopsy` = ifelse(is.na(`Resection or Biopsy`),"NA",`Resection or Biopsy`)) |> 
  dplyr::mutate(`Resection or Biopsy` = factor(`Resection or Biopsy`, levels = c("Resection", "Biopsy", "NA"))) |> 
  dplyr::mutate(type = factor(type, levels = c("primary", "recurrence", "combined")))



ggplot(plt, aes(x = `Resection or Biopsy`, y = rna.signature.C1.collagen.2022, col=`Resection or Biopsy`)) +
  facet_grid(cols = vars(type)) +
  ggbeeswarm::geom_quasirandom() +
  ggsignif::geom_signif(
    comparisons = list(c("Resection", "Biopsy")),
    test = "wilcox.test",
    col = "black",
    tip_length = 0
  ) +
  theme_bw() +
  scale_color_manual(values=c('Resection'='black', 'Biopsy'='black','NA'='gray60'), guide="none") +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(caption = paste0("G-SAM: n=", nrow(plt_both), " samples"), y = "C1/col signature") +
  scale_y_continuous(expand = expansion(mult = .075))



ggsave("output/figures/2022_figure_S13i.pdf", width=8.3 / 2,height=8.3/4.5, scale=2)




# figure S13j: C1 x MGMT ----




plt_r1 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "primary") |>
  dplyr::select(sid, pid, `MGMT`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "primary")


plt_r2 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "recurrence") |>
  dplyr::select(sid, pid, `MGMT`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "recurrence")


plt_both <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::select(sid, pid, `MGMT`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "combined")


plt <- rbind(plt_r1, plt_r2, plt_both) |>
  #dplyr::filter(pid %in% tmp.metadata.paired$pid) |> 
  dplyr::mutate(`MGMT` = ifelse(is.na(`MGMT`),"NA",`MGMT`)) |> 
  dplyr::mutate(`MGMT` = factor(`MGMT`, levels = c("Unmethylated", "Methylated", "NA"))) |> 
  dplyr::mutate(type = factor(type, levels = c("primary", "recurrence", "combined")))



ggplot(plt, aes(x = `MGMT`, y = rna.signature.C1.collagen.2022, col=MGMT)) +
  facet_grid(cols = vars(type)) +
  ggbeeswarm::geom_quasirandom() +
  ggsignif::geom_signif(
    comparisons = list(
      c("Methylated", "Unmethylated")
    ),
    test = "wilcox.test",
    col = "black",
    tip_length = 0
  ) +
  theme_bw() +
  scale_color_manual(values=c('Methylated'='black', 'Unmethylated'='black','NA'='gray60'),guide="none") +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(caption = paste0("G-SAM: n=", nrow(plt_both), " samples"), y = "C1/col signature", x="MGMT methylation status") +
  scale_y_continuous(expand = expansion(mult = .075))



ggsave("output/figures/2022_figure_S13j.pdf", width=8.3 / 2,height=8.3/4.5, scale=2)




# figure S13k: C1 x TumorLocation ----



plt_r1 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "primary") |>
  dplyr::select(sid, pid, `Tumor location`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "primary")


plt_r2 <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::filter(resection == "recurrence") |>
  dplyr::select(sid, pid, `Tumor location`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "recurrence")


plt_both <- tmp.metadata |>
  dplyr::filter(!is.na(rna.signature.C1.collagen.2022)) |>
  dplyr::select(sid, pid, `Tumor location`, rna.signature.C1.collagen.2022) |>
  dplyr::mutate(type = "combined")


plt <- rbind(plt_r1, plt_r2, plt_both) |>
  #dplyr::filter(pid %in% tmp.metadata.paired$pid) |> 
  dplyr::mutate(`Tumor location` = as.character(`Tumor location`)) |> 
  #dplyr::mutate(`Tumor location` = ifelse(`Tumor location` == "Fossa posterior","Fossa\nposterior",`Tumor location`)) |> 
  dplyr::mutate(`Tumor location` = factor(`Tumor location`, levels=c("Temporal","Occipital","Frontal","Parietal","Fossa posterior", "NA"))) |> 
  dplyr::mutate(type = factor(type, levels = c("primary", "recurrence", "combined")))



ggplot(plt, aes(x = `Tumor location`, y = rna.signature.C1.collagen.2022, col=`Tumor location`)) +
  facet_grid(cols = vars(type)) +
  ggbeeswarm::geom_quasirandom() +
  ggpubr::stat_compare_means(method = "anova", label.x = 4.25) +
  # ggsignif::geom_signif( # becomes a visual mess
  #   comparisons = list(
  #     c("Temporal", "Occipital"),
  #     c("Temporal", "Frontal"),
  #     c("Temporal", "Parietal"),
  #     c("Temporal", "Fossa\nposterior"),
  #     
  #     c("Occipital","Frontal"),
  #     c("Occipital","Parietal"),
  #     c("Occipital","Fossa\nposterior"),
  #     
  #     
  #     c("Frontal", "Fossa\nposterior"),
  #     c("Parietal", "Fossa\nposterior")
  #     
  #   ),
  #   y_position=c(14,16,18,20,22,24,26,14,16),
  #   test = "wilcox.test",
  #   col = "black",
  #   tip_length = 0
  # ) +
  theme_bw() +
  scale_color_manual(values=c(
    "Temporal"="black",
    "Occipital"="black",
    "Frontal"="black",
    "Parietal"="black",
    "Fossa posterior"="black",
    'NA'='gray60'), guide="none") +
  theme(
    # text = element_text(family = 'Arial'), seems to require a postscript equivalent
    # strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold", size = rel(1)),
    # axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.25)
  ) +
  labs(caption = paste0("G-SAM: n=", nrow(plt_both), " samples"), y = "C1/col signature") +
  scale_y_continuous(expand = expansion(mult = .075))


ggsave("output/figures/2022_figure_S13k.pdf", width=8.3 / 1,height=8.3/4.5, scale=2)


