#!/usr/bin/env R 

# load libs ----


library(patchwork) # avoid calling infix functions pls


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
  dplyr::filter(.data$blacklist.pca == F)  |> 
  dplyr::filter(.data$pat.with.IDH == F) |> 
  dplyr::filter(.data$sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::filter(.data$batch != 'old') |> 
  dplyr::filter(.data$tumour.percentage.dna >= 15) |> 
  dplyr::select(dplyr::contains("rna.signature") | `sid` | `pid` | `GITS.150.svm.2022.subtype` | `resection`) |> 
  dplyr::mutate(resection = ifelse(resection == "r1","primary","recurrence")) |> 
  dplyr::filter(!is.na(.data$rna.signature.C1.collagen.2022))


tmp.metadata.paired <- tmp.metadata |> 
  tidyr::pivot_wider(id_cols =  pid,
                     names_from = resection, 
                     values_from = -c(pid, resection)) |> 
  as.data.frame() |> 
  dplyr::filter(!is.na(sid_primary) & !is.na(sid_recurrence)) |> # only complete pairs for these stats

  dplyr::mutate(delta.rna.signature.C0.fuzzy.2022 = rna.signature.C0.fuzzy.2022_recurrence - rna.signature.C0.fuzzy.2022_primary) |> 
  dplyr::mutate(delta.rna.signature.C1.collagen.2022 = rna.signature.C1.collagen.2022_recurrence - rna.signature.C1.collagen.2022_primary) |> 
  dplyr::mutate(delta.rna.signature.C2.endothelial.2022 = rna.signature.C2.endothelial.2022_recurrence - rna.signature.C2.endothelial.2022_primary) |> 
  dplyr::mutate(delta.rna.signature.C3.oligodendrocyte.2022 = rna.signature.C3.oligodendrocyte.2022_recurrence - rna.signature.C3.oligodendrocyte.2022_primary) |> 
  dplyr::mutate(delta.rna.signature.C4.neuron.2022 = rna.signature.C4.neuron.2022_recurrence - rna.signature.C4.neuron.2022_primary) |> 

  dplyr::left_join(
    gsam.patient.metadata |> 
      dplyr::select(
        studyID,
        
        survivalDays, 
        survivalFromSecondSurgeryDays,
        status,
        
        age,
        gender,
        performanceAtSecondSurgery,
        
        treatedWithTMZ,
        treatedWithRT,
        bevacizumab.before.recurrence,
        
        mgmtStability,HM,
        
        AXIN2,APC,JAK2,
        RB1,MSH2,BRCA1,
        BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,
        NF1,ERBB3,
        EGFR,TP53BP1,TP53,PIK3R1,PIK3CA,TSC2,
        SETD2,PDGFRA,
        
        cnStatusCDKN2ABs,
        cnStatusEGFRs,
        cnStatusRB1s,
        cnStatusNF1s,
        cnStatusCDK4s,
        cnStatusMDM2s
        
        ) |> 
      dplyr::rename(event = status)
      , by=c('pid'='studyID'), suffix=c('','')
  ) |> 
  dplyr::mutate(event = ifelse(.data$event == "Deceased",1,0)) |> 
  dplyr::mutate(Deceased = dplyr::recode(event, "1" = "Yes", "0" = "No" )) |> 
  dplyr::mutate(rank = order(order(delta.rna.signature.C1.collagen.2022, delta.rna.signature.C1.collagen.2022, pid))) |> 
  dplyr::mutate(em.pc.status = ifelse(.data$`rna.signature.C1.collagen.2022_recurrence` > .data$`rna.signature.C1.collagen.2022_primary`, "increase", "decrease")) |> 
  dplyr::mutate(`MGMT meth` = dplyr::recode( mgmtStability,
                                             "Stable methylated" = "Stable",
                                             "Stable unmethylated" = "Wildtype" ), mgmtStability = NULL ) |> 
  dplyr::mutate(`Treatment: Beva` = ifelse(bevacizumab.before.recurrence, "Yes","No"), bevacizumab.before.recurrence=NULL) |> 
  dplyr::rename(`Treatment: TMZ` = treatedWithTMZ) |> 
  dplyr::rename(`Treatment: RT` = treatedWithRT) |> 
  dplyr::rename(`GITS subtype R1` = .data$GITS.150.svm.2022.subtype_primary) |> 
  dplyr::rename(`GITS subtype R2` = .data$GITS.150.svm.2022.subtype_recurrence) |> 
  dplyr::mutate(`Age above 50` = ifelse(age > 50,"Yes","No")) |> 
  dplyr::mutate(gender = as.character(gender)) |> 
  dplyr::rename(Sex = gender) |> 
  dplyr::mutate(`KPS 70 or above` = factor(ifelse(is.na(performanceAtSecondSurgery) | performanceAtSecondSurgery >= 70, "Yes","No"), levels=c('Yes','No'))) |> 
  dplyr::mutate(performanceAtSecondSurgery = as.factor(as.character(performanceAtSecondSurgery))) |> 
  dplyr::mutate(daysToProgression = (survivalDays - survivalFromSecondSurgeryDays)) |> 
  dplyr::mutate(progression.event = 1)





tmp.metadata <- tmp.metadata |>
  dplyr::left_join(
    tmp.metadata.paired |> dplyr::select(pid, rank, em.pc.status ),# stats only accessible through pairing
    by=c('pid'='pid'),suffix=c('','')
  ) |> 
  dplyr::left_join(
    gsam.patient.metadata |> 
      dplyr::select(
        studyID,
        
        survivalDays, 
        survivalFromSecondSurgeryDays,
        status,
        
        age,
        gender
      ) |> 
      dplyr::rename(event = status)
    , by=c('pid'='studyID'), suffix=c('','')
  ) |> 
  dplyr::mutate(event = ifelse(.data$event == "Deceased",1,0)) |> 
  dplyr::mutate(`age.above.50` = ifelse(age > 50,"Yes","No"))




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



ggsave("output/figures/2022_figure_S13a.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
ggsave("output/figures/2022_figure_S13a.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


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



ggsave("output/figures/2022_figure_S13b.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
ggsave("output/figures/2022_figure_S13b.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


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



ggsave("output/figures/2022_figure_S13c.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
ggsave("output/figures/2022_figure_S13c.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


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




ggsave("output/figures/2022_figure_S13d.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
ggsave("output/figures/2022_figure_S13d.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


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




ggsave("output/figures/2022_figure_S13e.pdf", width=8.3 / 5,height=8.3/5*1.17, scale=2)
ggsave("output/figures/2022_figure_S13e.svg", width=8.3 / 5,height=8.3/5*1.17, scale=2)


rm(p1, p2, plt.1, plt.2, k, d.prim, d.rec)




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
  dplyr::mutate(value = gsub('^Yes|Stable$','Yes / stable',value)) %>% 
  dplyr::mutate(value = gsub('^No|Wildtype$','No / wildtype',value)) %>% 
  dplyr::mutate(panel = case_when(
    grepl("Deceased|Age|Sex|KPS", variable) ~ "A",
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
                                "No / wildtype"="white",
                                
                                "Gained / increased / female"="#bb5f6c", # rood #bb5f6c
                                "Lost / decreased / male"="#79b1b1", # lichtlauw #79b1b1
                                "Yes / stable"="#2e415e", # donker blauw #2e415e
                                #"Yes" = "#2e415e", # zelfde als stable #2e415e
                                
                                "NA"="grey")) + 
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom') +
  labs(caption="G-SAM: n=122 pairs")



#### export ----

p1 / p2 / p3 + patchwork::plot_layout(heights = c(0.8,0.8,1.95))


ggsave("output/figures/2022_figure_6abc.pdf", width=8.3 / 2,height=8.3/2 * 0.75, scale=2.2)




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
                  )
    #case_when
    #is.na(`MGMT meth`) ~ as.character(NA),
    #`MGMT meth` %in% c("Gained", "Stable") ~ "Methylated",
    #T ~ "Unmethylated"
  ,levels=c('Unmethylated','Methylated')))


### figure 6d: signature x svvl ----



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





## R2 ----


### figure 6e: KM collagen ----

surv_object <- survival::Surv(time = tmp.metadata.paired$survivalFromSecondSurgeryDays, event=tmp.metadata.paired$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired)
p1 <- survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      palette = c(
                        'C1/col signature: high'=alpha('#CB75A4',0.7),
                        'C1/col signature: low'=alpha('#009E74',0.7)
                      ),
                      legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                    'C1.col.signature=low'='C1/col signature: low'),
                      xlab="Survival time from recurrence")

ggsave("output/figures/2022_figure_6e.pdf", width=8.3 / 2 * 0.45,height=8.3/5 * 0.7, scale=3,  plot = p1)


### figure s6f COX ----

#plot(tmp.metadata.paired$performanceAtSecondSurgery, tmp.metadata.paired$delta.rna.signature.C1.collagen.2022)


# fit.cox <- coxph(surv_object ~ `C6.ECM.signature.R2.status`  , data = tmp.metadata.paired)
# ggforest(fit.cox)
# 
# fit.cox <- coxph(surv_object ~ `C5.signature.R2.status`  , data = tmp.metadata.paired)
# ggforest(fit.cox)
# 
# fit.cox <- coxph(surv_object ~ `C4.signature.R2.status`  , data = tmp.metadata.paired)
# ggforest(fit.cox)
# 
# fit.cox <- coxph(surv_object ~ `C3.endothelial.signature.R2.status`  , data = tmp.metadata.paired)
# ggforest(fit.cox)
# 
# fit.cox <- coxph(surv_object ~ `C2.oligodendrocyte.signature.R2.status` , data = tmp.metadata.paired)
# ggforest(fit.cox)
# 
# fit.cox <- coxph(surv_object ~ `C2.oligodendrocyte.signature.R2.status` , data = tmp.metadata.paired)
# ggforest(fit.cox)



### figure 13x: KM collagen full ----


surv_object <- survival::Surv(time = tmp.metadata.paired$survivalFromSecondSurgeryDays, event=tmp.metadata.paired$event)
fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             `Sex` +
                             `KPS 70 or above` +
                             `Treatment: Beva` +
                             `Treatment: TMZ` +
                             `C0/fuzzy signature at Rec.` +
                             `C1/col signature at Rec.` +
                             `C2/endo signature at Rec.` +
                             `C3/olig signature at Rec.` +
                             `C4/neu signature at Rec.`
                           ,
                           data = tmp.metadata.paired)
survminer::ggforest(fit.cox, data = tmp.metadata.paired)
ggsave("output/figures/2022_figure_S13f.pdf", width=8.3 / 2,height=8.3/3.4, scale=2)




#### MGMT at Rec. ----

tmp.metadata.paired.mgmt <- tmp.metadata.paired |> 
  dplyr::filter(!is.na(`MGMT at Rec.`))

surv_object <- survival::Surv(time = tmp.metadata.paired.mgmt$survivalFromSecondSurgeryDays, event=tmp.metadata.paired.mgmt$event)
fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             `Sex` +
                             `KPS 70 or above` +
                             `Treatment: Beva` +
                             `Treatment: TMZ` +
                             `MGMT at Rec.` +
                             `C0/fuzzy signature at Rec.` +
                             `C1/col signature at Rec.` +
                             `C2/endo signature at Rec.` +
                             `C3/olig signature at Rec.` +
                             `C4/neu signature at Rec.`
                           ,
                           data = tmp.metadata.paired.mgmt)
survminer::ggforest(fit.cox, data = tmp.metadata.paired.mgmt)
ggsave("output/figures/2022_figure_S13g.pdf", width=8.3 / 2,height=8.3/3.4, scale=2)


plt <- tmp.metadata.paired |> 
  dplyr::select(`MGMT at Rec.`, `C1/col signature at Rec.`) |> 
  table() |> 
  as.data.frame() |> 
  dplyr::rename(`MGMT at Rec.` = `MGMT.at.Rec.`) |> 
  dplyr::rename(`C1/col signature at Rec.` = `C1.col.signature.at.Rec.`)
ggplot(plt, aes(x=`C1/col signature at Rec.`,y=`Freq`,fill=`MGMT at Rec.`)) +
  geom_bar(stat="identity") +
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
  )


ggsave("output/figures/2022_figure_S13h.pdf", width=8.3 / 4,height=8.3/4, scale=2)


plot(tmp.metadata.paired |>  dplyr::select(C0.fuzzy.signature, `C1.col.signature`) |> table())
plot(tmp.metadata.paired |>  dplyr::select(C3.olig.signature, `C4.neu.signature`) |> table())
plot(tmp.metadata.paired |>  dplyr::select(`C1.col.signature`,`MGMT at Rec.`) |> table())





#tmp.metadata.paired |>  dplyr::select(C0.fuzzy.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C2.endo.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C3.olig.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C4.neu.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C4.neu.signature, C3.olig.signature) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C4.neu.signature, C3.olig.signature) |> table() |> fisher.test()



tmp.metadata.paired |>
  dplyr::filter(!is.na(`MGMT at Rec.`)) |> 
  dplyr::select(`MGMT at Rec.`, `C1/col signature at Rec.`) |>  
  table()


tmp.metadata.paired |>
  dplyr::filter(!is.na(`MGMT at Rec.`)) |> 
  dplyr::select(`MGMT at Rec.`, `C1/col signature at Rec.`) |>  
  table() |>  
  fisher.test()


#ggsave("output/figures/R2_survival_signatures.pdf",height=5.5,width=12)






## R1 ----



surv_object <- survival::Surv(time = tmp.metadata.paired$survivalDays, event=tmp.metadata.paired$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.prim` , data = tmp.metadata.paired)
survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      palette = c(
                        'C1/col signature [Prim.]: high'=alpha('#CB75A4',0.7),
                        'C1/col signature [Prim.]: low'=alpha('#009E74',0.7)
                      ),
                      legend.labs=c('C1.col.signature=high'='C1/col signature [Prim.]: high',
                                    'C1.col.signature=low'='C1/col signature [Prim.]: low'),
                      xlab="Survival time from primary")




## Time to progress ----

# RS Question: those patients that acquired high COL signature, did they have shorter time to recurrence?


surv_object <- survival::Surv(time = tmp.metadata.paired$daysToProgression, event=tmp.metadata.paired$progression.event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature.rec` , data = tmp.metadata.paired)
survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      palette = c(
                        'C1/col signature [Rec.]: high'=alpha('#CB75A4',0.7),
                        'C1/col signature [Rec.]: low'=alpha('#009E74',0.7)
                      ),
                      legend.labs=c('C1.col.signature=high'='C1/col signature [Rec.]: high',
                                    'C1.col.signature=low'='C1/col signature [Rec.]: low'),
                      xlab="Time to progression/recurrence from primary")




fit.cox <- survival::coxph(surv_object ~
                             `MGMT at Rec.` +
                             `C1/col signature at Rec.` ,
                           
                           data = tmp.metadata.paired)
survminer::ggforest(fit.cox, data = tmp.metadata.paired)






