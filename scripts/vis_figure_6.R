#!/usr/bin/env R 

# load libs ----


# library(tidyverse)
# library(reshape2)
# library(RColorBrewer)
# library(cowplot)
# library(patchwork)
# library(survival)
# library(survminer)
# 


#  load data ----


source('scripts/R/palette.R')

source('scripts/R/gsam_rna-seq_expression.R') # recursively calls metadata
source('scripts/R/gsam_metadata.R') # recursively calls metadata



delta <- function(x) {
  n <- length(x)
  
  return ( x[2:n] - x[1:(n-1)] )
}


## prep data table ----


tmp.metadata <- gsam.rna.metadata |> 
  dplyr::filter(.data$blacklist.pca == F)  |> 
  dplyr::filter(.data$pat.with.IDH == F) |> 
  dplyr::filter(.data$sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
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

  dplyr::mutate( delta.rna.signature.C0.fuzzy.2022 = rna.signature.C0.fuzzy.2022_recurrence - rna.signature.C0.fuzzy.2022_primary) |> 
  dplyr::mutate( delta.rna.signature.C1.collagen.2022 = rna.signature.C1.collagen.2022_recurrence - rna.signature.C1.collagen.2022_primary) |> 
  dplyr::mutate( delta.rna.signature.C2.endothelial.2022 = rna.signature.C2.endothelial.2022_recurrence - rna.signature.C2.endothelial.2022_primary) |> 
  dplyr::mutate( delta.rna.signature.C3.oligodendrocyte.2022 = rna.signature.C3.oligodendrocyte.2022_recurrence - rna.signature.C3.oligodendrocyte.2022_primary) |> 
  dplyr::mutate( delta.rna.signature.C4.neuron.2022 = rna.signature.C4.neuron.2022_recurrence - rna.signature.C4.neuron.2022_primary) |> 

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
        
        cnStatusCDKN2ABs,
        AXIN2,APC,JAK2,
        RB1,MSH2,BRCA1,
        BRCA2,ATM,SETD2,ARID2,KMT2C,KMT2D,
        NF1,ERBB3,
        EGFR,TP53BP1,TP53,PIK3R1,PIK3CA,TSC2,
        
        cnStatusRB1s,
        cnStatusNF1s,
        cnStatusCDK4s,
        cnStatusMDM2s,
        
        SETD2,PDGFRA
        
        ) |> 
      dplyr::rename(event = status)
      , by=c('pid'='studyID'), suffix=c('','')
  ) |> 
  dplyr::mutate(event = ifelse(.data$event == "Deceased",1,0)) |> 
  dplyr::mutate(deceased = dplyr::recode(event, "1" = "Yes", "0" = "No" )) |> 
  dplyr::mutate(rank = order(order(delta.rna.signature.C1.collagen.2022, delta.rna.signature.C1.collagen.2022, pid))) |> 
  dplyr::mutate(em.pc.status = ifelse(.data$`rna.signature.C1.collagen.2022_recurrence` > .data$`rna.signature.C1.collagen.2022_primary`, "increase", "decrease")) |> 
  dplyr::mutate(`MGMT meth` = dplyr::recode( mgmtStability, "Stable methylated" = "Stable", "Stable unmethylated" = "Wildtype" ), mgmtStability = NULL ) |> 
  dplyr::mutate(`treatment: Beva` = ifelse(bevacizumab.before.recurrence, "Yes","No"), bevacizumab.before.recurrence=NULL) |> 
  dplyr::rename(`treatment: TMZ` = treatedWithTMZ) |> 
  dplyr::rename(`treatment: RT` = treatedWithRT) |> 
  dplyr::rename(`GITS subtype R1` = .data$GITS.150.svm.2022.subtype_primary) |> 
  dplyr::rename(`GITS subtype R2` = .data$GITS.150.svm.2022.subtype_recurrence) |> 
  dplyr::mutate(`Age above 50` = ifelse(age > 50,"Yes","No")) |> 
  dplyr::mutate(gender = as.character(gender)) |> 
  dplyr::rename(Sex = gender) |> 
  dplyr::mutate(`KPS 70 or above` = ifelse(is.na(performanceAtSecondSurgery) | performanceAtSecondSurgery >= 70, "Yes","No")) |> 
  dplyr::mutate(performanceAtSecondSurgery = as.factor(as.character(performanceAtSecondSurgery)))





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




## determine cut-off C0 / fuz ----
# ### from all R2
# 
# d <- tmp.metadata |> 
#   dplyr::filter(resection == "recurrence") |> 
#   dplyr::arrange(rna.signature.C0.fuzzy.2022) |> 
#   dplyr::pull(rna.signature.C0.fuzzy.2022)
# 
# 
# plot(delta(d))
# 
# 
# k <- 128
# c0.fuz.cutoff.s <- d[(k-1):k] |> 
#   sum()  |> 
#   (function(x){return(x/2)})()
# 
# 
# plot(d)
# abline(h=c0.fuz.cutoff.s)


### from paired ----

d <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C0.fuzzy.2022_recurrence) |>
  dplyr::pull(rna.signature.C0.fuzzy.2022_recurrence)


plot(delta(d))

rm(k)
#c0.fuz.cutoff.p <- d[(k-1):k] |> sum() |> (function(x){return(x/2)})()
c0.fuz.cutoff.p <- d |> median()


plt.1 <- data.frame(y = d) |> 
  dplyr::rename(`C0/fuzzy signature at Rec.` =y) |> 
  dplyr::mutate(x = order(order(1:n())))

plt.2 <- data.frame(y = delta(d)) |> 
  dplyr::rename(`delta C0/fuzzy signature at Rec.` = y ) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=`C0/fuzzy signature at Rec.`, x=x)) +
  geom_hline(yintercept=c0.fuz.cutoff.p, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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

p2 <- ggplot(plt.2, aes(y=`delta C0/fuzzy signature at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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


p1 / p2




## determine cut-off C1 / col ----
### from all R2
# 
# d <- tmp.metadata |> 
#   dplyr::filter(resection == "recurrence") |> 
#   dplyr::arrange(rna.signature.C1.collagen.2022) |> 
#   dplyr::pull(rna.signature.C1.collagen.2022)
# 
# 
# plot(delta(d))
# 
# 
# k <- 91
# c1.em.cutoff.s <- d[(k-1):k] |> 
#   sum()  |> 
#   (function(x){return(x/2)})()
# 
# 
# plot(d)
# abline(h=c1.em.cutoff.s, col="red")



### from paired ----

d <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C1.collagen.2022_recurrence) |>
  dplyr::pull(rna.signature.C1.collagen.2022_recurrence)



k <- 82+1
c1.em.cutoff.p <- d[(k-1):k] |> sum() |> (function(x){return(x/2)})()


plt.1 <- data.frame(y = d) |> 
  dplyr::rename(`C1/col signature at Rec.` =y) |> 
  dplyr::mutate(x = order(order(1:n())))

plt.2 <- data.frame(y = delta(d)) |> 
  dplyr::rename(`delta C1/col signature at Rec.` = y ) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=`C1/col signature at Rec.`, x=x)) +
  geom_hline(yintercept=c1.em.cutoff.p, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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

p2 <- ggplot(plt.2, aes(y=`delta C1/col signature at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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


p1 / p2




## determine cut-off C2 / endo ----
# ### from all R2 
# 
# d <- tmp.metadata |> 
#   dplyr::filter(resection == "recurrence") |> 
#   dplyr::arrange(rna.signature.C2.endothelial.2022) |> 
#   dplyr::pull(rna.signature.C2.endothelial.2022)
# 
# 
# plot(delta(d))
# 
# 
# # k <- 82
# # c2.endo.cutoff.s <- d[(k-1):k] |> 
# #   sum()  |> 
# #   (function(x){return(x/2)})()
# c2.endo.cutoff.s <- d[(k-1):k] |> 
#   median()
# 
# 
# plot(d)
# abline(h=c2.endo.cutoff.s)


### from paired ----


d <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C2.endothelial.2022_recurrence) |>
  dplyr::pull(rna.signature.C2.endothelial.2022_recurrence)


k <- 52
c2.endo.cutoff.p <- d[(k-1):k] |> sum() |> (function(x){return(x/2)})()


plt.1 <- data.frame(y = d) |> 
  dplyr::rename(`C2/endo signature at Rec.` =y) |> 
  dplyr::mutate(x = order(order(1:n())))

plt.2 <- data.frame(y = delta(d)) |> 
  dplyr::rename(`delta C2/endo signature at Rec.` = y ) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=`C2/endo signature at Rec.`, x=x)) +
  geom_hline(yintercept=c2.endo.cutoff.p, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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

p2 <- ggplot(plt.2, aes(y=`delta C2/endo signature at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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


p1 / p2



## determine cut-off C3 / olig ----
# ### from all R2
# 
# d <- tmp.metadata |> 
#   dplyr::filter(resection == "recurrence") |> 
#   dplyr::arrange(rna.signature.C3.oligodendrocyte.2022) |> 
#   dplyr::pull(rna.signature.C3.oligodendrocyte.2022)
# 
# 
# plot(delta(d))
# 
# 
# k <- 82
# c3.olig.em.cutoff.s <- d[(k-1):k] |> 
#   sum()  |> 
#   (function(x){return(x/2)})()
# 
# 
# plot(d)
# abline(h=c3.olig.cutoff.s)
# abline(v=k-1,col="red")
# abline(v=k,col="red")


### from paired ----

d <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C3.oligodendrocyte.2022_recurrence) |>
  dplyr::pull(rna.signature.C3.oligodendrocyte.2022_recurrence)


#k <- 76
#c3.olig.cutoff.p <- d[(k-1):k] |> sum() |> (function(x){return(x/2)})()
rm(k)
c3.olig.cutoff.p <- d |>  median()


plt.1 <- data.frame(y = d) |> 
  dplyr::rename(`C3/olig signature at Rec.` =y) |> 
  dplyr::mutate(x = order(order(1:n())))

plt.2 <- data.frame(y = delta(d)) |> 
  dplyr::rename(`delta C3/olig signature at Rec.` = y ) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=`C3/olig signature at Rec.`, x=x)) +
  geom_hline(yintercept=c3.olig.cutoff.p, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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

p2 <- ggplot(plt.2, aes(y=`delta C3/olig signature at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  #geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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


p1 / p2




## determine cut-off C4 / neur ----
# ### from all R2
# 
# d <- tmp.metadata |> 
#   dplyr::filter(resection == "recurrence") |> 
#   dplyr::arrange(rna.signature.C4.neuron.2022) |> 
#   dplyr::pull(rna.signature.C4.neuron.2022)
# 
# 
# plot(delta(d))
# 
# 
# k <- 90
# c4.neur.cutoff.s <- d[(k-1):k] |> 
#   sum()  |> 
#   (function(x){return(x/2)})()
# 
# 
# plot(d)
# abline(h=c4.neur.cutoff.s)


### from paired ----


d <- tmp.metadata.paired |>
  dplyr::arrange(rna.signature.C4.neuron.2022_recurrence) |>
  dplyr::pull(rna.signature.C4.neuron.2022_recurrence)


plot(delta(d))

k <- 83
c4.neur.cutoff.p <- d[(k-1):k] |> sum() |> (function(x){return(x/2)})()



plt.1 <- data.frame(y = d) |> 
  dplyr::rename(`C4/neu signature at Rec.` =y) |> 
  dplyr::mutate(x = order(order(1:n())))

plt.2 <- data.frame(y = delta(d)) |> 
  dplyr::rename(`delta C4/neu signature at Rec.` = y ) |> 
  dplyr::mutate(x = order(order(1:n())) + 0.5)


p1 <- ggplot(plt.1, aes(y=`C4/neu signature at Rec.`, x=x)) +
  geom_hline(yintercept=c4.neur.cutoff.p, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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

p2 <- ggplot(plt.2, aes(y=`delta C4/neu signature at Rec.`, x=x)) +
  geom_hline(yintercept= 0, lwd=1.5, color = "gray50") +
  #geom_vline(xintercept=k - 0.5, lwd=1.5, color = "red") +
  geom_line() +
  geom_point() +
  xlim(1,nrow(plt.1)) +
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


p1 / p2






## plot(s) ----

### panel A: svvl ----

plt <- tmp.metadata.paired |> 
  dplyr::mutate(ttp = -1 * (survivalDays - survivalFromSecondSurgeryDays)) |> 
  dplyr::rename(tfp = survivalFromSecondSurgeryDays) |> 
  dplyr::select(pid, rank, event, ttp, tfp) |> 
  tidyr::pivot_longer(cols=c(ttp,tfp))


ggplot(plt, aes(x=reorder(pid, rank),y=value, group=pid)) +
  geom_line(lwd=2, fill="gray60") +
  geom_hline(yintercept=0) +
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
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.45, size=5.5)  ) +
  labs(x = NULL, y="Survival")


### panel B: signature arrows ----


plt <- tmp.metadata |>  # only resections from complete pairs
  dplyr::group_by(pid) |> 
  dplyr::filter(n() == 2) |> 
  dplyr::ungroup()


ggplot(plt, aes(x = reorder(pid, rank), y=`rna.signature.C1.collagen.2022`, col=em.pc.status, group=pid))  +
  geom_hline(yintercept=2.5, lty=1, color = "#FFFFFF44",lwd=3) +
  geom_point(data = subset(plt, resection == "primary"), pch=19, cex=1.2, alpha=0.25) +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.75 )  + 
  labs(x=NULL, col=NULL, y="ECM signature",
       caption="G-SAM: n=122 pairs"
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
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  ) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c(
    "increase"="#bb5f6c", # rood #bb5f6c
    "decrease"="#79b1b1" # lichtlauw #79b1b
  )) +
  geom_hline(yintercept=c1.em.cutoff.p, lty=3, color = "red",lwd=0.5)


### panel C: squares ----

#'@todo x-check @ fig2 fgh -- mgmt


plt <- tmp.metadata.paired |>
  dplyr::select(-em.pc.status, -sid_primary, -sid_recurrence, -survivalDays, -survivalFromSecondSurgeryDays, -event, -age, -mgmt.status.rec) |> 
  dplyr::select(-contains(".signature")) |> 
  reshape2::melt(id.vars=c("pid", "rank")) |> 
  dplyr::mutate(value = gsub('^Gained$','Gained/increased/female',value)) %>% 
  dplyr::mutate(value = gsub('^Lost$','Lost/decreased/male',value)) %>% 
  dplyr::mutate(value = gsub('^Female$','Gained/increased/female',value)) %>% 
  dplyr::mutate(value = gsub('^Male$','Lost/decreased/male',value)) %>% 
  dplyr::mutate(value = gsub('^Yes|Stable$','Yes/Stable',value)) %>% 
  dplyr::mutate(value = gsub('^No|Wildtype$','No/Wildtype',value)) %>% 
  dplyr::mutate(panel = case_when(
    grepl("deceased|age|sex|kps", variable) ~ "A",
    grepl("treatment", variable) ~ "B",
    grepl("subtype", variable) ~ "C",
    variable %in% c("HM",  "MGMT meth") ~ "E",
    T ~ "D"
  ))


ggplot(plt, aes(x = reorder(pid, rank), y = variable, fill = value)) +
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
                                "No/Wildtype"="white",
                                
                                "Gained/increased/female"="#bb5f6c", # rood #bb5f6c
                                "Lost/decreased/male"="#79b1b1", # lichtlauw #79b1b1
                                "Yes/Stable"="#2e415e", # donker blauw #2e415e
                                #"Yes" = "#2e415e", # zelfde als stable #2e415e
                                
                                "NA"="grey")) + 
  #coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = 'bottom')







### stats ? ----





## re-newed SVVL analysis [low/high] ----



tmp.metadata.paired <- tmp.metadata.paired |> 
  dplyr::mutate(`C0.fuzzy.signature` = ifelse(rna.signature.C0.fuzzy.2022_recurrence >           c0.fuz.cutoff.p,  "high", "low")) |> 
  dplyr::mutate(`C1.col.signature`   = ifelse(rna.signature.C1.collagen.2022_recurrence >        c1.em.cutoff.p,   "high", "low")) |> 
  dplyr::mutate(`C2.endo.signature`  = ifelse(rna.signature.C2.endothelial.2022_recurrence >     c2.endo.cutoff.p ,"high", "low")) |> 
  dplyr::mutate(`C3.olig.signature`  = ifelse(rna.signature.C3.oligodendrocyte.2022_recurrence > c3.olig.cutoff.p ,"high", "low")) |> 
  dplyr::mutate(`C4.neu.signature`   = ifelse(rna.signature.C4.neuron.2022_recurrence >          c4.neur.cutoff.p ,"high", "low")) |> 

  dplyr::mutate(`C0.fuzzy.signature` = factor(`C0.fuzzy.signature`, levels=c('low','high')))  |> 
  dplyr::mutate(`C1.col.signature` = factor(`C1.col.signature`, levels=c('low','high'))) |> 
  dplyr::mutate(`C2.endo.signature` = factor(`C2.endo.signature`, levels=c('low','high'))) |> 
  dplyr::mutate(`C3.olig.signature` = factor(`C3.olig.signature`, levels=c('low','high'))) |> 
  dplyr::mutate(`C4.neu.signature` = factor(`C4.neu.signature`, levels=c('low','high'))) |> 
  
  dplyr::mutate(`C0/fuzzy signature at Rec.` = `C0.fuzzy.signature`)  |> 
  dplyr::mutate(`C1/col signature at Rec.` = `C1.col.signature`) |> 
  dplyr::mutate(`C2/endo signature at Rec.` = `C2.endo.signature`) |> 
  dplyr::mutate(`C3/olig signature at Rec.` = `C3.olig.signature`) |> 
  dplyr::mutate(`C4/neu signature at Rec.` = `C4.neu.signature`) |> 
  
  dplyr::mutate(`MGMT at Rec.` = case_when(
    is.na(`MGMT meth`) ~ as.character(NA),
    `MGMT meth` %in% c("Gained", "Stable") ~ "Methylated",
    T ~ "Unmethylated"
  )) 



### panel E: KM collagen ----

surv_object <- survival::Surv(time = tmp.metadata.paired$survivalFromSecondSurgeryDays, event=tmp.metadata.paired$event)
fit1 <- survival::survfit(surv_object ~  `C1.col.signature` , data = tmp.metadata.paired)
survminer::ggsurvplot(fit1, data = tmp.metadata.paired, pval = TRUE, risk.table=T, tables.y.text = FALSE,
                      palette = c(
                        'C1/col signature: high'=alpha('#CB75A4',0.7),
                        'C1/col signature: low'=alpha('#009E74',0.7)
                      ),
                      legend.labs=c('C1.col.signature=high'='C1/col signature: high',
                                    'C1.col.signature=low'='C1/col signature: low'),
                      xlab="Survival time from recurrence")



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




fit.cox <- survival::coxph(surv_object ~
                             `Age above 50` +
                             Sex +
                           `KPS 70 or above` +
                             `treatment: Beva` +
                             `treatment: TMZ` +
#                             mgmt.status.rec +

#`C0/fuzzy signature at Rec.` +
`C1/col signature at Rec.` +
`C2/endo signature at Rec.` +
`C3/olig signature at Rec.` +
`C4/neu signature at Rec.`
,
                 data = tmp.metadata.paired)
survminer::ggforest(fit.cox
                    
                    #legend.labs=c('Sex'='sex',
#                                  'Age > 50'='age.above.50'),
                    )


#tmp.metadata.paired |>  dplyr::select(C0.fuzzy.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C2.endo.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C3.olig.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C4.neu.signature, `C1.col.signature`) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C4.neu.signature, C3.olig.signature) |> table() |> fisher.test()
#tmp.metadata.paired |>  dplyr::select(C4.neu.signature, C3.olig.signature) |> table() |> fisher.test()



fit.cox <- survival::coxph(surv_object ~
                             #age.above.50 +
                             #sex +
                             #kps.70.or.above +
                             `MGMT at Rec.` +
                             `C1/col signature at Rec.`
                             #`C2.endo.signature` +
                             #`C3.olig.signature` +
                             #`C4.neu.signature` +
                             #`C0.fuzzy.signature`
                           ,
                           data = tmp.metadata.paired)
survminer::ggforest(fit.cox)


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




plot( -components.paired$survivalFromSecondSurgeryDays , components.paired$C5.signature.R2)


plot( components.paired$C6.ECM.signature.R2 , components.paired$C5.signature.R2)
plot( components.paired$C6.ECM.signature.R2 , components.paired$C4.signature.R2)


plot( -components.paired$survivalFromSecondSurgeryDays , components.paired$C5.signature.R2)




