#!/usr/bin/env R 

# load libs ----


# library(tidyverse)
# library(reshape2)
# library(RColorBrewer)
# library(cowplot)
# library(patchwork)
# library(survival)
# library(survminer)


#  load data ----


source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')
source('scripts/R/palette.R')

#source('scripts/R/gsam_rna-seq_expression.R') # recursively calls metadata
source('scripts/R/gsam_metadata.R') # recursively calls metadata


# longitudinal components overview ----


## neuron ----


plt <- gsam.rna.metadata |> 
  dplyr::select(contains(".component") | contains(".component") | 'sid' | 'pid' | resection ) |> 
  dplyr::filter(!is.na(neuron.component)) |> 
  dplyr::arrange(pid, resection, sid) |> 
  dplyr::group_by(pid) |> 
  dplyr::mutate(paired.status = factor(case_when(
    dplyr::n() > 1 ~ "paired",
    resection == "r1" ~ "R1 / primary",
    resection == "r2" ~ "R2 / recurrence",
    T ~ "??"),levels=c("paired","R1 / primary","R2 / recurrence"))) |> 
  dplyr::ungroup() |> 
  dplyr::left_join(
    gsam.rna.metadata |>
      dplyr::filter(!is.na(neuron.component)) |> 
      dplyr::select(pid,resection,neuron.component) |> 
      tidyr::pivot_wider(names_from = resection, values_from = neuron.component) |> 
      dplyr::mutate(neuron.component.delta = r2 - r1) |> 
      dplyr::mutate(neuron.component.delta = ifelse(is.na(neuron.component.delta),0, neuron.component.delta)) |> 
      dplyr::mutate(r1 = NULL) |> 
      dplyr::mutate(r2 = NULL),
    by=c('pid'='pid')
  ) |>
  dplyr::mutate(neuron.component.delta = -1* neuron.component.delta) |> 
  dplyr::mutate(neuron.component = -1* neuron.component) |> 
  dplyr::mutate(order = order(order(as.numeric(paired.status), neuron.component.delta, neuron.component))) 
  
  

ggplot(plt, aes(x=reorder(pid, order), y=neuron.component)) +
  facet_grid(cols = vars(paired.status), scales = "free",space="free_x") +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.6 )  +
  geom_point(data = subset(plt, paired.status != "paired" | (paired.status == "paired" & resection == "r1"))) +
  theme_bw() +
  theme(    axis.text.x = element_text(angle = 90, vjust = 0.45, size=5.5)  ) +
  labs(y = "Neuron RNA signature / contribution", x=NULL)


ggsave("/tmp/neuron.pdf", width=11,height=4.5)



## ECM ----





plt <- gsam.rna.metadata |> 
  dplyr::select(contains(".component") | contains(".component") | 'sid' | 'pid' | resection ) |> 
  dplyr::filter(!is.na(extracellular.matrix.component)) |> 
  dplyr::arrange(pid, resection, sid) |> 
  dplyr::group_by(pid) |> 
  dplyr::mutate(paired.status = factor(case_when(
    dplyr::n() > 1 ~ "paired",
    resection == "r1" ~ "R1 / primary",
    resection == "r2" ~ "R2 / recurrence",
    T ~ "??"),levels=c("paired","R1 / primary","R2 / recurrence"))) |> 
  dplyr::ungroup() |> 
  dplyr::left_join(
    gsam.rna.metadata |>
      dplyr::filter(!is.na(extracellular.matrix.component)) |> 
      dplyr::select(pid,resection,extracellular.matrix.component) |> 
      tidyr::pivot_wider(names_from = resection, values_from = extracellular.matrix.component) |> 
      dplyr::mutate(extracellular.matrix.component.delta = r2 - r1) |> 
      dplyr::mutate(r1 = NULL) |> 
      dplyr::mutate(r2 = NULL),
    by=c('pid'='pid')
  ) |>
  dplyr::mutate(order = order(order(as.numeric(paired.status), extracellular.matrix.component.delta, extracellular.matrix.component))) 



ggplot(plt, aes(x=reorder(pid, order), y=extracellular.matrix.component)) +
  facet_grid(cols = vars(paired.status), scales = "free",space="free_x") +
  geom_path(arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.125, "inches")) , alpha = 0.6 )  +
  geom_point(data = subset(plt, paired.status != "paired" | (paired.status == "paired" & resection == "r1"))) +
  theme_bw() +
  theme(    axis.text.x = element_text(angle = 90, vjust = 0.45, size=5.5)  ) +
  labs(y = "ECM/collagen signature / contribution", x=NULL)


ggsave("/tmp/ECM.pdf", width=11,height=4.5)




