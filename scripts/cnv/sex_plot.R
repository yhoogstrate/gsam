#!/usr/bin/env R


setwd("~/projects/gsam")

# ---- load libs ----

library(ggplot2)
library(limma)


# ---- load functions ----


# ---- load data----

#source("scripts/R/gsam_metadata.R") - auto loaded in script below
source("scripts/R/cnv_matrix.R")

# there is a difference in size of results files of CNVKit
# they either have 27493 lines (all in ExtraSequencing_CN_data_29Mar2018/)
# or have 27287 lines (195/374) = 50%


job_gg_theme <- theme(
  text = element_text(family = 'Helvetica'),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.title.y = element_text(angle=90,vjust =2),
  axis.title.x = element_text(vjust = -0.2),
  axis.text = element_text(),
  axis.line = element_line(colour="black"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
  panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted')
)


# ---- Sex Plot Batch-1 ----

sel <- T# gsam.cnv.metadata$batch == "b1"
cnv_matrix.b1 <- cnv_matrix[,sel]
cnv_matrix.b1.X <- cnv_matrix.b1[cnv_matrix_genes$chr %in% c("chrX"),]
cnv_matrix.b1.Y <- cnv_matrix.b1[cnv_matrix_genes$chr %in% c("chrY"),]
rm(cnv_matrix.b1)

cnv_matrix.b1.X <- colSums(cnv_matrix.b1.X)
cnv_matrix.b1.Y <- colSums(cnv_matrix.b1.Y)

#plot(cnv_matrix.b1.X, cnv_matrix.b1.Y)

tmp <- data.frame(
  x = cnv_matrix.b1.X,
  y = cnv_matrix.b1.Y,
  batch = gsam.cnv.metadata[sel,]$batch ,
  gender = as.factor(gsam.cnv.metadata[sel,]$gender) , 
  sid = paste0(colnames(cnv_matrix)[sel],".",gsam.cnv.metadata[sel,]$batch )
)
rownames(tmp ) <- colnames(cnv_matrix[sel])


mislabels <- subset(tmp, 
              (y > -750) & (gender == "Female")
                |
              (y <= -750) & (gender == "Male"))


gg <- ggplot(tmp, aes(x=x, y=y, label=sid)) + 
  geom_point(aes(col=gender, shape=batch)) +
  geom_encircle(aes(x=x, y=y), 
                data=subset(tmp, y > -750),
                color=1, 
                size=1, 
                expand=0.08) +
  geom_encircle(aes(x=x, y=y), 
              data=subset(tmp, y <= -750),
              color=1, 
              size=1, 
              expand=0.08) +
  geom_text_repel(
    nudge_y       = 200 - mislabels$y,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x",
    data = mislabels,
    ) +
  job_gg_theme +
  ylim(NA, 150)

plot(gg)







