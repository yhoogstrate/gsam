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



# ---- Sex Plot Batch-1 ----

sel <- T# gsam.cnv.metadata$batch == "b1"
cnv_matrix.b1 <- cnv_matrix[,sel]
cnv_matrix.b1.X <- cnv_matrix.b1[cnv_matrix_genes$chr %in% c("chrX"),]
cnv_matrix.b1.Y <- cnv_matrix.b1[cnv_matrix_genes$chr %in% c("chrY"),]
rm(cnv_matrix.b1)

cnv_matrix.b1.X <- colSums(cnv_matrix.b1.X)
cnv_matrix.b1.Y <- colSums(cnv_matrix.b1.Y)

plot(cnv_matrix.b1.X, cnv_matrix.b1.Y)

tmp <- data.frame(
  x = cnv_matrix.b1.X,
  y = cnv_matrix.b1.Y,
  batch = gsam.cnv.metadata[sel,]$batch ,
  gender = as.factor(gsam.cnv.metadata[sel,]$gender)
)


tmp2 <- tmp[tmp$y > -750, ]


gg <- ggplot(tmp, aes(x=x, y=y)) + 
  geom_point(aes(col=gender, shape=batch)) +
  geom_encircle(aes(x=x, y=y), 
                data=tmp2, 
                color=1, 
                size=2, 
                expand=0.08) +
  geom_encircle(aes(x=x, y=y), 
              data=tmp2, 
              color=1, 
              size=2, 
              expand=0.08) 
plot(gg)







