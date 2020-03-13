#!/usr/bin/env R

"
Goal of this script is to determine the tumour percentage somehow

This is technically possible with CNVKit + bam files + loads of other analysis, but may most likely take weeks



https://cnvkit.readthedocs.io/en/stable/heterogeneity.html

Requires VCF Files, which require running mutect, which require the bam files, of which ~25% are missing
"

setwd("~/projects/gsam")

# ---- load libs ----

library(ggplot2)
library(limma)


# ---- load functions ----


# ---- load data----

#source("scripts/R/gsam_metadata.R") - auto loaded in script below
source("scripts/R/cnv_matrix.R")


# ---- one sample example ----
# 5,6 = pair without further gains and losses

i = 1
j = 2

off <- median(cnv_matrix [, i] - cnv_matrix[,j])
s1 <- cnv_matrix [, i ]
s2 <- (cnv_matrix[,j] + off)

plot(c(-2.2, 2.2),c(-2.2, 2.2), type="n",xlab=colnames(cnv_matrix)[i],ylab=colnames(cnv_matrix)[j])
points(s1, s2, pch=19,cex=0.2,col=rgb(0,0,0,0.05))

#ss1 = s1[abs(s1) < 0.25]
#ss2 = s2[abs(s1) < 0.25]
#d <- lm(ss2 ~ ss1)
#abline(d,col="red")


#st1 = s1[abs(s1) >= 0.34 & abs(s1) <= 1]
#st2 = s2[abs(s1) >= 0.34 & abs(s1) <= 1]
#d <- lm(st2 ~ st1)
#abline(d,col="red")


#data <- data.frame(s1=s1,s2=s2)
#d <- loess(s2 ~ s1, data=data) 

#j <- order(data$s1)
#lines(data$s1[j],d$fitted[j],col="red",lwd=1)


d <- s2-s1
d <- d[abs(d) < 1.6]
hist(d,breaks=1000)



