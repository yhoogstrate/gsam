#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load libs ----

library(ggplot2)

# ---- load data ----

source('scripts/R/gsam_metadata.R')

# ---- DV200 (all) ----


plt <- gsam.rna.metadata[order(gsam.rna.metadata$giga.DV200),]

plt$blacklist <- as.character(plt$blacklist.pca)
plt$blacklist[plt$blacklist.heavy.dna.contamination] <- paste0(plt$blacklist[plt$blacklist.heavy.dna.contamination], ", DNA-Cont.")


ggplot(plt, aes(x = giga.DV200,y =  pct.rRNA.by.chrUn.gl000220, col=blacklist)) +
  geom_point()

ggplot(plt, aes(x = giga.DV200,y =  plt$`giga.ng/µL.Bio-.analyzer`, col=blacklist)) +
  geom_point()


ggplot(plt, aes(y = `giga.DV200`, x=order(order(`giga.RNA.quantity.à.atteindre.(selon.DV200)`)), group=blacklist, col=blacklist)) +
  geom_point(cex=0.6)



ggplot(plt, aes(y = `giga.DV200`, x=`giga.nM.on.qPCR`, group=blacklist, col=blacklist)) +
  geom_point(cex=0.6) + 
  facet_grid(. ~ giga.Plate)  + 
  labs(x="nM on qPCR",y = "DV200", faced="plate")






