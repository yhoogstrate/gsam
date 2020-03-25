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
ggsave("output/figures/qc/giga_stats_dv200x_nm-qpcr_plate.pdf")



ggplot(plt, aes(y = `giga.DV200`, x=`giga.nM.on.qPCR`, group=blacklist, col=blacklist)) +
  geom_point(cex=0.6) + 
  labs(x="nM on qPCR",y = "DV200", faced="plate")
ggsave("output/figures/qc/giga_stats_dv200x_nm-qpcr.pdf")



ggplot(plt, aes(y = `giga.DV200`, x=`giga.ng/µL.Bio-.analyzer`, group=blacklist, col=blacklist)) +
  geom_point(cex=0.6) + 
  labs(x="ng/µL.Bio-.analyzer",y = "DV200", faced="plate")
ggsave("output/figures/qc/giga_stats_dv200x_ng_ul_bioanalyser_plate.pdf")



ggplot(plt, aes(y = `giga.DV200`, x=`giga.ng/µL.Bio-.analyzer`, group=blacklist, col=blacklist)) +
  geom_point(cex=0.6) + 
  facet_grid(. ~ giga.Plate)  + 
  labs(x="ng/µL.Bio-.analyzer",y = "DV200", faced="plate")
ggsave("output/figures/qc/giga_stats_dv200x_ng_ul_bioanalyser.pdf")



gsam.rna.metadata$wt.reads.v3 <- NULL
gsam.rna.metadata$vIII.reads.v3 <- NULL
gsam.rna.metadata$vIII.percentage <- NULL
gsam.rna.metadata$qPCR.ct.EGFR.wt <- NULL
gsam.rna.metadata$qPCR.ct.EGFR.vIII <- NULL
gsam.rna.metadata$qPCR.percent.EGFR.vIII <- NULL
gsam.rna.metadata$blacklist.gender.mislabeling <- NULL

write.csv(gsam.rna.metadata,"/tmp/g-sam.rna-seq.metadata.csv")
write.csv2(gsam.rna.metadata,"/tmp/g-sam.rna-seq.metadata.csv2")
write.table(gsam.rna.metadata,"/tmp/g-sam.rna-seq.metadata.tab.txt")


