#!/usr/bin/env R


# ---- initialization setup ----

wd <- "/home/youri/projects/gsam"
wd <- "/home/yhoogstrate/projects/gsam"

setwd(wd)

rm(wd)


# ---- load: libs ----

library(ggplot2)
library(pheatmap)



# ---- load: functions ----

all_mutations <- read.delim("data/DNA/GBM_allSamples_PassVariants_1736_GBM_TERT_05Jan2018_flags_VAFs.txt", stringsAsFactors = F)
all_mutations <- all_mutations[ all_mutations$Gene %in% c("IDH1", "IDH2"),]

is_idh1 <- all_mutations$Gene == "IDH1" &  gsub("^([A-Z][0-9]+).+$","\\1",all_mutations$Protein) == "R132"
is_idh2 <- dna_idh$Gene == "IDH2" &  gsub("^([A-Z][0-9]+).+$","\\1",dna_idh$Protein) %in% c("R140", "R172")

