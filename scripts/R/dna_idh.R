#!/usr/bin/env R

dna_idh <- read.delim("data/DNA/GBM_allSamples_PassVariants_1736_GBM_TERT_05Jan2018_flags_VAFs.txt",stringsAsFactors=F,sep="\t")
dna_idh$Normal <- NULL
dna_idh$Sample <- as.factor(dna_idh$Sample)
dna_idh <- dna_idh[ dna_idh$Gene %in% c("IDH1", "IDH2"),]

is_idh1 <- dna_idh$Gene == "IDH1" &  gsub("^([A-Z][0-9]+).+$","\\1",dna_idh$Protein) == "R132"
is_idh2 <- dna_idh$Gene == "IDH2" &  gsub("^([A-Z][0-9]+).+$","\\1",dna_idh$Protein) %in% c("R140", "R172")
#R140 and 20% at R172. 

dna_idh <- dna_idh[is_idh1 | is_idh2 ,]

rm(is_idh1, is_idh2)

is_idh <- as.character(dna_idh$Sample)


tmp <- read.delim("data/DNA/sample codes sanger gsam.txt",stringsAsFactors=FALSE)
tmp <- tmp[,match(c("donor_ID", "PD_ID"), colnames(tmp))]


match(is_idh  , tmp$PD_ID)
is_idh <- tmp[match(is_idh  , tmp$PD_ID),]$donor_ID

rm(tmp, dna_idh)


