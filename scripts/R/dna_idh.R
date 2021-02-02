#!/usr/bin/env R

print("SCRIPT DEPRECATED - IDH STATUS IS PART OF gsam.cnv.metadata (scripts/R/gsam_metadata.R)")
stopifnot(FALSE)



# 
# 
# dna_idh <- read.delim("data/gsam/DNA/GBM_allSamples_PassVariants_1736_GBM_TERT_05Jan2018_flags_VAFs.txt",stringsAsFactors=F,sep="\t")
# 
# dna_idh_tested <- data.frame(Sample=unique(dna_idh$Sample))
# dna_idh_tested$Sample <- as.character(dna_idh_tested$Sample)
# 
# dna_idh$Normal <- NULL
# #dna_idh$Sample <- dna_idh$Sample
# dna_idh <- dna_idh[ dna_idh$Gene %in% c("IDH1", "IDH2") ,]
# 
# is_idh1 <- dna_idh$Gene == "IDH1" &  gsub("^([A-Z][0-9]+).+$","\\1",dna_idh$Protein) == "R132"
# is_idh2 <- dna_idh$Gene == "IDH2" &  gsub("^([A-Z][0-9]+).+$","\\1",dna_idh$Protein) %in% c("R140", "R172")
# #R140 and 20% at R172. 
# 
# 
# #is_idh1b <- dna_idh$Sample %in%  c("PD29228a","PD29228c","PD29263a","PD29264a2","PD29264c","PD30239c","PD30242a3","PD30242c3","PD36768a","PD36768c","PD36770a","PD36770c","PD36772c","PD36783a","PD36783c")
# 
# 
# 
# #dna_idh <- dna_idh[is_idh1 | is_idh2 | is_idh1b ,] # n=25
# dna_idh <- dna_idh[is_idh1 | is_idh2  ,] # n=18
# 
# rm(is_idh1, is_idh2)
# 
# 
# dna_idh <- dna_idh[,colnames(dna_idh) %in% c("Sample", "ID", "ShriramVAF")]
# #print(dim(dna_idh))
# 
# dna_idh <- merge(dna_idh_tested, dna_idh, by.x = "Sample", by.y = "Sample" , all.x=T)
# rm(dna_idh_tested)
# #print(dim(dna_idh_tested))
# #print(dim(dna_idh))
# 
# dna_idh[is.na(dna_idh$ID),]$ID <- "-"
# dna_idh[is.na(dna_idh$ShriramVAF),]$ShriramVAF <- "-"
# colnames(dna_idh) <- c("Sample","IDH.mut", "IDH.mut.VAF")
# 
# tmp <- gsam.cnv.metadata[,colnames(gsam.cnv.metadata) %in% c("sid","donor_ID")]
# dna_idh <- merge(dna_idh, tmp, by.x="Sample", by.y = "sid", all.x=TRUE)
# rm(tmp)
# 
# 
# # show all duplicated entries:
# #dna_idh[dna_idh$donor_ID  %in% dna_idh$donor_ID[!is.na(dna_idh$donor_ID) & duplicated(dna_idh$donor_ID)],]
# # all are idh negative so removal can be done safely
# dna_idh$duplicate <- !is.na(dna_idh$donor_ID) & duplicated(dna_idh$donor_ID)




