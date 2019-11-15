#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load data ----

source("scripts/R/gsam_metadata.R")

# ---- ----

counts_stranded <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt",stringsAsFactors = F,comment="#")
colnames(counts_stranded) <- gsub("^[^\\.]+\\.([^\\]+)\\.Aligned.sorted.+$","\\1",colnames(counts_stranded),fixed=F)

rownames(counts_stranded) <- counts_stranded$Geneid
counts_stranded$Geneid <- NULL
counts_stranded$Chr <- NULL
counts_stranded$Start <- NULL
counts_stranded$End <- NULL
counts_stranded$Strand <- NULL
counts_stranded$Length <- NULL


# 


counts_unstranded <- read.delim("output/tables/gsam_featureCounts_readcounts.unstranded.txt",stringsAsFactors = F,comment="#")
colnames(counts_unstranded) <- gsub("^[^\\.]+\\.RNA.alignments.([^\\]+)\\.Aligned.sorted.+$","\\1",colnames(counts_unstranded),fixed=F)

rownames(counts_unstranded) <- counts_unstranded$Geneid
counts_unstranded$Geneid <- NULL
counts_unstranded$Chr <- NULL
counts_unstranded$Start <- NULL
counts_unstranded$End <- NULL
counts_unstranded$Strand <- NULL
counts_unstranded$Length <- NULL




# 

counts_stranded[counts_stranded < 1 ] <- 1

counts_antistranded <- counts_unstranded - counts_stranded
counts_antistranded[counts_antistranded < 1 ] <- 1


# ---- selecet only genes that have some reads ----

sel <- rowSums(counts_stranded) > 372 * 15
counts_stranded <- counts_stranded[sel,]
counts_antistranded <- counts_antistranded[sel,]
counts_unstranded <- counts_unstranded[sel,]
rm(sel)

# 

#stranded_antistranded_ratio <- log( (counts_stranded + 0.0001) / (counts_antistranded + 0.0001) )
stranded_antistranded_ratio <- log( (counts_stranded ) / (counts_antistranded ) )
no_dna_c <- colMeans(stranded_antistranded_ratio) >= 2.6
dna_c <- colMeans(stranded_antistranded_ratio) < 2.6 & colMeans(stranded_antistranded_ratio) >= 1.8
heavy_dna_c <- colMeans(stranded_antistranded_ratio) < 1.8

plot(c(1, nrow(stranded_antistranded_ratio)) , c(-8,12) , type='n',ylab="log(counts stranded / count anti-stranded)")
abline(h=0, col="gray",lty=2)
for(i in which(no_dna_c)) {
  lines(sort(stranded_antistranded_ratio[,i], decreasing=T))
}
for(i in which(dna_c)) {
  lines(sort(stranded_antistranded_ratio[,i], decreasing=T),col="orange")
}
for(i in which(heavy_dna_c)) {
  lines(sort(stranded_antistranded_ratio[,i], decreasing=T),col="red")
}
rm(i)
rm(counts_antistranded, counts_stranded, counts_unstranded)

df <- data.frame(sid = colnames(stranded_antistranded_ratio),
                 ratio.stranded.antistranded = colMeans(stranded_antistranded_ratio),
                 ratio.stranded.antistranded.lin = exp(colMeans(stranded_antistranded_ratio)), 
                 ratio.stranded.antistranded.dna = rep(NA, ncol(stranded_antistranded_ratio)) )
df$ratio.stranded.antistranded.dna [no_dna_c] <- "No-DNA-contamination"
df$ratio.stranded.antistranded.dna [dna_c] <- "Moderate-DNA-contamination"
df$ratio.stranded.antistranded.dna [heavy_dna_c] <- "Heavy-DNA-contamination"
df$ratio.stranded.antistranded.dna <- as.factor(df$ratio.stranded.antistranded.dna)
df$sid <- gsub(".","-",as.character(df$sid), fixed=T)

rm(dna_c, heavy_dna_c, no_dna_c, stranded_antistranded_ratio)
gsam.rna.metadata <- merge(gsam.rna.metadata, df, by.x="sid",by.y="sid")


