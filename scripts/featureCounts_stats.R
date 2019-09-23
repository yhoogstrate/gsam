#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)

# ---- load data ----
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")

d <- read.delim("output/tables/featureCounts_gsam_1st96.exon-level.txt",stringsAsFactors = F,comment="#")

colnames(d) <- gsub("X.data.users.youri.mnt.neurogen.ro.gsam.RNA.alignments.","",colnames(d),fixed=T)
colnames(d) <- gsub(".sambamba.dedup.bam","",colnames(d),fixed=T)

rownames(d) <- d$Geneid
d$Geneid <- NULL
d$Chr <- NULL
d$Start <- NULL
d$End <- NULL
d$Strand <- NULL
d$Length <- NULL

#plot(sort(log(colSums(d))))
#abline(h=log(2000000))
# min 2mljn reads to take the low-Q samples out
d <- d[,colSums(d) > 2000000]
#plot(sort(log(colSums(d))))
#abline(h=log(2000000))


png("output/figures/featureCounts_stats.png",width=3*480,height=2*480,res=2*72)

sel <- order(colSums(d))
plot(  c(1, ncol(d)) , c(0, max(colSums(d)[sel] / 1000000)) ,xlab = "Sample nr. / order", ylab = "Million reads counted in exons", type="n",main="Exon counted reads")

#max
y <- max( colSums(d) ) / 1000000
lines(c(1, nrow(d)), c(y, y), col="darkgray", lty=2)
points(colSums(d)[sel] / 1000000,pch=19,cex=0.8)
text(ncol(d)/2, y, paste0("max=",round(y,2),"M"),pos=1,col="black")

#min
y <- min( colSums(d) ) / 1000000
lines(c(1, nrow(d)), c(y, y), col="darkgray", lty=2)
points(colSums(d)[sel] / 1000000,pch=19,cex=0.8)
text(ncol(d)/2, y, paste0("min=",round(y,2),"M"),pos=3,col="black")

text(    1:ncol(d),    (colSums(d)[sel] / 1000000) - 0.6,    gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(d)) , srt=90 ,pch = 19 , cex=0.6, pos=1   )

dev.off()



