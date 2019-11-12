#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(ggplot)
library(tidyr)
library(scales)

# ---- load data ----
#source("scripts/R/ligands.R")
#source("scripts/R/gsam_metadata.R")

data.qc <- read.table("data/output/tables/gsam_featureCounts_readcounts.txt.summary",stringsAsFactors = F,comment="#",header=T,row.names = 1)
colnames(data.qc) <- gsub("processed.(.+).Aligned.sortedByCoord.out.markdup.bam","\\1",colnames(data.qc))
#data.qc <- data.qc[rowSums(data.qc) > 0,]

t <- data.frame(t(data.qc))
t <- t[order(t$Assigned,decreasing=F, rownames(t)),]
t <- t[,order(colSums(t),decreasing=T)]
colnames(t)[colnames(t) == "Assigned"] <- "ZAssigned"
t$x <- 1:nrow(t)
t <- t[,colSums(t) > 0]
t$Sample <- as.factor(rownames(t))
t <- gather(t, type, count, -Sample, -x)
t$count <- t$count / 1000000


ggplot(t, aes(x = x ,y = count, fill=type, label=Sample)) +
  coord_flip() + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_y_continuous(labels = unit_format(unit = "M"))





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



