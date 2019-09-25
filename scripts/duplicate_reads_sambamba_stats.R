#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.delim("output/tables/duplication-stats.txt",stringsAsFactors = F)

e <- d[order(d[,2]),]

plot(c(1,nrow(e)), c(0,100), type="n",xlab = "sample", ylab = "percentage duplicated reads (gray)", main=paste0("duplication ratio GSAM [mean=",round(100.0*mean(e$percentage.duplicated.reads),1),"% in ",nrow(e)," samples]"))
for(i in 1:nrow(e)) {
  lines(c(i,i), c(e[i,2] * 100,100), col="green", lwd=1.5)
  lines(c(i,i), c(0,e[i,2] * 100), col="darkgray", lwd=1.5)
}


