#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.delim("output/tables/v3_extract_readcounts_rotated.txt",stringsAsFactors = F)
e <- d[!is.na(d$resection.1.wt) & !is.na(d$resection.2.wt),]


png("output/figures/v3_plot_expression_per_resection.png", width=480*3,height=480*3, res=72*3)

off <- 0.05
plot(c(1.0 - off,2.0 + off), c(0,100),type="n", ylab="%vIII",xlab="resection",main=paste0("vIII over resection/time in ",nrow(e)," matching samples"))

for(i in 1:nrow(e)) {
  p1 <- 100.0 *   e[i,4] / (e[i,3] + e[i,4])
  p2 <- 100.0 *   e[i,7] / (e[i,6] + e[i,7])
  
  points(1, p1, pch=19,cex=0.6)
  points(2, p2, pch=19,cex=0.6)
  lines(c(1,2), c(p1, p2))
  
  text(1,p1, e[i,1],cex=0.65,pos=2,col="darkgray")
  text(2,p2, e[i,1],cex=0.65,pos=4,col="darkgray")
}

dev.off()




