#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.delim("output/tables/duplicate_reads_sambamba_stats.txt",stringsAsFactors = F)

e <- d[order(d[,2]),]
e$ampli_factor <- paste0(round (1 / (1 - (e$percentage.duplicate.reads/100)),1),'x' )



png("output/figures/duplicate_reads_sambamba_stats.png",width=3*480,height=2*480,res=1.8*72)
plot(c(1,nrow(e)), c(-4,100),
     type="n",xlab = "sample",
     ylab = "percentage duplicate reads (red)",
     main=paste0("duplicatie reads ratio GSAM [mean=",
      round(mean(e$percentage.duplicate.reads),1),"% in ",nrow(e)," samples]"))

for(i in 1:nrow(e)) {
  lines(c(i,i), c(0, 100-e[i,2]), col="chartreuse4", lwd=1.5)
  lines(c(i,i), c(100 - e[i,2] , 100 ), col="red", lwd=1.5)
  
  text(i + 1, -1,  e[i,]$sample.id, cex=0.6, srt=90,pos=2)
}

dy <- 3.5
rect(1-1,79-dy,nrow(e)+1,79+dy,border=NA,col=rgb(1,1,1,0.65))

for(i in 1:nrow(e)) {
  text(i  + 1, 80 + 2, e[i,]$ampli_factor, cex=0.7, srt=90,pos=2)
}

text(0  , 80 + 2, "d=1/(1-p)", cex=0.7, srt=90,pos=2,col=rgb(0.2,0.2,0.2))


dev.off()




