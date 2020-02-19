#!/usr/bin/env R

setwd("~/projects/gsam")

d <- read.delim("output/tables/star_percentage_mapped.txt",stringsAsFactors = F)
d$percentage.uniquely.mapped <- as.numeric(gsub("%","",d$percentage.uniquely.mapped,fixed=T))
d$percentage.multimap <- as.numeric(gsub("%","",d$percentage.multimap,fixed=T))


e <- d[order(d$percentage.uniquely.mapped,decreasing=TRUE),]



png("output/figures/star_percentage_mapped.png", width=480*3,height=480*2, res=72*2)

plot(c(1,nrow(e)), c(-10,100), type="n", xlab="Sample nr. / order",
     ylab="Percentage aligned reads unique (black) / multi-map (red)", main="STAR mapping stats")
for(i in 1:nrow(e)) {
  lines(c(i,i), c(0,e[i,]$percentage.uniquely.mapped), lwd=1.5)
  lines(c(i,i), c(e[i,]$percentage.uniquely.mapped,e[i,]$percentage.uniquely.mapped + e[i,]$percentage.multimap), col="red", lwd=1.5)
  text(i-0.35, -1, pos=1, e[i,]$sample , srt=90, cex=0.6) 
}


dev.off()




