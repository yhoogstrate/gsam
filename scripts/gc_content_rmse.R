#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.delim("output/tables/gc_content_rmse.txt",stringsAsFactors = F)
e <- d[order(d$RMSE, d$sample.id),]

# mediane waarden in niet outlier samples:
# a     c     t     g
# 25.11 24.61 25.61 24.67


png("output/figures/gc_content_rmse.png", width=480*4,heigh=480*2,res=72*2)
plot(c(1,nrow(e)), c(0,50), type="n", ylab="",xlab="Sample nr. / order",main="GC deviation")

for(i in 1:nrow(e)) {
  if(i %% 2 == 0) {
    rect(i - 0.5 , 0 , i + 0.5 ,100 , col=rgb(0.95,0.95,0.95,0.55) , border=NA)
  }
}

abline(h=25,col=rgb(0.85, 0.85, 0.85), lty=2)

points(1:nrow(e), e$RMSE , pch=19 , cex=0.36)


points(1:nrow(e) , e$percentage.A, col=rgb(0,0.5,0.0),pch="-", cex=1.5)
points(1:nrow(e) , e$percentage.C, col=rgb(0.5,0.5,0),pch="-", cex=1.5)
points(1:nrow(e) , e$percentage.T, col=rgb(0.5,0,0.5),pch="-", cex=1.5)
points(1:nrow(e) , e$percentage.G, col="blue",pch="-", cex=1.5)

text(1:nrow(e) - 0.8,  45, e$sample.id , pos=4, srt=90, cex=0.75)

legend(0, 20, pch=c("-","-","-","-","o"), col=c(rgb(0,0.5,0.0),rgb(0.5,0.5,0),rgb(0.5,0,0.5),"blue","black"), c("%A","%C","%T","%G","RMSE from 25%") )

dev.off()
