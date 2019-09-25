#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.delim("output/tables/low_complexity_reads.txt", header=T, stringsAsFactors = F)
d$percentage_numerical <- as.numeric(gsub("%","",d$percentage.low.complexity.reads,fixed=T))
e <- d[order(d$percentage_numerical),]

png("output/figures/low_complexity_reads.png",width=480*3,height=1.7*480,res=72*1.8)
plot(c(1, max(nrow(e))) , c(0, max(e$percentage_numerical) + 1), main="Percentage low-complexity discarded reads",type="n", xlab="Sample nr. / order", ylab = "Percentage reads discarted due to low complexity")
for(i in 1:nrow(e)) {
  print(i)
  lines(c(i,i), c(0,e[i,]$percentage_numerical),col="darkgray",lty=1)
  text(i,e[i,]$percentage_numerical + 0.5, e[i,]$sample, srt=90, pos=3,cex=0.6)
}

points(1:nrow(e) , e$percentage_numerical, pch=19,cex=0.6)
dev.off()




#g <- data.frame(
#  x = 1:nrow(e),
#  y_max = e$percentage_numerical,
#  y_min = rep(0, nrow(e))
#  )

#library(ggplot2)
#ggplot(g, aes(x=x, y=y_max)) + geom_point(shape=1) 




