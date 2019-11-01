#!/usr/bin/env R



wd <- "/home/youri/projects/gsam"
wd <- "/home/yhoogstrate/projects/gsam"

setwd(wd)

rm(wd)


# ---- load: libs ----

# ---- load: functions ----

# ---- load: data ----
# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")

#gsam.metadata$v3.rnaseq.v3.percentage


# ----  ----

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

# ---- Final_qPCR_EGFR_GSAM.csv -----

tmp <- read.delim("data/DNA/Final_qPCR_EGFR_GSAM.csv",sep=",",stringsAsFactors = F)



tmp2 <- gsam.metadata$pid[duplicated(gsam.metadata$pid)] # ids of which two time points have been RNA-sequenced
tmp2 <- gsam.metadata[gsam.metadata$pid %in% tmp2,]
tmp2.1 <- tmp2[tmp2$resection == 1,]
tmp2.2 <- tmp2[tmp2$resection == 2,]

stopifnot(sum(tmp2.1$pid != tmp2.2$pid) == 0)

tmp2 <- data.frame(
  pid  = tmp2.1$pid , 
  primary.v3.rnaseq.v3.percentage = tmp2.1$v3.rnaseq.v3.percentage ,
  primary.v3.rnaseq.wt.reads.v3 = tmp2.1$v3.rnaseq.wt.reads.v3 ,
  primary.v3.rnaseq.vIII.reads.v3 = tmp2.1$v3.rnaseq.vIII.reads.v3 ,
  
  secondary.v3.rnaseq.v3.percentage = tmp2.2$v3.rnaseq.v3.percentage ,
  secondary.v3.rnaseq.wt.reads.v3 = tmp2.2$v3.rnaseq.wt.reads.v3 ,
  secondary.v3.rnaseq.vIII.reads.v3 = tmp2.2$v3.rnaseq.vIII.reads.v3 )
rm(tmp2.1, tmp2.2)

#tmp2 <- tmp2[!is.na(tmp2$primary.v3.rnaseq.v3.percentage) & !is.na(tmp2$secondary.v3.rnaseq.v3.percentage),]
tmp2$logfc.rna.seq <- log2((tmp2$primary.v3.rnaseq.v3.percentage + 0.00001) / (tmp2$secondary.v3.rnaseq.v3.percentage + 0.00001) )
print(dim(tmp2))

tmp3 <- merge(tmp2, tmp, by.x = "pid", by.y = "patientID")
print(dim(tmp3))



#intersect(tmp2$pid , tmp$patientID)
plot(tmp$percentageEGFRvIII , tmp$recurrent_percentageEGFRvIII, cex=0.5, pch=19)
text(tmp$percentageEGFRvIII , tmp$recurrent_percentageEGFRvIII, tmp$patientID, cex=0.75,pos=2)

high_mut <- tmp$patientID %in% c("FAN","HAC","AQA","AOA","BAB","JAG","GAR","AAT","BAK","JAC","FAG","JAF","BAD","EAZ","CAO")
points(tmp[high_mut,]$percentageEGFRvIII , tmp[high_mut,]$recurrent_percentageEGFRvIII, cex=0.5, pch=19,col="red")

plt.data <- data.frame(x= tmp$percentageEGFRvIII ,  y=tmp$recurrent_percentageEGFRvIII, pid=tmp$patientID)

#lines(c(0,100),c(0,100))
#ggplot(plt.data, aes(x=x , y = y),pid=tmp$patientID)
#ggplot(plt.data, aes(x=x , y = y)) + geom_point() + geom_smooth()
#ggplot(plt.data, aes(x=x , y = y)) + geom_point() + geom_smooth() + geom_abline(intercept = 0,slope=1)
ggplot(plt.data, aes(x=x , y = y, label=pid)) + geom_point() +  geom_abline(intercept = 0,slope=1) + ggrepel::geom_label_repel()

filter = abs(plt.data$x - plt.data$y) > 4
ggplot(plt.data, aes(x=x , y = y, label=pid)) + geom_point() +  geom_abline(intercept = 0,slope=1) + ggrepel::geom_label_repel(data=plt.data[filter,]) + stat_binhex(bins = 50)




