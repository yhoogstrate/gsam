#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.table("ignore/plate-layout.txt",header=T,stringsAsFactors = F)
d[,4] <- gsub("A",1,d[,4])
d[,4] <- gsub("B",2,d[,4])
d[,4] <- gsub("C",3,d[,4])
d[,4] <- gsub("D",4,d[,4])
d[,4] <- gsub("E",5,d[,4])
d[,4] <- gsub("F",6,d[,4])
d[,4] <- gsub("G",7,d[,4])
d[,4] <- gsub("H",8,d[,4])

e <- read.delim("ignore/qc-stats.txt",header=T,sep="\t",stringsAsFactors = F)
e$Sample <- gsub("^Pfrench_([^_]+)_.+$","\\1", e$Sample )
e[,2] <- as.numeric(gsub("%","",e[,2]))
e[,3] <- as.numeric(gsub("%","",e[,3]))
e[,4] <- as.numeric(gsub(",","",e[,4]))
colnames(e) <- c("Sample", "duplicate.percentage", "exon.counted.percentage","exon.counted.reads")


png("output/figures/qc_stats_1st_96.png",width=1200,height=800)
plot(e$duplicate.percentage, e$exon.counted.percentage, pch=19, cex=1)
text(e$duplicate.percentage, e$exon.counted.percentage, e$Sample,pos=2, cex=1.2)
dev.off()


png("output/figures/percentage_duplicates_by_plate_layout.png",width=1200,height=800)
plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",main=paste0("Duplicates by plate layout (red=",max(e$duplicate.percentage),"%, black=",min(e$duplicate.percentage),"%)"),xlab="Plate layout: 1:12",ylab="Plate layout: A:H")
for(x in 1:12) {
  for(y in 1:8) {
    sid <- d[d$A.H == y & d$X1.12 == x,]$seqID
    
    qual = e[e$Sample == sid,]
    col_dup <- (qual$duplicate.percentage - min(e$duplicate.percentage)) / (max(e$duplicate.percentage) - min(e$duplicate.percentage))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_dup ,0.0,0.0) )
    
    text(x, y, sid,cex=1.4,col="white")
  }
}

abline(h=(1:(12+1)) - 0.5)
abline(v=(1:(12+1)) - 0.5)
dev.off()







png("output/figures/percentage_exon_counts_by_plate_layout.png",width=1200,height=800)
plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",main=paste0("Low percentage exon counted reads (red=",min(e$exon.counted.percentage),"%, black=",max(e$exon.counted.percentage),"%)"),xlab="Plate layout: 1:12",ylab="Plate layout: A:H")
for(x in 1:12) {
  for(y in 1:8) {
    sid <- d[d$A.H == y & d$X1.12 == x,]$seqID
    
    qual = e[e$Sample == sid,]
    col_discarded <- 1 - ( (qual$exon.counted.percentage - min(e$exon.counted.percentage)) / max(e$exon.counted.percentage - min(e$exon.counted.percentage)))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_discarded ,0.0,0.0) )
    
    text(x, y, sid,cex=1.4,col="white")
  }
}

abline(h=(1:(12+1)) - 0.5)
abline(v=(1:(12+1)) - 0.5)

dev.off()






f <- read.delim('ignore/rna-volumes-s1-96.txt',stringsAsFactors = F,header=T,sep="\t")

png("output/figures/duplicates_x_ng_rna_input.png",width=1200,height=800)
plot(f[match(e$Sample, f$seqID),]$ng.DNA.ul ,  e$duplicate.percentage, xlab= "ng RNA input", ylab="% duplicates" )
dev.off()

png("output/figures/duplicates_x_Vsample.png",width=1200,height=800)
plot(f[match(e$Sample, f$seqID),]$Vsample..uL. ,  e$duplicate.percentage, xlab= "Vsample (uL)", ylab="% duplicates" )
dev.off()

png("output/figures/duplicates_x_Vh2o.png",width=1200,height=800)
plot(f[match(e$Sample, f$seqID),]$Vh2o..uL. ,  e$duplicate.percentage, xlab= "Vh2o (uL)", ylab="% duplicates" )
dev.off()





