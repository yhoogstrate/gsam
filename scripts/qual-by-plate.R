#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- load plate layout ----
p <- read.table("ignore/plate-layout.txt",header=T,stringsAsFactors = F)

p[,4] <- gsub("A",1,p[,4])
p[,4] <- gsub("B",2,p[,4])
p[,4] <- gsub("C",3,p[,4])
p[,4] <- gsub("D",4,p[,4])
p[,4] <- gsub("E",5,p[,4])
p[,4] <- gsub("F",6,p[,4])
p[,4] <- gsub("G",7,p[,4])
p[,4] <- gsub("H",8,p[,4])


# ---- load stats ----
source("scripts/R/gsam_metadata.R")
d <- gsam.metadata

# ---- plot: ----

png("output/figures/qc_stats_1st_96.png",width=1200,height=800)

plot(d$percentage.duplicate.reads , d$percentage.exon.counts, pch=19, cex=1, xlab="Percentage duplicate reads" , ylab = "Percentage reads mapped to exons")
text(d$percentage.duplicate.reads , d$percentage.exon.counts, d$sample,pos=2, cex=1.2)

dev.off()


# ---- correlation low complexity & GC RMSE ----

png("output/figures/low_complexity_x_gc_error.png", width=2*480, height = 2 * 480 , res =72  * 2)

plot(d$percentage.low.complexity.reads , d$RMSE, pch=19, cex=1, log="xy", xlab="Percentage low complexity reads", ylab="RMSE (error/deviation) from 25% ACUG distribution", main = "Correlation low complexity & GC")
text(d$percentage.low.complexity.reads , d$RMSE, d$sample.id, pch=19, pos=2, cex=0.8, col="darkgray")

dev.off()



# ---- plot duplicate reads x plate indeling ----

png("output/figures/percentage_duplicates_by_plate_layout.png",width=1200,height=800, res=72*1.4)

plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",
     main=paste0("Percentage Duplicates by plate layout (red=",
                  max(d$percentage.duplicate.reads),
                 "%, black=",
                 min(d$percentage.duplicate.reads),"%)"),
     xlab="Plate layout: 1:12",ylab="Plate layout: A:H")

for(x in 1:12) {
  for(y in 1:8) {
    sid <- p[p$A.H == y & p$X1.12 == x,]$seqID
    
    qual = d[d$sample.id == sid,]
    col_dup <- (qual$percentage.duplicate.reads - min(d$percentage.duplicate.reads)) / (max(d$percentage.duplicate.reads) - min(d$percentage.duplicate.reads))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_dup ,0.0,0.0) )
    
    text(x, y, sid,cex=1,col="white")
  }
}

abline(h=(1:(12+1)) - 0.5)
abline(v=(1:(12+1)) - 0.5)


dev.off()

# ---- low complexity x plate layout

png("output/figures/percentage_low_complexity_reads_by_plate_layout.png",width=1200,height=800, res=72*1.4)

plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",
     main=paste0("Percentage Low Complexity reads by plate layout (red=",
                  max(d$percentage.low.complexity.reads),
                 "%, black=",
                 min(d$percentage.low.complexity.reads),"%)"),
     xlab="Plate layout: 1:12",ylab="Plate layout: A:H")

for(x in 1:12) {
  for(y in 1:8) {
    sid <- p[p$A.H == y & p$X1.12 == x,]$seqID
    
    qual = d[d$sample.id == sid,]
    
    # log transfor to make gradient better
    col_dup <- (log(qual$percentage.low.complexity.reads) - log(min(d$percentage.low.complexity.reads))) / (log(max(d$percentage.low.complexity.reads)) - log(min(d$percentage.low.complexity.reads)))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_dup ,0.0,0.0) )
    
    text(x, y, sid,cex=1,col="white")
  }
}

abline(h=(1:(12+1)) - 0.5)
abline(v=(1:(12+1)) - 0.5)


dev.off()


# ---- plot plate layout x volumes ----

png("output/figures/percentage_exon_counts_by_plate_layout.png",width=1200,height=800,res=72*1.4)


plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",
     main=paste0("Low percentage exon counted reads (red=",round(min(d$percentage.exon.counts),2),
                 "%, black=",round(max(d$percentage.exon.counts),2),"%)"),
     xlab="Plate layout: 1:12",ylab="Plate layout: A:H")
for(x in 1:12) {
  for(y in 1:8) {
    sid <- p[p$A.H == y & p$X1.12 == x,]$seqID
    qual <- d[d$sample.id == sid,]
    
    col_discarded <- 1 - ( (qual$percentage.exon.counts - min(d$percentage.exon.counts)) / max(d$percentage.exon.counts - min(d$percentage.exon.counts)))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_discarded ,0.0,0.0) )
    
    text(x, y, sid,cex=1,col="white")
  }
}

abline(h=(1:(12+1)) - 0.5)
abline(v=(1:(12+1)) - 0.5)


dev.off()






f <- read.delim('ignore/rna-volumes-s1-96.txt',stringsAsFactors = F,header=T,sep="\t")

png("output/figures/duplicates_x_ng_rna_input.png",width=1200,height=800)
plot(f[match(e$Sample, f$seqID),]$ng.DNA.ul ,  e$percentage.duplicate.reads, xlab= "ng RNA input", ylab="% duplicates" )
dev.off()

png("output/figures/duplicates_x_Vsample.png",width=1200,height=800)
plot(f[match(e$Sample, f$seqID),]$Vsample..uL. ,  e$percentage.duplicate.reads, xlab= "Vsample (uL)", ylab="% duplicates" )
dev.off()

png("output/figures/duplicates_x_Vh2o.png",width=1200,height=800)
plot(f[match(e$Sample, f$seqID),]$Vh2o..uL. ,  e$percentage.duplicate.reads, xlab= "Vh2o (uL)", ylab="% duplicates" )
dev.off()





