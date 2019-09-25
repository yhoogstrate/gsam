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
# duplicate reads
d <- read.delim("output/tables/duplicate_reads_sambamba_stats.txt",header=T,sep="\t",stringsAsFactors = F)

# low complexity reads
e <- read.delim("output/tables/low_complexity_reads.txt",header=T,sep="\t",stringsAsFactors = F)
e$percentage.low.complexity.reads <- gsub("%","",e$percentage.low.complexity.reads,fixed=T)
d <- merge(d, e,by.x="sample.id", by.y = "sample")
dim(d)
rm(e)

# mapping stats
e <- read.delim("output/tables/star_percentage_mapped.txt",header=T,sep="\t",stringsAsFactors = F)
e$percentage.uniquely.mapped <- gsub("%","",e$percentage.uniquely.mapped,fixed=T)
e$percentage.multimap <- gsub("%","",e$percentage.multimap,fixed=T)
d <- merge(d, e,by.x="sample.id", by.y = "sample")
dim(d)
rm(e)

# featureCounts stats
e <- read.delim("output/tables/featureCounts_gsam_1st96.exon-level.txt.summary",header=T,sep="\t",stringsAsFactors = F,row.names=1)
colnames(e) <- gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(e))
e <- t(e) # transpose
e <- data.frame(e)
e$read.count <- as.numeric(rowSums(e))
e$sample <-  rownames(e)
d <- merge(d, e,by.x="sample.id", by.y = "sample")
dim(d)
rm(e)


# ---- plot: ----

png("output/figures/qc_stats_1st_96.png",width=1200,height=800)

plot(d$percentage.duplicate.reads , d$Assigned/ d$read.count, pch=19, cex=1)
text(d$percentage.duplicate.reads , d$Assigned/ d$read.count, d$sample,pos=2, cex=1.2)

dev.off()


# ---- plot ----

png("output/figures/percentage_duplicates_by_plate_layout.png",width=1200,height=800)
plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",main=paste0("Duplicates by plate layout (red=",max(e$percentage.duplicate.reads),"%, black=",min(e$percentage.duplicate.reads),"%)"),xlab="Plate layout: 1:12",ylab="Plate layout: A:H")
for(x in 1:12) {
  for(y in 1:8) {
    sid <- d[d$A.H == y & d$X1.12 == x,]$seqID
    
    qual = e[e$Sample == sid,]
    col_dup <- (qual$percentage.duplicate.reads - min(e$percentage.duplicate.reads)) / (max(e$percentage.duplicate.reads) - min(e$percentage.duplicate.reads))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_dup ,0.0,0.0) )
    
    text(x, y, sid,cex=1.4,col="white")
  }
}

abline(h=(1:(12+1)) - 0.5)
abline(v=(1:(12+1)) - 0.5)
dev.off()







png("output/figures/percentage_exon_counts_by_plate_layout.png",width=1200,height=800)
plot(c(1-0.5,12+0.5), c(1-0.5,8+0.5), type="n",main=paste0("Low percentage exon counted reads (red=",min(e$Assigned),"%, black=",max(e$Assigned),"%)"),xlab="Plate layout: 1:12",ylab="Plate layout: A:H")
for(x in 1:12) {
  for(y in 1:8) {
    sid <- d[d$A.H == y & d$X1.12 == x,]$seqID
    
    qual = e[e$Sample == sid,]
    col_discarded <- 1 - ( (qual$Assigned - min(e$Assigned)) / max(e$Assigned - min(e$Assigned)))
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = rgb( col_discarded ,0.0,0.0) )
    
    text(x, y, sid,cex=1.4,col="white")
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





