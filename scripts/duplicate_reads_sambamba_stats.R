#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load libraries ----

library(ggplot2)

# ---- load functions ----

# ---- load themes ----

source('scripts/R/youri_gg_theme.R')
source('scripts/R/job_gg_theme.R')

# ---- load data ----

# not in the mood to also run sambamba
#d <- read.delim("data/output/tables/duplicate_reads_sambamba_stats.txt",stringsAsFactors = F)
gsam.fastp.duplication.stats <- read.delim("output/tables/fastp_duplication_statistics.txt" ,stringsAsFactors = F) %>%
  dplyr::arrange(duplication.rate, decreasing=T)


gsam.fastp.duplication.stats$run <- gsub("^.+_([^_]{7,15})_S.+$","\\1", gsam.fastp.duplication.stats$filename)
gsam.fastp.duplication.stats$lane <- gsub("^.+_(L[0-9]+)_.+$","\\1",gsam.fastp.duplication.stats$filename)
gsam.fastp.duplication.stats$run.lane <- paste0(gsam.fastp.duplication.stats$run,"_",gsam.fastp.duplication.stats$lane)

gsam.fastp.duplication.stats$sid <- gsub("^(.+)-[0-9]+_.+$","\\1",gsam.fastp.duplication.stats$filename)

gsam.fastp.duplication.stats <- gsam.fastp.duplication.stats %>%
  dplyr::mutate(new.batch = )


# ---- ----

# ,col=type,fill=type

plt <- gsam.fastp.duplication.stats


ggplot(plt , aes(x = input.reads ,y = duplication.rate, col=run, label=sid)) + 
  geom_point() + 
  geom_text(data=subset(plt, run=="AHKK5WDSXX" & lane=="L002"), col="black", size=2.5) +
  youri_gg_theme

ggsave('output/figures/qc/duplicate_read_stats_01.png',width=7,height=6)


ggplot(plt , aes(x = input.reads ,y = duplication.rate, col=run, label=sid)) + 
  geom_point() + 
  geom_text(data=subset(plt, run=="AHKK5WDSXX" & lane=="L002"), col="black", size=2.5) +
  youri_gg_theme

ggsave('output/figures/qc/duplicate_read_stats_02.png',width=7,height=6)



# ------


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




