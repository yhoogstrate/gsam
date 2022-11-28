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
  dplyr::rename(sid = sample.id) %>%
  dplyr::arrange(duplication.rate, decreasing=T) %>%
  dplyr::mutate(run = gsub("^.+_([^_]{7,15})_S.+$","\\1", filename)) %>%
  dplyr::mutate(lane = gsub("^.+_(L[0-9]+)_.+$","\\1",filename) ) %>%
  dplyr::mutate(run.lane <- paste0(run,"_", lane) ) %>%
  dplyr::mutate(new.batch = grepl("new", sid)) %>%
  dplyr::mutate(name.strip = ifelse(new.batch == T ,gsub("-new", "", .$sid, fixed=T),"") ) %>%
  dplyr::mutate(old.batch = sid %in% name.strip) %>%
  dplyr::mutate(name.strip = NULL) %>%
  dplyr::mutate(batch = dplyr::case_when ((old.batch == T & new.batch == F) ~ "old" , (old.batch == F & new.batch == T) ~ "new" , (old.batch == F & new.batch == F) ~ "single" )) %>%
  dplyr::mutate(batch = as.factor(batch)) %>%
  dplyr::mutate(ampli_factor = paste0(round (1 / (1 - (duplication.rate )),1),'x' )) # fastp underestimates the percentage duplicates?



gsam.markdup.stats <- read.delim("output/tables/qc/flagstat/flagstats.matrix.txt" ,stringsAsFactors = F) %>%
  dplyr::mutate(pct.duplicates = duplicates / in.total..QC.passed.reads...QC.failed.reads. * 100) %>%
  dplyr::arrange(pct.duplicates , decreasing = T) %>%
  dplyr::mutate(ampli_factor = paste0(round (1 / (1 - (0.01 * pct.duplicates )),1),'x' )) %>% # fastp underestimates the percentage duplicates?
  dplyr::mutate(new.batch = grepl("new", sid)) %>%
  dplyr::mutate(name.strip = ifelse(new.batch == T ,gsub("-new", "", .$sid, fixed=T),"") ) %>%
  dplyr::mutate(old.batch = sid %in% name.strip) %>%
  dplyr::mutate(name.strip = NULL) %>%
  dplyr::mutate(batch = dplyr::case_when ((old.batch == T & new.batch == F) ~ "old" , (old.batch == F & new.batch == T) ~ "new" , (old.batch == F & new.batch == F) ~ "single" )) %>%
  dplyr::mutate(batch = as.factor(batch)) %>%
  dplyr::mutate(sid = gsub(".txt","", .$sid, fixed=T))




# ---- ----

# ,col=type,fill=type
plt <- gsam.fastp.duplication.stats

ggplot(plt , aes(x = input.reads ,y = duplication.rate, col=batch, label=sid)) + 
  geom_point() + 
  #ggrepel::geom_text_repel(data=subset(plt, new.batch == T), col="black", size=2.5) +
  #ggrepel::geom_text_repel(data=subset(plt, old.batch == T), col="black", size=2.5) +
  youri_gg_theme +
  scale_color_manual(name = NULL, values = c('old'='#59B2E6', 'new'='#009E74', 'single'='#CB75A4'))  # , 




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
plot(c(1,nrow(plt)), c(-4,100),
     type="n",xlab = "sample",
     ylab = "percentage duplicate reads (red)",
     main=paste0("duplicatie reads ratio GSAM [mean=",
      round(mean( plt$duplication.rate * 100),1),"% in ",nrow(plt)," FASTQ samples]"))

for(i in 1:nrow(plt)) {
  lines(c(i,i), c(0, 100 - (plt$duplication.rate[i] * 100)), col="chartreuse4", lwd=1.5)
  lines(c(i,i), c(100 - (plt$duplication.rate[i] * 100) , 100 ), col="red", lwd=1.5)
  
  text(i + 1, -1,  plt$sample.id[i] , cex=0.6, srt=90,pos=2)
}

dy <- 3.5
rect(1-1,79-dy,nrow(e)+1,79+dy,border=NA,col=rgb(1,1,1,0.65))

for(i in 1:nrow(e)) {
  text(i  + 1, 80 + 2, e[i,]$ampli_factor, cex=0.7, srt=90,pos=2)
}

text(0  , 80 + 2, "d=1/(1-p)", cex=0.7, srt=90,pos=2,col=rgb(0.2,0.2,0.2))


dev.off()

# ---- plot per bam ----

plt <- gsam.markdup.stats %>%
  dplyr::mutate(order = reorder(sid, pct.duplicates) )


ggplot(plt , aes(x = order,y = 100 - pct.duplicates, label=sid )) + 
  geom_segment(aes(y = 0, yend = 100 - pct.duplicates, xend=reorder(sid, pct.duplicates)), col="darkgreen", lwd=1.5, alpha=0.3) + 
  geom_segment(aes(yend = 100, xend=reorder(sid, pct.duplicates)), col="red" , lwd=1.5, alpha=0.3) + 
  geom_point(data=subset(plt, batch=="new")) + 
  geom_text(aes(y = 100 - pct.duplicates + 7), data=subset(plt, batch=="new" ), col="black", size=2.5, srt=90) + 
  geom_text(aes(y = 100 - pct.duplicates - 5), data=subset(plt, batch=="old" ), col="black", size=2.5, srt=90) + 
  geom_point(data=subset(plt, batch=="old" ) , pch=2 )

ggsave("output/figures/qc/duplicates_reads_sambamba.pdf")



