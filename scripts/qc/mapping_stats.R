#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load libraries ----

library(ggplot2)

# ---- load functions ----

# ---- load themes ----

source('scripts/R/youri_gg_theme.R')
source('scripts/R/job_gg_theme.R')

# ---- load data ----

source('scripts/R/gsam_metadata.R')

gsam.multimap.stats <- read.delim("output/tables/star_percentage_mapped.txt",stringsAsFactors = F)
gsam.multimap.stats$percentage.uniquely.mapped <- as.numeric(gsub("%","",gsam.multimap.stats$percentage.uniquely.mapped,fixed=T))
gsam.multimap.stats$percentage.multimap <- as.numeric(gsub("%","",gsam.multimap.stats$percentage.multimap,fixed=T))

gsam.fastp.duplication.stats <- read.delim("output/tables/fastp_duplication_statistics.txt" ,stringsAsFactors = F)
gsam.fastp.duplication.stats <- gsam.fastp.duplication.stats[order(gsam.fastp.duplication.stats$duplication.rate,decreasing=T),]
gsam.fastp.duplication.stats$run <- gsub("^.+_([^_]{7,15})_S.+$","\\1", gsam.fastp.duplication.stats$filename)
gsam.fastp.duplication.stats$lane <- gsub("^.+_(L[0-9]+)_.+$","\\1",gsam.fastp.duplication.stats$filename)
gsam.fastp.duplication.stats$run.lane <- paste0(gsam.fastp.duplication.stats$run,"_",gsam.fastp.duplication.stats$lane)
gsam.fastp.duplication.stats$sid <- gsub("^(.+)-[0-9]+_.+$","\\1",gsam.fastp.duplication.stats$filename)
gsam.fastp.duplication.stats$sid <- gsub("-repl.+$","-replicate",gsam.fastp.duplication.stats$sid)

gsam.rna.seqruns <- data.frame(sid=gsam.rna.metadata$sid, runs = NA, runlanes= NA, mean.r1.len= NA, mean.r2.len=NA)
for(i in 1:nrow(gsam.rna.seqruns)) {
  sid <- gsam.rna.seqruns$sid[i]
  #sid='FAQ1'
  e <- gsam.fastp.duplication.stats[gsam.fastp.duplication.stats$sid == sid,]
  if(nrow(e) == 0) {
    print(sid)
  }
  stopifnot(nrow(e) >=1)
  e$pct.run <- e$input.reads / sum(e$input.reads) * 100
  e <- e[e$pct.run > 10,]
  gsam.rna.seqruns$runs[i] <- paste(sort(e$run), collapse=",")
  gsam.rna.seqruns$runlanes[i] <- paste(sort(e$run.lane), collapse=",")
  
  e$pct.run <- e$input.reads / sum(e$input.reads)
  gsam.rna.seqruns$mean.r1.len[i] <- sum(e$pct.run * e$read1_mean_length)
  gsam.rna.seqruns$mean.r2.len[i] <- sum(e$pct.run * e$read2_mean_length)
  
  rm(sid)
}


# pct mapped, pct counted
gsam.count.stats <- read.delim('output/tables/gsam_featureCounts_readcounts.txt.summary',stringsAsFactors = F,row.names=1)
colnames(gsam.count.stats) <- gsub(".Aligned.sortedByCoord.out.markdup.bam","",colnames(gsam.count.stats),fixed=T)
colnames(gsam.count.stats) <- gsub("processed.","",colnames(gsam.count.stats),fixed=T)
colnames(gsam.count.stats) <- gsub(".","-",colnames(gsam.count.stats),fixed=T)
gsam.count.stats <- data.frame(t(gsam.count.stats))
gsam.count.stats$pct.Assigned <- gsam.count.stats$Assigned / rowSums(gsam.count.stats) * 100
gsam.count.stats$sid <- rownames(gsam.count.stats)

gsam.count.stats <- gsam.count.stats[match(gsam.rna.metadata$sid, rownames(gsam.count.stats)),]
stopifnot(sum(gsam.rna.metadata$sid != rownames(gsam.count.stats)) == 0)


# ---- fig: unique/dup stats ----

tmp.1 <- gsam.multimap.stats[,1:2]
tmp.1$type <- 'Uniquely mapped'
colnames(tmp.1)[2] <- 'Percentage'
tmp.1$uniquely.mapped <-  gsam.multimap.stats$percentage.uniquely.mapped

tmp.2 <- gsam.multimap.stats[,c(1,3)]
tmp.2$type <- 'Multi-mapped'
colnames(tmp.2)[2] <- 'Percentage'
tmp.2$uniquely.mapped <-  gsam.multimap.stats$percentage.uniquely.mapped

plt <- rbind(tmp.1, tmp.2)
print(dim(plt))
rm(tmp.1, tmp.2)



ggplot(plt , aes(x = reorder(sample, -uniquely.mapped) ,y = Percentage,col=type,fill=type)) +
  geom_bar(stat="identity") +
  labs(x = NULL, y = "Percentage of input reads") +
  ylim(0, 100) + youri_gg_theme


ggsave("output/figures/qc/mapping_stats.png", width=8, height=6.5)

# ---- duplicate reads stats ----


rl <- sort(unique(gsam.rna.seqruns$runlanes))



plt <- data.frame(sid = as.character(gsam.rna.metadata$sid))
print(dim(plt))

plt <- merge(plt,
             gsam.rna.seqruns,
             by.x = 'sid',
             by.y = 'sid')
print(dim(plt))

plt <- merge(plt,
             gsam.count.stats,
             by.x = 'sid',
             by.y = 'sid')
print(dim(plt))


plt <- merge(plt,
             gsam.rna.metadata,
             by.x = 'sid',
             by.y = 'sid')
print(dim(plt))



# rl[4] = slecht
# rl[7] & rl[9] = best slecht
# 5 & 10 & 11 = goed
plt$runlanes2 <- as.character(plt$runlanes)
plt[(plt$runlanes2 %in% c(rl[4], rl[7],  rl[9],rl[12])) == FALSE,]$runlanes2 <- "other"
#plt[(plt$runlanes2 %in% c(rl[12])) == FALSE,]$runlanes2 <- "other"

ggplot(plt, aes(x = pct.duplicate.STAR, y=pct.Assigned , shape = blacklist.pca, col=storage.box)) +
  geom_point()


ggplot(subset(plt, runs == "BHM33WDSXX" ), aes(x = pct.duplicate.STAR, y=pct.Assigned , shape = blacklist.pca, col=runlanes)) +
  geom_point()

ggplot(subset(plt, runs == "AHKK5WDSXX,BHJTTYDSXX" ), aes(x = pct.duplicate.STAR, y=pct.Assigned , shape = blacklist.pca, col=runlanes)) +
  geom_point()

ggplot(subset(plt, runs == "AHLGTCDSXX" ), aes(x = pct.duplicate.STAR, y=pct.Assigned , shape = blacklist.pca, col=runlanes)) +
  geom_point()

ggplot(subset(plt, runs == "AHKK5WDSXX" ), aes(x = pct.duplicate.STAR, y=pct.Assigned , shape = blacklist.pca, col=runlanes)) +
  geom_point()

ggplot(subset(plt, runs == "AHM2NJDSXX" ), aes(x = pct.duplicate.STAR, y=pct.Assigned , shape = blacklist.pca, col=runlanes)) +
  geom_point()



ggplot(plt, aes(x = mean.r1.len, y=pct.Assigned , shape = blacklist.pca, col=runlanes)) +
  geom_point()

ggplot(plt, aes(x = mean.r1.len, y=RMSE , shape = blacklist.pca, shape=blacklist.pca, label=sid)) +
  geom_point() + 
  geom_text_repel(data = subset(plt, mean.r1.len < 105 | RMSE > 4.5))

ggplot(plt, aes(x = mean.r1.len, y=RMSE , shape = blacklist.pca, col=runlanes2)) +
  geom_point()



# ---- barplot ----


plt <- data.frame(sid = as.character(gsam.rna.metadata$sid))
print(dim(plt))

plt <- merge(plt,
             gsam.rna.seqruns,
             by.x = 'sid',
             by.y = 'sid')
print(dim(plt))

plt <- merge(plt,
             gsam.count.stats,
             by.x = 'sid',
             by.y = 'sid')
print(dim(plt))


plt <- merge(plt,
             gsam.rna.metadata,
             by.x = 'sid',
             by.y = 'sid')
print(dim(plt))



plt.1 <- aggregate(blacklist.pca == T ~ runlanes, plt, sum)
plt.1$type <- "is blacklisted"
colnames(plt.1)[2] <- "frequency"
plt.2 <- aggregate(blacklist.pca == F ~ runlanes, plt, sum)
plt.2$type <- "is not blacklisted"
colnames(plt.2)[2] <- "frequency"
plt <- rbind(plt.1, plt.2)


ggplot(plt, aes(x = runlanes, y = frequency, col = type)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





