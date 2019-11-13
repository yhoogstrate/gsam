#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(ggplot)
library(tidyr)
library(scales)
library(cowplot)

# ---- load data ----
#source("scripts/R/ligands.R")
#source("scripts/R/gsam_metadata.R")

data.qc <- read.table("data/output/tables/gsam_featureCounts_readcounts.txt.summary",stringsAsFactors = F,comment="#",header=T,row.names = 1)
colnames(data.qc) <- gsub("processed.(.+).Aligned.sortedByCoord.out.markdup.bam","\\1",colnames(data.qc))

barplot_theme <- theme(
  text = element_text(family = 'Helvetica'),
  axis.text.x = element_text(angle = 45,size=10),
  axis.text.y = element_text(size=5),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.text = element_text(),
  axis.line = element_line(colour="black"),
  panel.grid.major.y = element_line(colour = 'grey90', linetype = 'dotted')
)

# ---- plot stats in millions, as is ----

t <- data.frame(t(data.qc))
t <- t[order(t$Assigned,t$Unassigned_Duplicate,decreasing=F),]
t <- t[,order(colSums(t),decreasing=T)]
t <- t[,colSums(t) > 0]
t$Sample <- as.factor(rownames(t))
t <- gather(t, type, count, -Sample)
t$count <- t$count / 1000000
t$Sample <- factor(as.character(t$Sample), 
          levels=unique(as.character(t$Sample)))

# define order of coloring
t$type <- factor(t$type, levels=c(
  "Unassigned_Chimera"    ,  "Unassigned_Ambiguity" ,    "Unassigned_MultiMapping", "Unassigned_Duplicate", "Unassigned_NoFeatures",   "Assigned"
   )
)


ggplot(t, aes(x = Sample ,y = count, fill=type, label=Sample)) +
  coord_flip() + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_y_continuous(labels = unit_format(unit = "M")) + 
  scale_fill_brewer(palette = "Set3") + 
  barplot_theme
ggsave("output/tables/qc/featureCounts.summary.mln.pdf",height=16*1.5,width=6*1.5)



# ---- plot stats as percentage ----

t <- data.frame(t(data.qc))
t <- t[,order(colSums(t),decreasing=T)]
t <- t[,colSums(t) > 0]
t <- t / rowSums(t) * 100
t <- t[order(t$Assigned,t$Unassigned_Duplicate,decreasing=F),]
t$Sample <- as.character(rownames(t))
t <- gather(t, type, count, -Sample)
t$Sample <- factor(as.character(t$Sample), levels=unique(as.character(t$Sample))) # relabel by order
t$type <- factor(t$type, levels=c("Unassigned_Chimera","Unassigned_Ambiguity","Unassigned_MultiMapping","Unassigned_Duplicate","Unassigned_NoFeatures","Assigned"))# define order of coloring

ggplot(t, aes(x = Sample ,y = count, fill=type, label=Sample)) +
  coord_flip() + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_y_continuous(labels = unit_format(unit = "%")) + 
  scale_fill_brewer(palette = "Set3") + 
  barplot_theme



ggsave("output/tables/qc/featureCounts.summary.per.pdf",height=16*1.5,width=6*1.5)


# ---- cowplot both ----

plotQC.featureCounts.summary.R <- function(featureCounts.DF) {
  t <- data.frame(t(featureCounts.DF))
  t <- t[order(t$Assigned,t$Unassigned_Duplicate,decreasing=F),]
  t <- t[,order(colSums(t),decreasing=T)]
  t <- t[,colSums(t) > 0]
  
  
  t1 <- t
  t1 <- t1 / 1000000
  t1$Sample <- as.factor(rownames(t1))
  t1 <- gather(t1, type, count, -Sample)
  t1$Sample <- factor(as.character(t1$Sample), levels=unique(as.character(t1$Sample))) # relabel by order
  t1$type <- factor(t1$type, levels=c("Unassigned_Chimera","Unassigned_Ambiguity","Unassigned_MultiMapping","Unassigned_Duplicate","Unassigned_NoFeatures","Assigned"))# define order of coloring
  
  t2 <- t
  t2 <- t2 / rowSums(t2) * 100
  t2$Sample <- as.factor(rownames(t2))
  t2 <- gather(t2, type, count, -Sample)
  t2$Sample <- factor(as.character(t2$Sample), levels=unique(as.character(t2$Sample))) # relabel by order
  t2$type <- factor(t2$type, levels=c("Unassigned_Chimera","Unassigned_Ambiguity","Unassigned_MultiMapping","Unassigned_Duplicate","Unassigned_NoFeatures","Assigned"))# define order of coloring
  
  
  plot_grid(
    ggplot(t1, aes(x = Sample ,y = count, fill=type, label=Sample)) +
      coord_flip() + 
      geom_bar(stat = "identity", position = "stack") + 
      scale_y_continuous(labels = unit_format(unit = "M")) + 
      scale_fill_brewer(palette = "Set3") + 
      barplot_theme
    ,
    ggplot(t2, aes(x = Sample ,y = count, fill=type, label=Sample)) +
      coord_flip() + 
      geom_bar(stat = "identity", position = "stack") + 
      scale_y_continuous(labels = unit_format(unit = "%")) + 
      scale_fill_brewer(palette = "Set3") +
      barplot_theme
    ,
    align="h"
  )
}

t <- data.qc
rownames(t) <- gsub(".replicate","r",rownames(t))
plotQC.featureCounts.summary.R(t)

ggsave("output/tables/qc/featureCounts.summary.pdf",height=16*1.5,width=6*1.5*2)

  

# ---- oude scripts ----

png("output/figures/featureCounts_stats.png",width=3*480,height=2*480,res=2*72)

sel <- order(colSums(d))
plot(  c(1, ncol(d)) , c(0, max(colSums(d)[sel] / 1000000)) ,xlab = "Sample nr. / order", ylab = "Million reads counted in exons", type="n",main="Exon counted reads")

#max
y <- max( colSums(d) ) / 1000000
lines(c(1, nrow(d)), c(y, y), col="darkgray", lty=2)
points(colSums(d)[sel] / 1000000,pch=19,cex=0.8)
text(ncol(d)/2, y, paste0("max=",round(y,2),"M"),pos=1,col="black")

#min
y <- min( colSums(d) ) / 1000000
lines(c(1, nrow(d)), c(y, y), col="darkgray", lty=2)
points(colSums(d)[sel] / 1000000,pch=19,cex=0.8)
text(ncol(d)/2, y, paste0("min=",round(y,2),"M"),pos=3,col="black")

text(    1:ncol(d),    (colSums(d)[sel] / 1000000) - 0.6,    gsub("^[^_]+_([^_]+)_.+$","\\1",colnames(d)) , srt=90 ,pch = 19 , cex=0.6, pos=1   )

dev.off()



