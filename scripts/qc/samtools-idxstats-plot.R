#!/usr/bin/env R

# ----config and settings ----

setwd("~/projects/gsam")

cutoff <- 0.1 # cutoff in percentage for a loci to include in plot - denote that there are quite some rRNA genes on alternate loci

# ---- load libs ----

library(ggplot2)
library(cowplot)

# ---- load functions ----

get_samtools_idxstats_rna <- function() {
  df <- read.table("output/tables/idxstats/samtools.indexstats.matrix.txt",header=T,stringsAsFactor=F)
  idx <- df[,1:2]
  df <- df[,-(1:2)]
  
  idx <- idx[rowSums(df) > 0,]
  df <- df[rowSums(df) > 0,]
  
  rownames(df) <- idx$ref
  colnames(df) <- gsub(".samtools.idxstats.txt","",colnames(df),fixed=T)
  
  return( list(idx= idx, data = df))
}

rowMax <- function(df) {
  return (apply(df, 1, FUN=max))
}


# ---- load data ----

rna.samtools.idxstats <- get_samtools_idxstats_rna()

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

source('scripts/R/gsam_metadata.R')
source('scripts/rna/dna-contamination.R')

# ---- make jamin-plot ----

# select those chromosomes that vary enough
rna.samtools.idxstats$data.pct <- t(t(rna.samtools.idxstats$data) / colSums(rna.samtools.idxstats$data) * 100.0)

colSums(rna.samtools.idxstats$data.pct[, 1:10])

m <- rowMax(rna.samtools.idxstats$data.pct)
o <- order(m)
#plot(c(0,max(o)+3),c(0,120),type="n")
#points(1:max(o), m[o],pch=19,col="red")
s = m[o] > 0.8
#points((1:max(o))[s], m[o][s],pch=19,col="blue")
#text(1:max(o) - 1, m[o] + 4, idx$chr.name[o],srt=90,pos=4,cex=0.6)


sel <- rna.samtools.idxstats$idx$ref[o][s]

idx2 <- rna.samtools.idxstats$idx[rna.samtools.idxstats$idx$ref %in% sel,]

df2 <- data.frame(rna.samtools.idxstats$data.pct[rna.samtools.idxstats$idx$ref %in% sel,])
df2 <- df2[,order(df2[rownames(df2) == "chrUn_gl000220",],decreasing=T)]
df2a <- data.frame(sid = gsub(".","-",as.character(colnames(df2)),fixed=T)) # to preserve order
df2$ref <- idx2$ref

# rotate / gather table
df2g <- gather(df2, sample, percentage, -ref)
df2g$sample <- factor(as.character(df2g$sample), levels=unique(as.character(df2g$sample))) # relabel by order


df2a$x <- 1:nrow(df2a)
df2a <- merge(df2a, gsam.rna.metadata, by.x="sid", by.y="sid")
df2a <- df2a[order(df2a$x),]
df2a$sid <- factor(as.character(df2a$sid[df2a$x]), levels = as.character(df2a$sid[df2a$x]))
#factor(as.character(df2a$sid[df2a$x]), levels = as.character(df2a$sid[df2a$x]))




plot_grid(
  ggplot(df2g, aes(x = sample ,y = percentage, fill=ref, label=sample)) +
  #coord_flip() + 
  geom_bar(stat = "identity", position = "stack",colour="black") + 
  scale_y_continuous(labels = unit_format(unit = "%")) + 
  theme_bw() + 
  barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 )
    )
,
ggplot(df2a, aes(x = sid ,y = ratio.stranded.antistranded.lin, fill=df2a$ratio.stranded.antistranded.dna, label=sid)) +
  #coord_flip() + 
  geom_bar(stat = "identity", position = "stack",colour="black") + 
  scale_y_continuous() + 
  theme_bw() + 
  barplot_theme + 
  theme( axis.title.y = element_text(size = 11) ,
         axis.text.x = element_text(angle = 90, size = 5 )
         )


,align="v", axis="tblr",ncol=1  , rel_heights = c(3, 1.6))

ggsave("output/figures/qc/samtools.idxstats.pdf",height=7*1.5,width=16*1.5)




