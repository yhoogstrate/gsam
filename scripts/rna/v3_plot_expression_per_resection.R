#!/usr/bin/env R


setwd("~/projects/gsam")


# ---- load: libs ----

library(ggplot2)
library(ggrepel)
library(cowplot)

source("scripts/R/job_gg_theme.R")

# ---- load: functions ----

# ---- load: data ----

source("scripts/R/vIII.R")
vIII.rot$sid <- rownames(vIII.rot)

tmp <- read.csv('data/RNA/Final_qPCR_EGFR_GSAM.csv',stringsAsFactors = F)
tmp$recurrent_patientSampleID <- NULL # intuitive
tmp$recurrent_sampleID <- NULL # intuitive
tmp$patientSampleID <- NULL # intuitive
tmp$X <- NULL
tmp$type <- NULL
tmp$recurrent_type <- NULL
colnames(tmp) <- paste0("qPCR.",colnames(tmp))
tmp$qPCR.delta_percentage <- tmp$qPCR.recurrent_percentageEGFRvIII - tmp$qPCR.percentageEGFRvIII

vIII.rot <- merge(x=vIII.rot, y = tmp, by.x="sid" , by.y = "qPCR.recurrent_patientID", all.x = T)
rm(tmp)


# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")

source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)



# ---- plot v3 x resection [linear] ----

png("output/figures/rna/v3_plot_expression_per_resection.png", width=480*3,height=480*3, res=72*3)


off <- 0.05
plot(c(1.0 - off,2.0 + off), c(1,100),type="n", ylab="%vIII",xlab="resection",main=paste0("vIII over resection/time in ",nrow(e)," matching samples"))

for(i in 1:nrow(e)) {
  p1 <- 100.0 *   e[i,4] / (e[i,3] + e[i,4])
  p2 <- 100.0 *   e[i,7] / (e[i,6] + e[i,7])
  
  #p1 <- p1 / 100 * 99 + 1
  #p2 <- p2 / 100 * 99 + 1
  
  points(1, p1, pch=19,cex=0.6)
  points(2, p2, pch=19,cex=0.6)
  lines(c(1,2), c(p1, p2))
  
  text(1,p1, e[i,1],cex=0.65,pos=2,col="darkgray")
  text(2,p2, e[i,1],cex=0.65,pos=4,col="darkgray")
}



dev.off()

# ---- plot v3 x resection [log2] ----

png("output/figures/rna/v3_plot_expression_per_resection_log.png", width=480*3,height=480*3, res=72*3)

off <- 0.05
plot(c(1.0 - off,2.0 + off), c(1,100),type="n", ylab="%vIII",xlab="resection",main=paste0("vIII over resection/time in ",nrow(e)," matching samples"),log="y")

for(i in 1:nrow(e)) {
  p1 <- 100.0 *   e[i,4] / (e[i,3] + e[i,4])
  p2 <- 100.0 *   e[i,7] / (e[i,6] + e[i,7])
  
  p1 <- p1 / 100 * 99 + 1
  p2 <- p2 / 100 * 99 + 1
  
  points(1, p1, pch=19,cex=0.6)
  points(2, p2, pch=19,cex=0.6)
  lines(c(1,2), c(p1, p2))
  
  text(1,p1, e[i,1],cex=0.65,pos=2,col="darkgray")
  text(2,p2, e[i,1],cex=0.65,pos=4,col="darkgray")
}

dev.off()


# ---- plot v3 x resection [log2 delta/ratio] ----



outliers1 <- subset(tmp, 
                    delta.percentage > 10)
outliers2 <- subset(tmp, 
                    delta.percentage < -11)
ggplot(tmp, aes(x=x, y=delta.percentage, label=sid)) + 
  geom_bar(aes(x=x, y=delta.percentage),stat="identity",width=0.7) + 
  geom_point(aes(col=v3.stat)) +
  ylim(-100, 100) + 

  labs(title = "GSAM: Change in percentage EGFR vIII reads (RNA-seq)",
       subtitle=paste0(format(Sys.time(), "[%d-%b-%y]")),
       y = "Change in percentage vIII",
       x = "GSAM patient (ordered on difference)") + 
  job_gg_theme

ggsave("output/figures/rna/GSAM_percentage_EGFRvIII_barplot.png")




outliers1 <- subset(tmp, 
                    delta.percentage > 10)
outliers2 <- subset(tmp, 
                    delta.percentage < -11)
ggplot(tmp, aes(x=x, y=delta.percentage, label=sid)) + 
  geom_bar(aes(x=x, y=delta.percentage),stat="identity",width=0.7) + 
  geom_point(aes(col=v3.stat)) +
  ylim(-100, 100) + 
    geom_text_repel(
    nudge_y       = 100 - outliers1$delta.percentage,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x",
    size=0.7*5,
    data = outliers1) + 
  geom_text_repel(
    nudge_y       = -100 - outliers2$delta.percentage,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x",
    size=0.7*5,
    data = outliers2) +
  labs(title = "GSAM: Change in percentage EGFR vIII reads (RNA-seq)",
       subtitle=paste0(format(Sys.time(), "[%d-%b-%y]")),
       y = "Change in percentage vIII",
       x = "GSAM patient (ordered on difference)") + 
  job_gg_theme

ggsave("output/figures/rna/GSAM_percentage_EGFRvIII_barplot_labels.png")



# ----- barplot style ----




tmp <- vIII.rot$resection.1.sum >= 15  & vIII.rot$resection.2.sum >= 15 & !is.na(vIII.rot$resection.1.sum) & !is.na(vIII.rot$resection.2.sum)
tmp <- vIII.rot[tmp,]
tmp <- tmp[order(tmp$delta.percentage, rownames(tmp)),]
tmp$x <- order(tmp$delta.percentage, rownames(tmp))
tmp$sid <- rownames(tmp)


outliers1 <- subset(tmp, 
                    delta.percentage > 10)
outliers2 <- subset(tmp, 
                    delta.percentage < -11)

 

plot_grid(

    ggplot(tmp, aes(x=x,y=delta.percentage, label=sid) ) + 
      geom_rect(
        fill = "gray90",
        aes(xmin=x - 0.33,
            xmax=x + 0.33 ,
            ymin=0,
            ymax=100),
        data=subset(tmp, is.na(tmp$qPCR.delta_percentage))
        ) +
    geom_curve( curvature = 0,
                arrow = arrow(length = unit(0.01, "npc")),
                aes(x = x , y = qPCR.percentageEGFRvIII , 
                    xend = x, yend = qPCR.recurrent_percentageEGFRvIII),
                data = subset(tmp, tmp$qPCR.percentageEGFRvIII != 0 | tmp$qPCR.recurrent_percentageEGFRvIII != 0)
                ) +
    ylim(0, 100) +
    geom_point(
      pch = '-',
      size = 4,
      data = subset(tmp, tmp$qPCR.percentageEGFRvIII != 0 | tmp$qPCR.recurrent_percentageEGFRvIII != 0),
      aes(y=qPCR.percentageEGFRvIII)
    ) +
    geom_point(
      #pch = '-',
      size = 1,
      aes(y=qPCR.percentageEGFRvIII),
      col="red",
      data = subset(tmp, tmp$qPCR.percentageEGFRvIII == 0 & tmp$qPCR.recurrent_percentageEGFRvIII == 0)
    ) + 
    labs(title = "GSAM: Percentage EGFR vIII (qPCR table)",
         y = "Percentage vIII (qPCR)",x=NULL) + 
    scale_colour_manual(name = 'Legend', 
                        guide = 'legend',
                        values = c('MA50' = 'red',
                                   'MA200' = 'blue'), 
                        labels = c('SMA(50)',
                                   'SMA(200)')) +
    job_gg_theme
  ,
  ggplot(tmp, aes(x=x,y=delta.percentage, label=sid) ) + 
    geom_curve( curvature = 0,
                arrow = arrow(length = unit(0.01, "npc")),
                aes(x = x , y = resection.1.p , xend = x, yend = resection.2.p),
                data = subset(tmp, tmp$resection.1.p != 0 | tmp$resection.2.p != 0 )
                ) +
    ylim(0, 100) +
    geom_point(
      pch = '-',
      size = 4,
      data = subset(tmp, tmp$resection.1.p != 0 | tmp$resection.2.p != 0),
      aes(y=resection.1.p)
    ) +
    geom_point(
      #pch = '-',
      size = 1,
      aes(y=resection.1.p),
      col="red",
      data = subset(tmp, tmp$resection.1.p == 0 & tmp$resection.2.p == 0  )
    ) + 
    labs(y = "Percentage vIII (RNA-seq)",
         x = "GSAM patient (ordered on difference in RNA-seq)") + 
    scale_colour_manual(name = 'Legend', 
                        guide = 'legend',
                        values = c('MA50' = 'red',
                                   'MA200' = 'blue'), 
                        labels = c('SMA(50)',
                                   'SMA(200)')) +
    job_gg_theme      ,
  align="hv", axis="tblr",ncol=1 )

ggsave("output/figures/rna/GSAM_percentage_EGFRvIII_arrowplot.png")


  #  geom_text_repel(
    #    aes(y = resection.2.p),
    #    nudge_y       = 100 - outliers1$resection.1.p,
    #    segment.size  = 0.2,
    #    segment.color = "grey50",
    #    #    direction     = "x",
    #    size=0.7*5,
    #    data = outliers1) + 
  #  geom_text_repel(
    #    aes(y = resection.2.p),
    #    nudge_y       = - outliers2$resection.1.p,
    #    segment.size  = 0.2,
    #segment.color = "grey50",
    #    direction     = "x",
    #    size=0.7*5,
#    data = outliers2)







# -----

# find LR of EGFR 
e <- expression_matrix_full
cond <- factor(paste0("r",gsub("^.+([0-9])$","\\1",colnames(e))))
dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds,blind=T))
e.vst <- e.vst[rownames(e.vst) %in% c("ENSG00000146648.18_5","ENSG00000124882.4_3","ENSG00000109321.11_4","ENSG00000138798.12_4"),]

e.vst.r1 <- e.vst[,gsub("^.+([0-9])$","\\1",colnames(e.vst)) == "1"]
e.vst.r2 <- e.vst[,gsub("^.+([0-9])$","\\1",colnames(e.vst)) == "2"]
shared <- intersect(gsub("[0-9]$","",colnames(e.vst.r1)), gsub("[0-9]$","",colnames(e.vst.r2)))

e.vst.r1 <- e.vst.r1[, match(shared, gsub("[0-9]$","",colnames(e.vst.r1)))]
e.vst.r2 <- e.vst.r2[, match(shared, gsub("[0-9]$","",colnames(e.vst.r2)))]

# ADD EGFR
lr_egfr <- log2(e.vst.r2[rownames(e.vst.r2) == "ENSG00000146648.18_5",]  / e.vst.r1[rownames(e.vst.r1) == "ENSG00000146648.18_5",])
lr_egfr <- data.frame(lr_egfr)
lr_egfr$sid <-  gsub("[0-9]$","",rownames(lr_egfr) )
tmp <- merge(tmp, lr_egfr, by.x ="sid", by.y = "sid", all.x=TRUE)

# ADD EREG
lr_ereg <- log2(e.vst.r2[rownames(e.vst.r2) == "ENSG00000124882.4_3",]  / e.vst.r1[rownames(e.vst.r1) == "ENSG00000124882.4_3",])
lr_ereg <- data.frame(lr_ereg)
lr_ereg$sid <-  gsub("[0-9]$","",rownames(lr_ereg) )
tmp <- merge(tmp, lr_ereg, by.x ="sid", by.y = "sid", all.x=TRUE)

# ADD egf
lr_egf <- log2(e.vst.r2[rownames(e.vst.r2) == "ENSG00000138798.12_4",]  / e.vst.r1[rownames(e.vst.r1) == "ENSG00000138798.12_4",])
lr_egf <- data.frame(lr_egf)
lr_egf$sid <-  gsub("[0-9]$","",rownames(lr_egf) )
tmp <- merge(tmp, lr_egf, by.x ="sid", by.y = "sid", all.x=TRUE)


#tmp <- tmp[order(tmp$lr_egfr, rownames(tmp)),]
#tmp$x <- order(tmp$lr_egfr, rownames(tmp))

tmp <- tmp[order(tmp$delta.percentage, rownames(tmp)),]
tmp$x <- order(tmp$delta.percentage, rownames(tmp))

plot_grid(
  ggplot(tmp, aes(x=x, y=lr_egfr, label=sid)) + 
    geom_bar(aes(x=x, y=lr_egfr),stat="identity",width=0.7) + 
    geom_point(aes(col=v3.stat)) +
    ylim(-1, 1)
  ,
  ggplot(tmp, aes(x=x, y=lr_ereg, label=sid)) + 
    geom_bar(aes(x=x, y=lr_ereg),stat="identity",width=0.7) + 
    geom_point(aes(col=v3.stat)) +
    ylim(-0.5, 0.5)
  ,
  ggplot(tmp, aes(x=x, y=lr_egf, label=sid)) + 
    geom_bar(aes(x=x, y=lr_egf),stat="identity",width=0.7) + 
    geom_point(aes(col=v3.stat)) +
    ylim(-0.5, 0.5)
  ,
  ggplot(tmp, aes(x=x, y=delta.percentage, label=sid)) + 
    geom_bar(aes(x=x, y=delta.percentage),stat="identity",width=0.7) + 
    geom_point(aes(col=v3.stat)) +
    ylim(-100, 100)
  , 
  align="hv", axis="tblr",ncol=1
)


#ggsave("output/figures/rna/vIII_changes.png")
ggsave("output/figures/rna/egfr_changes.png")




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




