#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load libs ----



# ---- load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/youri_gg_theme.R')



gsam.gc.content <- read.delim("output/tables/gc_content_rmse.txt",stringsAsFactors = F)
gsam.gc.content <- gsam.gc.content[order(gsam.gc.content$RMSE, tmp$sample.id),]
gsam.gc.content$filename <- factor(gsam.gc.content$filename , levels=gsam.gc.content$filename)

# mediane waarden in niet outlier samples:
# a     c     t     g
# 25.11 24.61 25.61 24.67



# ---- make plot ----

tmp.A <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.A","RMSE")]
tmp.A$type <- "A"
colnames(tmp.A)[3] <- "Percentage"

tmp.C <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.C","RMSE")]
tmp.C$type <- "C"
colnames(tmp.C)[3] <- "Percentage"

tmp.T <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.T","RMSE")]
tmp.T$type <- "T"
colnames(tmp.T)[3] <- "Percentage"

tmp.G <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.G","RMSE")]
tmp.G$type <- "G"
colnames(tmp.G)[3] <- "Percentage"

tmp.RMSE <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","RMSE")]
tmp.RMSE$type <- "RMSE"
colnames(tmp.RMSE)[3] <- "Percentage"
tmp.RMSE$RMSE <-tmp.RMSE$Percentage

plt <- rbind(tmp.A, tmp.C, tmp.T, tmp.G, tmp.RMSE)

ggplot(plt , aes(x = reorder(sample.id, RMSE) ,y = Percentage, group=filename, col=type)) +
  geom_point(pch="-",size=3) +
  labs(x="GSAM RNA-Seq sample", y="Percentage",col="Percentage of:") + 
  ylim(0, 50) + theme( axis.text.x = element_text(angle = 90, size = 5 )) + job_gg_theme


ggsave("output/figures/qc/gc_conent_rmse.pdf", width=12,height=6)
ggsave("output/figures/qc/gc_conent_rmse.png", width=12 * 0.75,height=6)


# ---- ----

e <- tmp
#png("output/figures/gc_content_rmse.png", width=480*4,heigh=480*2,res=72*2)
plot(c(1,nrow(e)), c(0,50), type="n", ylab="",xlab="Sample nr. / order",main="GC deviation")

for(i in 1:nrow(e)) {
  if(i %% 2 == 0) {
    rect(i - 0.5 , 0 , i + 0.5 ,100 , col=rgb(0.95,0.95,0.95,0.55) , border=NA)
  }
}

abline(h=25,col=rgb(0.85, 0.85, 0.85), lty=2)

points(1:nrow(e), e$RMSE , pch=19 , cex=0.36)


points(1:nrow(e) , e$percentage.A, col=rgb(0,0.5,0.0),pch="-", cex=1.5)
points(1:nrow(e) , e$percentage.C, col=rgb(0.5,0.5,0),pch="-", cex=1.5)
points(1:nrow(e) , e$percentage.T, col=rgb(0.5,0,0.5),pch="-", cex=1.5)
points(1:nrow(e) , e$percentage.G, col="blue",pch="-", cex=1.5)

text(1:nrow(e) - 0.8,  45, e$sample.id , pos=4, srt=90, cex=0.75)

legend(0, 20, pch=c("-","-","-","-","o"), col=c(rgb(0,0.5,0.0),rgb(0.5,0.5,0),rgb(0.5,0,0.5),"blue","black"), c("%A","%C","%T","%G","RMSE from 25%") )

#dev.off()

