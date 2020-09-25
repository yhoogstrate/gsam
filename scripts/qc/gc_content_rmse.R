#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load libs ----

library(tidyverse)

# ---- load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')



gsam.gc.content <- read.delim("output/tables/qc/gc_content_rmse.txt",stringsAsFactors = F) %>%
  dplyr::arrange(RMSE, sample.id) %>%
  dplyr::mutate(filename = factor(filename , levels=filename)) %>% 
  dplyr::mutate(new.batch = grepl("new", sample.id)) %>%
  dplyr::mutate(name.strip = ifelse(new.batch == T ,gsub("-new", "", .$sample.id, fixed=T),"") ) %>%
  dplyr::mutate(old.batch = sample.id %in% name.strip) %>%
  dplyr::mutate(name.strip = NULL) %>%
  dplyr::mutate(batch = dplyr::case_when ((old.batch == T & new.batch == F) ~ "old" , (old.batch == F & new.batch == T) ~ "new" , (old.batch == F & new.batch == F) ~ "single" )) %>%
  dplyr::mutate(batch = as.factor(batch))
  


# mediane waarden in niet outlier samples:
# a     c     t     g
# 25.11 24.61 25.61 24.67



# ---- make plot (ggplot) ----

tmp.A <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.A","RMSE", "batch")]
tmp.A$type <- "A"
colnames(tmp.A)[3] <- "Percentage"

tmp.C <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.C","RMSE", "batch")]
tmp.C$type <- "C"
colnames(tmp.C)[3] <- "Percentage"

tmp.T <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.T","RMSE", "batch")]
tmp.T$type <- "T"
colnames(tmp.T)[3] <- "Percentage"

tmp.G <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","percentage.G","RMSE", "batch")]
tmp.G$type <- "G"
colnames(tmp.G)[3] <- "Percentage"

tmp.RMSE <- gsam.gc.content[,colnames(gsam.gc.content) %in% c("sample.id","filename","RMSE", "batch")]
tmp.RMSE$type <- "RMSE"
colnames(tmp.RMSE)[3] <- "Percentage"
tmp.RMSE$RMSE <-tmp.RMSE$Percentage

plt <- rbind(tmp.A, tmp.C, tmp.T, tmp.G, tmp.RMSE) %>%
  dplyr::mutate(order = reorder(sample.id, RMSE))

# , alpha=batch
ggplot(plt , aes(x = order ,y = Percentage, group=sample.id, col=type)) +
  #geom_segment(aes(y=0, yend=50, x = order, xend=order), data=plt %>% dplyr::filter(batch == "new") %>% dplyr::mutate(order = factor(order, levels = levels(plt$order)))) + 
  geom_point(pch="-",size=3) +
  #geom_vline( aes(xintercept = order), data=plt %>% dplyr::filter(batch == "new" ) , alpha=0.1 ) + 
  geom_vline( aes(xintercept = order), data=plt %>% dplyr::filter(batch == "old" ) , alpha=0.1 , col="red") + 
  labs(x="GSAM RNA-Seq sample", y="Percentage",col="Percentage of:") + 
  ylim(0, 50) +
  theme( axis.text.x = element_text(angle = 90, size = 5 )) +
  job_gg_theme +
  scale_alpha_manual(NULL, values=c('new'=1.0, 'old'=0.45, 'single'=0.1))



ggsave("output/figures/qc/gc_conent_rmse.pdf", width=12,height=6)
ggsave("output/figures/qc/gc_conent_rmse.png", width=12 * 0.75,height=6)




# ---- plot (base R) ----

e <- gsam.gc.content
png("output/figures/gc_content_rmse.png", width=480*4,heigh=480*2,res=72*2)
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

dev.off()

