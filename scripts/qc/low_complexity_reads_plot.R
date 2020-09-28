#!/usr/bin/env R

# ---- load cfg -----

setwd("~/projects/gsam")

# ---- load libs ----

library(tidyverse)

source('scripts/R/youri_gg_theme.R')
source('scripts/R/job_gg_theme.R')

# ---- load data ----

plt <- read.delim("data/output/tables/qc/low_complexity_reads.txt", header=T, stringsAsFactors = F) %>%
  dplyr::mutate(percentage_numerical = as.numeric(gsub("%","", .$percentage.low.complexity.reads,fixed=T))) %>%
  dplyr::arrange(percentage_numerical) %>%
  dplyr::mutate(new.batch = grepl("new", sample)) %>%
  dplyr::mutate(name.strip = ifelse(new.batch == T ,gsub("-new", "", .$sample, fixed=T),"") ) %>%
  dplyr::mutate(old.batch = sample %in% name.strip) %>%
  dplyr::mutate(name.strip = NULL) %>%
  dplyr::mutate(batch = dplyr::case_when ((old.batch == T & new.batch == F) ~ "old" , (old.batch == F & new.batch == T) ~ "new" , (old.batch == F & new.batch == F) ~ "single" )) %>%
  dplyr::mutate(batch = as.factor(batch))

# ---- plot 1 base R ----

#png("output/figures/low_complexity_reads.png",width=480*3,height=1.7*480,res=72*1.8)
plot(c(1, max(nrow(e))) , c(0, max(e$percentage_numerical) + 1), main="Percentage low-complexity discarded reads",type="n", xlab="Sample nr. / order", ylab = "Percentage reads discarted due to low complexity",
     col = as.numeric(grepl("new", e$sample)) + 1)
for(i in 1:nrow(e)) {
  print(i)
  lines(c(i,i), c(0,e[i,]$percentage_numerical),col="darkgray",lty=1)
  text(i - 3,e[i,]$percentage_numerical + 0.25, e[i,]$sample, srt=90, pos=4,cex=0.6,
      col = as.numeric(grepl("new", e$sample)) + 1 )
}
points(1:nrow(e) , e$percentage_numerical, pch=19,cex=0.6, col = as.numeric(grepl("new", e$sample)) + 1 )
#dev.off()



# ---- plot 1 gg ----

ggplot(plt, aes(x=reorder(sample, percentage_numerical), fill = batch, y=percentage_numerical, label=sample)) +
  geom_point(col = "black",pch=21) +
  ggrepel::geom_text_repel(aes(col=batch),srt = 90, data=subset(plt, batch !="single"))
  #+ job_gg_theme


ggsave("output/figures/qc/low_complexity_reads_plot_gg.pdf",width=20,height=8)








