#!/usr/bin/env R

setwd("~/projects/gsam")

# there is a difference in size of results files of CNVKit
# they either have 27493 lines (all in ExtraSequencing_CN_data_29Mar2018/)
# or have 27287 lines (195/374) = 50%

# @todo: synch with neuro-genomic and change data directory

d <- read.table("output/tables/cnv_copynumber-ratio.cnr_all.txt",stringsAsFactors = F,header=T)
naPerRow <- apply(d, 1, function(x) sum(is.na(x)))
x <- nrow(d)
# there are four samples of which the CNV data is missing at this moment, and are always NA
d <- d[naPerRow <= 4,]
print(paste0("Removal of ",x - nrow(d) , " CNV regions due to batch effects in files" ))

genes <- d[,1:4]
d <- d[,-(1:4)]

dd <-

