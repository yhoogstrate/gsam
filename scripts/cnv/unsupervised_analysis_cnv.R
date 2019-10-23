#!/usr/bin/env R

# there is a difference in size of results files of CNVKit
# they either have 27493 lines (all in ExtraSequencing_CN_data_29Mar2018/)
# or have 27287 lines (195/374) = 50%

# @todo: synch with neuro-genomic and change data directory

d <- read.table("output/tables/cnv_copynumber-ratio.cnr_all.txt")

#test_27493 <- read.table("data/DNA/cnvkit/first_series/PD36773c__tumour.cnr", header=T)
#test_27287 <- read.table("/home/yhoogstrate/mnt/neuro-genomic-ro/gsam/DNA/cnvkit/first_series/PD29177a2__tumour.cnr" ,header=T)
#t1 <- test_27493
#t2 <- test_27287
#t1$loc <- paste0(t1$chromosome,":",t1$start,"-",t1$end)
#t2$loc <- paste0(t2$chromosome,":",t2$start,"-",t2$end)
#isct <- intersect(t1$loc, t2$loc)
#length(isct)





