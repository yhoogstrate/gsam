#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

d <- read.delim("output/tables/duplication-stats.txt",stringsAsFactors = F)

e <- d[order(d[,2]),]