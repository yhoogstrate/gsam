#!/usr/bin/env R

#ensembl_to_geneid <- function() {
#  d <- read.table("/home/youri/bio/gxf/gencode.v29.annotation.gff3",sep="\t",comment.char="#",stringsAsFactors = F)
#  d$V2 <- NULL
#  d$V6 <- NULL
#  d$V7 <- NULL
#  d$V8 <- NULL
#  d <- d[d$V3 == "gene",]
#  d$gene.id <- gsub(';.+','',gsub('.+gene_name=','',d$V9))
#  d$ens.id <- gsub(';.+','',gsub('.+gene_id=','',d$V9))
#  d$V9 <- NULL
#  d$gene.id <- paste(d$gene.id,"|",d$V1,":",round(((d$V4 + d$V5) / 2)/1000000),"M",sep="")
#  
#  d$V1 <- NULL
#  d$V3 <- NULL
#  d$V4 <- NULL
#  d$V5 <- NULL
#  rownames(d) <- NULL
#  
#  ens_ids <- d
#  rm(d)
#  save(ens_ids, file="tmp/ens_ids.RData")
#}

load("tmp/ens_ids.RData")

