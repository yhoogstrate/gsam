#!/usr/bin/env R

cutoff <- 0.1 # cutoff in percentage for a loci to include in plot - denote that there are quite some rRNA genes on alternate loci

idxstats <- list.files('output/tables/idxstats', pattern = '.txt', full.names = T)


idx <- NA
df <- data.frame()
for(idxstat in idxstats) {
  t <- read.delim(idxstat,header=F,stringsAsFactors = F)
  colnames(t)[3]  <- gsub(".+/([^/]+).flagstats.txt","\\1",idxstat)
  t$V4 <- NULL
  colnames(t)[1]  <- "chr.name"
  colnames(t)[2]  <- "chr.len"
  
  if(is.null(idx) | ncol(df) == 0) {
    idx <- t
    idx[,3] <- NULL
    
    df <- t
  } else {
    t$chr.name <- NULL
    t$chr.len <- NULL
    df <- cbind(df,t)
  }
  
  #rm(t)
}
df$chr.name <- NULL
df$chr.len <- NULL


df <- t(df)
df <- df / rowSums(df) * 100.0
df <- t(df)


# select those that vary enough
m <- rowMax(df)
o <- order(m)
#plot(c(0,max(o)+3),c(0,120),type="n")
#points(1:max(o), m[o],pch=19,col="red")
s = m[o] > 0.8
#points((1:max(o))[s], m[o][s],pch=19,col="blue")
#text(1:max(o) - 1, m[o] + 4, idx$chr.name[o],srt=90,pos=4,cex=0.6)


sel <- idx$chr.name[o][s]

idx <- idx[idx$chr.name %in% sel,]
df <- df[idx$chr.name %in% sel,]


