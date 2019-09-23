#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)

# ---- load data ----
source("scripts/R/ligands.R")

d <- read.delim("output/tables/featureCounts_gsam_1st96.exon-level.txt",stringsAsFactors = F,comment="#")

colnames(d) <- gsub("X.data.users.youri.mnt.neurogen.ro.gsam.RNA.alignments.","",colnames(d),fixed=T)
colnames(d) <- gsub(".sambamba.dedup.bam","",colnames(d),fixed=T)

rownames(d) <- d$Geneid
d$Geneid <- NULL
d$Chr <- NULL
d$Start <- NULL
d$End <- NULL
d$Strand <- NULL
d$Length <- NULL

#plot(sort(log(colSums(d))))
#abline(h=log(2000000))
# min 2mljn reads to take the low-Q samples out
d <- d[,colSums(d) > 2000000]
#plot(sort(log(colSums(d))))
#abline(h=log(2000000))


# separate resections
sid <- gsub("^[^_]+_([A-Z]+[0-9]).+$","\\1",colnames(d))
pid <- gsub("^[^_]+_([A-Z]+)[0-9]+.+$","\\1",colnames(d))
res <- gsub("^[^_]+_[^0-9]+","",colnames(d))
res <- gsub("_.+$","",res )

shared <- intersect( pid[res == "1"] , pid[res == "2"] )
shared <- sort(shared)


e <- d[, match(c(paste0(shared, "1"), paste0(shared, "2")), sid) ]
cond <- as.factor( c(rep("1", length(shared)) , rep("2", length(shared))))
rm(pid)




f <- e[ rowSums(e) >= ncol(e) * 2.5 ,]
f.vst <- DESeqDataSetFromMatrix(f, DataFrame(cond), ~cond)
f.vst <- assay(vst(f.vst))
colnames(f.vst) <-  gsub("^[^_]+_([A-Z]+[0-9]).+$","\\1",colnames(f.vst))

f.vst.1 <- f.vst[,1:sum(cond == "1")]
f.vst.2 <- f.vst[,(sum(cond == "1") + 1):ncol(f.vst)]
rm(f.vst)


#order_ligands <- match(c("ENSG00000146648", "ENSG00000109321", "ENSG00000124882", "ENSG00000182585", "ENSG00000163235", "ENSG00000138798", "ENSG00000113070", "ENSG00000174808")   ,   gsub("\\..+$","",rownames(d)))


off <- 0.05

for(ligand in names(tt2)) {
  gene <- tt2[ligand]
  ens <- gsub("\\..+$","",ligand)
  
  print(as.character(gene))
  
  dd <- match(ens , gsub("\\..+$","",rownames(d)) )
  g <- match(ens , gsub("\\..+$","",rownames(f.vst.1)) )

  if(is.na(g)) {
    print(paste0("IS NA: ", gene))
  }
  else {
    png(filename=paste0("output/figures/expression_",gene,"_over_time.png"), width = 480*2,height=480*2, res=72*2)
    
    ymin <- min( f.vst.1[g,] , f.vst.2[g,]   )
    ymax <- max( f.vst.1[g,] , f.vst.2[g,]   )
    plot(c(1.0 - off,2.0 + off),c(ymin,ymax),type="n", ylab=paste0("Expression ",as.character(gene)), xlab="Resection")
    
    for(i in 1:ncol(f.vst.1)) {
      lines(c(1,2),c( f.vst.1[g,i] , f.vst.2[g,i] ) )
      
      text(1, f.vst.1[g,i], shared[i],cex=0.65,pos=2,col="darkgray")
      text(2, f.vst.2[g,i], shared[i],cex=0.65,pos=4,col="darkgray")
    }
    
    dev.off()
  }
}




library(pheatmap)
pheatmap(cor(f.vst))












