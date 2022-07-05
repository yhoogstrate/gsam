#!/usr/bin/env R

setwd("/home/youri/projects/gsam")

# ---- libs ----
library(DESeq2)

# ---- load data ----
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")

# ---- analysis ----


d <- expression_matrix

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

tmp.metadata <- gsam.metadata[ match(colnames(f.vst), gsam.metadata$sample), ]
tmp.metadata$v3 <- '?'
tmp.metadata[ tmp.metadata$vIII.reads.v3 / (tmp.metadata$vIII.reads.v3 + tmp.metadata$wt.reads.v3) > 0.01 ,]$v3 <- 'v3+'
tmp.metadata[ tmp.metadata$vIII.reads.v3 / (tmp.metadata$vIII.reads.v3 + tmp.metadata$wt.reads.v3) <= 0.01 ,]$v3 <- 'v3-'


colnames(f.vst) <- paste0 ( colnames(f.vst) , "." , tmp.metadata$v3 )
cond <- as.factor(tmp.metadata$v3)

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



#  PCA
library(vegan)
library(rgl)
library(ape)

# actual PCoA analysis
dds.pcoa=pcoa(vegdist(t(f.vst),method="manhattan")/1000)
scores=dds.pcoa$vectors

# plotting
plot(scores[,1], scores[,2],col=as.numeric(as.factor(cond)))
ordispider(scores,cond,label=T)
ordiellipse(scores,cond)

# interactive 3d plot - can rotate it by dragging mouse
radiusScale=2 # change it if the spheres come out too small or too big in the next one
plot3d(scores[,1], scores[,2],scores[,3],col=as.numeric(cond),type="s",radius=radiusScale*as.numeric(cond))

# formal permutation-based analysis of variance 
#adonis(t(vsd)~factor1*factor2,data=conditions,method="manhattan")  

rownames(tmp.metadata) <- paste0(tmp.metadata$sample, "." , tmp.metadata$v3)
#library(ggfortify)
#pc <- prcomp(t(f.vst))
#autoplot(pc, label=T,pos=3,data=tmp.metadata, colour="v3")

compx <- 1
compy <- 2
plot(pc$x[,compx], pc$x[,compy],type="n")

for(patient in sort(unique(tmp.metadata$pid))) {
  print(patient)
  r1 <- paste0(patient, "1")
  r1 <- pc$x[ match(r1, gsub("\\..+","",rownames(pc$x))) ,]

  r2 <- paste0(patient, "2")
  r2 <- pc$x[ match(r2, gsub("\\..+","",rownames(pc$x))) ,]
  
  x1 <- r1[compx]
  y1 <- r1[compy]
    
  x2 <- r2[compx]
  y2 <- r2[compy]
  
  lines(c(x1, x2), c(y1, y2))

}

points(pc$x[,compx], pc$x[,compy],col=as.factor(tmp.metadata$v3),pch=19,cex=0.8)
text(pc$x[,compx], pc$x[,compy],tmp.metadata$sample,cex=0.7,pos=3)




# pca 500 most DE genes, v3 samples cluster together, the rest is different?

ntop <- 500
variances <- rowVars(f.vst)
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- f.vst[select,]
pc <- prcomp(t(high_variance_genes), retx=T)




source('scripts/R/ensembl_to_geneid.R')


# DE vIII - non vIII
cond <- as.factor(tmp.metadata$v3)
cond <- gsub("+",".pos",as.character(cond),fixed=T)
cond <- gsub("-",".neg",as.character(cond),fixed=T)
cond <- as.factor(cond)
dds2 <- DESeqDataSetFromMatrix(f, DataFrame(cond), ~cond)
dds2 <- DESeq(dds2)
res <- results(dds2)
res <- res[order(res$padj),]
sig <- res[res$padj < 0.01,]
sig[rownames(sig) %in% tt1,]

# lncRNA EGFR = DE
sig$ensid <- ens_ids[ match(gsub("\\..+$","",rownames(sig)) , gsub("\\..+$","",ens_ids$ens.id)) ,]$gene.id



# plot v3 percentage against:
# TBX5 ~ ENSG00000089225
sel <- gsub("\\..+$","",rownames(f.vst)) == "ENSG00000089225"
plot(f.vst[sel,] , tmp.metadata$vIII.reads.v3 / (tmp.metadata$vIII.reads.v3 + tmp.metadata$wt.reads.v3) , xlab="TBX5 expr", ylab="% vIII")

# CDK4 ~ ENSG00000135446
sel <- gsub("\\..+$","",rownames(f.vst)) == "ENSG00000135446"
plot(f.vst[sel,] , tmp.metadata$vIII.reads.v3 / (tmp.metadata$vIII.reads.v3 + tmp.metadata$wt.reads.v3) , xlab="CDK4 expr", ylab="% vIII")
text(f.vst[sel,] , tmp.metadata$vIII.reads.v3 / (tmp.metadata$vIII.reads.v3 + tmp.metadata$wt.reads.v3) , tmp.metadata$sample, cex=0.85,srt=90,pos=3)


# SEC61G ~ ENSG00000132432
sel <- gsub("\\..+$","",rownames(f.vst)) == "ENSG00000132432"
plot(f.vst[sel,] , tmp.metadata$vIII.reads.v3 / (tmp.metadata$vIII.reads.v3 + tmp.metadata$wt.reads.v3) , xlab="SEC61G [chr7] expr", ylab="% vIII")





# ----
cond <- as.factor(tmp.metadata$resection)
cond <- gsub("+",".pos",as.character(cond),fixed=T)
cond <- gsub("-",".neg",as.character(cond),fixed=T)
cond <- as.factor(cond)
dds2 <- DESeqDataSetFromMatrix(f, DataFrame(cond), ~cond)
dds2 <- DESeq(dds2)
res <- results(dds2)
res <- res[order(res$padj),]
sig <- res[res$padj < 0.05,]
sig[rownames(sig) %in% tt1,]

# lncRNA EGFR = DE
sig$ensid <- ens_ids[ match(gsub("\\..+$","",rownames(sig)) , gsub("\\..+$","",ens_ids$ens.id)) ,]$gene.id


# plot v3 percentage against:
# TBX5 ~ ENSG00000089225
sel <- gsub("\\..+$","",rownames(f.vst)) == "ENSG00000238217"
plot(as.factor(tmp.metadata$resection) , f.vst[sel,]  , xlab="resection", ylab="LINC01877|chr2:200M")


