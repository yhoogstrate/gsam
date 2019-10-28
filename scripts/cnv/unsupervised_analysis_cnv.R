#!/usr/bin/env R


setwd("~/projects/gsam")

# ---- load libs ----

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(limma)
library(corrplot)


# ---- load functions ----


# ---- load data----

source("scripts/R/gsam_metadata.R")

source("scripts/R/cnv_matrix.R")
cnv_matrix_autosomes <- cnv_matrix[!cnv_matrix_genes$chr %in% c("chrX", "chrY"),]
cnv_matrix_genes_autosomes <- cnv_matrix_genes[!cnv_matrix_genes$chr %in% c("chrX", "chrY"),]
print(paste0("Removal of ",nrow(cnv_matrix) - nrow(cnv_matrix_autosomes) , " CNV regions at chrX & chrY" ))

source("scripts/R/gsam.dna_seq.metadata.R")

source("scripts/R/chrom_sizes.R")

# there is a difference in size of results files of CNVKit
# they either have 27493 lines (all in ExtraSequencing_CN_data_29Mar2018/)
# or have 27287 lines (195/374) = 50%



# ---- PCA ----

"
We would like to exclude chrX and chrY regions unless we would like to sample swaps (with we also need to do)
"

#dd <- dist(d)

ntop <- 5000
variances <- rowVars(as.matrix(f))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- f[select,]

dim(high_variance_genes)



pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.45,
     main=paste0("G-SAM DNA CNV PCA [chrX+chrY] (pc",pc1," & pc",pc2,")"),col=as.integer(cond), pch=(as.numeric(cond2) - 1) * 18  + 1 )
legend("bottomright", c(unique(as.character(cond)), "batch1", "batch2"),col=c(as.integer(unique(cond)), "darkgray", "darkgray"),pch=c(15,15,15,1,19))



# ---- MDS ----

dd <- dist(high_variance_genes)
dd2 <- dist(t(high_variance_genes))

fit <- cmdscale(dd,eig=TRUE, k=2)
fit2 <- cmdscale(dd2,eig=TRUE, k=2)


x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="p",pch=19,cex=0.5)




c <- cor(t(f))
plot1 <- corrplot(cor(t(f)), 
                  method="square",
                  order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black")


#chrs_hg19_s[genes$chr]


# ---- t-SNE ----


# ---- hierarchical clustering ---- 


e <- cnv_matrix_autosomes[,1:50]

# striking difference in this sample:
plot(c(-5,5),c(-5,5),type="n")
points(e[,1], e[,2], pch=19, cex=0.2)

cor(e[,1], e[,2])

f <- e
for(i in 1:ncol(f)) {
  f[,i] <- f[,i] - median( f[,i][ f[,i] <= 1 & f[,i] >= -1 ] )
}
#f[f > 10] <- 10 # extreme  outlier reduction
#f[f < -10] <- -10 # extreme  outlier reduction

dist(t(e[,1:2]))
dist(t(f[,1:2]))

# normalise? median should be 0?
# colMedians(as.matrix(d))

#ntop <- 2500
#variances <- rowVars(as.matrix(e))
#select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
#high_variance_genes <- e[select,]


dist_mat <- dist(t(f), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)






# random correlation

plot(c(-3, 7), c(-3, 7), type="n")
points(d[,3] , d[,4], pch=19,cex=0.15)

lines(c(-2,10),c(-2,10),col="red")

# ---- corrplot ----
f <- e
for(i in 1:ncol(f)) {
  #f[,i] <- f[,i] - median( f[,i][ f[,i] <= 1 & f[,i] >= -1 ] ) 
}

cc <- cor(cnv_matrix_autosomes, method="spearman")
library(corrplot)
corrplot(cc,
         method="square",
         order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black")

dd <- dist(cc)
sim.by.hclust <- hclust(dd)

pdf("output/figures/cnv/hclust_on_spearmancor_on_log2_cnv_values.pdf", width=32)
plot(sim.by.hclust,cex=0.5)
dev.off()


# measure distance between pairs and median with any other




# ---- pca ----


e <- d[,]#1:250]

# normalise? median should be 0?
# colMedians(as.matrix(d))

ntop <- 2500
variances <- rowVars(as.matrix(e))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e[select,]

# scree plot = pca variance plot
# screeplot(pc)

pc <- prcomp(t(high_variance_genes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.7,pch=19, main=paste0("G-SAM RNA PCA (pc",pc1," & pc",pc2,")"))#,col=as.numeric(cond)+1)

#((pc$sdev / sum(pc$sdev) * 100)


# ---- Pheatmap ----

# ---- Sex Plot ----

cnv_matrix_allosomes <- cnv_matrix[cnv_matrix_genes$chr %in% c("chrX", "chrY"),]
cnv_matrix_genes_allosomes <- cnv_matrix_genes[cnv_matrix_genes$chr %in% c("chrX", "chrY"),]

# batch as number of lines in cnvkit file
batch1 <- c("PD29174a2","PD29174c2","PD29177a2","PD29177c2","PD29178a2","PD29178c2","PD29179a2","PD29179c","PD29181c2","PD29183a","PD29183c","PD29184a2","PD29184c2","PD29185a","PD29185c2","PD29186a","PD29186c","PD29188a","PD29188c2","PD29189a","PD29189c2","PD29190a","PD29190c2","PD29192a","PD29192c2","PD29193a","PD29193c","PD29194a","PD29194c","PD29195a2","PD29195c2","PD29196a","PD29196c","PD29197a2","PD29197c2","PD29198a2","PD29198c","PD29199a2","PD29199c2","PD29200a2","PD29200c2","PD29201a2","PD29201c2","PD29202a2","PD29202c","PD29203a2","PD29203c2","PD29204a2","PD29204c2","PD29205a2","PD29205c2","PD29206a2","PD29206c2","PD29207a2","PD29207c2","PD29208a2","PD29208c","PD29209a2","PD29209c2","PD29212a2","PD29212c2","PD29214a2","PD29214c2","PD29216a2","PD29216c2","PD29217a2","PD29217c","PD29218a2","PD29218c2","PD29220a2","PD29220c","PD29221c2","PD29221c","PD29222a2","PD29222c","PD29223a2","PD29223c","PD29224a","PD29224c2","PD29225a2","PD29225c","PD29226a","PD29226c","PD29227a2","PD29227c2","PD29228a","PD29228c","PD29229a","PD29229c2","PD29230a2","PD29230c2","PD29232a2","PD29232c","PD29233a","PD29233c2","PD29235a2","PD29235c2","PD29236a2","PD29236c2","PD29239a2","PD29239c2","PD29240a2","PD29240c","PD29241a","PD29241c2","PD29242a","PD29242c2","PD29243a","PD29243c2","PD29244a2","PD29246a","PD29246c2","PD29249a2","PD29249c2","PD29250a2","PD29250c2","PD29255a","PD29255c2","PD29256a2","PD29256c","PD29257a2","PD29257c2","PD29258a","PD29258c2","PD29262a2","PD29262c","PD29263a","PD29263c2","PD29264a2","PD29264c","PD29265a","PD29265c2","PD29266a2","PD29266c2","PD30207a","PD30207c","PD30208a","PD30208c","PD30209a","PD30209c","PD30210a","PD30210c","PD30211a","PD30211c","PD30216a","PD30216c","PD30217a","PD30217c","PD30218a","PD30218c","PD30219a","PD30219c","PD30221a","PD30221c","PD30224a","PD30224c","PD30225a","PD30225c","PD30229a","PD30229c","PD30236a","PD30236c","PD30240a","PD30240c","PD30243a","PD30243c","PD30244a","PD30244c","PD30246a","PD30246c","PD30247a","PD30248a","PD30248c","PD30249a","PD30249c","PD30250a","PD30250c","PD30252a","PD30252c","PD30253c","PD30255a","PD30255c","PD30256a","PD30256c","PD30257a","PD30258a","PD30258c","PD30260c","PD30261a","PD30261c","PD30262a","PD30262c","PD30263c","PD30264a","PD30264c")
batch <- as.factor(colnames(cnv_matrix_allosomes) %in% batch1)
cnv_matrix_allosomes <- removeBatchEffect(cnv_matrix_allosomes, batch)


# sex
cond <- gsam.dna_seq.metadata[match(colnames(cnv_matrix_allosomes), gsam.dna_seq.metadata$PD_ID ),]$donor_sex

cond <- as.factor(gsam.dna_seq.metadata[match(colnames(cnv_matrix_allosomes), gsam.dna_seq.metadata$PD_ID ),]$COSMIC_HIST)


# primary(a) / recurrent(c)
#cond <- as.factor(gsub("^PD[0-9]+","",gsub("[0-9]+$","",colnames(cnv_matrix_allosomes))))

# date processed / batch
#cond <- as.factor(gsam.dna_seq.metadata[match(colnames(cnv_matrix_allosomes), gsam.dna_seq.metadata$PD_ID ),]$Date.Received)


#png("output/figures/cnv/sex-plot_dna_cnv.png",width=2*480, height=2*480, res=2*72)



pc <- prcomp(t(cnv_matrix_allosomes))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.45,
     main=paste0("G-SAM DNA CNV PCA [chrX+chrY] (pc",pc1," & pc",pc2,")"),col=as.integer(cond), pch=(as.numeric(cond2) - 1) * 18  + 1 )
legend("bottomright", c(unique(as.character(cond)), "batch1", "batch2"),col=c(as.integer(unique(cond)), "darkgray", "darkgray"),pch=c(15,15,15,1,19))



#dev.off()



rm(cnv_matrix_allosomes, cnv_matrix_genes_allosomes)







# ---- batch effect plot ----

"
Met name chrY heeft grote verschillen, maar ook chr1,
chr2, chr11, chr16, chr20, chr21 en chrX laten een
duidelijk verschil zien in batches.
"

png("output/figures/cnv/batch-effect_pca_cnvkit_cnr.png",width=2*480,height=2*480,res=2*72)

batch <- as.factor(gsub("^.+\\.","",colnames(cnv_matrix)))
sids <- gsub("\\.b[0-9]$","",colnames(cnv_matrix))
cond <- gsam.dna_seq.metadata[match(sids, gsam.dna_seq.metadata$PD_ID ),]$donor_sex

cond <- as.factor(paste0(as.character(cond), ".", gsub("b","batch",as.character(batch),fixed=T)))

tmp <- cnv_matrix[cnv_matrix_genes$chr %in% c("chr22"),]
pc <- prcomp(t(tmp))
pc1 <- 1
pc2 <- 2
plot(pc$x[,c(pc1,pc2)],  cex=0.50,
     main=paste0("G-SAM DNA CNV PCA (pc",pc1," & pc",pc2,")"),col=as.integer(cond), pch=19)

#legend("bottomright", 
#        unique(as.character(cond)),
#       col=as.integer(unique(cond)),
#       pch=19)


dev.off()



"
Show from how many of the samples the batch is different per resection

batch affected samples:

PD29181a	PD29181c2
PD29244a2	PD29244c2
PD30253a	PD30253c
PD30260a	PD30260c
PD30263a	PD30263c
"

cnv_matrix_allosomes <- cnv_matrix#[cnv_matrix_genes$chr %in% c("chrX","chrY"),]
#cnv_matrix_genes_allosomes <- cnv_matrix_genes[cnv_matrix_genes$chr %in% c("chrX", "chrY"),]

sids <-colnames(cnv_matrix)
batch1 <- c("PD29174a2","PD29174c2","PD29177a2","PD29177c2","PD29178a2","PD29178c2","PD29179a2","PD29179c","PD29181c2","PD29183a","PD29183c","PD29184a2","PD29184c2","PD29185a","PD29185c2","PD29186a","PD29186c","PD29188a","PD29188c2","PD29189a","PD29189c2","PD29190a","PD29190c2","PD29192a","PD29192c2","PD29193a","PD29193c","PD29194a","PD29194c","PD29195a2","PD29195c2","PD29196a","PD29196c","PD29197a2","PD29197c2","PD29198a2","PD29198c","PD29199a2","PD29199c2","PD29200a2","PD29200c2","PD29201a2","PD29201c2","PD29202a2","PD29202c","PD29203a2","PD29203c2","PD29204a2","PD29204c2","PD29205a2","PD29205c2","PD29206a2","PD29206c2","PD29207a2","PD29207c2","PD29208a2","PD29208c","PD29209a2","PD29209c2","PD29212a2","PD29212c2","PD29214a2","PD29214c2","PD29216a2","PD29216c2","PD29217a2","PD29217c","PD29218a2","PD29218c2","PD29220a2","PD29220c","PD29221c2","PD29221c","PD29222a2","PD29222c","PD29223a2","PD29223c","PD29224a","PD29224c2","PD29225a2","PD29225c","PD29226a","PD29226c","PD29227a2","PD29227c2","PD29228a","PD29228c","PD29229a","PD29229c2","PD29230a2","PD29230c2","PD29232a2","PD29232c","PD29233a","PD29233c2","PD29235a2","PD29235c2","PD29236a2","PD29236c2","PD29239a2","PD29239c2","PD29240a2","PD29240c","PD29241a","PD29241c2","PD29242a","PD29242c2","PD29243a","PD29243c2","PD29244a2","PD29246a","PD29246c2","PD29249a2","PD29249c2","PD29250a2","PD29250c2","PD29255a","PD29255c2","PD29256a2","PD29256c","PD29257a2","PD29257c2","PD29258a","PD29258c2","PD29262a2","PD29262c","PD29263a","PD29263c2","PD29264a2","PD29264c","PD29265a","PD29265c2","PD29266a2","PD29266c2","PD30207a","PD30207c","PD30208a","PD30208c","PD30209a","PD30209c","PD30210a","PD30210c","PD30211a","PD30211c","PD30216a","PD30216c","PD30217a","PD30217c","PD30218a","PD30218c","PD30219a","PD30219c","PD30221a","PD30221c","PD30224a","PD30224c","PD30225a","PD30225c","PD30229a","PD30229c","PD30236a","PD30236c","PD30240a","PD30240c","PD30243a","PD30243c","PD30244a","PD30244c","PD30246a","PD30246c","PD30247a","PD30248a","PD30248c","PD30249a","PD30249c","PD30250a","PD30250c","PD30252a","PD30252c","PD30253c","PD30255a","PD30255c","PD30256a","PD30256c","PD30257a","PD30258a","PD30258c","PD30260c","PD30261a","PD30261c","PD30262a","PD30262c","PD30263c","PD30264a","PD30264c")
batch <- as.character(sids %in% batch1)

pids <- gsub("^(PD[0-9]+).+","\\1",sids)

for(pid in  unique(sort(pids))) {
  sel <- pids == pid
  #print(sum(sel))
  #print(batch[sel])
  
  if ( length(unique(batch[sel])) > 1) {
    print(sids[sel])
  }
}


# plot diff PD29181a ~ PD29181c2 to see probes

plot(c(1,nrow(cnv_matrix_allosomes)),c(-10,10),type="n")
points(cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD29181a"] - cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD29181c2"], pch=19,cex=0.5)

plot(c(1,nrow(cnv_matrix_allosomes)),c(-10,10),type="n")
points(cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD29244a2"] - cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD29244c2"], pch=19,cex=0.5)

plot(c(1,nrow(cnv_matrix_allosomes)),c(-10,10),type="n")
points(cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD30253a"] - cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD30253c"], pch=19,cex=0.5)

plot(c(1,nrow(cnv_matrix_allosomes)),c(-10,10),type="n")
points(cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD30260a"] - cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD30260c"], pch=19,cex=0.5)

plot(c(1,nrow(cnv_matrix_allosomes)),c(-10,10),type="n")
points(cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD30263a"] - cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD30263c"], pch=19,cex=0.5)


plot(c(1,nrow(cnv_matrix_allosomes)),c(-10,10),type="n")
points(cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD29174a2"] - cnv_matrix_allosomes[,colnames(cnv_matrix_allosomes) == "PD29174c2"], pch=19,cex=0.5)


# ---- all probe change plot ----

d <- cnv_matrix[,!colnames(cnv_matrix) %in% c("PD29181a","PD29181c2","PD29244a2","PD29244c2","PD30253a","PD30253c","PD30260a","PD30260c","PD30263a","PD30263c")]
sid <- gsub("[0-9]$","",colnames(d))
colnames(d) <- sid
c <- gsub("^PD[0-9]+","",sid)
da <- d[, c=="a"]
dc <- d[, c=="c"]

shared <- intersect( gsub("a","",colnames(da)) , gsub("c","",colnames(dc)) )

da <- da[,match(shared, gsub("a","",colnames(da)))]
dc <- dc[,match(shared, gsub("c","",colnames(dc)))]

# gsub("a","",colnames(da)) == gsub("c","",colnames(dc))
stopifnot ( length( gsub("a","",colnames(da)) == gsub("c","",colnames(dc)) ) == 178 )


delta <- da - dc





