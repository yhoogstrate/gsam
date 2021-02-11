#!/usr/bin/env R

# load libs ----

library(clusterRepro)
library(matrixStats)


library(ggplot2)
library(pheatmap)
library(rgl)
library(gridExtra)


source('scripts/R/youri_gg_theme.R')
source('scripts/R/gsam_metadata.R')


# 01: load raw Array data ----
 
glio_dat <- read.table("data/gravendeel/RDamClust_RMA_Dai.txt" , header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

rownames(glio_dat) <- glio_dat[,3]

glio_dat$var <- apply(glio_dat[,5:296],1,var)
glio_dat <-tail(glio_dat[order(glio_dat$var), ], 5000)
glio_dat <- glio_dat[,5:296]
glio_dat <- glio_dat %>%
  dplyr::mutate_all( function(x) as.numeric(as.character(x)))
glio_gnames <- rownames(glio_dat)
sum(duplicated(glio_gnames)) # no duplicates




# 02: load the labels ----
#glio_lab <- read.csv(file="D:/my documents/arrays/array data summary files/TCGA/juiste_labels_Olivier.csv", sep=';', header=F)
glio_lab <- read.table("data/gravendeel/juiste_labels_Olivier_merge.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)



# 03: centroid function ----

# create a loop that calculates the centroids per cluster
# and puts these in a matrix
centroid_fun <- function(cl, dt, annt) {

	clust_lab <- unique(cl)
	CL <- matrix(clust_lab, 1, length(clust_lab))
	Centroids <- matrix(0, nrow(dt), ncol(CL))

	for (i in 1:ncol(CL)) {
		Centroids[,i] <- rowMeans(subset(dt,
		select=cl==CL[i]))
	}

	colnames(Centroids) <-CL
	rownames(Centroids) <-annt
	list(Centroids=Centroids)
}

fun1 <- centroid_fun(glio_lab$V2, glio_dat, glio_gnames)

dim(fun1$Centroids)
fun1$Centroids[0,]

glio_cen <- fun1$Centroids
colnames(glio_cen)


# 04: match centroids and set2 by rownames function ----
match_fun <- function(cent, set2)	{
	cent.N <- rownames(cent)
	set2.N <- rownames(set2)
	common <- intersect(cent.N, set2.N)
	cent.I <- which(cent.N %in% common)
	set2.I <- which(set2.N %in% common)
	cent.S <- cent[cent.I,]
	set2.S <- set2[set2.I,]
	cent.O <- cent.S[order(rownames(cent.S)),]
	set2.O <- set2.S[order(rownames(set2.S)),]
	list(cent.O=cent.O, set2.O=set2.O)
}



# # load dataset to be clustered [GSE4271] ----
# 
# phil<- read.delim('GSE4271-GPL96_series_matrix.txt', sep='\t', header=T,  skip=42)
# phil<- phil[-c(1:37),] # exclude first 37 columns?
# #phil2<-phil
# for(i in 2:101) { # 2:101 zijn numerical, dus expressie
#   phil[,i]<- as.numeric(as.character(phil[,i]))
#   phil[,i]<- log2(phil[,i]) # transform to ~ normal data?
# }
# 
# # this estimates the variance?
# 
# phil$var <- apply(phil[,2:100],1, var) # komt achteraan
# phil <- phil[order(-phil$var), ] # is dit om zo dubbele per-gen probe met hoogste variantie te pakken?
# phil <- phil[!duplicated(phil[,1]),]
# row.names(phil)<-phil[,1] # col 1 = gene.id
# 
# 
# ## reduce to 1 probe per gene?
# annot<-read.delim('D:/My Documents/arrays/array data summary files/HU133A TCGA rma.summary.txt',sep='\t', header=T,row.names=1)
# annot<-annot[,-c(1:236,239)]
# head(annot)
# phil <- merge(annot,phil, by='row.names')
# phil <- phil[!duplicated(phil[,3]),]
# row.names(phil)<-phil[,3]
# 
# phil <- phil[,-c(1:4,105)]
# 
# colnames(phil)<-gsub("..HG.U133A.","",colnames(phil))
# 
# 
# 
# 
# 
# 
# # select only most relevant clusters
# myvars <- c('0','9', '16', '17', '18', '22', '23')
# glio_cs <- glio_cen[,myvars]
#       dim(glio_cs) 
#       glio_cs[0,]
#       rownames(glio_cs)
# fun2 <- match_fun(glio_cs, phil)
# nrow(fun2$cent.O) # check the number of genes in common
# 
# 
# ## The clusterRepro package
# 
# 
# 
# Result <- IGP.clusterRepro(phil, glio_cs)
# cl.assign<-as.data.frame(Result$Class)
# row.names(cl.assign)<- gsub('.AVG_Signal', '', row.names(cl.assign))
# row.names(cl.assign)<- gsub('X', '', row.names(cl.assign))
#   cl.assign<- as.matrix(Result$Class)  # make a matrix in which all the clusters are renamed
#   cl.assign<- gsub('1', 'IGS-0', cl.assign)
#   cl.assign<- gsub('2', 'IGS-9', cl.assign)
#   cl.assign<- gsub('5', 'IGS-18', cl.assign)
#   cl.assign<- gsub('6', 'IGS-22', cl.assign)
#   cl.assign<- gsub('3', 'IGS-16', cl.assign)
#   cl.assign<- gsub('7', 'IGS-23', cl.assign)    
#   cl.assign<- gsub('4', 'IGS-17', cl.assign)
#   cl.assign<- cbind(Result$Class, cl.assign)
# 
# ## make a report
# 
# write.table(cl.assign, 'output.txt'  ,sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")


#rownames(fun3$set2.O) [rownames(fun3$set2.O) != rownames(fun3$cent.O) ]
#rownames(fun3$cent.O) [rownames(fun3$cent.O) %in%  rownames(fun3$set2.O) == F]

#sum(duplicated(rownames(fun3$set2.O)))
#  "GGT1" "SOD2



# G-SAM ----
## Load data ----
source('scripts/R/gsam_rna-seq_expression.R')
gsam.vst <- gsam.rnaseq.expression.vst %>%
 `rownames<-`( gsub("^[^~|]+[~|](.+)[~|].+$", "\\1" , rownames(.) ,fixed=F) )

# remove duplicates, first order by sd (from high to low)
sd <- rowSds(gsam.vst)
gsam.vst <- gsam.vst[order(- sd), ] # is dit om zo dubbele per-gen probe met hoogste variantie te pakken?
gsam.vst <- gsam.vst[!duplicated(rownames(gsam.vst)),]
rm(sd)


## Normalise ----

# towards same mean and same sd as original array data

fun2 <- match_fun(glio_cs, gsam.vst)
stopifnot(nrow(fun2$cent.O) == nrow(fun2$set2.O))

# Run withoutadditional normalisation
# cl.assign.1 <- IGP.clusterRepro(fun2$set2.O, fun2$cent.O) %>%
#   purrr::pluck('Class') %>%
#   as.data.frame()  %>%
#   `colnames<-`(c("Result.Class")) %>%
#   dplyr::mutate(Gravendeel.Centroid.Class = case_when(
#     Result.Class == 1 ~  'IGS-0',
#     Result.Class == 2 ~  'IGS-9',
#     Result.Class == 3 ~  'IGS-16',
#     Result.Class == 4 ~  'IGS-17',
#     Result.Class == 5 ~  'IGS-18',
#     Result.Class == 6 ~  'IGS-22',
#     Result.Class == 7 ~  'IGS-23'
#   )) %>%
#   dplyr::mutate(Gravendeel.Centroid.Class = case_when( # c('9', '17', '18', '22')
#     Result.Class == 1 ~  'IGS-9',
#     Result.Class == 2 ~  'IGS-17',
#     Result.Class == 3 ~  'IGS-18',
#     Result.Class == 4 ~  'IGS-22'
#   ))  %>%
#   dplyr::mutate(Gravendeel.Centroid.Class = as.factor(Gravendeel.Centroid.Class))
# 
# # make a matrix in which all the clusters are renamed
# plot(cl.assign.1$Class)
# rm(cl.assign.1$Gravendeel.Centroid.Class) # NOT NORMALISED against the array data (!)



# perform median scaling and sd correction to fit it w/ array data
glio_sub <- fun2$set2.O[,1:2] %>% # subset of original data
  as.data.frame() %>%
  `colnames<-`(c("c1","c2")) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::left_join(glio_dat %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column('gid')
                   , by=c('gid' = 'gid')) %>%
  dplyr::select(- c("c1", "c2") ) %>%
  tibble::column_to_rownames('gid') %>%
  as.matrix()

sd = rowSds(glio_sub) # sd in original array data
med = rowMedians(glio_sub) # median in original array data


fun3 <- fun2
fun3$set2.O <-  fun3$set2.O - rowMedians(fun3$set2.O) # median/mean centering
fun3$set2.O <- fun3$set2.O / rowSds(fun3$set2.O) # scale back to sd of 1

fun3$set2.O <-  fun3$set2.O * sd # scale back to sd in original array data
fun3$set2.O <- fun3$set2.O + med # add up the mean from the original array data



## Full classification: 7 classes  ----

myvars <- c('0','9', '16', '17', '18', '22', '23')
cl.assign.normalised.full <- IGP.clusterRepro(fun3$set2.O, fun3$cent.O[,myvars]) %>%
  purrr::pluck('Class') %>%
  as.data.frame()  %>%
  `colnames<-`(c("Result.Class")) %>%
  dplyr::mutate(Gravendeel.Centroid.Class.Full = case_when(
    Result.Class == 1 ~  'IGS-0',
    Result.Class == 2 ~  'IGS-9',
    Result.Class == 3 ~  'IGS-16',
    Result.Class == 4 ~  'IGS-17',
    Result.Class == 5 ~  'IGS-18',
    Result.Class == 6 ~  'IGS-22',
    Result.Class == 7 ~  'IGS-23'
  )) %>%
  dplyr::mutate(Gravendeel.Centroid.Class.Full = as.factor(Gravendeel.Centroid.Class.Full)) %>%
  dplyr::mutate(Result.Class = NULL) %>%
  tibble::rownames_to_column('sid')

plot(cl.assign.normalised.full$Gravendeel.Centroid.Class.Full) # no probabilities; no error values


## Subset classification: 4 major classes  ----


myvars <- c('9', '17', '18', '22')
cl.assign.normalised.sub <- IGP.clusterRepro(fun3$set2.O, fun3$cent.O[,myvars]) %>%
  purrr::pluck('Class') %>%
  as.data.frame()  %>%
  `colnames<-`(c("Result.Class")) %>%
  dplyr::mutate(Gravendeel.Centroid.Class.Subset = case_when( # c('9', '17', '18', '22')
    Result.Class == 1 ~  'IGS-9',
    Result.Class == 2 ~  'IGS-17',
    Result.Class == 3 ~  'IGS-18',
    Result.Class == 4 ~  'IGS-22'
  ))  %>%
  dplyr::mutate(Gravendeel.Centroid.Class.Subset = as.factor(Gravendeel.Centroid.Class.Subset)) %>%
  dplyr::mutate(Result.Class = NULL) %>%
  tibble::rownames_to_column('sid')

plot(cl.assign.normalised.sub$Gravendeel.Centroid.Class.Subset) # no probabilities; no error values


write.table(dplyr::full_join(cl.assign.normalised.full,
                             cl.assign.normalised.sub,
                             by = c('sid' = 'sid')),
            "output/tables/gravendeel_centroid_classification_gsam.txt")



# TODO: consider KNN approach?



# visualise & inspect ----


pc <- prcomp(t(fun2$set2.O))
#screeplot(pc) # at least 3 good coordinates

plt <- pc$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(cl.assign.normalised , by=c('sid' = 'sid'))


ggplot(plt, aes(x=PC1, y=PC3, col=Gravendeel.Centroid.Class) ) +
  geom_point() +
  youri_gg_theme


d <- plt %>%
  dplyr::select(c('sid','PC1','PC2','PC3','PC4')) %>%
  tibble::column_to_rownames('sid')
d <- dist(d) # euclidean distances between the rows
d <- dist(t(fun2$set2.O))
fit <- cmdscale(d, eig=TRUE, k=2)


plt.mds <- fit$points %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sid') %>%
  dplyr::left_join(cl.assign.normalised , by=c('sid' = 'sid')) %>%
  dplyr::left_join(gsam.rna.metadata %>% dplyr::select(c('sid','gliovis.majority_call','NMF.123456.PCA.LDA.class')) , by=c('sid'='sid'))


p1 <- ggplot(plt.mds, aes(x=V1, y=V2, col=Gravendeel.Centroid.Class)) + 
  geom_point() +
  labs(x=paste0("MDS1 on ",nrow(fun2$set2.O)," G-genes"),
       y=paste0("MDS1 on ",nrow(fun2$set2.O)," G-genes"),
       col="Gravendeel centroid fit") +
  youri_gg_theme

p2 <- ggplot(plt.mds, aes(x=V1, y=V2, col=NMF.123456.PCA.LDA.class)) + 
  geom_point() +
  labs(x=paste0("MDS1 on ",nrow(fun2$set2.O)," G-genes"),
       y=paste0("MDS1 on ",nrow(fun2$set2.O)," G-genes"),
       col="Wang subtypes NMF + LDA re-fit") +
  youri_gg_theme

p3 <- ggplot(plt.mds, aes(x=V1, y=V2, col=gliovis.majority_call)) + 
  geom_point() +
  youri_gg_theme


library(patchwork)
p1 + p2
ggsave("output/figures/Gravendeel_Verhaak_MDS.pdf", width=10.9, height=4.94)


# kunnen 0, 16, 23 weg?
# 0: check tumour percentage(!)
# 16: lijkt wel te kunnen
# 23: lijkt weg te kunnen

# eventueel 17 en 9?
# 9: nee, mooi cluster




#dist <- dist(pc$x[,1:3])
#plot(hclust(dist, method = "complete", leaf_labels=plt$Gravendeel.Centroid.Class) )

d <- dist(t(scale(t(pc$x[,1:3]))))
hc <- hclust(d, method = "complete")
dhc <- as.dendrogram(hc,hang=0.1)
##ddata <- dendro_data(dhc, type="rectangle")
#ggdendrogram(hc, leaf_labels=plt$Gravendeel.Centroid.Class)

df <- data.frame(cluster = cutree(hc,345)) %>%
  dplyr::mutate(sid = rownames(.)) %>%
  dplyr::left_join(cl.assign.normalised, by=c('sid'='sid')) %>%
  dplyr::left_join( gsam.rna.metadata %>% dplyr::select(c('sid','gliovis.majority_call')) , by=c('sid'='sid'))

df <- df[hc$order,] %>%
  dplyr::mutate(cluster = 1:nrow(.))

#ggplot(segment(ddata)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

p1 <- ggdendrogram(hc)

p2 <- ggplot(df, aes(cluster , y=1, fill=Gravendeel.Centroid.Class)) +
  geom_tile()  +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")

p3 <- ggplot(df, aes(cluster , y=1, fill=gliovis.majority_call)) +
  geom_tile()  +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")


gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
gp3<-ggplotGrob(p3)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5], gp3$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

grid.arrange(gp1, gp2, gp3, ncol=1,heights=c(4/5,1/10,1/10))




library(pheatmap)
pheatmap(fun3$set2.O, show_row_names = FALSE)


