#!/usr/bin/env R

setwd("~/projects/gsam")

# ---- load data ----

source("scripts/R/gsam_metadata.R")

barplot_theme <- theme(
  text = element_text(family = 'Helvetica'),
  axis.text.x = element_text(angle = 45,size=10),
  axis.text.y = element_text(size=5),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.text = element_text(),
  axis.line = element_line(colour="black"),
  panel.grid.major.y = element_line(colour = 'grey90', linetype = 'dotted')
)
source('scripts/R/job_gg_theme.R')


# ---- ----

counts_stranded <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt",stringsAsFactors = F,comment="#")
colnames(counts_stranded) <- gsub("^[^\\.]+\\.([^\\]+)\\.Aligned.sorted.+$","\\1",colnames(counts_stranded),fixed=F)

rownames(counts_stranded) <- counts_stranded$Geneid
counts_stranded$Geneid <- NULL
counts_stranded$Chr <- NULL
counts_stranded$Start <- NULL
counts_stranded$End <- NULL
counts_stranded$Strand <- NULL
counts_stranded$Length <- NULL


# 


counts_unstranded <- read.delim("output/tables/gsam_featureCounts_readcounts.unstranded.txt",stringsAsFactors = F,comment="#")
colnames(counts_unstranded) <- gsub("^[^\\.]+\\.RNA.alignments.([^\\]+)\\.Aligned.sorted.+$","\\1",colnames(counts_unstranded),fixed=F)

rownames(counts_unstranded) <- counts_unstranded$Geneid
counts_unstranded$Geneid <- NULL
counts_unstranded$Chr <- NULL
counts_unstranded$Start <- NULL
counts_unstranded$End <- NULL
counts_unstranded$Strand <- NULL
counts_unstranded$Length <- NULL




# 

counts_stranded[counts_stranded < 1 ] <- 1

counts_antistranded <- counts_unstranded - counts_stranded
counts_antistranded[counts_antistranded < 1 ] <- 1


# ---- selecet only genes that have some reads ----

sel <- rowSums(counts_stranded) > 372 * 15
counts_stranded <- counts_stranded[sel,]
counts_antistranded <- counts_antistranded[sel,]
counts_unstranded <- counts_unstranded[sel,]
rm(sel)

# 

#stranded_antistranded_ratio <- log( (counts_stranded + 0.0001) / (counts_antistranded + 0.0001) )
stranded_antistranded_ratio <- log( (counts_stranded ) / (counts_antistranded ) )
no_dna_c <- colMeans(stranded_antistranded_ratio) >= 2.6
dna_c <- colMeans(stranded_antistranded_ratio) < 2.6 & colMeans(stranded_antistranded_ratio) >= 1.8
heavy_dna_c <- colMeans(stranded_antistranded_ratio) < 1.8

plot(c(1, nrow(stranded_antistranded_ratio)) , c(-8,12) , type='n',ylab="log(counts stranded / count anti-stranded)")
abline(h=0, col="gray",lty=2)
for(i in which(no_dna_c)) {
  lines(sort(stranded_antistranded_ratio[,i], decreasing=T))
}
for(i in which(dna_c)) {
  lines(sort(stranded_antistranded_ratio[,i], decreasing=T),col="orange")
}
for(i in which(heavy_dna_c)) {
  lines(sort(stranded_antistranded_ratio[,i], decreasing=T),col="red")
}
rm(i)
rm(counts_antistranded, counts_stranded, counts_unstranded)

df <- data.frame(sid = colnames(stranded_antistranded_ratio),
                 ratio.stranded.antistranded = colMeans(stranded_antistranded_ratio),
                 ratio.stranded.antistranded.lin = exp(colMeans(stranded_antistranded_ratio)), 
                 ratio.stranded.antistranded.dna = rep(NA, ncol(stranded_antistranded_ratio)) )
df$ratio.stranded.antistranded.dna [no_dna_c] <- "No-DNA-contamination"
df$ratio.stranded.antistranded.dna [dna_c] <- "Moderate-DNA-contamination"
df$ratio.stranded.antistranded.dna [heavy_dna_c] <- "Heavy-DNA-contamination"
df$ratio.stranded.antistranded.dna <- as.factor(df$ratio.stranded.antistranded.dna)
df$sid <- gsub(".","-",as.character(df$sid), fixed=T)

rm(dna_c, heavy_dna_c, no_dna_c, stranded_antistranded_ratio)
gsam.rna.metadata <- merge(gsam.rna.metadata, df, by.x="sid",by.y="sid")


# ---- stddev *pm ----

"
We nemen aan dat wanneer er hoge mate van DNA contaminatie is, dat de read-count per gen linear gecorreleerd is met de grootte van ieder gen

De variantie van FPKM/TPM achtige waardes over alle genen zou dan extreem laag moeten zijn
"

counts_stranded <- read.delim("data/output/tables/gsam_featureCounts_readcounts.txt",stringsAsFactors = F,comment="#")
colnames(counts_stranded) <- gsub("^[^\\.]+\\.([^\\]+)\\.Aligned.sorted.+$","\\1",colnames(counts_stranded),fixed=F)
len <- counts_stranded$Length
rownames(counts_stranded) <- counts_stranded$Geneid
counts_stranded$Geneid <- NULL
counts_stranded$Chr <- NULL
counts_stranded$Start <- NULL
counts_stranded$End <- NULL
counts_stranded$Strand <- NULL
counts_stranded$Length <- NULL

sel <- rowSums(counts_stranded) > 372 * 15
counts_stranded <- counts_stranded[sel,]


tt <- counts_stranded
tt <- (tt / len) 
tt <- (tt / colSums(tt))

colvar1 <- colVars( as.matrix ( tt ))
names(colvar1) <- colnames(tt)
plot(log(sort(colvar1)))


tt <- counts_stranded
tt <- (tt / len) 
#tt <- (tt / colSums(tt))

colvar2 <- colVars( as.matrix ( tt ))
names(colvar2) <- colnames(tt)
plot(log(sort(colvar2)))


# ---- plot -----

tmp <-  gsam.rna.metadata
order <- order(  scale(tmp$ratio.stranded.antistranded.lin + scale(tmp$pct.nofeature.STAR)) * (2* scale(tmp$pct.rRNA.by.chrUn.gl000220)) )
tmp <- tmp[order,]
tmp$x <- 1:nrow(tmp)
tmp$sid <- factor(as.character(tmp$sid[tmp$x]), levels = as.character(tmp$sid[tmp$x]))


#names(colvar1)[match(gsub("-",".",tmp$sid), names(colvar1))] == tmp$sid
tmp$colvar1 <- log(  colvar1[match(gsub("-",".",tmp$sid), names(colvar1))]  )
tmp$colvar2 <- log(  colvar2[match(gsub("-",".",tmp$sid), names(colvar2))]  )


plot_grid(
  ggplot(tmp, aes(x = sid ,y = ratio.stranded.antistranded.lin, fill=ratio.stranded.antistranded.dna, label=sid)) +
    #coord_flip() + 
    geom_bar(stat = "identity", position = "stack",colour="black") + 
    scale_y_continuous() + 
    theme_bw() + 
    barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 ))
  ,
  ggplot(tmp, aes(x = sid ,y =  pct.nofeature.STAR , label=sid)) +
    #coord_flip() + 
    geom_bar(stat = "identity", position = "stack",colour="black") + 
    scale_y_continuous() + 
    theme_bw() + 
    barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 ))
  ,
  ggplot(tmp, aes(x = sid ,y =  tmp$pct.rRNA.by.chrUn.gl000220 , label=sid)) +
    #coord_flip() + 
    geom_bar(stat = "identity", position = "stack",colour="black") + 
    scale_y_continuous() + 
    theme_bw() + 
    barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 )),
  ggplot(tmp, aes(x = sid ,y =  STAR.Assigned , label=sid)) +
    #coord_flip() + 
    geom_bar(stat = "identity", position = "stack",colour="black") + 
    scale_y_continuous() + 
    theme_bw() + 
    barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 )),
  ggplot(tmp, aes(x = sid ,y =  colvar1 , label=sid)) +
    #coord_flip() + 
    geom_bar(stat = "identity", position = "stack",colour="black") + 
    scale_y_continuous() + 
    theme_bw() + 
    barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 )),
  ggplot(tmp, aes(x = sid ,y =  colvar2 , label=sid)) +
    #coord_flip() + 
    geom_bar(stat = "identity", position = "stack",colour="black") + 
    scale_y_continuous() + 
    theme_bw() + 
    barplot_theme + 
    theme( axis.title.y = element_text(size = 11) ,
           axis.text.x = element_text(angle = 90, size = 5 ))
  ,align="v", axis="tblr",ncol=1  )


ggsave("output/figures/qc/dna-contamination_or_rRNA.pdf",height=20,width=16*1.5)


# ---- in/exclusion fig ----

tmp <- gsam.rna.metadata
tmp$isolation.order <- as.numeric(gsub("[\\.\\-][0-9]+$","",tmp$isolation.id))

order <- order(tmp$isolation.order , tmp$blacklist.heavy.dna.contamination,
               tmp$storage.box, tmp$plate, gsam.rna.metadata$sid)
tmp <- tmp[order,]
rm(order)
tmp$x <- 1:nrow(tmp)
tmp$sid <- factor(as.character(tmp$sid[tmp$x]), levels = as.character(tmp$sid[tmp$x]))



cols_low.assigned <- c("FALSE" = "white", "TRUE" = "black",
                       "24D3"="#000000",
                       "25A1"="#111111",
                       "25B1"="darkgreen",#222222",
                       "25B3"="#333333",
                       "26A2"="blue",#444444",
                       "26B1"="#555555",
                       "27A1"="#666666",
                       "27A2"="#777777",
                       "27A4"="yellow",#"#888888",
                       "27B2"="#999999",
                       "27B3"="#AAAAAA",
                       "27B4"="red",#BBBBBB",
                       "27C2"="#CCCCCC",
                       "27C3"="#DDDDDD",
                       "28C3"="#EEEEEE",
                       "28C4"="#FFFFFF",
                       "plate1" = "#444444",
                       "plate2" = "#888888",
                       "plate3" = "#BBBBBB",
                       "plate4" = "#FFFFFF"
                       )
ggplot( tmp , aes(x = sid ,y = 1 ,  fill = factor(plate)) )   +
  scale_fill_manual(values = cols_low.assigned) + 
  geom_tile(color = "black", aes(fill=factor(blacklist.too.low.assigned))) + 
  geom_tile(color = "black", aes(fill = factor(blacklist.heavy.dna.contamination) , y = 2)) + 
  geom_tile(color = "black", aes(fill = factor(storage.box) , y = 3)) + 
  geom_tile(color = "black", aes(fill = factor(plate) , y = 4)) + 
  job_gg_theme 




