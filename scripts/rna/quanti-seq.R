#!/urs/bin/env R

# ---- load libs ----
library(plyr)
library(dplyr)
library(ggplot2)
library(Biobase)

source('scripts/R/job_gg_theme.R')

# ---- load data ----



dir <- 'output/tables/rsem/'
files.RSEM <- list.files(dir, pattern = '.genes.results', full.names = T)

data.TPM.NGSProtocol <- tibble::as_tibble(do.call(cbind, lapply(files.RSEM, function(x){ y <- readr::read_delim(x, delim = '\t'); return(y$TPM)})))
colnames(data.TPM.NGSProtocol) <- gsub('(.*).genes.results', '\\1', basename(files.RSEM))

data.TPM.NGSProtocol <- as.matrix(data.TPM.NGSProtocol)
data.RSEM <- readr::read_delim(files.RSEM[[1]], delim = '\t')[,1]

gtf <- "/mnt/data/ccbc_environment/project/gsam-neurology/ref/star-hg19/gencode.v31lift37.annotation.gtf"
geneInfo <- rtracklayer::import.gff(gtf)
geneInfo <- tibble::as_tibble(GenomicRanges::mcols(geneInfo))
geneInfo <- geneInfo %>% dplyr::distinct(gene_id, gene_name)

data.RSEM <- data.RSEM %>% dplyr::inner_join(geneInfo, by = c('gene_id' = 'gene_id'))

rownames(data.TPM.NGSProtocol) <- data.RSEM$gene_name


data.quantiseq <- immunedeconv::deconvolute(data.TPM.NGSProtocol, method = 'quantiseq', tumor = T, scale_mrna = T)

save(data.quantiseq, file = 'output/tables/quantiseq/data.quantiseq.RData')


# ----  ----

load('output/tables/quantiseq/data.quantiseq.RData')

# clusteren op waarden kan misschien beter op een transformatie, voglens job de cosine similary en daar de dist van

data.quantiseq.melt <- reshape2::melt(data.quantiseq)

#data.quantiseq.melt <- data.quantiseq.melt %>% dplyr::inner_join(sampleInfo, c('variable' = 'Accession'))

data.quantiseq.melt$resection <- as.factor(gsub("^.+([0-9])$","\\1",
               gsub("-replicate","",as.character(data.quantiseq.melt$Sample),fixed=T)   
               ))

clusterData <- scale(data.quantiseq[2:ncol(data.quantiseq)]  , center = T, scale = T)

#a = clusterData[10,]
#b = as.numeric(as.matrix(data.quantiseq)[10,-1])
#plot(a, b)






d <- dist(t(clusterData), method = 'euclidean')
hc <- hclust(d, method="ward.D2")

clusteringLabels <- hc$labels[hc$order]

data.quantiseq.melt$Sample <- factor(data.quantiseq.melt$variable, levels = clusteringLabels)

plot.dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, size = 1, labels = F)

# Order of cell-types
data.quantiseq.melt$cell_type <- factor(data.quantiseq.melt$cell_type, levels = rev(unique(c("Dendritic cell", data.quantiseq[order(rowSums(data.quantiseq[2:ncol(data.quantiseq)])),]$cell_type))))

labels <- data.quantiseq.melt[data.quantiseq.melt$cell_type == levels(data.quantiseq.melt$cell_type)[1],]

#plot.quantiseq <- 
ggplot(data.quantiseq.melt, aes(x = Sample, y = value, fill = cell_type)) + 
  geom_bar(stat = 'identity', width = 1, color = 'black', size = .1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(type = 'seq', palette = 'Paired', direction = 1, guide = guide_legend(title = 'Cell Populations', title.position = 'top', title.hjust = 0.5, nrow = 2, keywidth = 0.5, keyheight = 0.5)) +
  # Axis titles.
  labs(x = paste0("G-SAM\n(n = ",ncol(data.quantiseq) - 1,")"), y = 'Relative frequency\n(TIL10)') + 
  theme(legend.position = 'bottom', 
        axis.title.y = element_text(size = 8), 
        text=element_text(size=10, family='Helvetica'),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = NA),
        panel.border = element_rect(fill = NA, colour = 'grey20')
  )  + 
  geom_point(aes(x = Sample, y = -0.025, col=resection),pch=22,data=labels,size=0.2)

ggsave("output/figures/rna/quanti-seq/profiles-overview.pdf",width=16*1.7,height=6*1.5)






tmp <- data.quantiseq.melt[data.quantiseq.melt$cell_type == levels(data.quantiseq.melt$cell_type)[6],]
tmp <- tmp[order(tmp$resection, tmp$value),]
tmp$x <- 1:nrow(tmp)
labels <- tmp[tmp$cell_type == as.character(tmp$cell_type)[1],]
ggplot(tmp, aes(x = x, y = value, fill = cell_type)) + 
  geom_bar(stat = 'identity', width = 1, color = 'black', size = .1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(type = 'seq', palette = 'Paired', direction = 1, guide = guide_legend(title = 'Cell Populations', title.position = 'top', title.hjust = 0.5, nrow = 2, keywidth = 0.5, keyheight = 0.5)) +
  # Axis titles.
  labs(x = paste0("G-SAM\n(n = ",ncol(data.quantiseq) - 1,")"), y = 'Relative frequency\n(TIL10)') + 
  theme(legend.position = 'bottom', 
        axis.title.y = element_text(size = 8), 
        text=element_text(size=10, family='Helvetica'),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = NA),
        panel.border = element_rect(fill = NA, colour = 'grey20')
  )  + 
  geom_point(aes(x = x, y = -0.025, col=resection),pch=22,data=labels,size=0.2)







