# Set wd

# ---- initialization setup ----

setwd("~/projects/gsam")

# ---- load: libs ----

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(fgsea)
library(limma)
library(EnhancedVolcano)


# ---- load: functions ----

source("scripts/R/ensembl_to_geneid.R") # obsolete? can be replaced with the get_ensembl function
ensembl_genes <- get_ensembl_hsapiens_gene_ids()


# ---- load: data ----
# Load the data (expression values; metadata) into data frames

source("scripts/R/ligands.R")# @ todo gsub suffixes in ensembl ids
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
gene_matrix$Chr <- gsub("^([^;]+);.+$","\\1",gene_matrix$Chr)
gene_matrix$Start <- gsub("^([^;]+);.+$","\\1",gene_matrix$Start)

source("scripts/R/dna_idh.R")

source("scripts/R/chrom_sizes.R")

source("scripts/R/youri_gg_theme.R")


# ---- EGFR ligand analysis ----

# subset to all genes with avg coun >= 3 or ligand


e <- expression_matrix_full
m <- gsam.rna.metadata[match(gsub(".","-",colnames(e),fixed=T), gsam.rna.metadata$sid),]
print(dim(e))

gsub(".","-",colnames(e),fixed=T) == m$sid # check ids

e <- e[,m$blacklist.pca == F]
m <- m[m$blacklist.pca == F,]

gsub(".","-",colnames(e),fixed=T) == m$sid # check ids

# remove lowcounts
print(dim(e))
sel <- rowSums( e )/nrow(e) > 2 | gsub("\\..+","",rownames(e)) %in% gsub("\\..+","",tt1)
e <- e[sel,]
print(dim(e))

#
cond <- as.factor(paste0('c',round(runif(ncol(e),1,2))))
dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond)
e.vst <- assay(vst(dds, blind=T))




#ENSG00000146648.17_3 ENSG00000109321.10_3  ENSG00000124882.3_2  ENSG00000182585.9_2 
#"EGFR"               "AREG"               "EREG"               "EPGN" 
#ENSG00000163235.15_3 ENSG00000138798.11_3  ENSG00000113070.7_2 ENSG00000174808.11_2 
#"TGFA"                "EGF"              "HBEGF"                "BTC



#plot(e.vst.egfr, e.vst.areg, col=as.numeric(as.factor(m$resection)))
#plot(e.vst.egfr, e.vst.ereg, col=as.numeric(as.factor(m$resection)))

plt <- data.frame(
  e.vst.egfr = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["EGFR"]),],
  e.vst.ereg = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["EREG"]),],
  e.vst.areg = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["AREG"]),],
  e.vst.hbegf = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["HBEGF"]),],
  e.vst.egf = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["EGF"]),],
  e.vst.btc = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["BTC"]),],
  e.vst.tgfa = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["TGFA"]),],
  e.vst.epgn = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["EPGN"]),],
  e.vst.ngf = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["NGF"]),],
  
  e.vst.erbb2 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["ERBB2"]),],
  e.vst.erbb3 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["ERBB3"]),],
  e.vst.erbb4 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["ERBB4"]),],

  e.vst.nrg1 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["NRG1"]),],
  e.vst.nrg2 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["NRG2"]),],
  e.vst.nrg3 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["NRG3"]),],
  
  e.vst.fgfr3 = e.vst[gsub("\\..+$","",rownames(e.vst)) == gsub("\\..+$","",tt1["FGFR3"]),],
  vIII.percentage = m$vIII.percentage >= 1,
  sid = m$sid,
  pid = m$pid,
  resection = m$resection
  )


plt$e.vst.ereg <- plt$e.vst.ereg - min(plt$e.vst.ereg)  + 0.5
plt$e.vst.ereg[plt$res == "r2"] <- -1 * plt$e.vst.ereg[plt$res == "r2"]

ggplot(plt, aes(x=e.vst.egfr, y=e.vst.ereg, group=pid)) +
  #geom_line(col="lightgray") + 
  #geom_point(aes(col = m$vIII.percentage > 1))
  geom_point(aes(col=resection)) + ylim(-max(abs(plt$e.vst.ereg)) , max(abs(plt$e.vst.ereg)) )


# ---- BTC ~ andere genen in locus ----
# AAU: extreem hoog BTC


ggplot(plt, aes(x=e.vst.erbb3, y=e.vst.nrg2, group=pid, label=sid)) +
  #geom_text(data=subset(plt, pid %in% c('GAA','BAO','FAD','GAG','EBX','ACA','AIA','CB2','FAF','AOA','FAI','ECD','CBV'))) +
  geom_point(aes(col=resection)) +   youri_gg_theme 
ggsave("/tmp/ss1.png")

ggplot(plt, aes(x=e.vst.erbb3, y=e.vst.nrg2, group=pid, label=sid)) +
  #geom_text(data=subset(plt, pid %in% c('GAA','BAO','FAD','GAG','EBX','ACA','AIA','CB2','FAF','AOA','FAI','ECD','CBV'))) +
  geom_point(aes(col=e.vst.egfr > 12.3))  +   youri_gg_theme 
ggsave("/tmp/ss2.png")

ggplot(plt, aes(x=e.vst.erbb3, y=e.vst.nrg2, group=pid, label=sid,shape=resection)) +
  #geom_text(data=subset(plt, pid %in% c('GAA','BAO','FAD','GAG','EBX','ACA','AIA','CB2','FAF','AOA','FAI','ECD','CBV'))) +
  geom_point(aes(col=e.vst.fgfr3 > 10.75))  +   youri_gg_theme 

ggplot(plt, aes(x=e.vst.erbb3, y=e.vst.fgfr3, group=pid, label=sid,shape=resection)) +
  #geom_text(data=subset(plt, pid %in% c('GAA','BAO','FAD','GAG','EBX','ACA','AIA','CB2','FAF','AOA','FAI','ECD','CBV'))) +
  geom_point(aes(col=resection))  +   youri_gg_theme 



ggplot(plt, aes(x=e.vst.egfr, y=e.vst.nrg2, group=pid, label=sid)) +
  geom_point(aes(col=resection)) +
  #geom_text(data=subset(plt, pid %in% c('AAU','AAS','GAQ','FAI'))) + 
  #geom_text(data=subset(plt, pid %in% c('CBV','AAG','HAB','GAQ','AAW','CAF'))) + 
  geom_text(data=subset(plt, pid %in% c('GAA','BAO','FAD','GAG','EBX','ACA','AIA','CB2','FAF','AOA','FAI','ECD','CBV'))) +
  youri_gg_theme 

ggsave('output/figures/rna/egfr-nrg2_labeled.pdf')


ggplot(plt, aes(x=e.vst.egfr, y=e.vst.ngf, label=sid)) +
  geom_point(aes(col=resection)) +
  #geom_smooth(method='lm', formula=y~x) + 
  labs(x = "EGFR expression (VST norm)", y = "NRG2 expression (VST norm)") +
  youri_gg_theme 

ggsave('output/figures/rna/egfr-nrg2_unlabeled.pdf')






ggplot(plt, aes(x=e.vst.egfr, y=e.vst.nrg2, label=sid)) +
  geom_point(aes(col=resection)) +
  #geom_smooth(method='lm', formula=y~x) + 
  labs(x = "EGFR expression (VST norm)", y = "NRG2 expression (VST norm)") +
  youri_gg_theme 






ggplot(subset(plt, pid %in% c('GAA','BAO','FAD','GAG','EBX','ACA','AIA','CB2','FAF','AOA','FAI','ECD','CBV')), aes(x=e.vst.egfr, y=e.vst.nrg2, label=sid)) +
  #geom_point(aes(col=resection)) +
  geom_point(aes(col=vIII.percentage >= 1, shape=resection)) +
  #geom_smooth(method='lm', formula=y~x) + 
  geom_text_repel() +
  labs(x = "EGFR expression (VST norm)", y = "NRG2 expression (VST norm)") +
  youri_gg_theme 






ggplot(plt, aes(x=e.vst.areg, y=e.vst.nrg2, label=sid)) +
  geom_point(aes(col=resection)) +
  #geom_smooth(method='lm', formula=y~x) + 
  labs(x = "AREG expression (VST norm)", y = "NRG2 expression (VST norm)") +
  youri_gg_theme 
ggsave('output/figures/rna/areg-nrg2_unlabeled.pdf')




  



#geom_point(aes(col = m$vIII.percentage > 1))