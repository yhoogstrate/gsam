#!/usr/bin/env R

# ---- basic loading ----

setwd("~/projects/gsam")

# ---- load libs ----

library("readxl")


# ---- load data ----

source("scripts/R/vIII.R")


tmp <- "data/docs/EGFR_and_EGFRvIII_status_of_Neuro-Oncology_study_PMC5762005/EGFR_amp_v_Expression_raw_data.xlsx"
tmp.1 <- read_excel(tmp,sheet=2) # RNA sheet
tmp.1$`PA number` <- gsub(" \\([12]\\)","",tmp.1$`PA number`)
tmp.2 <- read_excel(tmp,sheet=3) # identifier sheet
tmp.2$Id <- NULL
tmp.3 <- tmp.2
tmp.2$`Union 2nd surg.PA number` <- NULL
tmp.3$`Union 1 st surg.PA number` <- NULL
tmp.2$sid <- paste0(tmp.2$patient_inst_code,"1")
tmp.3$sid <- paste0(tmp.3$patient_inst_code,"2")
tmp.2$patient_inst_code <- NULL
tmp.3$patient_inst_code <- NULL

tmp.1 <- merge(tmp.1, tmp.2 , by.x = "PA number", by.y = "Union 1 st surg.PA number", all.x=T)
tmp.1 <- merge(tmp.1, tmp.3 , by.x = "PA number", by.y = "Union 2nd surg.PA number", all.x=T)

tmp.1$sid <- NA
sel <- !is.na(tmp.1$sid.x) & is.na(tmp.1$sid.y)
tmp.1$sid[sel] <- tmp.1$sid.x[sel]
sel <- is.na(tmp.1$sid.x) & !is.na(tmp.1$sid.y)
tmp.1$sid[sel] <- tmp.1$sid.y[sel]

tmp.1$sid.x <- NULL
tmp.1$sid.y <- NULL

tmp.1$`POP4 1` <- NULL
tmp.1$`POP4 2` <- NULL
tmp.1$`RPL30 1` <- NULL
tmp.1$`RPL30 2` <- NULL



sel <- tmp.1$`EGFR 1` == "Undetermined"
tmp.1$`EGFR 1` [sel] <- NA
tmp.1$`EGFR 1`<- as.double(tmp.1$`EGFR 1`)

sel <- tmp.1$`EGFR 2` == "Undetermined"
tmp.1$`EGFR 2` [sel] <- NA
tmp.1$`EGFR 2`<- as.double(tmp.1$`EGFR 2`)



tmp.1$vIII.pct.11 <- 100 - (1/(1 + 2 ^ (tmp.1$`EGFR 1` - tmp.1$`EGFRvIII 1`)) * 100)
tmp.1$vIII.pct.11[!is.na(tmp.1$vIII.pct.11) & tmp.1$`EGFRvIII 1` >= 40] <- 0

tmp.1$vIII.pct.22 <- 100 - (1/(1 + 2 ^ (tmp.1$`EGFR 2` - tmp.1$`EGFRvIII 2`)) * 100)
tmp.1$vIII.pct.22[!is.na(tmp.1$vIII.pct.22) & tmp.1$`EGFRvIII 2` >= 40] <- 0


tmp.1$vIII.pct.avg <- NA

sel <- !is.na(tmp.1$vIII.pct.11) & !is.na(tmp.1$vIII.pct.22)
tmp.1$vIII.pct.avg[sel] <- (tmp.1$vIII.pct.11[sel] + tmp.1$vIII.pct.22[sel]) / 2

sel <- !is.na(tmp.1$vIII.pct.11) & is.na(tmp.1$vIII.pct.22)
tmp.1$vIII.pct.avg[sel] <- tmp.1$vIII.pct.11[sel] 

sel <- is.na(tmp.1$vIII.pct.11) & !is.na(tmp.1$vIII.pct.22)
tmp.1$vIII.pct.avg[sel] <-  tmp.1$vIII.pct.22[sel]


nrow(tmp.1)
tmp.1 <- tmp.1[!is.na(tmp.1$vIII.pct.avg),]
nrow(tmp.1)
tmp.1 <- tmp.1[order(tmp.1$vIII.pct.avg),]


pdf("output/figures/qc/vIII.qPCR.neuro-onco.pdf",height=5,width=8)
plot(c(1, nrow(tmp.1)), c(-10,100),type="n",xlab="Sample; ordered by EGFRvIII percentage", ylab="Percentage EGFRvIII/EGFRwt determined by RT-qPCR", main="qPCR variance by 2x2 RT-qPCRs ~ Neuro-Onco paper 2015")
points(1:nrow(tmp.1), tmp.1$vIII.pct.11, pch=19, cex=0.45)
points(1:nrow(tmp.1), tmp.1$vIII.pct.22, pch=19, cex=0.45)
#text(1, 80, expression("p = 100 - (100 x (1 / (1 + 2^(Ct[wt] - Ct[vIII])))", pos=4)
text(1, 90, expression(y == 100 - (100 * frac(1, 1 + 2^{Ct^{Wt} - Ct^{vIII}}))   ), pos=4,cex=1.1)

text(1:nrow(tmp.1), -7, tmp.1$sid,srt=90,cex=0.3)



dev.off()


# ---- ----


g1 <- tmp.1[,c(1,7,8)]
colnames(g1) <- c("PA number",   "sid", "vIII.pct")
g2 <- tmp.1[,c(1,7,9)]
colnames(g2) <- c("PA number",   "sid","vIII.pct")
g1$col <- "black"
g1$col[duplicated(g1$`PA number`)] <- "red"
g2$col <- "black"
g2$col[duplicated(g2$`PA number`)] <- "red"
g <- rbind(g1, g2)
g <- g[!is.na(g$vIII.pct),]


order <- aggregate(g$vIII.pct, by=list(g$`PA number`), FUN=mean)
order <- order[order(order$x, order$Group.1),]
order$x <- 1:nrow(order)
g <- merge(g, order, by.x="PA number", by.y = "Group.1")
g <- g[order(g$x),]

gmin <- aggregate(g$vIII.pct, by=list(g$`PA number`), FUN=min)
gmin <- gmin[order(match(gmin$Group.1 , g$`PA number`),decreasing=F),]

gmax <- aggregate(g$vIII.pct, by=list(g$`PA number`), FUN=max)
gmax <- gmax[order(match(gmax$Group.1 , g$`PA number`),decreasing=F),]


plot(c(1, max(g$x)), c(-10,100),type="n",xlab="Sample; ordered by EGFRvIII percentage", ylab="Percentage EGFRvIII/EGFRwt determined by RT-qPCR", main="qPCR variance by 2x2 RT-qPCRs ~ Neuro-Onco paper 2015")
for(n in 1:length(gmin$x)) {
  lines(c(n, n), c(gmin$x[n], gmax$x[n]), col="gray")
}
points(g$x ,g$vIII.pct, col=g$col, pch=19, cex=0.45)
text(g$x, -7, g$sid,srt=90,cex=0.6)



is <- intersect(g$sid , gsam.rna.metadata$sid)
qpcr.gsam <- gsam.rna.metadata[is %in% is,colnames(gsam.rna.metadata) %in% c("sid","qPCR.percent.EGFR.vIII")]
qpcr.gsam <- qpcr.gsam[!is.na(qpcr.gsam$qPCR.percent.EGFR.vIII),]
#g[match(qpcr.gsam$sid, g$sid),]$sid == qpcr.gsam$sid
points( g[match(qpcr.gsam$sid, g$sid),]$x , qpcr.gsam$qPCR.percent.EGFR.vIII, col="purple")



# CAB2
# CAS1



# ---- calc correlation ----


