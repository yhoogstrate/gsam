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

tmp.1$vIII.pct.12 <- 100 - (1/(1 + 2 ^ (tmp.1$`EGFR 1` - tmp.1$`EGFRvIII 2`)) * 100)
tmp.1$vIII.pct.12[!is.na(tmp.1$vIII.pct.12) & tmp.1$`EGFRvIII 2` >= 40] <- 0

tmp.1$vIII.pct.21 <- 100 - (1/(1 + 2 ^ (tmp.1$`EGFR 2` - tmp.1$`EGFRvIII 1`)) * 100)
tmp.1$vIII.pct.21[!is.na(tmp.1$vIII.pct.21) & tmp.1$`EGFRvIII 1` >= 40] <- 0

tmp.1$vIII.pct.22 <- 100 - (1/(1 + 2 ^ (tmp.1$`EGFR 2` - tmp.1$`EGFRvIII 2`)) * 100)
tmp.1$vIII.pct.22[!is.na(tmp.1$vIII.pct.22) & tmp.1$`EGFRvIII 2` >= 40] <- 0


tmp.1$EGFR.avg <- NA
sel <- !is.na(tmp.1$`EGFR 1`) & !is.na(tmp.1$`EGFR 2`)
tmp.1$EGFR.avg[sel] <- (tmp.1$`EGFR 1`[sel] + tmp.1$`EGFR 2`[sel]) / 2
sel <- !is.na(tmp.1$`EGFR 1`) & is.na(tmp.1$`EGFR 2`)
tmp.1$EGFR.avg[sel] <- tmp.1$`EGFR 1`[sel]
sel <- is.na(tmp.1$`EGFR 1`) & !is.na(tmp.1$`EGFR 2`)
tmp.1$EGFR.avg[sel] <- tmp.1$`EGFR 2`[sel]


tmp.1$EGFRvIII.avg <- NA
sel <- !is.na(tmp.1$`EGFRvIII 1`) & !is.na(tmp.1$`EGFRvIII 2`)
tmp.1$EGFRvIII.avg[sel] <- (tmp.1$`EGFRvIII 1`[sel] + tmp.1$`EGFRvIII 2`[sel]) / 2
sel <- !is.na(tmp.1$`EGFRvIII 1`) & is.na(tmp.1$`EGFRvIII 2`)
tmp.1$EGFRvIII.avg[sel] <- tmp.1$`EGFRvIII 1`[sel]
sel <- is.na(tmp.1$`EGFRvIII 1`) & !is.na(tmp.1$`EGFRvIII 2`)
tmp.1$EGFRvIII.avg[sel] <- tmp.1$`EGFRvIII 2`[sel]

# het is deze:
tmp.1$vIII.pct.avg <-  100 - (1/(1 + 2 ^ (tmp.1$EGFR.avg - tmp.1$EGFRvIII.avg)) * 100)
tmp.1$vIII.pct.avg[tmp.1$EGFRvIII.avg >= 40 ] <-  0

nrow(tmp.1)
tmp.1 <- tmp.1[!is.na(tmp.1$vIII.pct.avg),]
nrow(tmp.1)
tmp.1 <- tmp.1[order(tmp.1$vIII.pct.avg),]


pdf("output/figures/qc/vIII.qPCR.neuro-onco.pdf",height=5,width=8)
plot(c(1, nrow(tmp.1)), c(-10,100),type="n",xlab="Sample; ordered by EGFRvIII percentage", ylab="Percentage EGFRvIII/EGFRwt determined by RT-qPCR", main="qPCR variance by 2x2 RT-qPCRs ~ Neuro-Onco paper 2015")
for(x in 1:nrow(tmp.1)) {
  #x = 125
  d <- c(tmp.1$vIII.pct.11[x], tmp.1$vIII.pct.12[x], tmp.1$vIII.pct.21[x], tmp.1$vIII.pct.22[x])
  d <- d[!is.na(d)]
  d <- sort(d)
  
  if(length(d) == 4) {
    lines(c(x, x), c(d[2], d[3]), col=rgb(0,0,0,0.20),lwd=3)
    d <- c(d[1], d[4])
  }
  
  if(length(d) == 2) {
    lines(c(x, x), c(d[1], d[2]), col=rgb(0,0,0,0.15),lwd=2.5)
    d <- c(d[1], d[4])
  }
  
  
}
points(1:nrow(tmp.1), tmp.1$vIII.pct.avg, pch=19, cex=0.45)
#text(1, 80, expression("p = 100 - (100 x (1 / (1 + 2^(Ct[wt] - Ct[vIII])))", pos=4)
text(1, 90, expression(y == 100 - (100 * frac(1, 1 + 2^{Ct^{Wt} - Ct^{vIII}}))   ), pos=4,cex=1.1)

text(1:nrow(tmp.1), -7, tmp.1$sid,srt=90,cex=0.3)

legend(0, 75, "%vIII by mean Ct Values",pch=19, cex=1)

dev.off()



# ---- calc correlation ----



