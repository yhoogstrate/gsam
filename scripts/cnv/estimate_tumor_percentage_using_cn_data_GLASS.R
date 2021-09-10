#!/usr/bin/env R 

# load libs ----

library(dplyr)

# load data ----

source("scripts/R/chrom_sizes.R")

# metadata, samples
sids <- read.table("data/gsam/data/GLASS_GBM_R1-R2/glass_transcript_counts.txt",sep="\t",header=T,stringsAsFactors = F) %>% 
  colnames() %>%
  gsub(".","-",. , fixed=T) %>%
  gsub("^([^-]+-[^-]+-[^-]+-[^-]+).+$","\\1", . ) %>%
  purrr::discard(. %in% c("target_id","length"))



# SELECT * FROM syn18477378 [variants_gatk_seg]
dat <- read.csv('data/glass_cnv/variants_gatk_seg_glass.csv')

dat <- dat %>%
  dplyr::mutate(ROW_ID=NULL, ROW_VERSION=NULL) %>%
  dplyr::mutate(chrom = paste0("chr",chrom) ) %>%
  dplyr::filter(chrom %in% c('chr23', 'chr24') == F) %>% # chr23 is most likely chrX
  dplyr::mutate(sample.short = gsub("","", gsub("^([^-]+-[^-]+-[^-]+-[^-]+).+$","\\1", aliquot_barcode ) ) ) %>%
  dplyr::filter(sample.short %in% sids) %>%
  dplyr::filter(grepl('WXS', aliquot_barcode)) %>% # WGS are not as clean, and only present for TCGA
  #dplyr::filter(aliquot_barcode != "TCGA-06-0190-R1-01D-WGS-P20F5P") %>%
  #dplyr::filter(aliquot_barcode != "TCGA-06-0190-TP-01D-WGS-BM9MHX")
  dplyr::filter(aliquot_barcode != "TCGA-06-0190-R1-02D-WXS-LE9JNK") %>%
  dplyr::filter(aliquot_barcode != 'GLSS-SM-R101-R1-02D-WXS-PUTOC9') %>%
  dplyr::filter(aliquot_barcode != 'GLSS-SM-R101-TP-02D-WXS-TKL136') %>%
  dplyr::filter(aliquot_barcode != 'GLSS-SM-R107-R1-02D-WXS-W7UY7B') %>%
  dplyr::filter(aliquot_barcode != 'GLSS-SM-R110-TP-02D-WXS-BQFEME') %>%
  dplyr::filter(aliquot_barcode != 'TCGA-06-0125-R1-02D-WXS-FB1H59') %>%
  dplyr::filter(aliquot_barcode != 'TCGA-06-0125-TP-02D-WXS-FQUL0U') %>%
  dplyr::filter(aliquot_barcode != 'TCGA-06-0210-R1-02D-WXS-5M0ICH') %>%
  dplyr::filter(aliquot_barcode != 'TCGA-06-0211-R1-03D-WXS-M22HJ0')

  #dplyr::filter(grepl("TCGA-14-1034",aliquot_barcode))

dat <- dat %>%
  dplyr::mutate(aliquot_barcode = as.factor(aliquot_barcode))






# calc distances ----



tpc.estimate = data.frame()
for(bc in levels(dat$aliquot_barcode)) {
  
  a = dat %>%
    dplyr::filter(aliquot_barcode == bc) %>%
    dplyr::mutate(chrom.offset = chrs_hg19_s[ chrom ] ) %>%
    dplyr::mutate(x1 = start + chrom.offset, x2 = end + chrom.offset) %>%
    dplyr::mutate(outlier = num_points <= 60 ) %>%
    dplyr::filter(outlier == F) %>%
    dplyr::filter(log2_copy_ratio < 1.1) %>%
    dplyr::filter(log2_copy_ratio > -2.0)

  #m <- median ( rep(a$log2_copy_ratio, a$num_points) )
  #a$log2_copy_ratio <- a$log2_copy_ratio - m
  

  #plot(c(min(a$x1), max(a$x2)) , c(-y.scale, y.scale) , type = 'n' ) 
  png(paste0("output/figures/cnv/glass/", bc ,".tpc.estimate.png"),width=11 * 100,height=6 * 100)
  
  plot(c(min(a$x1), max(a$x2)) , c(-2, 2) , type = 'n' )
  abline(v = chrs_hg19_e, col="gray")
  for(i in 1:nrow(a)) {
    b = a[i,]
    lines(c(b$x1, b$x2), c(b$log2_copy_ratio, b$log2_copy_ratio) , col = (b$num_points > 60) + 1 , lwd=2 )
  }
  
  if(bc == "GLSS-HF-2548-R1-01D-WXS-DN2HD4") { a <- a %>% dplyr::filter(chrom %in% c('chr2','chr7','chr10') == T) }
  else if(bc == "GLSS-HF-2548-TP-01D-WXS-CPWUK8") { a <- a %>% dplyr::filter(chrom %in% c('chr2','chr7','chr10') == T) }
  else if(bc == "GLSS-HF-2829-R1-01D-WXS-W19NLP") { a <- a %>% dplyr::filter(chrom %in% c('chr10') == T) }
  else if(bc == "GLSS-HF-2829-TP-01D-WXS-1OMY6S") { a <- a %>% dplyr::filter(chrom %in% c('chr10') == T) }
  else if(bc == "GLSS-HF-2869-TP-01D-WXS-PNBQFR") { a <- a %>% dplyr::filter(chrom %in% c('chr10') == T) }
  else if(bc == "GLSS-HF-2998-R1-01D-WXS-LCQISE") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10') == T) }
  else if(bc == "GLSS-HF-2998-TP-01D-WXS-0R9KZC") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10') == T) }
  else if(bc == "GLSS-HF-3081-R2-01D-WXS-4LBL2G") { a <- a %>% dplyr::filter(chrom %in% c('chr10') & start >= 43087061 ) }
  else if(bc == "GLSS-HF-3081-TP-01D-WXS-V0EK51") { a <- a %>% dplyr::filter(chrom %in% c('chr10') & start >= 43087061 ) }
  else if(bc == "GLSS-HF-3162-R1-01D-WXS-3HVQJ6") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10')) }
  else if(bc == "GLSS-HF-3162-TP-01D-WXS-HTQZ6B") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10')) }
  else if(bc == "GLSS-SM-R056-R2-01D-WXS-A7U2KL") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R056-TP-01D-WXS-1BJVV2") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R060-R1-01D-WXS-MHAES3") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr7','chr10')) }
  else if(bc == "GLSS-SM-R060-TP-01D-WXS-WWG8SV") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R061-R1-01D-WXS-T05F82") { a <- a %>% dplyr::filter(chrom %in% c('chr2','chr3','chr10')) }
  else if(bc == "GLSS-SM-R061-TP-01D-WXS-UZBY1Q") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R063-R1-01D-WXS-SF2PHH") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) }
  else if(bc == "GLSS-SM-R063-TP-01D-WXS-UQC9Y9") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) }
  else if(bc == "GLSS-SM-R064-R1-01D-WXS-GFA2BL") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R064-TP-01D-WXS-16OKWU") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R067-R1-01D-WXS-WG6N7X") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr8','chr9','chr10')) }
  else if(bc == "GLSS-SM-R067-TP-01D-WXS-LTGC1O") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr8','chr9','chr10')) }
  else if(bc == "GLSS-SM-R068-R1-01D-WXS-KBROK6") { a <- a %>% dplyr::filter(chrom %in% c('chr8','chr10')) }
  else if(bc == "GLSS-SM-R068-TP-01D-WXS-5G9ZR6") { a <- a %>% dplyr::filter(chrom %in% c('chr8','chr10')) }
  else if(bc == "GLSS-SM-R070-R1-01D-WXS-83NH8Q") { a <- a %>% dplyr::filter(chrom %in% c('chr1','chr9')) }
  else if(bc == "GLSS-SM-R070-TP-01D-WXS-VLPMYS") { a <- a %>% dplyr::filter(chrom %in% c('chr1','chr9')) }
  else if(bc == "GLSS-SM-R071-R1-01D-WXS-LZGI7S") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R071-TP-01D-WXS-45QCP2") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R072-R1-01D-WXS-XOT8LM") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R072-TP-01D-WXS-QZZEMI") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R099-R1-01D-WXS-I52D70") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R099-TP-01D-WXS-5Z45O8") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R100-R1-01D-WXS-O61DWB") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr13','chr14')) }
  else if(bc == "GLSS-SM-R100-TP-01D-WXS-9D7AET") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr13','chr14')) }
  else if(bc == "GLSS-SM-R101-R1-01D-WXS-04ZGTS") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R101-TP-01D-WXS-QXZ8TD") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R103-R1-01D-WXS-0YX7OK") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R103-TP-01D-WXS-TBB9NO") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R106-R1-01D-WXS-B6X7IN") { a <- a %>% dplyr::filter(chrom %in% c('chr7') & log2_copy_ratio > 0) }
  else if(bc == "GLSS-SM-R106-TP-01D-WXS-8H4HYR") { a <- a %>% dplyr::filter(chrom %in% c('chr7') & log2_copy_ratio > 0) }
  else if(bc == "GLSS-SM-R107-R1-01D-WXS-66033B") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R107-TP-01D-WXS-NIBO3N") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "GLSS-SM-R108-R1-01D-WXS-YRYRWN") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R108-TP-01D-WXS-KC4OXZ") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "GLSS-SM-R110-R1-01D-WXS-I5HNRJ") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) }
  else if(bc == "GLSS-SM-R110-TP-01D-WXS-5XSI28") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) }
  else if(bc == "GLSS-SM-R107-TP-01D-WXS-NIBO3N") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "TCGA-06-0190-R1-01D-WXS-ODJETQ") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "TCGA-06-0190-TP-01D-WXS-6BX2M1") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "TCGA-06-0210-R1-01D-WXS-4MJVD2") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }
  else if(bc == "TCGA-06-0210-TP-01D-WXS-MGZTBD") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) }

  else if(bc == "TCGA-06-0211-R1-02D-WXS-X1GNKU") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  else if(bc == "TCGA-06-0211-TP-01D-WXS-AS2XJL") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) }
  
  
  else if(bc == "TCGA-14-1034-R1-01D-WXS-NUKYNI") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr13','chr22')) }
  else if(bc == "TCGA-14-1034-TP-01D-WXS-B0AXT9") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr13','chr22')) }
  
  
  

  # code
  out <- data.frame()
  for(frac in  1:100 / 100) {
    fc.p.4 <- ((1 - frac) * 2 + frac * 4) / 2
    fc.p.3 <- ((1 - frac) * 2 + frac * 3) / 2
    fc.n.1 <- ((1 - frac) * 2 + frac * 1) / 2
    
    lfc.p.4 <- log2(fc.p.4)
    lfc.p.3 <- log2(fc.p.3)
    lfc.n.1 <- log2(fc.n.1)
    
    dists <- c()
    dist <- 0
    for(i in 1:nrow(a)) {
      
      e <- a[i,]
      
      d <- c(
        (((e$log2_copy_ratio - 0) * e$num_points)^2),
        (((e$log2_copy_ratio - lfc.p.4) * e$num_points)^2) * 1.1 * 999, # penalize, do never prefer
        (((e$log2_copy_ratio - lfc.p.3) * e$num_points)^2),
        (((e$log2_copy_ratio - lfc.n.1) * e$num_points)^2)
      )
      d <- min(d)
      dists <- c(dists, d)
      
      dist <- dist + d
    }
    
    out <- rbind(out, data.frame(pct = frac * 100,
                                 lfc.3p = lfc.p.3,
                                 lfc.4p = lfc.p.4,
                                 lfc.n = lfc.n.1,
                                 dist = dist) )
  }
  
  
  r = out %>%
    dplyr::arrange(dist) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(sample = bc,
                  sample.short = gsub("^([^-]+-[^-]+-[^-]+-[^-]+).+$","\\1",bc))
  tpc.estimate <- rbind(tpc.estimate, r)
  
  abline(h=r$lfc.3p, col="blue")
  abline(h=r$lfc.4p, col="blue")
  abline(h=r$lfc.n, col="blue")
  
  text(0, 2, paste0("Estimated tumor percentage [",bc,"]: ",r$pct , "%" )  ,pos=4)
  
  dev.off()
  
}

write.table(tpc.estimate , "output/tables/cnv/tumor.percentage.estimate_glass.txt")





