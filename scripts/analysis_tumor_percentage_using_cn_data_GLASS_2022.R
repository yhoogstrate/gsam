#!/usr/bin/env R 

# load libs ----


# load data ----

source("scripts/R/chrom_sizes.R")

if(!exists('glass.gbm.rnaseq.metadata.all.samples')) {
  source('scripts/load_glass_expression_data.R')
}


# analysis ----


sel <- glass.gbm.rnaseq.metadata.all.samples |>
  dplyr::filter(lapply(excluded, length) == 0) 



#'@warning there are quite some with both WGS and WXS

# data from synapse portal
dat <- 'data/gsam/data/GLASS_GBM_R1-R2/variants_gatk_seg.syn31121137.all_samlpes.tsv' |> 
  read.table(sep="\t",header=T,stringsAsFactors = F) |> 
  dplyr::mutate(ROW_ID=NULL, ROW_VERSION=NULL) |> 
  dplyr::mutate(portion_barcode = gsub("^([^\\-]+)-([^\\-]+)-([^\\-]+)-([^\\-]+)-([0-9]+).+$$","\\1-\\2-\\3-\\4-\\5",aliquot_barcode)) |> 
  dplyr::filter(portion_barcode %in% sel$portion_barcode) |> 
  dplyr::mutate(chrom = paste0("chr",chrom) )  |> 
  #dplyr::filter(grepl('WXS', aliquot_barcode)) %>% # WGS are not as clean, and only present for TCGA  
  dplyr::filter(chrom %in% c('chr23', 'chr24') == F) |>  # chr23 is most likely chrX
  dplyr::filter(aliquot_barcode != "GLSS-CU-R017-R1-01D-WXS-C0CBCW") |>  # too noisy data
  dplyr::filter(aliquot_barcode != "GLSS-HF-3118-R1-01D-WXS-QRF6VZ") |>  # too noisy data
  dplyr::filter(aliquot_barcode != "GLSS-HF-3162-R1-01D-WXS-3HVQJ6") |> 
  dplyr::filter(aliquot_barcode != "GLSS-HF-3162-TP-01D-WXS-HTQZ6B") |> 
  dplyr::filter(aliquot_barcode != "GLSS-HF-EE74-R1-01D-WXS-5AS51O") |> 
  dplyr::mutate(aliquot_barcode = as.factor(aliquot_barcode)) 


chr17.err <- dat |> 
  dplyr::filter(chrom == "chr17") |> 
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::filter(end == max(end))  |>  dplyr::filter(end > 140000000) |>  dplyr::pull(aliquot_barcode)


plt <- dat |> 
  dplyr::filter(chrom == "chr17") |> 
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::filter(end == max(end)) 
plot(sort(plt$end))


plt <- dat |> 
  dplyr::filter(chrom == "chr1") |> 
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::filter(end == max(end)) 
plot(sort(plt$end))


plt <- dat |> 
  dplyr::filter(chrom == "chr2") |> 
  dplyr::group_by(aliquot_barcode) |> 
  dplyr::filter(end == max(end)) |> 
  dplyr::mutate(chr17.err = aliquot_barcode %in% chr17.err)
ggplot(plt, aes(x=reorder(aliquot_barcode,end), y=end, col=chr17.err)) + 
  geom_point() +
  #geom_text(angle=90,size=2, aes(y=(max(plt$end) - min(plt$end))/2)) +
  theme_bw()



# calc distances ----



tpc.estimate = data.frame()
for(bc in levels(dat$aliquot_barcode)) {
  #bc = "TCGA-14-1402-R1-01D-WGS-2EHMQ2" # plot ex
  print(bc)
  
  m <- glass.gbm.rnaseq.metadata.all.samples |> 
    dplyr::filter(portion_barcode == gsub("^([^\\-]+)-([^\\-]+)-([^\\-]+)-([^\\-]+)-([0-9]+).+$$","\\1-\\2-\\3-\\4-\\5",bc))
  
  p.tpc.estimate <- m$purity.synapse.rna
  p.tpc.hoogstrate.2021 <- m$tumour.percentage.dna.imputed.rf.2021 / 100
  
  center <-  dat %>%
    dplyr::filter(aliquot_barcode == bc)
  center <- rep(center$log2_copy_ratio, center$num_points)
  
  a = dat %>%
    dplyr::filter(aliquot_barcode == bc) %>%
    dplyr::mutate(chrom.offset = chrs_hg19_s[ chrom ] ) %>%
    dplyr::mutate(x1 = start + chrom.offset, x2 = end + chrom.offset) %>%
    dplyr::mutate(outlier = num_points <= 60 ) %>%
    dplyr::filter(outlier == F) %>%
    dplyr::mutate(log2_copy_ratio = log2_copy_ratio - median(center)) %>%
    dplyr::filter(log2_copy_ratio < 1.1) %>%
    dplyr::filter(log2_copy_ratio > -2.0)


  

  png(paste0("output/figures/cnv/glass/2022/", bc ,".tpc.estimate.png"),width=11 * 100,height=6 * 100)
  
  
  plot(c(min(a$x1), max(a$x2)) , c(-2, 2) ,  type = 'n' )
  abline(v = chrs_hg19_e, col="gray")
  abline(h = 0, col="black",lty=2)
  for(i in 1:nrow(a)) {
    b = a[i,]
    lines(c(b$x1, b$x2), c(b$log2_copy_ratio, b$log2_copy_ratio) , col = (b$num_points > 60) + 1 , lwd=2.5 )
    #lines(c(b$x1, b$x2), c(b$log2_copy_ratio, b$log2_copy_ratio) , col = as.numeric(b$chrom == "chr17") + 1, lwd=2 )
  }
  
  
  # chr17 84622 190.863.195
  if(F) {
  } else if(bc == "GLSS-19-0266-R1-01D-WXS-0AUVTN") { a <- a %>% dplyr::filter(chrom %in% c('chr3', 'chr7','chr10'))
  } else if(bc == "GLSS-19-0273-R1-01D-WXS-ZEKWL8") { a <- a %>% dplyr::filter(chrom %in% c('chr3', 'chr7','chr10',"chr13","chr16"))
  } else if(bc == "GLSS-CU-P003-TP-01D-WXS-RDWJOX") { a <- a %>% dplyr::filter(chrom %in% c('chr4', 'chr5','chr6',"chr7","chr8","chr9","chr10","chr11"))
  } else if(bc == "GLSS-CU-P021-TP-01D-WXS-Y82P1D") { a <- a %>% dplyr::filter(chrom %in% c('chr2', 'chr5','chr15'))
  } else if(bc == "GLSS-CU-P046-TP-01D-WXS-T8ZFA7") { a <- a %>% dplyr::filter(chrom %in% c('chr2', 'chr6','chr17'))
  } else if(bc == "GLSS-CU-P055-TP-01D-WXS-FS5YF2") { a <- a %>% dplyr::filter(chrom %in% c('chr2', 'chr3','chr4',"chr5","chr10"))
  } else if(bc == "GLSS-CU-P056-R2-01D-WXS-MWYTTF") { a <- a %>% dplyr::filter(chrom %in% c('chr2', 'chr7','chr13'))
  } else if(bc == "GLSS-CU-P103-R1-01D-WXS-PFTRVM") { a <- a %>% dplyr::filter(chrom %in% c('chr10', 'chr13'))
  } else if(bc == "GLSS-CU-R005-R1-01D-WXS-78REC1") { a <- a %>% dplyr::filter(chrom %in% c("chr7","chr9",'chr10'))
  } else if(bc == "GLSS-CU-R006-R1-01D-WXS-OQXLD2") { a <- a %>% dplyr::filter(chrom %in% c("chr4","chr7",'chr10','chr13','chr17'))
  } else if(bc == "GLSS-CU-R010-R1-01D-WXS-LC99UJ") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-CU-R010-TP-01D-WXS-MWMVOK") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-CU-R018-R1-01D-WXS-PKDTNA") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-CU-R019-R1-01D-WXS-NAXQ9K") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-HF-2548-R1-01D-WXS-DN2HD4") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10') == T) 
  } else if(bc == "GLSS-HF-2548-TP-01D-WXS-CPWUK8") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-2829-R1-01D-WXS-W19NLP") { a <- a %>% dplyr::filter(chrom %in% c('chr10') == T) 
  } else if(bc == "GLSS-HF-2829-TP-01D-WXS-1OMY6S") { a <- a %>% dplyr::filter(chrom %in% c('chr10') == T) 
  } else if(bc == "GLSS-HF-2869-TP-01D-WXS-PNBQFR") { a <- a %>% dplyr::filter(chrom %in% c('chr10') == T) 
  } else if(bc == "GLSS-HF-2919-R1-01D-WXS-CG67L5") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-HF-2919-TP-01D-WXS-MIYI2Q") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-HF-2934-TP-01D-WXS-OIDP86") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-2998-R1-01D-WXS-LCQISE") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10') == T) 
  } else if(bc == "GLSS-HF-2998-TP-01D-WXS-0R9KZC") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10') == T) 
  } else if(bc == "GLSS-HF-3050-R1-01D-WXS-6ELHQ3") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-3050-TP-01D-WXS-ECTQ63") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-3081-R2-01D-WXS-4LBL2G") { a <- a %>% dplyr::filter(chrom %in% c('chr10') & start >= 43087061 ) 
  } else if(bc == "GLSS-HF-3081-TP-01D-WXS-V0EK51") { a <- a %>% dplyr::filter(chrom %in% c('chr10') & start >= 43087061 ) 
  } else if(bc == "GLSS-HF-3162-R1-01D-WXS-3HVQJ6") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10')) 
  } else if(bc == "GLSS-HF-3162-TP-01D-WXS-HTQZ6B") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10')) 
  } else if(bc == "GLSS-HF-50F3-R2-01D-WXS-HBGS0K") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-9A7A-TP-01D-WXS-W98DVW") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-DE05-TP-01D-WXS-44QP3Q") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-DF35-R1-01D-WXS-B4D43W") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HF-EE74-TP-01D-WXS-E5C4YL") { a <- a %>% dplyr::filter(chrom %in% c('chr7'))
  } else if(bc == "GLSS-HF-EE77-R3-01D-WXS-ZRNNYZ") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-HK-0003-R1-01D-WGS-R7P485") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr8','chr10','chr13','chr15'))
  } else if(bc == "GLSS-LU-00B9-R1-01D-WXS-LM12XS") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr19'))
  } else if(bc == "GLSS-LU-00C2-R1-01D-WXS-YYR15P") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr13','chr15'))
  } else if(bc == "GLSS-LU-00C4-R1-01D-WXS-KYPFOQ") { a <- a %>% dplyr::filter(chrom %in% c('chr10','chr13'))
  } else if(bc == "GLSS-LU-00C4-TP-01D-WXS-S4OOP1") { a <- a %>% dplyr::filter(chrom %in% c('chr10','chr14'))
  } else if(bc == "GLSS-LU-00C7-TP-01D-WXS-WFSIQ7") { a <- a %>% dplyr::filter(chrom %in% c('chr13'))
  } else if(bc == "GLSS-LU-0B12-R1-01D-WXS-NXWOWA") { a <- a %>% dplyr::filter(chrom %in% c('chr7'))
  } else if(bc == "GLSS-LU-0B13-R1-01D-WXS-F05JM5") { a <- a %>% dplyr::filter(chrom %in% c('chr9','chr12','chr13','chr17'))
  } else if(bc == "GLSS-MD-0003-R3-01D-WGS-I78LK2") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-MD-0014-R1-01D-WGS-GLH38H") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-MD-0022-R1-01D-WGS-F5KYZ6") { a <- a %>% dplyr::filter(chrom %in% c('chr7'))
  } else if(bc == "GLSS-MD-0022-TP-01D-WGS-ZAET2D") { a <- a %>% dplyr::filter(chrom %in% c('chr7'))
  } else if(bc == "GLSS-MD-0023-R1-01D-WGS-DQH2SF") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr7','chr9','chr10'))
  } else if(bc == "GLSS-MD-0026-TP-01D-WGS-2PPGFK") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr14'))
  } else if(bc == "GLSS-MD-0042-TP-01D-WGS-R5UYLI") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr13'))
  } else if(bc == "GLSS-MD-LP04-R1-01D-WXS-V7XG09") { a <- a %>% dplyr::filter(chrom %in% c('chr1','chr7','chr10','chr13'))
  } else if(bc == "GLSS-MD-LP04-TP-01D-WXS-05SV1D") { a <- a %>% dplyr::filter(chrom %in% c('chr1','chr7','chr10','chr13'))
  } else if(bc == "GLSS-SM-R056-R2-01D-WXS-A7U2KL") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R056-TP-01D-WXS-1BJVV2") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R060-R1-01D-WXS-MHAES3") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr7','chr10')) 
  } else if(bc == "GLSS-SM-R060-R3-01D-WXS-UXGCSX") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-SM-R060-TP-01D-WXS-WWG8SV") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R061-R1-01D-WXS-T05F82") { a <- a %>% dplyr::filter(chrom %in% c('chr2','chr3','chr10')) 
  } else if(bc == "GLSS-SM-R061-TP-01D-WXS-UZBY1Q") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R063-R1-01D-WXS-SF2PHH") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) 
  } else if(bc == "GLSS-SM-R063-TP-01D-WXS-UQC9Y9") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) 
  } else if(bc == "GLSS-SM-R064-R1-01D-WXS-GFA2BL") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R064-R2-01D-WXS-606Z2Z") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-SM-R064-TP-01D-WXS-16OKWU") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R066-R1-01D-WXS-2R3RRY") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-SM-R066-TP-01D-WXS-YJN4OE") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-SM-R067-R1-01D-WXS-WG6N7X") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr8','chr9','chr10')) 
  } else if(bc == "GLSS-SM-R067-TP-01D-WXS-LTGC1O") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr8','chr9','chr10')) 
  } else if(bc == "GLSS-SM-R068-R1-01D-WXS-KBROK6") { a <- a %>% dplyr::filter(chrom %in% c('chr8','chr10')) 
  } else if(bc == "GLSS-SM-R068-TP-01D-WXS-5G9ZR6") { a <- a %>% dplyr::filter(chrom %in% c('chr8','chr10')) 
  } else if(bc == "GLSS-SM-R070-R1-01D-WXS-83NH8Q") { a <- a %>% dplyr::filter(chrom %in% c('chr1','chr9')) 
  } else if(bc == "GLSS-SM-R070-TP-01D-WXS-VLPMYS") { a <- a %>% dplyr::filter(chrom %in% c('chr1','chr9')) 
  } else if(bc == "GLSS-SM-R071-R1-01D-WXS-LZGI7S") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R071-TP-01D-WXS-45QCP2") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R072-R1-01D-WXS-XOT8LM") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R072-TP-01D-WXS-QZZEMI") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R099-R1-01D-WXS-I52D70") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) 
  } else if(bc == "GLSS-SM-R099-TP-01D-WXS-5Z45O8") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) 
  } else if(bc == "GLSS-SM-R100-R1-01D-WXS-O61DWB") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr13','chr14')) 
  } else if(bc == "GLSS-SM-R100-TP-01D-WXS-9D7AET") { a <- a %>% dplyr::filter(chrom %in% c('chr6','chr13','chr14')) 
  } else if(bc == "GLSS-SM-R101-R1-01D-WXS-04ZGTS") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R101-R1-02D-WXS-PUTOC9") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-SM-R101-TP-01D-WXS-QXZ8TD") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R101-TP-02D-WXS-TKL136") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "GLSS-SM-R103-R1-01D-WXS-0YX7OK") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R103-TP-01D-WXS-TBB9NO") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R106-R1-01D-WXS-B6X7IN") { a <- a %>% dplyr::filter(chrom %in% c('chr10') ) 
  } else if(bc == "GLSS-SM-R106-TP-01D-WXS-8H4HYR") { a <- a %>% dplyr::filter(chrom %in% c('chr10') ) 
  } else if(bc == "GLSS-SM-R107-R1-01D-WXS-66033B") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R107-R1-02D-WXS-W7UY7B") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-SM-R107-TP-01D-WXS-NIBO3N") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "GLSS-SM-R108-R1-01D-WXS-YRYRWN") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R108-TP-01D-WXS-KC4OXZ") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "GLSS-SM-R110-R1-01D-WXS-I5HNRJ") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) 
  } else if(bc == "GLSS-SM-R110-TP-01D-WXS-5XSI28") { a <- a %>% dplyr::filter(chrom %in% c('chr7')) 
  } else if(bc == "GLSS-SN-0001-TP-01D-WGS-J7WGL3") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "GLSS-SN-0002-R1-01D-WGS-OMWX5F") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr14','chr18'))
  } else if(bc == "GLSS-SN-0004-R1-01D-WGS-FVYBOQ") { a <- a %>% dplyr::filter(chrom %in% c('chr4','chr5','chr10'))
  } else if(bc == "GLSS-SN-0009-TP-01D-WGS-Z9WB19") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9'))
  } else if(bc == "GLSS-SN-0009-R1-01D-WGS-X1YPCJ") { a <- a %>% dplyr::filter(chrom %in% c('chr7'))
  } else if(bc == "TCGA-06-0152-R1-01D-WGS-EK2VYI") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10'))
  } else if(bc == "TCGA-06-0190-R1-01D-WXS-ODJETQ") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "TCGA-06-0190-TP-01D-WXS-6BX2M1") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "TCGA-06-0210-R1-01D-WGS-PS45ZE") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10'))
  } else if(bc == "TCGA-06-0210-R1-01D-WXS-4MJVD2") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "TCGA-06-0210-TP-01D-WXS-MGZTBD") { a <- a %>% dplyr::filter(chrom %in% c('chr10')) 
  } else if(bc == "TCGA-06-0211-R1-02D-WGS-32OOWA") { a <- a %>% dplyr::filter(chrom %in% c('chr9','chr10'))
  } else if(bc == "TCGA-06-0211-R1-02D-WXS-X1GNKU") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10')) 
  } else if(bc == "TCGA-06-0211-TP-01D-WGS-JJPE5C") { a <- a %>% dplyr::filter(chrom %in% c('chr9','chr10'))
  } else if(bc == "TCGA-06-0211-TP-01D-WXS-AS2XJL") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "TCGA-06-0211-TP-01D-WXS-AS2XJL") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr9','chr10'))
  } else if(bc == "TCGA-14-1034-R1-01D-WGS-L9V6H0") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "TCGA-14-1034-R1-01D-WXS-NUKYNI") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "TCGA-14-1034-R1-01D-WXS-NUKYNI") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr13','chr22')) 
  } else if(bc == "TCGA-14-1034-TP-01D-WGS-2521IS") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "TCGA-14-1034-TP-01D-WXS-B0AXT9") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "TCGA-14-1034-R1-01D-WGS-L9V6H0") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10','chr13'))
  } else if(bc == "TCGA-14-1402-R1-01D-WGS-2EHMQ2") { a <- a %>% dplyr::filter(chrom %in% c('chr10'))
  } else if(bc == "TCGA-14-1402-R1-01D-WXS-QGNTYA") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  } else if(bc == "TCGA-19-0957-TP-01D-WXS-H59SRY") { a <- a %>% dplyr::filter(chrom %in% c('chr7','chr10'))
  
  }
  
  

  print(paste0(bc, " -> ", str_c(unique(a$chrom), collapse = ",")))
  

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
    dplyr::mutate(aliquot_barcode.wgs = bc,
                  #sample.short = gsub("^([^-]+-[^-]+-[^-]+-[^-]+).+$","\\1",bc),
                  portion_barcode == gsub("^([^\\-]+)-([^\\-]+)-([^\\-]+)-([^\\-]+)-([0-9]+).+$$","\\1-\\2-\\3-\\4-\\5",bc),
                  estimate.purity = p.tpc.estimate,
                  hoogstrate.rf.purity.2021 = p.tpc.hoogstrate.2021
                  )
  tpc.estimate <- rbind(tpc.estimate, r)
  
  abline(h=r$lfc.3p, col="blue")
  abline(h=r$lfc.4p, col="blue")
  abline(h=r$lfc.n, col="blue")
  
  
  abline(h= log2(((1 - p.tpc.estimate) * 2 + p.tpc.estimate * 4) / 2)  , col="gray80")
  abline(h= log2(((1 - p.tpc.estimate) * 2 + p.tpc.estimate * 3) / 2)  , col="gray80")
  abline(h= log2(((1 - p.tpc.estimate) * 2 + p.tpc.estimate * 1) / 2)  , col="gray80")
  

  
  abline(h= log2(((1 - p.tpc.hoogstrate.2021) * 2 + p.tpc.hoogstrate.2021 * 4) / 2)  , col=alpha("purple",0.3))
  abline(h= log2(((1 - p.tpc.hoogstrate.2021) * 2 + p.tpc.hoogstrate.2021 * 3) / 2)  , col=alpha("purple",0.3))
  abline(h= log2(((1 - p.tpc.hoogstrate.2021) * 2 + p.tpc.hoogstrate.2021 * 1) / 2)  ,  col=alpha("purple",0.3))
  
  
  text(0, 2, paste0("Estimated tumor-% [",bc,"]: ",r$pct , "%    (synapse: ",round(p.tpc.estimate * 100),",     hoogstrate RF 2021: ",round(p.tpc.hoogstrate.2021 * 100),")" )  ,pos=4)
  
  
  
  dev.off()
}


# write.table(tpc.estimate |> 
#               dplyr::select(portion_barcode, aliquot_barcode.wgs, pct, lfc.3p, lfc.4p, lfc.n, dist) |> 
#               dplyr::rename(tumor.purity.cnv.pct.2022 = pct)|> 
#               dplyr::rename(tumor.purity.cnv.lfc.3p.2022 = lfc.3p)|> 
#               dplyr::rename(tumor.purity.cnv.lfc.4p.2022 = lfc.4p)|> 
#               dplyr::rename(tumor.purity.cnv.lfc.n.2022 = lfc.n)|> 
#               dplyr::rename(tumor.purity.cnv.dist.2022 = dist)
#             , "output/tables/cnv/tumor.percentage.estimate_glass.2022.all_samples.txt")



plt <- tpc.estimate |> 
  dplyr::mutate(chr17.err = sample %in% chr17.err) |> 
  dplyr::left_join(
    glass.gbm.rnaseq.metadata.all.samples |> 
      dplyr::select(portion_barcode, aliquot_batch_synapse, predicted.GLASS.batch)
    ,by=c('portion_barcode'='portion_barcode')
  )


cor(plt$pct, plt$estimate.purity)

ggplot(plt, aes(x=pct, y=estimate.purity, col=aliquot_batch_synapse)) + # predicted.GLASS.batch, aliquot_batch_synapse
  geom_point() +
  theme_bw() +
  labs(x= "Tumor Purity (CNV detected)",y="Tumor Purity (GLASS/Estimate RNA call)")


ggplot(plt, aes(x=pct, y=hoostrate.rf.purity.2021, col=as.character(aliquot_batch_synapse))) + # predicted.GLASS.batch, aliquot_batch_synapse
  geom_point() +
  theme_bw() +
  labs(x= "Tumor Purity (CNV detected)",y="Tumor Purity (GLASS/Estimate RNA call)")


ggplot(plt, aes(x=pct, y=estimate.purity, col=chr17.err)) +
  geom_point() +
  facet_grid(cols = vars(chr17.err)) +
  theme_bw() +
  labs(x= "Tumor Purity (CNV detected)",y="Tumor Purity (GLASS/Estimate RNA call)")

ggplot(plt, aes(x=pct, y=estimate.purity, col=chr17.err)) +
  geom_point() +
  facet_grid(cols = vars(predicted.GLASS.batch)) +
  theme_bw() +
  labs(x= "Tumor Purity (CNV detected)",y="Tumor Purity (GLASS/Estimate RNA call)")


ggplot(plt, aes(x=pct, y=estimate.purity, col=chr17.err)) +
  geom_point() +
  facet_grid(cols = vars(aliquot_batch_synapse)) +
  theme_bw() +
  labs(x= "Tumor Purity (CNV detected)",y="Tumor Purity (GLASS/Estimate RNA call)")


ggplot(plt, aes(x=pct, y=hoostrate.rf.purity.2021, col=chr17.err)) +
  geom_point() +
  facet_grid(cols = vars(aliquot_batch_synapse)) +
  theme_bw() +
  labs(x= "Tumor Purity (CNV detected)",y="Tumor Purity (RF/hoogstrate/2021 RNA)")



ggplot(plt, aes(x=hoostrate.rf.purity.2021, y=estimate.purity, col=chr17.err)) +
  geom_point() +
  facet_grid(cols = vars(aliquot_batch_synapse)) +
  theme_bw() +
  labs(x= "Tumor Purity (RF/hoogstrate/2021 RNA)",y="Tumor Purity (GLASS/Estimate RNA call)")




a = read.csv('data/gsam/data/GLASS_GBM_R1-R2/variants_passgeno_20220531.csv')

aa <- a |> 
  dplyr::filter(grepl("CU-P102-R1",aliquot_barcode, fixed=T))


plot(sort(aa$af)) # VAFs also seem low, in concordance with DNA


ab <- a |> 
  dplyr::filter(grepl("CU-P102-TP",aliquot_barcode, fixed=T))


plot(sort(ab$af))



aa <- a |> 
  dplyr::filter(grepl("CU-R003-R1",aliquot_barcode, fixed=T))
plot(sort(aa$af)) # VAFs also seem low, in concordance with DNA





