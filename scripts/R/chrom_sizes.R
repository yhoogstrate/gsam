#!/usr/bin/env R

chrs_hg19 <- {{}}
chrs_hg19['chr1'] <- 249250621
chrs_hg19['chr2'] <- 243199373
chrs_hg19['chr3'] <- 198022430
chrs_hg19['chr4'] <- 191154276
chrs_hg19['chr5'] <- 180915260
chrs_hg19['chr6'] <- 171115067
chrs_hg19['chr7'] <- 159138663
chrs_hg19['chr8'] <- 146364022
chrs_hg19['chr9'] <- 141213431
chrs_hg19['chr10'] <- 135534747
chrs_hg19['chr11'] <- 135006516
chrs_hg19['chr12'] <- 133851895
chrs_hg19['chr13'] <- 115169878
chrs_hg19['chr14'] <- 107349540
chrs_hg19['chr15'] <- 102531392
chrs_hg19['chr16'] <- 90354753
chrs_hg19['chr17'] <- 81195210
chrs_hg19['chr18'] <- 78077248
chrs_hg19['chr19'] <- 59128983
chrs_hg19['chr20'] <- 63025520
chrs_hg19['chr21'] <- 48129895
chrs_hg19['chr22'] <- 51304566
chrs_hg19['chrX'] <- 155270560
chrs_hg19['chrY'] <- 59373566
n_hg19 = sum(chrs_hg19)


# cumulative_start:
chrs_hg19_s <- {{}}
chrs_hg19_s['chr1'] <- 0
chrs_hg19_s['chr2'] <- chrs_hg19['chr1'] + chrs_hg19_s['chr1']
chrs_hg19_s['chr3'] <- chrs_hg19['chr2'] + chrs_hg19_s['chr2']
chrs_hg19_s['chr4'] <- chrs_hg19['chr3'] + chrs_hg19_s['chr3']
chrs_hg19_s['chr5'] <- chrs_hg19['chr4'] + chrs_hg19_s['chr4']
chrs_hg19_s['chr6'] <- chrs_hg19['chr5'] + chrs_hg19_s['chr5']
chrs_hg19_s['chr7'] <- chrs_hg19['chr6'] + chrs_hg19_s['chr6']
chrs_hg19_s['chr8'] <- chrs_hg19['chr7'] + chrs_hg19_s['chr7']
chrs_hg19_s['chr9'] <- chrs_hg19['chr8'] + chrs_hg19_s['chr8']
chrs_hg19_s['chr10'] <- chrs_hg19['chr9'] + chrs_hg19_s['chr9']
chrs_hg19_s['chr11'] <- chrs_hg19['chr10'] + chrs_hg19_s['chr10']
chrs_hg19_s['chr12'] <- chrs_hg19['chr11'] + chrs_hg19_s['chr11']
chrs_hg19_s['chr13'] <- chrs_hg19['chr12'] + chrs_hg19_s['chr12']
chrs_hg19_s['chr14'] <- chrs_hg19['chr13'] + chrs_hg19_s['chr13']
chrs_hg19_s['chr15'] <- chrs_hg19['chr14'] + chrs_hg19_s['chr14']
chrs_hg19_s['chr16'] <- chrs_hg19['chr15'] + chrs_hg19_s['chr15']
chrs_hg19_s['chr17'] <- chrs_hg19['chr16'] + chrs_hg19_s['chr16']
chrs_hg19_s['chr18'] <- chrs_hg19['chr17'] + chrs_hg19_s['chr17']
chrs_hg19_s['chr19'] <- chrs_hg19['chr18'] + chrs_hg19_s['chr18']
chrs_hg19_s['chr20'] <- chrs_hg19['chr19'] + chrs_hg19_s['chr19']
chrs_hg19_s['chr21'] <- chrs_hg19['chr20'] + chrs_hg19_s['chr20']
chrs_hg19_s['chr22'] <- chrs_hg19['chr21'] + chrs_hg19_s['chr21']
chrs_hg19_s['chrX'] <- chrs_hg19['chr22'] + chrs_hg19_s['chr22']
chrs_hg19_s['chrY'] <- chrs_hg19['chrX'] + chrs_hg19_s['chrX']
#chrs_hg19_s['chrM'] <- chrs_hg19['chr'] + chrs_hg19_s['chr']


# cumulative:
chrs_hg19_e <- {{}}
chrs_hg19_e['chr1'] <- chrs_hg19['chr1']
chrs_hg19_e['chr2'] <- chrs_hg19['chr2'] + chrs_hg19_e['chr1']
chrs_hg19_e['chr3'] <- chrs_hg19['chr3'] + chrs_hg19_e['chr2']
chrs_hg19_e['chr4'] <- chrs_hg19['chr4'] + chrs_hg19_e['chr3']
chrs_hg19_e['chr5'] <- chrs_hg19['chr5'] + chrs_hg19_e['chr4']
chrs_hg19_e['chr6'] <- chrs_hg19['chr6'] + chrs_hg19_e['chr5']
chrs_hg19_e['chr7'] <- chrs_hg19['chr7'] + chrs_hg19_e['chr6']
chrs_hg19_e['chr8'] <- chrs_hg19['chr8'] + chrs_hg19_e['chr7']
chrs_hg19_e['chr9'] <- chrs_hg19['chr9'] + chrs_hg19_e['chr8']
chrs_hg19_e['chr10'] <- chrs_hg19['chr10'] + chrs_hg19_e['chr9']
chrs_hg19_e['chr11'] <- chrs_hg19['chr11'] + chrs_hg19_e['chr10']
chrs_hg19_e['chr12'] <- chrs_hg19['chr12'] + chrs_hg19_e['chr11']
chrs_hg19_e['chr13'] <- chrs_hg19['chr13'] + chrs_hg19_e['chr12']
chrs_hg19_e['chr14'] <- chrs_hg19['chr14'] + chrs_hg19_e['chr13']
chrs_hg19_e['chr15'] <- chrs_hg19['chr15'] + chrs_hg19_e['chr14']
chrs_hg19_e['chr16'] <- chrs_hg19['chr16'] + chrs_hg19_e['chr15']
chrs_hg19_e['chr17'] <- chrs_hg19['chr17'] + chrs_hg19_e['chr16']
chrs_hg19_e['chr18'] <- chrs_hg19['chr18'] + chrs_hg19_e['chr17']
chrs_hg19_e['chr19'] <- chrs_hg19['chr19'] + chrs_hg19_e['chr18']
chrs_hg19_e['chr20'] <- chrs_hg19['chr20'] + chrs_hg19_e['chr19']
chrs_hg19_e['chr21'] <- chrs_hg19['chr21'] + chrs_hg19_e['chr20']
chrs_hg19_e['chr22'] <- chrs_hg19['chr22'] + chrs_hg19_e['chr21']
chrs_hg19_e['chrX'] <- chrs_hg19['chrY'] + chrs_hg19_e['chr22']
chrs_hg19_e['chrY'] <- chrs_hg19['chrX'] + chrs_hg19_e['chrX']
#chrs_hg19_e['chrM'] <- 16569 + chrs_hg19_e['chr']



#chrs_hg38 <- {{}}
#chrs_hg38['chr1'] <- 248956422
#chrs_hg38['chr2'] <- 242193529
#chrs_hg38['chr3'] <- 198295559
#chrs_hg38['chr4'] <- 190214555
#chrs_hg38['chr5'] <- 181538259
#chrs_hg38['chr6'] <- 170805979
#chrs_hg38['chr7'] <- 159345973
#chrs_hg38['chr8'] <- 145138636
#chrs_hg38['chr9'] <- 138394717
#chrs_hg38['chr10'] <- 133797422
#chrs_hg38['chr11'] <- 135086622
#chrs_hg38['chr12'] <- 133275309
#chrs_hg38['chr13'] <- 114364328
#chrs_hg38['chr14'] <- 107043718
#chrs_hg38['chr15'] <- 101991189
#chrs_hg38['chr16'] <- 90338345
#chrs_hg38['chr17'] <- 83257441
#chrs_hg38['chr18'] <- 80373285
#chrs_hg38['chr19'] <- 58617616
#chrs_hg38['chr20'] <- 64444167
#chrs_hg38['chr21'] <- 46709983
#chrs_hg38['chr22'] <- 50818468
#chrs_hg38['chrX'] <- 156040895
#chrs_hg38['chrY'] <- 57227415
##chrs_hg38['chrM'] <- 16569
#n_hg38 = sum(chrs_hg38)

## cumulative:
#chrs_hg38_c <- {{}}
#chrs_hg38_c['chr1'] <- 248956422
#chrs_hg38_c['chr2'] <- 242193529 + chrs_hg38_c['chr1']
#chrs_hg38_c['chr3'] <- 198295559 + chrs_hg38_c['chr2']
#chrs_hg38_c['chr4'] <- 190214555 + chrs_hg38_c['chr3']
#chrs_hg38_c['chr5'] <- 181538259 + chrs_hg38_c['chr4']
#chrs_hg38_c['chr6'] <- 170805979 + chrs_hg38_c['chr5']
#chrs_hg38_c['chr7'] <- 159345973 + chrs_hg38_c['chr6']
#chrs_hg38_c['chr8'] <- 145138636 + chrs_hg38_c['chr7']
#chrs_hg38_c['chr9'] <- 138394717 + chrs_hg38_c['chr8']
#chrs_hg38_c['chr10'] <- 133797422 + chrs_hg38_c['chr9']
#chrs_hg38_c['chr11'] <- 135086622 + chrs_hg38_c['chr10']
#chrs_hg38_c['chr12'] <- 133275309 + chrs_hg38_c['chr11']
#chrs_hg38_c['chr13'] <- 114364328 + chrs_hg38_c['chr12']
#chrs_hg38_c['chr14'] <- 107043718 + chrs_hg38_c['chr13']
#chrs_hg38_c['chr15'] <- 101991189 + chrs_hg38_c['chr14']
#chrs_hg38_c['chr16'] <- 90338345 + chrs_hg38_c['chr15']
#chrs_hg38_c['chr17'] <- 83257441 + chrs_hg38_c['chr16']
#chrs_hg38_c['chr18'] <- 80373285 + chrs_hg38_c['chr17']
#chrs_hg38_c['chr19'] <- 58617616 + chrs_hg38_c['chr18']
#chrs_hg38_c['chr20'] <- 64444167 + chrs_hg38_c['chr19']
#chrs_hg38_c['chr21'] <- 46709983 + chrs_hg38_c['chr20']
#chrs_hg38_c['chr22'] <- 50818468 + chrs_hg38_c['chr21']
#chrs_hg38_c['chrX'] <- 156040895 + chrs_hg38_c['chr22']
#chrs_hg38_c['chrY'] <- 57227415 + chrs_hg38_c['chrX']
##chrs_hg38_c['chrM'] <- 16569 + chrs_hg38_c['chr']

## cumulative_start:
#chrs_hg38_s <- {{}}
#chrs_hg38_s['chr1'] <- 0
#chrs_hg38_s['chr2'] <- chrs_hg38['chr1'] + chrs_hg38_s['chr1']
#chrs_hg38_s['chr3'] <- chrs_hg38['chr2'] + chrs_hg38_s['chr2']
#chrs_hg38_s['chr4'] <- chrs_hg38['chr3'] + chrs_hg38_s['chr3']
#chrs_hg38_s['chr5'] <- chrs_hg38['chr4'] + chrs_hg38_s['chr4']
#chrs_hg38_s['chr6'] <- chrs_hg38['chr5'] + chrs_hg38_s['chr5']
#chrs_hg38_s['chr7'] <- chrs_hg38['chr6'] + chrs_hg38_s['chr6']
#chrs_hg38_s['chr8'] <- chrs_hg38['chr7'] + chrs_hg38_s['chr7']
#chrs_hg38_s['chr9'] <- chrs_hg38['chr8'] + chrs_hg38_s['chr8']
#chrs_hg38_s['chr10'] <- chrs_hg38['chr9'] + chrs_hg38_s['chr9']
#chrs_hg38_s['chr11'] <- chrs_hg38['chr10'] + chrs_hg38_s['chr10']
#chrs_hg38_s['chr12'] <- chrs_hg38['chr11'] + chrs_hg38_s['chr11']
#chrs_hg38_s['chr13'] <- chrs_hg38['chr12'] + chrs_hg38_s['chr12']
#chrs_hg38_s['chr14'] <- chrs_hg38['chr13'] + chrs_hg38_s['chr13']
#chrs_hg38_s['chr15'] <- chrs_hg38['chr14'] + chrs_hg38_s['chr14']
#chrs_hg38_s['chr16'] <- chrs_hg38['chr15'] + chrs_hg38_s['chr15']
#chrs_hg38_s['chr17'] <- chrs_hg38['chr16'] + chrs_hg38_s['chr16']
#chrs_hg38_s['chr18'] <- chrs_hg38['chr17'] + chrs_hg38_s['chr17']
#chrs_hg38_s['chr19'] <- chrs_hg38['chr18'] + chrs_hg38_s['chr18']
#chrs_hg38_s['chr20'] <- chrs_hg38['chr19'] + chrs_hg38_s['chr19']
#chrs_hg38_s['chr21'] <- chrs_hg38['chr20'] + chrs_hg38_s['chr20']
#chrs_hg38_s['chr22'] <- chrs_hg38['chr21'] + chrs_hg38_s['chr21']
#chrs_hg38_s['chrX'] <- chrs_hg38['chr22'] + chrs_hg38_s['chr22']
#chrs_hg38_s['chrY'] <- chrs_hg38['chrX'] + chrs_hg38_s['chrX']
##chrs_hg38_s['chrM'] <- chrs_hg38['chr'] + chrs_hg38_s['chr']
