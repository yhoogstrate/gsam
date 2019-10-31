

#md5sum /data/users/youri/neurogen-tmp/star_2_6__hg19/gencode.v29lift37.annotation.gff3
#md5sum /data/users/youri/neurogen-tmp/align-gsam/star-hg19/gencode.v29.annotation.gff3




featureCounts \
    -a star-hg19/gencode.v31lift37.annotation.gtf \
    -o "featureCounts_gsam_1st96.exon-level.txt" \
    --tmpDir "/dev/shm/" \
    -p \
    -s 2 \
    -T 64 \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAB1_1917_Plate1_E4_S35_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAC1_1920_Plate1_A4_S31_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAF2_1927_Plate1_D9_S21_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAG2_1930_Plate1_B10_S27_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAM1_1934_Plate1_B11_S35_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAM2_1942_Plate1_H11_S41_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAN2_1944_Plate1_F10_S31_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAP1_2658_Plate1_C7_S105_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAS1_1947_Plate1_E1_S11_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAS2_1949_Plate1_D6_S98_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAT1_1939_Plate1_F4_S36_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAT2_1950_Plate1_C9_S20_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAU1_1945_Plate1_F5_S92_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAV1_1946_Plate1_E6_S99_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAW1_1953_Plate1_F12_S47_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAW2_1956_Plate1_F11_S39_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAX1_1858_Plate1_A5_S87_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAX2_1957_Plate1_E10_S30_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AAY2_1959_Plate1_B4_S32_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_ACA2_1962_Plate1_G11_S40_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_ADA2_1964_Plate1_G10_S32_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AFA1_2545_Plate1_A12_S42_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AIA1_1904_Plate1_B9_S19_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AIA2_1974_Plate1_D1_S10_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_ALA2_1979_Plate1_H5_S94_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AMA2_1976_Plate1_A6_S95_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AOA2_1978_Plate1_A10_S26_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AQA2_1980_Plate1_G5_S93_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZA2_2669_Plate1_C1_S9_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZB1_2664_Plate1_E2_S19_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZB2_2670_Plate1_D2_S18_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZC2_2671_Plate1_E3_S27_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZD2_2673_Plate1_G3_S29_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZE1_2660_Plate1_B1_S8_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZE2_2672_Plate1_A1_S7_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZF1_2662_Plate1_H2_S22_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZG1_2661_Plate1_D3_S26_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZG2_2675_Plate1_C3_S25_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZH1_2666_Plate1_G2_S21_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_AZH2_2676_Plate1_F2_S20_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAC1_2097_Plate1_G9_S24_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAD1_2099_Plate1_A11_S34_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAD2_2100_Plate1_D7_S106_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAE1_2101_Plate1_G7_S109_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAE2_2102_Plate1_F7_S108_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAK2_2025_Plate1_E12_S46_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAO1_2190_Plate1_C4_S33_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAR1_2036_Plate1_A8_S111_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAT1_2038_Plate1_B6_S96_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAU2_2562_Plate1_E7_S107_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_BAW1_2119_Plate1_H9_S25_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CAF1_2755_Plate1_C6_S97_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CAF2_1907_Plate1_E5_S91_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CAO2_1912_Plate1_A9_S18_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CAV1_18641_Plate1_C11_S36_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CAV2_1872_Plate1_H10_S33_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBA1_1870_Plate1_D10_S29_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBE2_1873_Plate1_C10_S28_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBG1_1838-1_Plate1_H6_S102_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBR1_1853-1_Plate1_B5_S88_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBR2_1856_Plate1_H1_S14_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBS1_2739_Plate1_D4_S34_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBT1_2565_Plate1_A7_S103_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CBT2_2685_Plate1_F6_S100_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CCD1_1911_Plate1_D5_S90_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CCD2_1913_Plate1_F1_S12_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CCW1_1901_Plate1_H7_S110_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_CDA1_1892_Plate1_G1_S13_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EAG1_2197_Plate1_F8_S116_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EAJ1_2235_Plate1_E9_S22_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EAP2_2391_Plate1_D8_S114_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EAZ1_2265_Plate1_F9_S23_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EBF2_2278_Plate1_E11_S38_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EBH1_2281_Plate1_H12_S49_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EBN1_2399_Plate1_C12_S44_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EBP1_2295_Plate1_G12_S48_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EBU1_2627_Plate1_G8_S117_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_EBX2_2632_Plate1_H8_S118_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_ECN1_2635_Plate1_D11_S37_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_FAQ1_2381_Plate1_C8_S113_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_GAE2_2426_Plate1_B8_S112_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_HAA1_2519_Plate1_E8_S115_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAB1_2701_Plate1_B7_S104_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAB2_2702_Plate1_C2_S17_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAC2_2704_Plate1_G6_S101_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAD1_2705_Plate1_B2_S16_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAE1_2707_Plate1_A2_S15_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAE2_2709_Plate1_H3_S30_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAF1_2710_Plate1_B3_S24_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAG1_2712_Plate1_B12_S43_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAG2_2713_Plate1_C5_S89_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAI1_2716_Plate1_G4_S37_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAJ1_2718_Plate1_D12_S45_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAJ2_2720_Plate1_F3_S28_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAK1_2721_Plate1_A3_S23_001_hg19.sambamba-dedup.bam \
    /data/users/youri/mnt/neurogen-ro/gsam/RNA/alignments/Pfrench_JAK2_2722_Plate1_H4_S38_001_hg19.sambamba-dedup.bam
