#!/bin/bash

mkdir -p "processed/AAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAB1-1917_NGS19-J295_AHKK5WDSXX_S35_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAB1-1917_NGS19-J295_BHJTTYDSXX_S29_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAB1-1917_NGS19-J295_AHKK5WDSXX_S35_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAB1-1917_NGS19-J295_BHJTTYDSXX_S29_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAB2-2479_NGS19-J468_AHM2NJDSXX_S16_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAB2-2479_NGS19-J468_AHM2NJDSXX_S16_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAC1-1920_NGS19-J291_AHKK5WDSXX_S31_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAC1-1920_NGS19-J291_BHJTTYDSXX_S25_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAC1-1920_NGS19-J291_AHKK5WDSXX_S31_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAC1-1920_NGS19-J291_BHJTTYDSXX_S25_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAC2-1922_NGS19-J487_AHM2NJDSXX_S35_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAC2-1922_NGS19-J487_AHM2NJDSXX_S35_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAD1-2358_NGS19-J482_AHM2NJDSXX_S30_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAD1-2358_NGS19-J482_AHM2NJDSXX_S30_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAD2-1924_NGS19-J418_AHLGTCDSXX_S62_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAD2-1924_NGS19-J418_AHLGTCDSXX_S62_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAF1-1923_NGS19-J420_AHLGTCDSXX_S64_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAF1-1923_NGS19-J420_AHLGTCDSXX_S64_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAF2-1927_NGS19-J334_AHKK5WDSXX_S122_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAF2-1927_NGS19-J334_AHKK5WDSXX_S122_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAG1-1928_NGS19-J371_AHLGTCDSXX_S15_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAG1-1928_NGS19-J371_AHLGTCDSXX_S15_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAG2-1930_NGS19-J340_AHKK5WDSXX_S128_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAG2-1930_NGS19-J340_AHKK5WDSXX_S128_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAJ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAJ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAJ1-2603_NGS19-J391_AHLGTCDSXX_S35_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAJ1-2603_NGS19-J391_AHLGTCDSXX_S35_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAJ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAJ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAJ2-2604_NGS19-J504_AHM2NJDSXX_S231_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAJ2-2604_NGS19-J504_AHM2NJDSXX_S231_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAL1-2605_NGS19-J547_AHM2NJDSXX_S274_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAL1-2605_NGS19-J547_AHM2NJDSXX_S274_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAL2-2606_NGS19-J399_AHLGTCDSXX_S43_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAL2-2606_NGS19-J399_AHLGTCDSXX_S43_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAM1-1934_NGS19-J348_AHKK5WDSXX_S136_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAM1-1934_NGS19-J348_AHKK5WDSXX_S136_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAM2-1942_NGS19-J354_AHKK5WDSXX_S142_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAM2-1942_NGS19-J354_AHKK5WDSXX_S142_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAN1-2360_NGS19-J381_AHLGTCDSXX_S25_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAN1-2360_NGS19-J381_AHLGTCDSXX_S25_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAN2-1944_NGS19-J344_AHKK5WDSXX_S132_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAN2-1944_NGS19-J344_AHKK5WDSXX_S132_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAP1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAP1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAP1-2658_NGS19-J317_AHKK5WDSXX_S105_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAP1-2658_NGS19-J317_BHJTTYDSXX_S51_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAP1-2658_NGS19-J317_AHKK5WDSXX_S105_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAP1-2658_NGS19-J317_BHJTTYDSXX_S51_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAP2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAP2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAP2-2659_NGS19-J496_AHM2NJDSXX_S223_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAP2-2659_NGS19-J496_AHM2NJDSXX_S223_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAS1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAS1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAS1-1947_NGS19-J271_AHKK5WDSXX_S11_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAS1-1947_NGS19-J271_BHJTTYDSXX_S5_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAS1-1947_NGS19-J271_AHKK5WDSXX_S11_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAS1-1947_NGS19-J271_BHJTTYDSXX_S5_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAS2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAS2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAS2-1949_NGS19-J310_AHKK5WDSXX_S98_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAS2-1949_NGS19-J310_BHJTTYDSXX_S44_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAS2-1949_NGS19-J310_AHKK5WDSXX_S98_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAS2-1949_NGS19-J310_BHJTTYDSXX_S44_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAT1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAT1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAT1-1939_NGS19-J296_AHKK5WDSXX_S36_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAT1-1939_NGS19-J296_BHJTTYDSXX_S30_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAT1-1939_NGS19-J296_AHKK5WDSXX_S36_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAT1-1939_NGS19-J296_BHJTTYDSXX_S30_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAT2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAT2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAT2-1950_NGS19-J333_AHKK5WDSXX_S121_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAT2-1950_NGS19-J333_AHKK5WDSXX_S121_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAU1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAU1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAU1-1945_NGS19-J304_AHKK5WDSXX_S92_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAU1-1945_NGS19-J304_BHJTTYDSXX_S38_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAU1-1945_NGS19-J304_AHKK5WDSXX_S92_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAU1-1945_NGS19-J304_BHJTTYDSXX_S38_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAU2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAU2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAU2-1951_NGS19-J599_BHM33WDSXX_S51_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAU2-1951_NGS19-J599_BHM33WDSXX_S51_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAV1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAV1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAV1-1946_NGS19-J311_AHKK5WDSXX_S99_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAV1-1946_NGS19-J311_BHJTTYDSXX_S45_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAV1-1946_NGS19-J311_AHKK5WDSXX_S99_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAV1-1946_NGS19-J311_BHJTTYDSXX_S45_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAV2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAV2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAV2-1955_NGS19-J450_AHLGTCDSXX_S94_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAV2-1955_NGS19-J450_AHLGTCDSXX_S94_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAW1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAW1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAW1-1953_NGS19-J360_AHKK5WDSXX_S148_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAW1-1953_NGS19-J360_AHKK5WDSXX_S148_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAW2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAW2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAW2-1956_NGS19-J352_AHKK5WDSXX_S140_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAW2-1956_NGS19-J352_AHKK5WDSXX_S140_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAX1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAX1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAX1-1858_NGS19-J299_AHKK5WDSXX_S87_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAX1-1858_NGS19-J299_BHJTTYDSXX_S33_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAX1-1858_NGS19-J299_AHKK5WDSXX_S87_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAX1-1858_NGS19-J299_BHJTTYDSXX_S33_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAX2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAX2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAX2-1957_NGS19-J343_AHKK5WDSXX_S131_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAX2-1957_NGS19-J343_AHKK5WDSXX_S131_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AAY1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAY1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAY1-2543_NGS19-J390_AHLGTCDSXX_S34_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAY1-2543_NGS19-J390_AHLGTCDSXX_S34_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AAY2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AAY2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AAY2-1959_NGS19-J292_AHKK5WDSXX_S32_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AAY2-1959_NGS19-J292_BHJTTYDSXX_S26_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AAY2-1959_NGS19-J292_AHKK5WDSXX_S32_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AAY2-1959_NGS19-J292_BHJTTYDSXX_S26_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ABA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ABA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ABA1-2544_NGS19-J611_BHM33WDSXX_S63_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ABA1-2544_NGS19-J611_BHM33WDSXX_S63_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ABA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ABA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ABA2-1961_NGS19-J428_AHLGTCDSXX_S72_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ABA2-1961_NGS19-J428_AHLGTCDSXX_S72_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ACA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ACA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ACA1-1958_NGS19-J519_AHM2NJDSXX_S246_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ACA1-1958_NGS19-J519_AHM2NJDSXX_S246_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ACA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ACA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ACA2-1962_NGS19-J353_AHKK5WDSXX_S141_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ACA2-1962_NGS19-J353_AHKK5WDSXX_S141_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ADA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ADA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ADA1-1884_NGS19-J436_AHLGTCDSXX_S80_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ADA1-1884_NGS19-J436_AHLGTCDSXX_S80_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ADA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ADA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ADA2-1964_NGS19-J345_AHKK5WDSXX_S133_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ADA2-1964_NGS19-J345_AHKK5WDSXX_S133_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AFA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AFA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AFA1-2545_NGS19-J355_AHKK5WDSXX_S143_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AFA1-2545_NGS19-J355_AHKK5WDSXX_S143_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AFA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AFA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AFA2-1966_NGS19-J506_AHM2NJDSXX_S233_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AFA2-1966_NGS19-J506_AHM2NJDSXX_S233_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AHA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AHA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AHA1-1967_NGS19-J598_BHM33WDSXX_S50_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AHA1-1967_NGS19-J598_BHM33WDSXX_S50_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AHA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AHA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AHA2-2548_NGS19-J383_AHLGTCDSXX_S27_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AHA2-2548_NGS19-J383_AHLGTCDSXX_S27_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AIA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AIA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AIA1-1904_NGS19-J332_AHKK5WDSXX_S120_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AIA1-1904_NGS19-J332_AHKK5WDSXX_S120_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AIA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AIA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AIA2-1974_NGS19-J270_AHKK5WDSXX_S10_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AIA2-1974_NGS19-J270_BHJTTYDSXX_S4_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AIA2-1974_NGS19-J270_AHKK5WDSXX_S10_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AIA2-1974_NGS19-J270_BHJTTYDSXX_S4_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AKA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AKA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AKA1-1970_NGS19-J512_AHM2NJDSXX_S239_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AKA1-1970_NGS19-J512_AHM2NJDSXX_S239_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AKA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AKA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AKA2-1975_NGS19-J367_AHLGTCDSXX_S11_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AKA2-1975_NGS19-J367_AHLGTCDSXX_S11_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ALA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ALA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ALA1-1972_NGS19-J485_AHM2NJDSXX_S33_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ALA1-1972_NGS19-J485_AHM2NJDSXX_S33_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ALA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ALA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ALA2-1979_NGS19-J306_AHKK5WDSXX_S94_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/ALA2-1979_NGS19-J306_BHJTTYDSXX_S40_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ALA2-1979_NGS19-J306_AHKK5WDSXX_S94_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/ALA2-1979_NGS19-J306_BHJTTYDSXX_S40_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AMA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AMA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AMA1-2417_NGS19-J444_AHLGTCDSXX_S88_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AMA1-2417_NGS19-J444_AHLGTCDSXX_S88_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AMA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AMA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AMA2-1976_NGS19-J307_AHKK5WDSXX_S95_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AMA2-1976_NGS19-J307_BHJTTYDSXX_S41_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AMA2-1976_NGS19-J307_AHKK5WDSXX_S95_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AMA2-1976_NGS19-J307_BHJTTYDSXX_S41_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AOA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AOA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AOA1-1968_NGS19-J600_BHM33WDSXX_S52_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AOA1-1968_NGS19-J600_BHM33WDSXX_S52_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AOA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AOA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AOA2-1978_NGS19-J339_AHKK5WDSXX_S127_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AOA2-1978_NGS19-J339_AHKK5WDSXX_S127_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AQA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AQA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AQA1-1971_NGS19-J556_BHM33WDSXX_S8_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AQA1-1971_NGS19-J556_BHM33WDSXX_S8_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AQA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AQA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AQA2-1980_NGS19-J305_AHKK5WDSXX_S93_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AQA2-1980_NGS19-J305_BHJTTYDSXX_S39_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AQA2-1980_NGS19-J305_AHKK5WDSXX_S93_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AQA2-1980_NGS19-J305_BHJTTYDSXX_S39_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZA1-2663_NGS19-J376_AHLGTCDSXX_S20_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZA1-2663_NGS19-J376_AHLGTCDSXX_S20_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZA2-2669_NGS19-J269_AHKK5WDSXX_S9_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZA2-2669_NGS19-J269_BHJTTYDSXX_S3_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZA2-2669_NGS19-J269_AHKK5WDSXX_S9_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZA2-2669_NGS19-J269_BHJTTYDSXX_S3_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZB1-2664_NGS19-J279_AHKK5WDSXX_S19_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZB1-2664_NGS19-J279_BHJTTYDSXX_S13_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZB1-2664_NGS19-J279_AHKK5WDSXX_S19_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZB1-2664_NGS19-J279_BHJTTYDSXX_S13_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZB2-2670_NGS19-J278_AHKK5WDSXX_S18_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZB2-2670_NGS19-J278_BHJTTYDSXX_S12_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZB2-2670_NGS19-J278_AHKK5WDSXX_S18_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZB2-2670_NGS19-J278_BHJTTYDSXX_S12_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZC1-2665_NGS19-J396_AHLGTCDSXX_S40_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZC1-2665_NGS19-J396_AHLGTCDSXX_S40_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZC2-2671_NGS19-J287_AHKK5WDSXX_S27_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZC2-2671_NGS19-J287_BHJTTYDSXX_S21_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZC2-2671_NGS19-J287_AHKK5WDSXX_S27_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZC2-2671_NGS19-J287_BHJTTYDSXX_S21_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZD1-2668_NGS19-J595_BHM33WDSXX_S47_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZD1-2668_NGS19-J595_BHM33WDSXX_S47_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZD2-2673_NGS19-J289_AHKK5WDSXX_S29_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZD2-2673_NGS19-J289_BHJTTYDSXX_S23_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZD2-2673_NGS19-J289_AHKK5WDSXX_S29_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZD2-2673_NGS19-J289_BHJTTYDSXX_S23_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZE1-2660_NGS19-J268_AHKK5WDSXX_S8_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZE1-2660_NGS19-J268_BHJTTYDSXX_S2_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZE1-2660_NGS19-J268_AHKK5WDSXX_S8_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZE1-2660_NGS19-J268_BHJTTYDSXX_S2_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZE2-2672_NGS19-J267_AHKK5WDSXX_S7_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZE2-2672_NGS19-J267_BHJTTYDSXX_S1_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZE2-2672_NGS19-J267_AHKK5WDSXX_S7_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZE2-2672_NGS19-J267_BHJTTYDSXX_S1_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZF1-2662_NGS19-J282_AHKK5WDSXX_S22_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZF1-2662_NGS19-J282_BHJTTYDSXX_S16_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZF1-2662_NGS19-J282_AHKK5WDSXX_S22_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZF1-2662_NGS19-J282_BHJTTYDSXX_S16_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZF2-2674_NGS19-J364_AHLGTCDSXX_S8_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZF2-2674_NGS19-J364_AHLGTCDSXX_S8_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZG1-2661_NGS19-J286_AHKK5WDSXX_S26_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZG1-2661_NGS19-J286_BHJTTYDSXX_S20_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZG1-2661_NGS19-J286_AHKK5WDSXX_S26_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZG1-2661_NGS19-J286_BHJTTYDSXX_S20_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZG2-2675_NGS19-J285_AHKK5WDSXX_S25_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZG2-2675_NGS19-J285_BHJTTYDSXX_S19_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZG2-2675_NGS19-J285_AHKK5WDSXX_S25_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZG2-2675_NGS19-J285_BHJTTYDSXX_S19_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/AZH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZH1-2666_NGS19-J281_AHKK5WDSXX_S21_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZH1-2666_NGS19-J281_BHJTTYDSXX_S15_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZH1-2666_NGS19-J281_AHKK5WDSXX_S21_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZH1-2666_NGS19-J281_BHJTTYDSXX_S15_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/AZH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/AZH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/AZH2-2676_NGS19-J280_AHKK5WDSXX_S20_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/AZH2-2676_NGS19-J280_BHJTTYDSXX_S14_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/AZH2-2676_NGS19-J280_AHKK5WDSXX_S20_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/AZH2-2676_NGS19-J280_BHJTTYDSXX_S14_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAA1-2020_NGS19-J502_AHM2NJDSXX_S229_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAA1-2020_NGS19-J502_AHM2NJDSXX_S229_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAA2-2021_NGS19-J632_BHM33WDSXX_S84_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAA2-2021_NGS19-J632_BHM33WDSXX_S84_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAB1-2022_NGS19-J541_AHM2NJDSXX_S268_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAB1-2022_NGS19-J541_AHM2NJDSXX_S268_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAB2-2023_NGS19-J407_AHLGTCDSXX_S51_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAB2-2023_NGS19-J407_AHLGTCDSXX_S51_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAC1-2097_NGS19-J337_AHKK5WDSXX_S125_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAC1-2097_NGS19-J337_AHKK5WDSXX_S125_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAC2-2098_NGS19-J542_AHM2NJDSXX_S269_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAC2-2098_NGS19-J542_AHM2NJDSXX_S269_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAD1-2099_NGS19-J347_AHKK5WDSXX_S135_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAD1-2099_NGS19-J347_AHKK5WDSXX_S135_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAD2-2100_NGS19-J318_AHKK5WDSXX_S106_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAD2-2100_NGS19-J318_BHJTTYDSXX_S52_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAD2-2100_NGS19-J318_AHKK5WDSXX_S106_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAD2-2100_NGS19-J318_BHJTTYDSXX_S52_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAE1-2101_NGS19-J321_AHKK5WDSXX_S109_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAE1-2101_NGS19-J321_BHJTTYDSXX_S55_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAE1-2101_NGS19-J321_AHKK5WDSXX_S109_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAE1-2101_NGS19-J321_BHJTTYDSXX_S55_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAE2-2102_NGS19-J320_AHKK5WDSXX_S108_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAE2-2102_NGS19-J320_BHJTTYDSXX_S54_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAE2-2102_NGS19-J320_AHKK5WDSXX_S108_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAE2-2102_NGS19-J320_BHJTTYDSXX_S54_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAH1-2690_NGS19-J400_AHLGTCDSXX_S44_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAH1-2690_NGS19-J400_AHLGTCDSXX_S44_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAH2-2116_NGS19-J466_AHM2NJDSXX_S14_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAH2-2116_NGS19-J466_AHM2NJDSXX_S14_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAI1-2105_NGS19-J406_AHLGTCDSXX_S50_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAI1-2105_NGS19-J406_AHLGTCDSXX_S50_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAI2-2555_NGS19-J601_BHM33WDSXX_S53_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAI2-2555_NGS19-J601_BHM33WDSXX_S53_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAJ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAJ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAJ1-2107_NGS19-J497_AHM2NJDSXX_S224_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAJ1-2107_NGS19-J497_AHM2NJDSXX_S224_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAJ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAJ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAJ2-2108_NGS19-J539_AHM2NJDSXX_S266_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAJ2-2108_NGS19-J539_AHM2NJDSXX_S266_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAK1-2024_NGS19-J495_AHM2NJDSXX_S222_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAK1-2024_NGS19-J495_AHM2NJDSXX_S222_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAK2-2025_NGS19-J359_AHKK5WDSXX_S147_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAK2-2025_NGS19-J359_AHKK5WDSXX_S147_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAL1-2185_NGS19-J513_AHM2NJDSXX_S240_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAL1-2185_NGS19-J513_AHM2NJDSXX_S240_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAL2-2186_NGS19-J548_AHM2NJDSXX_S275_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAL2-2186_NGS19-J548_AHM2NJDSXX_S275_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAM1-2113_NGS19-J393_AHLGTCDSXX_S37_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAM1-2113_NGS19-J393_AHLGTCDSXX_S37_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAM2-2114_NGS19-J373_AHLGTCDSXX_S17_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAM2-2114_NGS19-J373_AHLGTCDSXX_S17_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAN1-2556_NGS19-J630_BHM33WDSXX_S82_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAN1-2556_NGS19-J630_BHM33WDSXX_S82_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAN2-2029_NGS19-J510_AHM2NJDSXX_S237_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAN2-2029_NGS19-J510_AHM2NJDSXX_S237_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAO1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAO1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAO1-2190_NGS19-J293_AHKK5WDSXX_S33_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAO1-2190_NGS19-J293_BHJTTYDSXX_S27_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAO1-2190_NGS19-J293_AHKK5WDSXX_S33_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAO1-2190_NGS19-J293_BHJTTYDSXX_S27_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAO2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAO2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAO2-2031_NGS19-J529_AHM2NJDSXX_S256_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAO2-2031_NGS19-J529_AHM2NJDSXX_S256_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAR1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAR1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAR1-2036_NGS19-J323_AHKK5WDSXX_S111_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAR1-2036_NGS19-J323_BHJTTYDSXX_S57_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAR1-2036_NGS19-J323_AHKK5WDSXX_S111_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAR1-2036_NGS19-J323_BHJTTYDSXX_S57_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAR2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAR2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAR2-2037_NGS19-J489_AHM2NJDSXX_S37_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAR2-2037_NGS19-J489_AHM2NJDSXX_S37_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAS1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAS1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAS1-2117_NGS19-J614_BHM33WDSXX_S66_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAS1-2117_NGS19-J614_BHM33WDSXX_S66_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAS2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAS2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAS2-2118_NGS19-J398_AHLGTCDSXX_S42_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAS2-2118_NGS19-J398_AHLGTCDSXX_S42_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAT1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAT1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAT1-2038_NGS19-J308_AHKK5WDSXX_S96_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAT1-2038_NGS19-J308_BHJTTYDSXX_S42_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAT1-2038_NGS19-J308_AHKK5WDSXX_S96_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAT1-2038_NGS19-J308_BHJTTYDSXX_S42_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAT2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAT2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAT2-2560_NGS19-J380_AHLGTCDSXX_S24_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAT2-2560_NGS19-J380_AHLGTCDSXX_S24_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAU1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAU1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAU1-2561_NGS19-J551_AHM2NJDSXX_S278_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAU1-2561_NGS19-J551_AHM2NJDSXX_S278_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAU2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAU2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAU2-2562_NGS19-J319_AHKK5WDSXX_S107_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/BAU2-2562_NGS19-J319_BHJTTYDSXX_S53_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAU2-2562_NGS19-J319_AHKK5WDSXX_S107_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/BAU2-2562_NGS19-J319_BHJTTYDSXX_S53_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAW1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAW1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAW1-2119_NGS19-J338_AHKK5WDSXX_S126_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAW1-2119_NGS19-J338_AHKK5WDSXX_S126_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAW2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAW2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAW2-2120_NGS19-J616_BHM33WDSXX_S68_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAW2-2120_NGS19-J616_BHM33WDSXX_S68_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAX1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAX1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAX1-2590_NGS19-J568_BHM33WDSXX_S20_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAX1-2590_NGS19-J568_BHM33WDSXX_S20_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAX2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAX2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAX2-2691_NGS19-J610_BHM33WDSXX_S62_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAX2-2691_NGS19-J610_BHM33WDSXX_S62_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/BAY1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAY1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAY1-2693_NGS19-J532_AHM2NJDSXX_S259_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAY1-2693_NGS19-J532_AHM2NJDSXX_S259_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/BAY2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/BAY2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/BAY2-2598_NGS19-J514_AHM2NJDSXX_S241_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/BAY2-2598_NGS19-J514_AHM2NJDSXX_S241_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAC1-2752_NGS19-J388_AHLGTCDSXX_S32_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAC1-2752_NGS19-J388_AHLGTCDSXX_S32_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAC2-2608_NGS19-J461_AHM2NJDSXX_S9_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAC2-2608_NGS19-J461_AHM2NJDSXX_S9_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAD1-2367_NGS19-J430_AHLGTCDSXX_S74_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAD1-2367_NGS19-J430_AHLGTCDSXX_S74_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAD2-2609_NGS19-J464_AHM2NJDSXX_S12_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAD2-2609_NGS19-J464_AHM2NJDSXX_S12_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CAF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAF1-2755_NGS19-J309_AHKK5WDSXX_S97_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CAF1-2755_NGS19-J309_BHJTTYDSXX_S43_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAF1-2755_NGS19-J309_AHKK5WDSXX_S97_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CAF1-2755_NGS19-J309_BHJTTYDSXX_S43_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CAF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAF2-1907_NGS19-J303_AHKK5WDSXX_S91_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CAF2-1907_NGS19-J303_BHJTTYDSXX_S37_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAF2-1907_NGS19-J303_AHKK5WDSXX_S91_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CAF2-1907_NGS19-J303_BHJTTYDSXX_S37_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CAO1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAO1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAO1-2480_NGS19-J462_AHM2NJDSXX_S10_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAO1-2480_NGS19-J462_AHM2NJDSXX_S10_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CAO1-replicate";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAO1-replicate" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAO1-replicate-2480_NGS19-J636_BHM33WDSXX_S88_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAO1-replicate-2480_NGS19-J636_BHM33WDSXX_S88_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CAO2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAO2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAO2-1912_NGS19-J331_AHKK5WDSXX_S119_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAO2-1912_NGS19-J331_AHKK5WDSXX_S119_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CAV1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAV1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAV1-1864_NGS19-J349_AHKK5WDSXX_S137_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAV1-1864_NGS19-J349_AHKK5WDSXX_S137_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CAV2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CAV2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CAV2-1872_NGS19-J346_AHKK5WDSXX_S134_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CAV2-1872_NGS19-J346_AHKK5WDSXX_S134_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBA1-1870_NGS19-J342_AHKK5WDSXX_S130_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBA1-1870_NGS19-J342_AHKK5WDSXX_S130_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBA2-1874_NGS19-J474_AHM2NJDSXX_S22_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBA2-1874_NGS19-J474_AHM2NJDSXX_S22_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBE1-1862_NGS19-J537_AHM2NJDSXX_S264_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBE1-1862_NGS19-J537_AHM2NJDSXX_S264_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBE2-1873_NGS19-J341_AHKK5WDSXX_S129_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBE2-1873_NGS19-J341_AHKK5WDSXX_S129_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBG1-1838_NGS19-J314_AHKK5WDSXX_S102_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CBG1-1838_NGS19-J314_BHJTTYDSXX_S48_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBG1-1838_NGS19-J314_AHKK5WDSXX_S102_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CBG1-1838_NGS19-J314_BHJTTYDSXX_S48_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBG2-1841_NGS19-J478_AHM2NJDSXX_S26_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBG2-1841_NGS19-J478_AHM2NJDSXX_S26_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBH1-1794_NGS19-J507_AHM2NJDSXX_S234_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBH1-1794_NGS19-J507_AHM2NJDSXX_S234_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBH2-1898_NGS19-J437_AHLGTCDSXX_S81_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBH2-1898_NGS19-J437_AHLGTCDSXX_S81_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBI1-2735_NGS19-J484_AHM2NJDSXX_S32_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBI1-2735_NGS19-J484_AHM2NJDSXX_S32_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBI2-1887_NGS19-J370_AHLGTCDSXX_S14_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBI2-1887_NGS19-J370_AHLGTCDSXX_S14_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBM1-1865_NGS19-J488_AHM2NJDSXX_S36_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBM1-1865_NGS19-J488_AHM2NJDSXX_S36_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBM2-1883_NGS19-J431_AHLGTCDSXX_S75_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBM2-1883_NGS19-J431_AHLGTCDSXX_S75_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBP1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBP1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBP1-2737_NGS19-J479_AHM2NJDSXX_S27_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBP1-2737_NGS19-J479_AHM2NJDSXX_S27_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBP2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBP2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBP2-2738_NGS19-J439_AHLGTCDSXX_S83_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBP2-2738_NGS19-J439_AHLGTCDSXX_S83_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBQ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBQ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBQ1-1847_NGS19-J382_AHLGTCDSXX_S26_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBQ1-1847_NGS19-J382_AHLGTCDSXX_S26_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBQ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBQ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBQ2-1863_NGS19-J394_AHLGTCDSXX_S38_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBQ2-1863_NGS19-J394_AHLGTCDSXX_S38_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBR1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBR1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBR1-1853_NGS19-J300_AHKK5WDSXX_S88_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CBR1-1853_NGS19-J300_BHJTTYDSXX_S34_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBR1-1853_NGS19-J300_AHKK5WDSXX_S88_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CBR1-1853_NGS19-J300_BHJTTYDSXX_S34_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBR2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBR2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBR2-1856_NGS19-J274_AHKK5WDSXX_S14_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CBR2-1856_NGS19-J274_BHJTTYDSXX_S8_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBR2-1856_NGS19-J274_AHKK5WDSXX_S14_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CBR2-1856_NGS19-J274_BHJTTYDSXX_S8_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBS1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBS1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBS1-2739_NGS19-J294_AHKK5WDSXX_S34_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CBS1-2739_NGS19-J294_BHJTTYDSXX_S28_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBS1-2739_NGS19-J294_AHKK5WDSXX_S34_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CBS1-2739_NGS19-J294_BHJTTYDSXX_S28_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBS2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBS2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBS2-2740_NGS19-J500_AHM2NJDSXX_S227_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBS2-2740_NGS19-J500_AHM2NJDSXX_S227_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBT1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBT1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBT1-2565_NGS19-J315_AHKK5WDSXX_S103_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CBT1-2565_NGS19-J315_BHJTTYDSXX_S49_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBT1-2565_NGS19-J315_AHKK5WDSXX_S103_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CBT1-2565_NGS19-J315_BHJTTYDSXX_S49_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBT2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBT2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBT2-2685_NGS19-J312_AHKK5WDSXX_S100_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CBT2-2685_NGS19-J312_BHJTTYDSXX_S46_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBT2-2685_NGS19-J312_AHKK5WDSXX_S100_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CBT2-2685_NGS19-J312_BHJTTYDSXX_S46_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CBV1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBV1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBV1-2741_NGS19-J540_AHM2NJDSXX_S267_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBV1-2741_NGS19-J540_AHM2NJDSXX_S267_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CBV2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CBV2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CBV2-2742_NGS19-J572_BHM33WDSXX_S24_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CBV2-2742_NGS19-J572_BHM33WDSXX_S24_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CCD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CCD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CCD1-1911_NGS19-J302_AHKK5WDSXX_S90_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CCD1-1911_NGS19-J302_BHJTTYDSXX_S36_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CCD1-1911_NGS19-J302_AHKK5WDSXX_S90_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CCD1-1911_NGS19-J302_BHJTTYDSXX_S36_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CCD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CCD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CCD2-1913_NGS19-J272_AHKK5WDSXX_S12_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CCD2-1913_NGS19-J272_BHJTTYDSXX_S6_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CCD2-1913_NGS19-J272_AHKK5WDSXX_S12_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CCD2-1913_NGS19-J272_BHJTTYDSXX_S6_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CCW1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CCW1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CCW1-1901_NGS19-J322_AHKK5WDSXX_S110_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CCW1-1901_NGS19-J322_BHJTTYDSXX_S56_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CCW1-1901_NGS19-J322_AHKK5WDSXX_S110_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CCW1-1901_NGS19-J322_BHJTTYDSXX_S56_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CCW2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CCW2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CCW2-2746_NGS19-J554_AHM2NJDSXX_S281_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CCW2-2746_NGS19-J554_AHM2NJDSXX_S281_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CCZ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CCZ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CCZ1-2457_NGS19-J553_AHM2NJDSXX_S280_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CCZ1-2457_NGS19-J553_AHM2NJDSXX_S280_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CCZ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CCZ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CCZ2-2458_NGS19-J470_AHM2NJDSXX_S18_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CCZ2-2458_NGS19-J470_AHM2NJDSXX_S18_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CDA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDA1-1892_NGS19-J273_AHKK5WDSXX_S13_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/CDA1-1892_NGS19-J273_BHJTTYDSXX_S7_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDA1-1892_NGS19-J273_AHKK5WDSXX_S13_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/CDA1-1892_NGS19-J273_BHJTTYDSXX_S7_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CDA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDA2-1908_NGS19-J408_AHLGTCDSXX_S52_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDA2-1908_NGS19-J408_AHLGTCDSXX_S52_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CDD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDD1-2615_NGS19-J424_AHLGTCDSXX_S68_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDD1-2615_NGS19-J424_AHLGTCDSXX_S68_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CDD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDD2-2616_NGS19-J561_BHM33WDSXX_S13_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDD2-2616_NGS19-J561_BHM33WDSXX_S13_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CDF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDF1-2472_NGS19-J443_AHLGTCDSXX_S87_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDF1-2472_NGS19-J443_AHLGTCDSXX_S87_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CDF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDF2-2456_NGS19-J458_AHLGTCDSXX_S102_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDF2-2456_NGS19-J458_AHLGTCDSXX_S102_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/CDH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDH1-2619_NGS19-J473_AHM2NJDSXX_S21_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDH1-2619_NGS19-J473_AHM2NJDSXX_S21_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/CDH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/CDH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/CDH2-2620_NGS19-J460_AHM2NJDSXX_S8_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/CDH2-2620_NGS19-J460_AHM2NJDSXX_S8_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/DAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/DAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/DAB1-2321_NGS19-J613_BHM33WDSXX_S65_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/DAB1-2321_NGS19-J613_BHM33WDSXX_S65_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/DAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/DAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/DAB2-2322_NGS19-J585_BHM33WDSXX_S37_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/DAB2-2322_NGS19-J585_BHM33WDSXX_S37_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/DAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/DAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/DAC1-2323_NGS19-J517_AHM2NJDSXX_S244_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/DAC1-2323_NGS19-J517_AHM2NJDSXX_S244_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/DAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/DAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/DAC2-2324_NGS19-J533_AHM2NJDSXX_S260_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/DAC2-2324_NGS19-J533_AHM2NJDSXX_S260_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/DAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/DAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/DAD1-2325_NGS19-J549_AHM2NJDSXX_S276_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/DAD1-2325_NGS19-J549_AHM2NJDSXX_S276_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/DAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/DAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/DAD2-2326_NGS19-J637_BHM33WDSXX_S89_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/DAD2-2326_NGS19-J637_BHM33WDSXX_S89_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAC1-2330_NGS19-J606_BHM33WDSXX_S58_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAC1-2330_NGS19-J606_BHM33WDSXX_S58_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAC2-2329_NGS19-J577_BHM33WDSXX_S29_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAC2-2329_NGS19-J577_BHM33WDSXX_S29_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAD1-2191_NGS19-J491_AHM2NJDSXX_S218_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAD1-2191_NGS19-J491_AHM2NJDSXX_S218_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAD2-2192_NGS19-J573_BHM33WDSXX_S25_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAD2-2192_NGS19-J573_BHM33WDSXX_S25_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAE1-2193_NGS19-J582_BHM33WDSXX_S34_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAE1-2193_NGS19-J582_BHM33WDSXX_S34_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAE2-2194_NGS19-J469_AHM2NJDSXX_S17_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAE2-2194_NGS19-J469_AHM2NJDSXX_S17_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAF1-2624_NGS19-J501_AHM2NJDSXX_S228_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAF1-2624_NGS19-J501_AHM2NJDSXX_S228_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAF2-2626_NGS19-J416_AHLGTCDSXX_S60_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAF2-2626_NGS19-J416_AHLGTCDSXX_S60_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAG1-2197_NGS19-J328_AHKK5WDSXX_S116_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/EAG1-2197_NGS19-J328_BHJTTYDSXX_S62_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAG1-2197_NGS19-J328_AHKK5WDSXX_S116_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/EAG1-2197_NGS19-J328_BHJTTYDSXX_S62_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAG2-2198_NGS19-J555_BHM33WDSXX_S7_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAG2-2198_NGS19-J555_BHM33WDSXX_S7_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAI1-2201_NGS19-J562_BHM33WDSXX_S14_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAI1-2201_NGS19-J562_BHM33WDSXX_S14_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAI2-2202_NGS19-J435_AHLGTCDSXX_S79_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAI2-2202_NGS19-J435_AHLGTCDSXX_S79_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAJ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAJ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAJ1-2235_NGS19-J335_AHKK5WDSXX_S123_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAJ1-2235_NGS19-J335_AHKK5WDSXX_S123_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAJ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAJ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAJ2-2387_NGS19-J423_AHLGTCDSXX_S67_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAJ2-2387_NGS19-J423_AHLGTCDSXX_S67_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAK1-2237_NGS19-J387_AHLGTCDSXX_S31_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAK1-2237_NGS19-J387_AHLGTCDSXX_S31_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAK2-2238_NGS19-J588_BHM33WDSXX_S40_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAK2-2238_NGS19-J588_BHM33WDSXX_S40_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAL1-2239_NGS19-J503_AHM2NJDSXX_S230_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAL1-2239_NGS19-J503_AHM2NJDSXX_S230_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAL2-2240_NGS19-J605_BHM33WDSXX_S57_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAL2-2240_NGS19-J605_BHM33WDSXX_S57_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAN1-2243_NGS19-J583_BHM33WDSXX_S35_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAN1-2243_NGS19-J583_BHM33WDSXX_S35_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAN2-2388_NGS19-J434_AHLGTCDSXX_S78_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAN2-2388_NGS19-J434_AHLGTCDSXX_S78_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAO1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAO1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAO1-2389_NGS19-J378_AHLGTCDSXX_S22_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAO1-2389_NGS19-J378_AHLGTCDSXX_S22_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAO2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAO2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAO2-2390_NGS19-J449_AHLGTCDSXX_S93_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAO2-2390_NGS19-J449_AHLGTCDSXX_S93_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAP1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAP1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAP1-2247_NGS19-J591_BHM33WDSXX_S43_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAP1-2247_NGS19-J591_BHM33WDSXX_S43_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAP2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAP2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAP2-2391_NGS19-J326_AHKK5WDSXX_S114_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/EAP2-2391_NGS19-J326_BHJTTYDSXX_S60_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAP2-2391_NGS19-J326_AHKK5WDSXX_S114_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/EAP2-2391_NGS19-J326_BHJTTYDSXX_S60_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAQ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAQ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAQ1-2249_NGS19-J525_AHM2NJDSXX_S252_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAQ1-2249_NGS19-J525_AHM2NJDSXX_S252_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAQ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAQ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAQ2-2250_NGS19-J535_AHM2NJDSXX_S262_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAQ2-2250_NGS19-J535_AHM2NJDSXX_S262_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAT1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAT1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAT1-2253_NGS19-J638_BHM33WDSXX_S90_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAT1-2253_NGS19-J638_BHM33WDSXX_S90_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAT2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAT2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAT2-2254_NGS19-J624_BHM33WDSXX_S76_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAT2-2254_NGS19-J624_BHM33WDSXX_S76_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAU1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAU1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAU1-2255_NGS19-J623_BHM33WDSXX_S75_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAU1-2255_NGS19-J623_BHM33WDSXX_S75_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAU2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAU2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAU2-2256_NGS19-J626_BHM33WDSXX_S78_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAU2-2256_NGS19-J626_BHM33WDSXX_S78_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAV1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAV1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAV1-2257_NGS19-J629_BHM33WDSXX_S81_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAV1-2257_NGS19-J629_BHM33WDSXX_S81_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAV2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAV2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAV2-2258_NGS19-J603_BHM33WDSXX_S55_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAV2-2258_NGS19-J603_BHM33WDSXX_S55_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAW1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAW1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAW1-2259_NGS19-J633_BHM33WDSXX_S85_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAW1-2259_NGS19-J633_BHM33WDSXX_S85_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAW2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAW2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAW2-2260_NGS19-J628_BHM33WDSXX_S80_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAW2-2260_NGS19-J628_BHM33WDSXX_S80_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAY1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAY1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAY1-2263_NGS19-J405_AHLGTCDSXX_S49_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAY1-2263_NGS19-J405_AHLGTCDSXX_S49_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAY2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAY2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAY2-2392_NGS19-J422_AHLGTCDSXX_S66_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAY2-2392_NGS19-J422_AHLGTCDSXX_S66_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EAZ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAZ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAZ1-2265_NGS19-J336_AHKK5WDSXX_S124_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAZ1-2265_NGS19-J336_AHKK5WDSXX_S124_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EAZ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EAZ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EAZ2-2266_NGS19-J476_AHM2NJDSXX_S24_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EAZ2-2266_NGS19-J476_AHM2NJDSXX_S24_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBB1-2269_NGS19-J528_AHM2NJDSXX_S255_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBB1-2269_NGS19-J528_AHM2NJDSXX_S255_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBB2-2270_NGS19-J508_AHM2NJDSXX_S235_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBB2-2270_NGS19-J508_AHM2NJDSXX_S235_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBC1-2271_NGS19-J534_AHM2NJDSXX_S261_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBC1-2271_NGS19-J534_AHM2NJDSXX_S261_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBC2-2272_NGS19-J518_AHM2NJDSXX_S245_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBC2-2272_NGS19-J518_AHM2NJDSXX_S245_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBD1-2273_NGS19-J480_AHM2NJDSXX_S28_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBD1-2273_NGS19-J480_AHM2NJDSXX_S28_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBD2-2274_NGS19-J579_BHM33WDSXX_S31_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBD2-2274_NGS19-J579_BHM33WDSXX_S31_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBF1-2277_NGS19-J494_AHM2NJDSXX_S221_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBF1-2277_NGS19-J494_AHM2NJDSXX_S221_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBF2-2278_NGS19-J351_AHKK5WDSXX_S139_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBF2-2278_NGS19-J351_AHKK5WDSXX_S139_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBG1-2279_NGS19-J459_AHM2NJDSXX_S7_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBG1-2279_NGS19-J459_AHM2NJDSXX_S7_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBG2-2280_NGS19-J527_AHM2NJDSXX_S254_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBG2-2280_NGS19-J527_AHM2NJDSXX_S254_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBH1-2281_NGS19-J362_AHKK5WDSXX_S150_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBH1-2281_NGS19-J362_AHKK5WDSXX_S150_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBH2-2282_NGS19-J520_AHM2NJDSXX_S247_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBH2-2282_NGS19-J520_AHM2NJDSXX_S247_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBK1-2285_NGS19-J438_AHLGTCDSXX_S82_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBK1-2285_NGS19-J438_AHLGTCDSXX_S82_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBK2-2286_NGS19-J413_AHLGTCDSXX_S57_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBK2-2286_NGS19-J413_AHLGTCDSXX_S57_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBL1-2287_NGS19-J392_AHLGTCDSXX_S36_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBL1-2287_NGS19-J392_AHLGTCDSXX_S36_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBL2-2288_NGS19-J631_BHM33WDSXX_S83_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBL2-2288_NGS19-J631_BHM33WDSXX_S83_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBM1-2289_NGS19-J492_AHM2NJDSXX_S219_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBM1-2289_NGS19-J492_AHM2NJDSXX_S219_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBM2-2398_NGS19-J369_AHLGTCDSXX_S13_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBM2-2398_NGS19-J369_AHLGTCDSXX_S13_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBN1-2399_NGS19-J357_AHKK5WDSXX_S145_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBN1-2399_NGS19-J357_AHKK5WDSXX_S145_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBN2-2400_NGS19-J395_AHLGTCDSXX_S39_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBN2-2400_NGS19-J395_AHLGTCDSXX_S39_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBO1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBO1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBO1-2401_NGS19-J509_AHM2NJDSXX_S236_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBO1-2401_NGS19-J509_AHM2NJDSXX_S236_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBO2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBO2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBO2-2294_NGS19-J471_AHM2NJDSXX_S19_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBO2-2294_NGS19-J471_AHM2NJDSXX_S19_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBP1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBP1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBP1-2295_NGS19-J361_AHKK5WDSXX_S149_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBP1-2295_NGS19-J361_AHKK5WDSXX_S149_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBR1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBR1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBR1-2297_NGS19-J490_AHM2NJDSXX_S38_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBR1-2297_NGS19-J490_AHM2NJDSXX_S38_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBR2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBR2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBR2-2298_NGS19-J410_AHLGTCDSXX_S54_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBR2-2298_NGS19-J410_AHLGTCDSXX_S54_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBU1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBU1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBU1-2627_NGS19-J329_AHKK5WDSXX_S117_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/EBU1-2627_NGS19-J329_BHJTTYDSXX_S63_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBU1-2627_NGS19-J329_AHKK5WDSXX_S117_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/EBU1-2627_NGS19-J329_BHJTTYDSXX_S63_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBU2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBU2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBU2-2628_NGS19-J543_AHM2NJDSXX_S270_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBU2-2628_NGS19-J543_AHM2NJDSXX_S270_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBV1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBV1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBV1-2629_NGS19-J538_AHM2NJDSXX_S265_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBV1-2629_NGS19-J538_AHM2NJDSXX_S265_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBV2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBV2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBV2-2630_NGS19-J523_AHM2NJDSXX_S250_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBV2-2630_NGS19-J523_AHM2NJDSXX_S250_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBW1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBW1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBW1-2403_NGS19-J446_AHLGTCDSXX_S90_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBW1-2403_NGS19-J446_AHLGTCDSXX_S90_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBW2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBW2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBW2-2404_NGS19-J445_AHLGTCDSXX_S89_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBW2-2404_NGS19-J445_AHLGTCDSXX_S89_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBX1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBX1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBX1-2631_NGS19-J419_AHLGTCDSXX_S63_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBX1-2631_NGS19-J419_AHLGTCDSXX_S63_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBX2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBX2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBX2-2632_NGS19-J330_AHKK5WDSXX_S118_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/EBX2-2632_NGS19-J330_BHJTTYDSXX_S64_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBX2-2632_NGS19-J330_AHKK5WDSXX_S118_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/EBX2-2632_NGS19-J330_BHJTTYDSXX_S64_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/EBY1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBY1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBY1-2205_NGS19-J627_BHM33WDSXX_S79_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBY1-2205_NGS19-J627_BHM33WDSXX_S79_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/EBY2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/EBY2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/EBY2-2206_NGS19-J617_BHM33WDSXX_S69_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/EBY2-2206_NGS19-J617_BHM33WDSXX_S69_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECA1-2309_NGS19-J368_AHLGTCDSXX_S12_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECA1-2309_NGS19-J368_AHLGTCDSXX_S12_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECA2-2310_NGS19-J511_AHM2NJDSXX_S238_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECA2-2310_NGS19-J511_AHM2NJDSXX_S238_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECB1-2311_NGS19-J417_AHLGTCDSXX_S61_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECB1-2311_NGS19-J417_AHLGTCDSXX_S61_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECB2-2312_NGS19-J575_BHM33WDSXX_S27_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECB2-2312_NGS19-J575_BHM33WDSXX_S27_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECD1-2313_NGS19-J522_AHM2NJDSXX_S249_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECD1-2313_NGS19-J522_AHM2NJDSXX_S249_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECD2-2314_NGS19-J545_AHM2NJDSXX_S272_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECD2-2314_NGS19-J545_AHM2NJDSXX_S272_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECE1-2406_NGS19-J452_AHLGTCDSXX_S96_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECE1-2406_NGS19-J452_AHLGTCDSXX_S96_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECE2-2407_NGS19-J451_AHLGTCDSXX_S95_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECE2-2407_NGS19-J451_AHLGTCDSXX_S95_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECF1-2408_NGS19-J448_AHLGTCDSXX_S92_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECF1-2408_NGS19-J448_AHLGTCDSXX_S92_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECF2-2409_NGS19-J447_AHLGTCDSXX_S91_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECF2-2409_NGS19-J447_AHLGTCDSXX_S91_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECG1-2410_NGS19-J457_AHLGTCDSXX_S101_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECG1-2410_NGS19-J457_AHLGTCDSXX_S101_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECG2-2411_NGS19-J456_AHLGTCDSXX_S100_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECG2-2411_NGS19-J456_AHLGTCDSXX_S100_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECH1-2412_NGS19-J455_AHLGTCDSXX_S99_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECH1-2412_NGS19-J455_AHLGTCDSXX_S99_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECH2-2413_NGS19-J590_BHM33WDSXX_S42_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECH2-2413_NGS19-J590_BHM33WDSXX_S42_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECI1-2414_NGS19-J454_AHLGTCDSXX_S98_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECI1-2414_NGS19-J454_AHLGTCDSXX_S98_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECI2-2415_NGS19-J453_AHLGTCDSXX_S97_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECI2-2415_NGS19-J453_AHLGTCDSXX_S97_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECK1-2633_NGS19-J483_AHM2NJDSXX_S31_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECK1-2633_NGS19-J483_AHM2NJDSXX_S31_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECK2-2634_NGS19-J397_AHLGTCDSXX_S41_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECK2-2634_NGS19-J397_AHLGTCDSXX_S41_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/ECN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECN1-2635_NGS19-J350_AHKK5WDSXX_S138_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECN1-2635_NGS19-J350_AHKK5WDSXX_S138_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/ECN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/ECN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/ECN2-2636_NGS19-J365_AHLGTCDSXX_S9_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/ECN2-2636_NGS19-J365_AHLGTCDSXX_S9_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAB1-2134_NGS19-J580_BHM33WDSXX_S32_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAB1-2134_NGS19-J580_BHM33WDSXX_S32_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAB2-2135_NGS19-J426_AHLGTCDSXX_S70_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAB2-2135_NGS19-J426_AHLGTCDSXX_S70_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAB2-replicate";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAB2-replicate" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAB2-replicate-2135_NGS19-J635_BHM33WDSXX_S87_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAB2-replicate-2135_NGS19-J635_BHM33WDSXX_S87_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAD1-2138_NGS19-J472_AHM2NJDSXX_S20_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAD1-2138_NGS19-J472_AHM2NJDSXX_S20_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAD2-2139_NGS19-J578_BHM33WDSXX_S30_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAD2-2139_NGS19-J578_BHM33WDSXX_S30_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAF1-2142_NGS19-J402_AHLGTCDSXX_S46_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAF1-2142_NGS19-J402_AHLGTCDSXX_S46_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAF2-2143_NGS19-J432_AHLGTCDSXX_S76_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAF2-2143_NGS19-J432_AHLGTCDSXX_S76_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAG1-2144_NGS19-J486_AHM2NJDSXX_S34_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAG1-2144_NGS19-J486_AHM2NJDSXX_S34_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAG2-2145_NGS19-J581_BHM33WDSXX_S33_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAG2-2145_NGS19-J581_BHM33WDSXX_S33_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAH2-2147_NGS19-J412_AHLGTCDSXX_S56_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAH2-2147_NGS19-J412_AHLGTCDSXX_S56_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAI1-2148_NGS19-J375_AHLGTCDSXX_S19_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAI1-2148_NGS19-J375_AHLGTCDSXX_S19_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAI2-2149_NGS19-J384_AHLGTCDSXX_S28_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAI2-2149_NGS19-J384_AHLGTCDSXX_S28_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAJ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAJ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAJ1-2150_NGS19-J429_AHLGTCDSXX_S73_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAJ1-2150_NGS19-J429_AHLGTCDSXX_S73_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAJ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAJ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAJ2-2151_NGS19-J363_AHLGTCDSXX_S7_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAJ2-2151_NGS19-J363_AHLGTCDSXX_S7_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAK1-2152_NGS19-J374_AHLGTCDSXX_S18_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAK1-2152_NGS19-J374_AHLGTCDSXX_S18_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAK2-2378_NGS19-J499_AHM2NJDSXX_S226_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAK2-2378_NGS19-J499_AHM2NJDSXX_S226_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAL1-2154_NGS19-J441_AHLGTCDSXX_S85_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAL1-2154_NGS19-J441_AHLGTCDSXX_S85_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAL2-2155_NGS19-J475_AHM2NJDSXX_S23_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAL2-2155_NGS19-J475_AHM2NJDSXX_S23_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAM1-2156_NGS19-J385_AHLGTCDSXX_S29_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAM1-2156_NGS19-J385_AHLGTCDSXX_S29_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAM2-2157_NGS19-J569_BHM33WDSXX_S21_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAM2-2157_NGS19-J569_BHM33WDSXX_S21_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAN1-2158_NGS19-J403_AHLGTCDSXX_S47_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAN1-2158_NGS19-J403_AHLGTCDSXX_S47_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAN2-2159_NGS19-J440_AHLGTCDSXX_S84_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAN2-2159_NGS19-J440_AHLGTCDSXX_S84_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAP1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAP1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAP1-2162_NGS19-J404_AHLGTCDSXX_S48_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAP1-2162_NGS19-J404_AHLGTCDSXX_S48_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAP2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAP2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAP2-2163_NGS19-J526_AHM2NJDSXX_S253_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAP2-2163_NGS19-J526_AHM2NJDSXX_S253_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/FAQ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAQ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAQ1-2381_NGS19-J325_AHKK5WDSXX_S113_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/FAQ1-2381_NGS19-J325_BHJTTYDSXX_S59_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAQ1-2381_NGS19-J325_AHKK5WDSXX_S113_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/FAQ1-2381_NGS19-J325_BHJTTYDSXX_S59_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/FAQ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/FAQ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/FAQ2-2382_NGS19-J524_AHM2NJDSXX_S251_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/FAQ2-2382_NGS19-J524_AHM2NJDSXX_S251_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAA1-2419_NGS19-J558_BHM33WDSXX_S10_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAA1-2419_NGS19-J558_BHM33WDSXX_S10_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAA2-2420_NGS19-J481_AHM2NJDSXX_S29_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAA2-2420_NGS19-J481_AHM2NJDSXX_S29_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAE1-2465_NGS19-J546_AHM2NJDSXX_S273_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAE1-2465_NGS19-J546_AHM2NJDSXX_S273_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAE2-2426_NGS19-J324_AHKK5WDSXX_S112_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/GAE2-2426_NGS19-J324_BHJTTYDSXX_S58_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAE2-2426_NGS19-J324_AHKK5WDSXX_S112_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/GAE2-2426_NGS19-J324_BHJTTYDSXX_S58_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAG1-2429_NGS19-J366_AHLGTCDSXX_S10_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAG1-2429_NGS19-J366_AHLGTCDSXX_S10_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAG2-2430_NGS19-J565_BHM33WDSXX_S17_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAG2-2430_NGS19-J565_BHM33WDSXX_S17_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAH1-2431_NGS19-J531_AHM2NJDSXX_S258_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAH1-2431_NGS19-J531_AHM2NJDSXX_S258_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAH2-2432_NGS19-J620_BHM33WDSXX_S72_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAH2-2432_NGS19-J620_BHM33WDSXX_S72_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAI1-2433_NGS19-J607_BHM33WDSXX_S59_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAI1-2433_NGS19-J607_BHM33WDSXX_S59_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAI2-2474_NGS19-J521_AHM2NJDSXX_S248_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAI2-2474_NGS19-J521_AHM2NJDSXX_S248_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAJ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAJ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAJ1-2435_NGS19-J593_BHM33WDSXX_S45_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAJ1-2435_NGS19-J593_BHM33WDSXX_S45_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAJ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAJ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAJ2-2436_NGS19-J621_BHM33WDSXX_S73_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAJ2-2436_NGS19-J621_BHM33WDSXX_S73_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAK1-2437_NGS19-J574_BHM33WDSXX_S26_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAK1-2437_NGS19-J574_BHM33WDSXX_S26_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAK2-2438_NGS19-J571_BHM33WDSXX_S23_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAK2-2438_NGS19-J571_BHM33WDSXX_S23_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAL1-2439_NGS19-J584_BHM33WDSXX_S36_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAL1-2439_NGS19-J584_BHM33WDSXX_S36_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAL2-2440_NGS19-J612_BHM33WDSXX_S64_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAL2-2440_NGS19-J612_BHM33WDSXX_S64_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAM1-2441_NGS19-J615_BHM33WDSXX_S67_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAM1-2441_NGS19-J615_BHM33WDSXX_S67_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAM2-2442_NGS19-J596_BHM33WDSXX_S48_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAM2-2442_NGS19-J596_BHM33WDSXX_S48_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAN1-2571_NGS19-J594_BHM33WDSXX_S46_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAN1-2571_NGS19-J594_BHM33WDSXX_S46_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAN2-2444_NGS19-J634_BHM33WDSXX_S86_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAN2-2444_NGS19-J634_BHM33WDSXX_S86_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAO1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAO1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAO1-2445_NGS19-J586_BHM33WDSXX_S38_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAO1-2445_NGS19-J586_BHM33WDSXX_S38_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAO2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAO2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAO2-2446_NGS19-J602_BHM33WDSXX_S54_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAO2-2446_NGS19-J602_BHM33WDSXX_S54_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAP1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAP1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAP1-2447_NGS19-J564_BHM33WDSXX_S16_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAP1-2447_NGS19-J564_BHM33WDSXX_S16_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAP2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAP2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAP2-2448_NGS19-J604_BHM33WDSXX_S56_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAP2-2448_NGS19-J604_BHM33WDSXX_S56_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAQ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAQ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAQ1-2449_NGS19-J608_BHM33WDSXX_S60_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAQ1-2449_NGS19-J608_BHM33WDSXX_S60_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAQ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAQ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAQ2-2450_NGS19-J619_BHM33WDSXX_S71_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAQ2-2450_NGS19-J619_BHM33WDSXX_S71_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAR1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAR1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAR1-2451_NGS19-J622_BHM33WDSXX_S74_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAR1-2451_NGS19-J622_BHM33WDSXX_S74_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAR2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAR2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAR2-2452_NGS19-J597_BHM33WDSXX_S49_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAR2-2452_NGS19-J597_BHM33WDSXX_S49_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/GAS1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAS1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAS1-2453_NGS19-J516_AHM2NJDSXX_S243_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAS1-2453_NGS19-J516_AHM2NJDSXX_S243_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAS2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAS2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAS2-2454_NGS19-J515_AHM2NJDSXX_S242_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAS2-2454_NGS19-J515_AHM2NJDSXX_S242_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/GAS2-replicate";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/GAS2-replicate" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/GAS2-repl_NGS19-J576_BHM33WDSXX_S28_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/GAS2-repl_NGS19-J576_BHM33WDSXX_S28_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAA1-2519_NGS19-J327_AHKK5WDSXX_S115_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/HAA1-2519_NGS19-J327_BHJTTYDSXX_S61_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAA1-2519_NGS19-J327_AHKK5WDSXX_S115_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/HAA1-2519_NGS19-J327_BHJTTYDSXX_S61_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAA2-2520_NGS19-J415_AHLGTCDSXX_S59_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAA2-2520_NGS19-J415_AHLGTCDSXX_S59_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAB1-2521_NGS19-J389_AHLGTCDSXX_S33_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAB1-2521_NGS19-J389_AHLGTCDSXX_S33_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAB2-2522_NGS19-J411_AHLGTCDSXX_S55_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAB2-2522_NGS19-J411_AHLGTCDSXX_S55_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAC1-2523_NGS19-J505_AHM2NJDSXX_S232_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAC1-2523_NGS19-J505_AHM2NJDSXX_S232_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAC2-2524_NGS19-J379_AHLGTCDSXX_S23_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAC2-2524_NGS19-J379_AHLGTCDSXX_S23_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAD1-2525_NGS19-J536_AHM2NJDSXX_S263_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAD1-2525_NGS19-J536_AHM2NJDSXX_S263_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAD2-2526_NGS19-J421_AHLGTCDSXX_S65_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAD2-2526_NGS19-J421_AHLGTCDSXX_S65_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAE1-2527_NGS19-J427_AHLGTCDSXX_S71_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAE1-2527_NGS19-J427_AHLGTCDSXX_S71_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAE2-2528_NGS19-J409_AHLGTCDSXX_S53_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAE2-2528_NGS19-J409_AHLGTCDSXX_S53_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAF1-2529_NGS19-J557_BHM33WDSXX_S9_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAF1-2529_NGS19-J557_BHM33WDSXX_S9_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAF2-2530_NGS19-J442_AHLGTCDSXX_S86_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAF2-2530_NGS19-J442_AHLGTCDSXX_S86_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAG1-2531_NGS19-J566_BHM33WDSXX_S18_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAG1-2531_NGS19-J566_BHM33WDSXX_S18_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAG2-2532_NGS19-J414_AHLGTCDSXX_S58_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAG2-2532_NGS19-J414_AHLGTCDSXX_S58_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAH1-2533_NGS19-J570_BHM33WDSXX_S22_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAH1-2533_NGS19-J570_BHM33WDSXX_S22_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAH2-2534_NGS19-J493_AHM2NJDSXX_S220_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAH2-2534_NGS19-J493_AHM2NJDSXX_S220_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAI1-2535_NGS19-J625_BHM33WDSXX_S77_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAI1-2535_NGS19-J625_BHM33WDSXX_S77_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAI2-2536_NGS19-J544_AHM2NJDSXX_S271_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAI2-2536_NGS19-J544_AHM2NJDSXX_S271_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/HAK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAK1-2584_NGS19-J589_BHM33WDSXX_S41_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAK1-2584_NGS19-J589_BHM33WDSXX_S41_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/HAK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/HAK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/HAK2-2540_NGS19-J530_AHM2NJDSXX_S257_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/HAK2-2540_NGS19-J530_AHM2NJDSXX_S257_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAB1-2701_NGS19-J316_AHKK5WDSXX_S104_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAB1-2701_NGS19-J316_BHJTTYDSXX_S50_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAB1-2701_NGS19-J316_AHKK5WDSXX_S104_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAB1-2701_NGS19-J316_BHJTTYDSXX_S50_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAB2-2702_NGS19-J277_AHKK5WDSXX_S17_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAB2-2702_NGS19-J277_BHJTTYDSXX_S11_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAB2-2702_NGS19-J277_AHKK5WDSXX_S17_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAB2-2702_NGS19-J277_BHJTTYDSXX_S11_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAC1-2703_NGS19-J592_BHM33WDSXX_S44_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAC1-2703_NGS19-J592_BHM33WDSXX_S44_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAC2-2704_NGS19-J313_AHKK5WDSXX_S101_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAC2-2704_NGS19-J313_BHJTTYDSXX_S47_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAC2-2704_NGS19-J313_AHKK5WDSXX_S101_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAC2-2704_NGS19-J313_BHJTTYDSXX_S47_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAD1-2705_NGS19-J276_AHKK5WDSXX_S16_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAD1-2705_NGS19-J276_BHJTTYDSXX_S10_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAD1-2705_NGS19-J276_AHKK5WDSXX_S16_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAD1-2705_NGS19-J276_BHJTTYDSXX_S10_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAD2-2706_NGS19-J559_BHM33WDSXX_S11_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAD2-2706_NGS19-J559_BHM33WDSXX_S11_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAE1-2707_NGS19-J275_AHKK5WDSXX_S15_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAE1-2707_NGS19-J275_BHJTTYDSXX_S9_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAE1-2707_NGS19-J275_AHKK5WDSXX_S15_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAE1-2707_NGS19-J275_BHJTTYDSXX_S9_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAE2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAE2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAE2-2709_NGS19-J290_AHKK5WDSXX_S30_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAE2-2709_NGS19-J290_BHJTTYDSXX_S24_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAE2-2709_NGS19-J290_AHKK5WDSXX_S30_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAE2-2709_NGS19-J290_BHJTTYDSXX_S24_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAF1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAF1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAF1-2710_NGS19-J284_AHKK5WDSXX_S24_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAF1-2710_NGS19-J284_BHJTTYDSXX_S18_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAF1-2710_NGS19-J284_AHKK5WDSXX_S24_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAF1-2710_NGS19-J284_BHJTTYDSXX_S18_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAF2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAF2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAF2-2711_NGS19-J587_BHM33WDSXX_S39_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAF2-2711_NGS19-J587_BHM33WDSXX_S39_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAG1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAG1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAG1-2712_NGS19-J356_AHKK5WDSXX_S144_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAG1-2712_NGS19-J356_AHKK5WDSXX_S144_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAG2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAG2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAG2-2713_NGS19-J301_AHKK5WDSXX_S89_L003_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAG2-2713_NGS19-J301_BHJTTYDSXX_S35_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAG2-2713_NGS19-J301_AHKK5WDSXX_S89_L003_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAG2-2713_NGS19-J301_BHJTTYDSXX_S35_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAH1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAH1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAH1-2714_NGS19-J401_AHLGTCDSXX_S45_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAH1-2714_NGS19-J401_AHLGTCDSXX_S45_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAH2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAH2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAH2-2715_NGS19-J498_AHM2NJDSXX_S225_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAH2-2715_NGS19-J498_AHM2NJDSXX_S225_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAI1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAI1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAI1-2716_NGS19-J297_AHKK5WDSXX_S37_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAI1-2716_NGS19-J297_BHJTTYDSXX_S31_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAI1-2716_NGS19-J297_AHKK5WDSXX_S37_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAI1-2716_NGS19-J297_BHJTTYDSXX_S31_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAI2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAI2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAI2-2717_NGS19-J618_BHM33WDSXX_S70_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAI2-2717_NGS19-J618_BHM33WDSXX_S70_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAJ1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAJ1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAJ1-2718_NGS19-J358_AHKK5WDSXX_S146_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAJ1-2718_NGS19-J358_AHKK5WDSXX_S146_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAJ2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAJ2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAJ2-2720_NGS19-J288_AHKK5WDSXX_S28_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAJ2-2720_NGS19-J288_BHJTTYDSXX_S22_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAJ2-2720_NGS19-J288_AHKK5WDSXX_S28_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAJ2-2720_NGS19-J288_BHJTTYDSXX_S22_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAK1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAK1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAK1-2721_NGS19-J283_AHKK5WDSXX_S23_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAK1-2721_NGS19-J283_BHJTTYDSXX_S17_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAK1-2721_NGS19-J283_AHKK5WDSXX_S23_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAK1-2721_NGS19-J283_BHJTTYDSXX_S17_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAK2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAK2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAK2-2722_NGS19-J298_AHKK5WDSXX_S38_L002_001.clean_yx_R1.fastq.gz,data2/fastq-clean/JAK2-2722_NGS19-J298_BHJTTYDSXX_S32_L001_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAK2-2722_NGS19-J298_AHKK5WDSXX_S38_L002_001.clean_yx_R2.fastq.gz,data2/fastq-clean/JAK2-2722_NGS19-J298_BHJTTYDSXX_S32_L001_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAL1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAL1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAL1-2723_NGS19-J550_AHM2NJDSXX_S277_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAL1-2723_NGS19-J550_AHM2NJDSXX_S277_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAL2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAL2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAL2-2724_NGS19-J425_AHLGTCDSXX_S69_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAL2-2724_NGS19-J425_AHLGTCDSXX_S69_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAM1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAM1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAM1-2725_NGS19-J386_AHLGTCDSXX_S30_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAM1-2725_NGS19-J386_AHLGTCDSXX_S30_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAM2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAM2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAM2-2726_NGS19-J377_AHLGTCDSXX_S21_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAM2-2726_NGS19-J377_AHLGTCDSXX_S21_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/JAN1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAN1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAN1-2727_NGS19-J372_AHLGTCDSXX_S16_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAN1-2727_NGS19-J372_AHLGTCDSXX_S16_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/JAN2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/JAN2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/JAN2-2728_NGS19-J433_AHLGTCDSXX_S77_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/JAN2-2728_NGS19-J433_AHLGTCDSXX_S77_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/KAA1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAA1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAA1-2572_NGS19-J563_BHM33WDSXX_S15_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAA1-2572_NGS19-J563_BHM33WDSXX_S15_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/KAA2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAA2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAA2-2573_NGS19-J463_AHM2NJDSXX_S11_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAA2-2573_NGS19-J463_AHM2NJDSXX_S11_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/KAB1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAB1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAB1-2697_NGS19-J477_AHM2NJDSXX_S25_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAB1-2697_NGS19-J477_AHM2NJDSXX_S25_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/KAB2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAB2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAB2-2575_NGS19-J609_BHM33WDSXX_S61_L003_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAB2-2575_NGS19-J609_BHM33WDSXX_S61_L003_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/KAC1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAC1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAC1-2576_NGS19-J567_BHM33WDSXX_S19_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAC1-2576_NGS19-J567_BHM33WDSXX_S19_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/KAC2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAC2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAC2-2577_NGS19-J465_AHM2NJDSXX_S13_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAC2-2577_NGS19-J465_AHM2NJDSXX_S13_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/KAD1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAD1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAD1-2578_NGS19-J467_AHM2NJDSXX_S15_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAD1-2578_NGS19-J467_AHM2NJDSXX_S15_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM


mkdir -p "processed/KAD2";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAD2" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAD2-2698_NGS19-J560_BHM33WDSXX_S12_L002_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAD2-2698_NGS19-J560_BHM33WDSXX_S12_L002_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

mkdir -p "processed/KAE1";
nice ./bin/STAR/bin/Linux_x86_64/STAR \
    --outFileNamePrefix "processed/KAE1" \
    --runThreadN 84 \
    --genomeDir ref/star-hg19 \
    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \
    --readFilesIn \
        data2/fastq-clean/KAE1-2580_NGS19-J552_AHM2NJDSXX_S279_L004_001.clean_yx_R1.fastq.gz \
        data2/fastq-clean/KAE1-2580_NGS19-J552_AHM2NJDSXX_S279_L004_001.clean_yx_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJoverhangMin 10 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --chimOutType WithinBAM SeparateSAMold \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --quantMode TranscriptomeSAM



# - - - - -

