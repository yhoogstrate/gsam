#!/usr/bin/env R

# add to gsam.metadata script
gsam.dna_seq.metadata <- read.delim("data/DNA/sample codes sanger gsam.txt")
#gsam.dna_seq.metadata <- gsam.dna_seq.metadata[,colnames(gsam.dna_seq.metadata) %in% c("donor_sex", "donor_ID", "PD_ID")]
gsam.dna_seq.metadata <- gsam.dna_seq.metadata[order(gsam.dna_seq.metadata$donor_ID, gsam.dna_seq.metadata$PD_ID),]
#gsam.dna_seq.metadata$dup <- duplicated(gsam.dna_seq.metadata$donor_ID)

