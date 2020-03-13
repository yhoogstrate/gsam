#!/usr/bin/env python

import os, pysam
from tqdm import tqdm

def concat(file_a, file_b, file_target):
    with open(file_target, 'w') as outfile:
        for fname in [file_a, file_b]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

with open("failed_disco_alignments.txt", "w") as fh_out:

    folder = 'output/bam-disco'
    for _ in tqdm(os.listdir(folder)):
        print("File: " + _)
        
        raw_sam = folder + '/' + _ + '/Chimeric.out.sam' 
        if not os.path.isfile(raw_sam) or (os.path.isfile(raw_sam) and os.path.getsize(raw_sam) == 0):
            fh_out.write("Missing or empty: " + raw_sam)
        else:
            new_sam = '/dev/shm/Chimeric.out.'+ _ +'.head.sam' 
            fixed_bam = folder + '/' + _ + '/Chimeric.out.'+ _ +'.fixed.bam' 

            if not os.path.isfile(fixed_bam) or (os.path.isfile(fixed_bam) and os.path.getsize(fixed_bam) == 0):
                concat('ref/star-hg19/sam-header.txt' ,raw_sam, new_sam)
                #print('\n')

                os.system('dr-disco fix "'+new_sam+'" "'+fixed_bam+'"')
                os.remove(new_sam)
                print("")

            else:
                print(" - skipping, file exists")
