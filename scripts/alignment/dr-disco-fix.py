#!/usr/bin/env python

import os, pysam
from tqdm import tqdm

def concat(file_a, file_b, file_target):
	with open(file_target, 'w') as outfile:
		for fname in [file_a, file_b]:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

folder = 'output/bam-disco'
for _ in tqdm(os.listdir(folder)):
	#print(_)
	raw_sam = folder + '/' + _ + '/Chimeric.out.sam' 
	new_sam = '/dev/shm/Chimeric.out.'+ _ +'.head.sam' 
	fixed_bam = folder + '/' + _ + '/Chimeric.out.'+ _ +'.fixed.bam' 
	#print(raw_sam)
	#print(new_sam)
	
	concat('ref/star-hg19/sam-header.txt' ,raw_sam, new_sam)
	#print('\n')

	
	os.system('dr-disco fix "'+new_sam+'" "'+fixed_bam+'"')
	os.remove(new_sam)

