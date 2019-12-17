#!/usr/bin/env python

import os, pysam
from tqdm import tqdm

p = ['data/RNA/alignments']
for _ in tqdm(os.listdir(p[0])):
	p = ['data/RNA/alignments']
	p.append(_)
	bam = "/".join(p) + "/Aligned.sortedByCoord.out.markdup.bam" 
	sam = "/".join(p) + "/Chimeric.out.sam"
	sam_full = "output/Chimeric.out." + _ + ".sam"
	with open(sam_full, "w") as fh_out:
		fh_out.write(pysam.view(bam,'-H')) # header
		with open(sam, 'r') as fh:
			for line in fh:
				fh_out.write(line)
	
	
	bam_fixed  = "output/Chimeric.out." + _ + ".fixed.bam"
	
	os.system('dr-disco fix "'+sam_full+'" "'+bam_fixed+'"')

