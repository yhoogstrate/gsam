#!/usr/bin/env python

import os
from tqdm import tqdm

path = "output/tables/idxstats"

files = sorted([_ for _ in os.listdir(path) if _[-14:] == ".flagstats.txt"])

chrs = []
chrs_idx = {}
first_sample = True

data = {}

for _ in tqdm(files):
	with open(os.path.join(path, _), 'r') as fh:
		for line in fh:
			line = line.strip().split("\t")
			#print(line)

			if first_sample:
				chrs.append(line[0])
				chrs_idx[line[0]] = line[1]
				data[line[0]] = []

			data[line[0]].append(line[2])
		
		if first_sample:
			first_sample = False



out = ["reference"] + files

with open("output/tables/idxstats/samtools.indexstats.matrix.txt","w") as fh:
	header = ["ref", "ref-len"] + files
	fh.write("\t".join(header) + "\n")
	for ref in chrs:
		fh.write("\t".join([ref, chrs_idx[ref]] + data[ref]) + "\n")


