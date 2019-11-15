#!/usr/bin/env python

import os
from tqdm import tqdm

path = 'output/tables/tin.py'
sids = [_ for _ in os.listdir(path) if os.path.isdir(path + '/'+ _)]

idx = {}
header = False

for sid in tqdm(sids):
	fname = path + "/" + sid + "/" + sid + ".summary.txt"
	with open(fname, 'r') as fh:
		for line in fh:
			line = line.strip().split("\t")
			if line[0] != 'Bam_file':
				idx[line[0]] = line
			else:
				if header == False:
					header = line


with open('output/tables/tin.py/tin.py.matrix.txt', 'w') as fh:
	fh.write("\t".join(header) + "\n")
	for _ in sorted(idx.values()):
		fh.write("\t".join(_))

