#!/usr/bin/env python

import os
from tqdm import tqdm

path = 'data/DNA/cnvkit/data'


"""
file sizes should be: 27493 or 27287   || first series always  27493

    20622 data/DNA/cnvkit/first_series/PD30245a__tumour.cnr
        0 data/DNA/cnvkit/first_series/PD30245c__tumour.cnr
        0 data/DNA/cnvkit/first_series/PD30247c__tumour.cnr
    18524 data/DNA/cnvkit/first_series/PD30257c__tumour.cnr
    

# Exclude
'PD36772a__tumour (1).cnr'
"""

failed = ["PD30245a__tumour.cnr", "PD30245c__tumour.cnr", "PD30247c__tumour.cnr", "PD30257c__tumour.cnr"]

idx = {}
gidx = {}
locs = []
samples = []

# 01. find all files
for _ in tqdm(sorted([_ for _ in os.listdir(path) if _[-4:] == ".cnr" and _.find(' (1)') == -1])):
	sid = _.split("__",1)[0]
	samples.append(sid)

	if _ not in failed:
		with open(path + "/" + _ , "r") as fh:
			for line in fh:
				line = line.strip()
				params = line.split("\t")
				if params[0] != "chromosome":
					#print(_, params)
					loc = params[0] + ":" + params[1] + "-" + params[2]

					if loc not in idx:
						locs.append(loc)
						idx[loc] = {}
						#print(_, params)
						gidx[loc] = [params[0], params[1], params[2], params[3]]
					
					if sid not in idx[loc]:
						idx[loc][sid] = {'depth': params[4], 'log2': params[5]}
	else:
		print("Skipping incomplete file: " + _)



with open("output/tables/cnv_copynumber-ratio.cnr_log2_all.txt", "w") as fh:
	out = "chromosome\tstart\tend\tgene"
	for sid in samples:
		out += "\t" + sid
	out += "\n"

	fh.write(out)
	
	for loc in tqdm(sorted(locs)):
		out = "\t".join(gidx[loc])
		for sid in samples:
			if sid in idx[loc]:
				out += "\t" + idx[loc][sid]['log2']
			else:
				out += "\t" + "NA"
		out += "\n"
		
		fh.write(out)
			

with open("output/tables/cnv_copynumber-ratio.cnr_depth_all.txt", "w") as fh:
	out = "chromosome\tstart\tend\tgene"
	for sid in samples:
		out += "\t" + sid
	out += "\n"

	fh.write(out)
	
	for loc in tqdm(sorted(locs)):
		out = "\t".join(gidx[loc])
		for sid in samples:
			if sid in idx[loc]:
				out += "\t" + idx[loc][sid]['depth']
			else:
				out += "\t" + "NA"
		out += "\n"
		
		fh.write(out)
			



