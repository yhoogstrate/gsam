#!/usr/bin/env python

import os
from tqdm import tqdm

paths = ['data/DNA/cnvkit/ExtraSequencing_CN_data_29Mar2018', 'data/DNA/cnvkit/first_series']

failed = ["PD21369b__normal.antitargetcoverage.cnn", "PD21369b__normal.targetcoverage.cnn", "PD30218c__tumour.targetcoverage.cnn", "PD30219a__tumour.antitargetcoverage.cnn", "PD30219a__tumour.targetcoverage.cnn", "PD30236a__tumour.antitargetcoverage.cnn", "PD30243a__tumour.antitargetcoverage.cnn", "PD30244c__tumour.targetcoverage.cnn", "PD30245c__tumour.antitargetcoverage.cnn", "PD30253c__tumour.antitargetcoverage.cnn", "PD30257a__tumour.targetcoverage.cnn", "PD30257c__tumour.antitargetcoverage.cnn", "PD30258c__tumour.targetcoverage.cnn", "PD30260c__tumour.targetcoverage.cnn"]

gidx = {}

idx = {}

locs = []
samples = []

# 01. find all files
for path in sorted(paths):
	for _ in tqdm(sorted([_ for _ in os.listdir(path) if _[-19:] == ".targetcoverage.cnn" ])):
		sid = _.split("__",1)[0]
		samples.append(sid)

		if _ not in failed:
			with open(path + "/" + _ , "r") as fh:
				for line in fh:
					line = line.strip()
					params = line.split("\t")
					if params[0] != "chromosome":
						loc = params[0] + ":" + params[1] + "-" + params[2]

						if loc not in idx:
							locs.append(loc)
							idx[loc] = {}
							gidx[loc] = [params[0], params[1], params[2], params[3]]
						
						if sid not in idx[loc]:
							idx[loc][sid] = params[4]




with open("output/tables/cnv_tagetcoverage.cnn_all.txt", "w") as fh:
	out = "chromosome\tstart\tend\tgene"
	for sid in samples:
		out += "\t" + sid
	out += "\n"

	fh.write(out)
	
	for loc in tqdm(locs):
		out = "\t".join(gidx[loc])
		for sid in samples:
			if sid in idx[loc]:
				out += "\t" + idx[loc][sid]
			else:
				out += "\t" + "NA"
		out += "\n"
		
		fh.write(out)
			
