#!/usr/bin/env python

import os
from tqdm import tqdm

path = 'data/DNA/cnvkit/data'




tt = {}
with open('data/DNA/sample codes sanger gsam.txt', 'r') as fh_in:
	for line in tqdm(fh_in):
		line = line.strip().split("\t")
		if line[0] != "Data Received":
			if line[13] in tt:
				print("ERROR DUPLICATE")
				import sys
				sys.exit(1)

			tt[line[13]] = line[5]
			#print( line[5]  + " -> " + line[13] )






printed_header = False

with open("output/tables/cnv_copynumber-ratio.cns_all.txt", "w") as fh_out:
	for _ in tqdm(sorted([_ for _ in os.listdir(path) if _[-4:] == ".cns" and _.find(' (1)') == -1])):
		sid = _.split("__",1)[0]

		#print(sid)

		with open(path + "/" + _, 'r') as fh_in:
			for line in fh_in:
				if line[0:10] != "chromosome":
					fh_out.write(tt[sid] + "\t" + sid + "\t" + line.strip() + "\n")
				else:
					if printed_header == False:
						fh_out.write("patient-id\tsample-id\t" + line.strip() + "\n")
						printed_header = True




