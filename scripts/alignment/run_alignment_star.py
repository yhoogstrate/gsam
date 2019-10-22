#!/usr/bin/env python

# 1. find all cleaned fastq pairs
# 2. merge sorted based on patient identifier
# 3. generate star run script

import os
import re

idx = {}

with open("scripts/alignment/run_alignment_star.sh", "w") as fh:
	for _ in sorted([_ for _ in os.listdir("data2/fastq-clean") if _[-9:] == ".fastq.gz"]):
		#sid = _.split("-",1)[0]
		sid = re.split("\-[0-9]+_NGS19", _)[0]
		sid = re.sub(r"-repl_.+$", "-replicate", sid)
		
		pid = sid[0:3]# + sid[4:]
		res = sid[3]
		srid = _.split(".clean",1)[0]
		mate = _.split("_yx_",1)[1].split(".",1)[0]

		if pid not in idx:
			idx[pid] = {}
		
		if sid not in idx[pid]:
			idx[pid][sid] = {}
		
		idx[pid][sid][srid] = {'R1': None, 'R2': None}
		
		
		#print(pid, sid, res, srid, mate)

for pid in sorted(idx):
	print(pid)
	for sid in sorted(idx[pid]):
		print("  " + sid + ":")
		for srid in sorted(idx[pid][sid]):
			print("    " + srid)
	print("")
