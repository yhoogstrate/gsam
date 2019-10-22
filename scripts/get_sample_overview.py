#!/usr/bin/env python

import json


failed = {}
err = "/mnt/ccbc/home/yhoogstrate/projects/gsam/verify.err.txt"
with open(err, "r") as fh:
	for line in fh:
		line = line.strip()
		if len(line) > 0:
			fn = line.split("/")[8].split(":")[0]

			if not fn in failed:
				failed[fn] = True






wget = "/mnt/neuro-genomic-1-ro/gsam/RNA/Liege_part_2/wget.txt"

idx = {}

with open(wget, "r") as fh:
	for line in fh:
		line = line.strip().lstrip('#')
		if len(line) > 0:
			line = line.split("files=")
			if len(line) > 1:
				fn = line[1].rstrip("'")
				sid = fn.split("-")[0] # sample id
				srun = fn[:-16] # sample run
				runid = fn.split("_")[2]
				pair = fn.split("_")[-2]

				# print(fn)
				# print(sid)
				# print(srun)
				# print(runid)
				# print(pair)
				# print("")
				
				if sid not in idx:
					idx[sid] = {'__failed' : False}

				if srun not in idx[sid]:
					idx[sid][srun] = {'R1':None,'R2':None}
				
				idx[sid][srun][pair] = (fn, runid)

				if fn in failed:
					idx[sid]['__failed'] = True
				#if srun == 'BHM33WDSXX':
				#	idx[sid]['__failed'] = True


for sid in idx:
	if idx[sid]['__failed']:
		print(sid + "   [FAILED SAMPLES]")
	else:
		print(sid + "   [SAMPLES COMPLETE]")
		
	for srun in [_ for _ in sorted(idx[sid]) if _ != "__failed"]:
		print("  " + srun)
		for pair in sorted(idx[sid][srun]):
			print("    " + pair + ": " + idx[sid][srun][pair][0] + " [" + idx[sid][srun][pair][1] + "]")

	print("")


#for _ in failed.keys():
#	print(_)



# # make first alignment code

# for sid in idx:
	# if not idx[sid]['__failed']:
		# print("print '" + sid + "   [SAMPLES COMPLETE]'")

		# # check if all cleaned files exist
		# for srun in [_ for _ in sorted(idx[sid]) if _ != "__failed"]:
			# #print("  " + srun)
			# for pair in sorted(idx[sid][srun]):
				# #AFA1-2545_NGS19-J355_AHKK5WDSXX_S143_L004_R1_001.fastq.gz
				# fn_clean = "/mnt/galaxian/data/users/youri/neurogen-tmp/align-gsam/fastq-clean/" # AAB1-1917_NGS19-J295_AHKK5WDSXX_S35_L002_001.clean_yx_R1.fastq.gz
				# fn_clean += idx[sid][srun][pair][0]
				# fn_clean = fn_clean.replace("_R1_001.fastq.gz","_001.clean_yx_R1.fastq.gz")
				# fn_clean = fn_clean.replace("_R2_001.fastq.gz","_001.clean_yx_R2.fastq.gz")
				# print("  " + fn_clean)
		# print('')

	# # for srun in [_ for _ in sorted(idx[sid]) if _ != "__failed"]:
		# # print("  " + srun)
		# # for pair in sorted(idx[sid][srun]):
			# # print("    " + pair + ": " + idx[sid][srun][pair][0] + " [" + idx[sid][srun][pair][1] + "]")

	# # print("")

