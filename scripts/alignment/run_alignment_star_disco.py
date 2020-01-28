#!/usr/bin/env python

# 1. find all cleaned fastq pairs
# 2. merge sorted based on patient identifier
# 3. generate star run script

import os
import re

idx = {}

with open("scripts/alignment/run_alignment_star_disco.sh", "w") as fh:
	fh.write("#!/bin/bash\n\n")
	fh.write("ulimit -n 35000 \n\n")
	
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
		
		if not srid in idx[pid][sid]:
			idx[pid][sid][srid] = {'R1': None, 'R2': None}

		idx[pid][sid][srid][mate] = _
		


	for pid in sorted(idx):
		print(pid)
		for sid in sorted(idx[pid]):
			print("  " + sid + ":")
			outdir = "output/bam-disco/" + sid
			
			fh.write('mkdir -p "' + outdir + '";\n')
			fh.write('nice ./bin/STAR/bin/Linux_x86_64/STAR \\\n')
			fh.write('    --outFileNamePrefix "' + outdir + '/" \\\n')
			fh.write('    --runThreadN 54 \\\n')
			fh.write('    --genomeDir ref/star-hg19 \\\n')
			fh.write('    --sjdbGTFfile ref/star-hg19/gencode.v31lift37.annotation.gtf \\\n')
			fh.write('    --readFilesIn \\\n')
			r1 = []
			r2 = []
			
			for srid in sorted(idx[pid][sid]):
				if idx[pid][sid][srid]['R1'] is not None and idx[pid][sid][srid]['R2'] is not None:
					r1.append('data2/fastq-clean/' + idx[pid][sid][srid]['R1'])
					r2.append('data2/fastq-clean/' + idx[pid][sid][srid]['R2'])
					#print("    " + srid + "\t\tCOMPLETE")
				else:
					print("    " + srid + "\t\tINCOMPLETE !!!!!!!!!!!!!!!")
					raise Exception("???")

			fh.write('        ' + ','.join(r1) + ' \\\n')
			fh.write('        ' + ','.join(r2) + ' \\\n')
			fh.write('    --readFilesCommand zcat \\\n')
			fh.write('    --outSAMtype BAM SortedByCoordinate \\\n')
			fh.write('    --alignSJoverhangMin 10 \\\n')
			fh.write('    --alignSJDBoverhangMin 1 \\\n')
			fh.write('    --alignIntronMin 20 \\\n')
			fh.write('    --alignMatesGapMax 10000 \\\n')
			fh.write('    --chimSegmentMin 12 \\\n')
			fh.write('    --chimJunctionOverhangMin 12 \\\n')
			fh.write('    --chimOutType WithinBAM SeparateSAMold \\\n')
			fh.write('    --outSAMstrandField intronMotif \\\n')
			fh.write('    --outFilterIntronMotifs None\n')
			fh.write('rm "'+outdir+'/Aligned.sortedByCoord.out.bam"\n')
			fh.write('\n\n')
		print("")
		fh.write('\n# - - - - -\n\n')

	"""

	mkdir -p "bam/Pfrench_EBU1_2627_Plate1_G8_S117_001_hg19/"
	nice ./bin/STAR/bin/Linux_x86_64/STAR \
		--runThreadN 84 \
		--genomeDir /data/users/youri/neurogen-tmp/align-gsam/star-hg19 \
		--sjdbGTFfile ref/star-hg19/star-hg19/gencode.v31lift37.annotation.gtf \
		--readFilesIn /dev/shm/Pfrench_EBU1_2627_Plate1_G8_S117_001.dedup.clean_R1.fastq.gz \
					  /dev/shm/Pfrench_EBU1_2627_Plate1_G8_S117_001.dedup.clean_R2.fastq.gz \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix "bam/Pfrench_EBU1_2627_Plate1_G8_S117_001_hg19/" \
		--alignSJoverhangMin 10 \
		--alignSJDBoverhangMin 1 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--chimSegmentMin 12 \
		--chimJunctionOverhangMin 12 \
		--chimOutType WithinBAM SeparateSAMold \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--quantMode TranscriptomeSAM



	With 

		--quantMode TranscriptomeSAM

	GeneCounts, and get both the Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab outputs.
	"""
