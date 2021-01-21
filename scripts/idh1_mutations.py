#!/usr/bin/env python

#path = 'data/DNA/snv_vcf'


from tqdm import tqdm
import os
import gzip


import glob


#for f in [_ for _ in os.listdir(path) if _[-7:] == ".vcf.gz"]:
    #print(f)
    #fn = os.path.join(path, f)
for fn in tqdm(glob.glob('data/gsam/DNA/dna_data_2020/*/*/decrypted/*/*/*cavem*.vcf.gz')):
    #print(fn)
    with gzip.open(fn, 'rb') as fh:

        idh = False
        sid = fn.split("_")[-2].split(".")[0]

        for line in fh:
            #if line[0:10] == "2\t20911311":# and  line.find("R132") != -1:
            if (line.find('IDH1') != -1 and line.find('R132') != -1):
                # or (line.find('IDH2') != -1 and (line.find('172') != -1 or line.find('140') != -1 ) ):
                mut = line.split(';VW=')[1].split("|")
                mut = mut[0] + mut[4][1:]
                
                vaf = str(round(100.0 * float(line.split("\t")[-1].split(":")[-1]))) + "%"
                
                coverage = line.split("DP=")[1].split(";")[0]
                
                print(sid + "\t" + mut + "\t" + line.split("\t")[6] + '\t' +  vaf + "\t" + coverage)
                idh = True

        if not idh:
            print(sid + "\t-\t-\tNA\tNA")

# there were no idh2's
#            elif line.find("IDH2") != -1:
#                if(line.find )

"""

IDH's from 313 VCF files


PD29173a ** large table     
PD29173c ** large table     
PD29199a2 ** large table    1559*
PD29199c2 ** large table    1559*
PD29216a2 ** large table    1559*
PD29216c2 ** large table    1559*



PD29228a	IDH1.R132G	PASS	38.0%	383
PD29228a ** large table

PD29228c	IDH1.R132G	PASS	44.0%	576
PD29228c ** large table

PD29263a	IDH1.R132H	PASS	32.0%	427
PD29263a ** large table

PD29264a2	IDH1.R132H	PASS	30.0%	425
PD29264a2 ** large table

PD29264c	IDH1.R132H	MNP  	12.0%	445

PD30239c	IDH1.R132H	PASS	47.0%	321
PD30239c ** large table

PD30242a3	IDH1.R132H	PASS	31.0%	410
PD30242a3 ** large table

PD30242c3	IDH1.R132H	MNP  	15.0%	417

PD36768a	IDH1.R132H	PASS	41.0%	298
PD36768a ** large table

PD36768c	IDH1.R132H	PASS	25.0%	321
PD36768c ** large table

PD36770a	IDH1.R132H	PASS	41.0%	384
PD36770a ** large table

PD36770c	IDH1.R132H	PASS	27.0%	414
PD36770c ** large table

PD36772c	IDH1.R132H	MNP  	15.0%	391

PD36783a	IDH1.R132H	PASS	41.0%	343
PD36783a ** large table

PD36783c	IDH1.R132H	PASS	37.0%	310
PD36783c ** large table

"""

