#!/usr/bin/env python

import os
from tqdm import tqdm
import json
import re

path = 'data/RNA/fastp_log'
splitregex = re.compile("[0-9]{3,}_NGS19")

with open("output/tables/fastp_duplication_statistics.txt", "w") as fh_w:
    fh_w.write("sample-id\tfilename\tinput-reads\tduplication-rate\tavg-duplication-per-read\n")
    
    for fn in tqdm(sorted([_ for _ in os.listdir(path) if _[-5:] == ".json"])):
        sid = splitregex.split(fn)[0].rstrip('-')
        #print(sid)
        with open(os.path.join(path, fn) , "r") as fh:
            j = json.loads(fh.read())

            #a = sum(j['read1_after_filtering']['content_curves']['A']) +  sum(j['read2_after_filtering']['content_curves']['A'])
            #print(j['summary']['before_filtering']['total_reads'])
            #print(j['duplication']['rate'])
            #
            #print('---------------\n')

            adpr = round(1.0 / (1.0 - float(j['duplication']['rate'])) , 4)
            fh_w.write(sid + "\t" + fn + "\t"+str(j['summary']['before_filtering']['total_reads'])+"\t"+str(j['duplication']['rate'])+"\t"+str(adpr)+"\n")

