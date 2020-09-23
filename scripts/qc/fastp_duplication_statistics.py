#!/usr/bin/env python

import os
from tqdm import tqdm
import json
import re

path = 'data/RNA/output/tables/fastp_log'
splitregex_19 = re.compile("[0-9]{3,}_NGS19")
#splitregex_20 = re.compile("[0-9]{3,}_NGS20")

with open("output/tables/fastp_duplication_statistics.txt", "w") as fh_w:
    fh_w.write("sample-id\tfilename\tinput-reads\tduplication-rate\tavg-duplication-per-read\tread1_mean_length\tread2_mean_length\n")
    
    
    for fn in tqdm(sorted([_ for _ in os.listdir(path) if _[-5:] == ".json"])):
        if fn.find("NGS19") != -1:
            sid = splitregex_19.split(fn)[0].rstrip('-')
        elif fn.find("NGS20") != -1:
            sid = fn.split("_NGS20")[0].split("-", 1)[1]+"-new"
        else:
            print(fn)
            import sys
            sys.exit(1)
            
        #print(sid)

        with open(os.path.join(path, fn) , "r") as fh:
            j = json.loads(fh.read())

            #a = sum(j['read1_after_filtering']['content_curves']['A']) +  sum(j['read2_after_filtering']['content_curves']['A'])
            #print(j['summary']['before_filtering']['total_reads'])
            #print(j['duplication']['rate'])
            #
            #print('---------------\n')

            adpr = round(1.0 / (1.0 - float(j['duplication']['rate'])) , 4)
            fh_w.write(
                sid + "\t" +
                fn + "\t" +
                str(j['summary']['before_filtering']['total_reads']) + "\t" +
                str(j['duplication']['rate']) + "\t" + 
                str(adpr) + "\t" +
                str(j['summary']['after_filtering']['read1_mean_length']) + "\t" +
                str(j['summary']['after_filtering']['read2_mean_length']) + 
                "\n")

