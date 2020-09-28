#!/usr/bin/env python

import os
from tqdm import tqdm
import json
import re


path = 'data/RNA/output/tables/fastp_log'


with open("output/tables/qc/low_complexity_reads.txt", "w") as fh_w:
    fh_w.write("sample\tpercentage low complexity reads\n")

    for fn in tqdm([_ for _ in os.listdir(path) if _[-5:] == ".json"]):
        if fn.find("NGS19") != -1:
            sid = fn.split("_NGS19")[0]
            sid = re.sub("\\-[0-9]{3,}","", sid)
        elif fn.find("NGS20") != -1:
            sid = fn.split("_NGS20")[0].split("-", 1)[1]+"-new"
        else:
            print(fn)
            import sys
            sys.exit(1)

        print(sid)

        with open(os.path.join(path, fn), 'r') as fh:
            c = fh.read()
            j = json.loads(c)
            
            lc = j["filtering_result"]["low_complexity_reads"]
            total = j["summary"]["before_filtering"]["total_reads"]
            lc_p = round( 100.0 * lc / total , 2)


            #print(str(lc) + "/" + str(total))
            fh_w.write(sid + "\t" + str(lc_p) + "%"+"\n")


