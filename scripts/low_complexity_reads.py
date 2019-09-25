#!/usr/bin/env python

import os
from tqdm import tqdm
import json


path = 'data/RNA/alignments_log'


with open("output/tables/low_complexity_reads.txt", "w") as fh_w:
    fh_w.write("sample\tpercentage low complexity reads\n")

    for fn in tqdm([_ for _ in os.listdir(path) if _[-5:] == ".json"]):
        with open(os.path.join(path, fn), 'r') as fh:
            c = fh.read()
            j = json.loads(c)
            
            lc = j["filtering_result"]["low_complexity_reads"]
            total = j["summary"]["before_filtering"]["total_reads"]
            lc_p = round( 100.0 * lc / total , 2)


            #print(str(lc) + "/" + str(total))
            fh_w.write(fn.split("_")[1] + "\t" + str(lc_p) + "%"+"\n")

