#!/usr/bin/env python

import os
from tqdm import tqdm
import json

path = 'data/RNA/alignments_log'

with open("output/tables/gc_content_rmse.txt", "w") as fh_w:
    fh_w.write("sample-id\tpercentage.A\tpercentage.C\tpercentage.T\tpercentage.G\tRMSE\n")
    
    for fn in tqdm(sorted([_ for _ in os.listdir(path) if _[-5:] == ".json"])):
        #print( os.path.join(path, fn))
        with open(os.path.join(path, fn) , "r") as fh:
            j = json.loads(fh.read())

            
            a = sum(j['read1_after_filtering']['content_curves']['A']) +  sum(j['read2_after_filtering']['content_curves']['A'])
            c = sum(j['read1_after_filtering']['content_curves']['C']) +  sum(j['read2_after_filtering']['content_curves']['C'])
            t = sum(j['read1_after_filtering']['content_curves']['T']) +  sum(j['read2_after_filtering']['content_curves']['T'])
            g = sum(j['read1_after_filtering']['content_curves']['G']) +  sum(j['read2_after_filtering']['content_curves']['G'])
            n = 1.0 * (a + c + t + g)

            pa = 100.0 * a / n
            pc = 100.0 * c / n
            pt = 100.0 * t / n
            pg = 100.0 * g / n

            ra = pow(25.0 - pa, 2)
            rc = pow(25.0 - pc, 2)
            rt = pow(25.0 - pt, 2)
            rg = pow(25.0 - pg, 2)

            rmse = (ra + rc + rt + rg) / 4.0
            rmse = pow(rmse, 0.5)

            fh_w.write(fn.split("_")[1] + "\t"+str(round(pa,2))+"\t"+str(round(pc,2))+"\t"+str(round(pt,2))+"\t"+str(round(pg,2))+"\t"+str(round(rmse, 3))+"\n")

