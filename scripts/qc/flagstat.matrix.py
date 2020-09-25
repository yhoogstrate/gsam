#!/usr/bin/env python

import os
from tqdm import tqdm


samples = []
idx = {}


p = 'output/tables/qc/flagstat'
for _ in tqdm(os.listdir(p)):
    fn = p + "/" + _
    if _[-4:] == ".txt" and _.find("matrix") == -1:
        samples.append(_)

        with open(fn, 'r') as fh_in:
            for line in fh_in:
                line = line.strip().split(" + 0 ")
                
                if line[1].find("singletons") != -1:
                    line[1] = "singletons"
                elif line[1].find("properly paired") != -1:
                    line[1] = "properly paired"
                elif line[1][0:6] == "mapped" != -1:
                    line[1] = "mapped"
            
                if line[1] not in idx:
                    idx[line[1]] = {}

                if _ in idx[line[1]]:
                    print("Error - duplicate sample entry: " + _)
                    import sys
                    sys.exit(1)
                else:
                    idx[line[1]][_] = line[0]
            

samples = sorted(samples)
keys = sorted(idx.keys())
with open("output/tables/qc/flagstat/flagstats.matrix.txt", "w") as fh_out:
    header = [ "sid" ] + keys
    fh_out.write("\t".join(header) + "\n")

    for sample in samples:
        line_out = [sample]
        for key in keys:
            line_out.append(idx[key][sample])

        fh_out.write("\t".join(line_out) + "\n")



