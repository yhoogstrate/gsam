#!/usr/bin/env python

import os
from tqdm import tqdm

with open("output/tables/duplicate_reads_sambamba_stats.txt", "w") as fh_w:
    fh_w.write("sample-id\tpercentage-duplicate-reads\n")

    path = "data/RNA/alignments/"
    for fn in tqdm([_ for _ in os.listdir(path) if _[-23:] == '.marked_dup_metrics.txt']):
        p = os.path.join(path, fn)
        #print(p)
        with open(p, 'r') as fh:
            std = None
            dups = None
            for line in fh:
                line = line.strip()
                if len(line) > 0 and line[0] != '[':
                    if line.find("sorted ") != -1:
                        std = line.strip().split(" ")[1]
                    elif line.find("found ") != -1:
                        dups = line.strip().split(" ")[1]
            
            #print(fn.split("_")[1] + "\t" +  str(round(100.0 * (0.5 * float(dups)) / float(std), 2)) + "%")
            fh_w.write(fn.split("_")[1] + "\t" +  str(round(100.0 * (0.5 * float(dups)) / float(std), 2)) + "\n")
