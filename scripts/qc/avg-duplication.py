#!/usr/bin/env python

path = "data/RNA/alignments/"

print("sample id\tpercentage-duplicated-reads\taverage-duplication-per-read = 1 / ( 1 - d )")

import os
for _ in os.listdir(path):
    if _[-15:] == "dup_metrics.txt":
        #print(path + _)
        with open(path + _) as fh:
            i = 0
            for line in fh:
                line = line.strip()
                if len(line) > 0 and line[0] != '#':
                    if i == 1:
                        line = line.split("\t")
                        ratio_dupl = float(line[8])
                        fact_dupl = 1.0 / (1.0 - ratio_dupl)
                        print(_.split("_hg19.")[0] + "\t" + str(ratio_dupl) + "\t" + str(round(fact_dupl, 3)) )
                    i += 1
