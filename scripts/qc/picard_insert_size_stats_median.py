#!/usr/bin/env python

import glob


idx = {}


for fn in glob.glob('output/tables/picard_insertsizemetrics/*.txt'):
    sid = fn.split("/")[-1].replace('.txt','')
    #print(sid)
    v = []

    with open(fn, 'r') as fh:
        for line in fh:
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                line = line.split("\t")
                if len(line) > 2 and line[0] != 'MEDIAN_INSERT_SIZE':
                    #print(line)
                    idx[sid] = int(line[0])


#print(idx)



fh = open("output/tables/picard_insertsizemetrics_combined.txt", "w")
fh.write("sample\tMEDIAN_INSERT_SIZE\n")
for k in sorted(idx.keys()):
    #print(k)
    fh.write(k + "\t" + str(idx[k]) + "\n")
fh.close()






