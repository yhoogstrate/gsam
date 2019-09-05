#!/usr/bin/env python

idx = {}

import os

print("sid\t\tresection-1-wt\tresection-1-v3\t\tresection-2-wt\tresection-2-v3")

path = "/mnt/galaxian/data/users/youri/neurogen-tmp/align-gsam/v3/"
for _ in os.listdir(path):
    with open(path + _ , "r") as fh:
        for line in fh:
            line = line.strip()
            if len(line) != "":
                line = line.split("\t")
                if line[0] != "sample":
                    sid = _.split("_")
                    resection = sid[1][3]
                    sid = sid[1][0:3]
                    
                    if sid not in idx:
                        idx[sid] = {'1': ["NA", "NA"], '2': ["NA", "NA"]}

                    idx[sid][resection] = line[1:]


for sid in sorted(idx.keys()):
    print (sid + "\t|\t" + idx[sid]["1"][0] + "\t" + idx[sid]["1"][1] + "\t|\t" +  idx[sid]["2"][0] + "\t" + idx[sid]["2"][1])

#print(ids[1][0:3] + "\t" + ids[1][3] + "\t" + "\t".join(line[1:]))
