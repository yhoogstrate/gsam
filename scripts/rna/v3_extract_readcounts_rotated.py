#!/usr/bin/env python

"""
Rotate "output/tables/v3_extract_readcounts.txt" into "output/tables/v3_extract_readcounts_rotated.txt"
"""

idx = {}

import os
from tqdm import tqdm

with open("output/tables/v3_extract_readcounts_rotated.txt", "w") as fh_out:
    fh_out.write("sid\t\tresection-1-wt\tresection-1-v3\t\tresection-2-wt\tresection-2-v3\n")

    with open("output/tables/v3_extract_readcounts.txt") as fh:
        for line in tqdm(fh):
            line = line.strip()
            if len(line) > 0:
                line = line.split("\t")
                if line[0] != "sample":
                    sid = line[0].split("alignments/")[1].split("/Aligned.s")[0]
                    resection = sid[3]
                    sid = sid[0:3]

                    if sid not in idx:
                        idx[sid] = {'1': ["NA", "NA"], '2': ["NA", "NA"]}

                    idx[sid][resection] = line[1:]


    for sid in sorted(idx.keys()):
        fh_out.write(sid + "\t|\t" + idx[sid]["1"][0] + "\t" + idx[sid]["1"][1] + "\t|\t" +  idx[sid]["2"][0] + "\t" + idx[sid]["2"][1] + "\n")

