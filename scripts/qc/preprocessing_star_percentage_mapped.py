#!/usr/bin/env python


import os
from tqdm import tqdm


with open("output/tables/star_percentage_mapped.txt", "w") as fh_w:
    fh_w.write("sample\tpercentage-uniquely-mapped\tpercentage-multimap\n")

    path = "data/RNA/alignments"
    for directory in tqdm([_ for _ in os.listdir(path)]):
        #print(directory)
        sid = directory
        #.split("_")[1]
        #print(os.listdir(  os.path.join(path, directory, "Log.final.out")  ))
        
        t = None
        um = None
        mm = None

        fn = os.path.join(path, directory, "Log.final.out")
        with open(fn, 'r') as fh:

            for line in fh:
                line = line.strip()
                
                if line.find("Number of input reads") != -1:
                    t =  int(line.split("|")[1].strip())
                elif line.find("Uniquely mapped reads number") != -1:
                    um = int(line.split("|")[1].strip())
                elif line.find("Number of reads mapped to multiple loci") != -1:
                    mm = int(line.split("|")[1].strip())

        fh_w.write(sid + "\t" + str(round(100.0 * um/t,2)) + "%\t" + str(round(100.0 * mm/t, 2)) + "%\n")
