#!/usr/bin/env python

import os
from tqdm import tqdm


encrypted_archives = {}
encryption_keys = {}

p = 'data/DNA/dna_data_2020/request_962/PDexport'
for _ in tqdm(os.listdir(p)):
    if _[-7:] == ".gz.gpg":
        fn = os.path.abspath(p + "/" + _)
        sid = _.replace('.tar.gz.gpg','')
        
        encrypted_archives[sid] = fn


p = 'data/DNA/dna_data_2020'
for _ in sorted(tqdm(os.listdir(p))):
    if _[-7:] == ".gpgkey":
        fn = os.path.abspath(p + "/" + _)
        sid = _.replace('.gpgkey','')
        
        with open(fn, 'r') as fh:
            encryption_keys[sid] = fh.read().strip().split()[1]
        
        


a = encrypted_archives.keys()
b = encryption_keys.keys()
c = set(a).intersection(set(b))

if len(a) != len(c):
    import sys
    print("error key / archive mismatch")
    sys.exit(1)
    



for _ in sorted(tqdm(encrypted_archives)):
    
    print('gpg -o \'/mnt/neuro-genomic-1-rw/gsam/'+sid+'.tar.gz\' --decrypt --passphrase '+encryption_keys[_]+' --batch --yes \''+encrypted_archives[_]+'\'')







