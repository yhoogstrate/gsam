#!/usr/bin/env python


from tqdm import tqdm
import os

# ---- find files ----

keys = {}
path = 'data/DNA/EncryptedZippedBAMs/PDexport/Gpgkey'
for _ in tqdm(os.listdir(path)):
  sid = _.replace('.gpgkey','')
  
  fn = os.path.abspath(path + '/' + _)
  
  with open(fn, 'r') as fh:
    c = fh.read()
    key = c.strip().split()[1]
  
  keys[sid] = {'keyfile': fn, 'key': key, 'archive': None}


#print(keys)

# ---- find matching archive -----



path = "/home/yhoogstrate/projects/gsam/data/DNA/dna_data_2020/request_962/PDexport"
for _ in os.listdir(path):
  if _[-4:] == ".gpg":
    sid = _.replace('.tar.gz.gpg','')
    #print(sid)
    
    if sid not in keys:
      keys[sid] = {'keyfile': None, 'key': None, 'archive': path + "/" + _}
    else:
      keys[sid]['archive'] = path + "/" + _


#import sys
#sys.exit(0)

# ---- print script ----
# out path
path = "/home/yhoogstrate/mnt/neuro-genomic-rw/gsam"
for sid in sorted(keys):
  
  # gpg -o 1736_PD29224a.tar.gz --decrypt --passphrase 167343F8F47511E7A27AD422B519FF62 --batch --yes 1736_PD29224a.tar.gz.gpg
  if(keys[sid]['keyfile'] != None):
    args =  ['gpg',
      '-o',
      path + "/" + sid + '.tar.gz',
      '--decrypt',
      '--passphrase',
      keys[sid]['key'],
      '--batch',
      '--yes',
      keys[sid]['archive']]
      
    print(" ".join(
     args
      ))
  else:
    print('# missing key: ' + sid + " -> " + keys[sid]['archive'])


