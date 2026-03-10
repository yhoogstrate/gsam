#!/usr/bin/env python

import subprocess
import os


def get_trailing_end_pos(fn):
   cmd = ['gunzip-split', '--list-only', fn]
   
   p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
   out, err = p.communicate()
   err = err.decode().strip()
   
   if len(err) == 0:
      return None
   else:
      txt = err.split("contains non-gzip data starting at byte")[1].strip("\r\n ./")
      return int(txt)



def find_new_data(fn, offset):
   with open(fn, "rb") as fh:
      fh.seek(offset - 1)
      
      rc = chr(0)
      while ord(rc) == 0:
         rc = fh.read(1)
      
      new_data = fh.tell()
   
      return new_data - 1



def split_file(fn_in, fn_out, local_end, new_start):
   print("Removing trailing data: " + str(local_end) + " - " + str(new_start))
   with open(fn_in, "rb") as fh_in:
      with open(fn_out, "wb") as fh_out:
         fh_out.write(fh_in.read(local_end))
         
         fh_in.seek(new_start)
         fh_out.write(fh_in.read())




def fix_gz_members(fn_old, fn_new):
   split = get_trailing_end_pos(fn_old)

   while split:
      new_start = find_new_data(fn_old, split)
      split_file(fn_old, fn_new, split, new_start)
      
      # make old new the new old
      os.rename(fn_new, fn_new + ".tmp")
      fn_old = fn_new + ".tmp"
      
      split = get_trailing_end_pos(fn_old)
   
   
   os.rename(fn_old, fn_new)




fix_gz_members("/home/youri/mnt/neuro-genomic-1-ro/gsam/RNA/fastq/Liege_part_2/ABA1-2544_NGS19-J611_BHM33WDSXX_S63_L003_R1_001.fastq.gz","/tmp/test.gz")




