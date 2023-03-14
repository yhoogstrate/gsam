#!/usr/bin/env python


import re

sids = ["AAB1","AAC1","AAD1","AAF1","AAF2","AAJ1","AAJ2","AAL1","AAM1","AAM2","AAN1","AAN2","AAP1-new","AAS1","AAS2","AAT1","AAT2","AAU1",
        "AAU2","AAV1","AAV2","AAW1","AAW2","AAX1","AAX2","AAY1","AAY2","ABA2","ACA1","ACA2","ADA1","ADA2","AFA1","AFA2","AHA1","AHA2","AIA1",
        "AIA2","AKA2","ALA1","ALA2","AMA1","AMA2","AOA1","AOA2","AQA1","AQA2","AZA1","AZA2","AZB1","AZB2","AZC1","AZC2","AZD1","AZD2","AZE1-new",
        "AZE2","AZF1","AZF2","AZG1","AZG2","AZH1","AZH2-new","BAA1","BAB1","BAB2","BAC1","BAC2","BAD1","BAD2","BAE1","BAE2","BAH1","BAH2","BAI1",
        "BAI2-new","BAJ1","BAJ2","BAK1","BAK2","BAL1","BAM1","BAM2","BAO1","BAR1","BAS1","BAS2","BAT1","BAT2-new","BAW1","BAW2","BAY1","BAY2-new",
        "CAC1-new","CAD1","CAD2","CAF2","CAO1-new","CAO2","CBA1","CBE2","CBH1","CBI1","CBI2","CBM1","CBM2","CBP1","CBP2","CBQ1-new","CBQ2","CBR1",
        "CBR2","CBS1-new","CBS2","CBT1","CBT2","CBV1","CCD1","CCD2","CCW1","CCW2","CCZ1-new","CCZ2","CDA1","CDA2","CDD1","CDD2","CDH1-new","DAC1",
        "DAC2","DAD1","DAD2","EAC1","EAC2","EAD1","EAD2","EAE1","EAE2","EAG1","EAG2","EAI1","EAI2","EAJ1","EAJ2","EAK1","EAK2","EAL1","EAN2","EAO1",
        "EAO2","EAP1","EAP2","EAQ1","EAQ2","EAT1","EAT2","EAU1","EAU2","EAV1","EAW1","EAY1","EAZ1","EAZ2","EBB1","EBB2","EBC1","EBC2","EBF1","EBF2",
        "EBG1","EBG2","EBH1","EBH2","EBK1","EBK2","EBL1","EBM2","EBO2","EBP1","EBR1","EBR2","EBU1","EBU2","EBV1","EBW2","EBX1","EBY1","EBY2","ECA1",
        "ECA2","ECD1","ECD2","ECE1","ECE2","ECG1","ECG2","ECH1","ECI1","ECI2","ECK1","ECK2","ECN1","ECN2","FAB1","FAB2-replicate","FAF1","FAF2",
        "FAG1","FAG2","FAH2","FAI1","FAI2","FAJ1","FAJ2","FAK1","FAK2","FAM1","FAN1","FAN2","FAP1","FAP2","FAQ1","FAQ2","GAA1","GAA2","GAE1-new",
        "GAE2","GAG1","GAG2","GAH1","GAH2","GAI1","GAI2","GAJ1","GAK1","GAK2","GAL1","GAM1-new","GAM2","GAN1","GAO1","GAO2-new","GAP1","GAP2","GAQ1",
        "GAQ2","GAR1","GAR2","GAS1-new","GAS2","HAA1","HAA2-new","HAB1","HAB2","HAC2","HAD1","HAE1-new","HAE2","HAF1-new","HAF2-new","HAG2","HAH1",
        "JAC1","JAC2","JAE1","JAE2","JAG1","JAG2","JAH1","JAH2","JAI1","JAI2","JAJ1","JAJ2","JAK1","JAK2","JAL1","JAL2","JAM1","JAM2","JAN1-new",
        "JAN2","KAA1-new","KAA2","KAB1","KAD1","KAE1-new"]

import glob

old = glob.glob("data/gsam/RNA/fastq/Liege_part_2/*.fastq.gz")
new = glob.glob("data/gsam/RNA/fastq/Liege_part_3/*.fastq.gz")

assert(len(old) + len(new) == 926)

# 
# idx_sid = {}
# for _ in sids:
#   if _.find("-new") != -1:
#     __ = _.replace("-new","")
#     path = "Liege_part_3"
#   else:
#     __ = _
#     path = "Liege_part_2"
#   
#   if __ in idx_sid:
#     print("Error")
#     print(sid)
#   else:
#     idx_sid[__] = path
#   


idx_old = {}
for _ in old:
  if _.find("_R1_") != -1:
    RR = "R1"
  else:
    RR = "R2"

  sid = _.split("/")[-1].split("_NGS19")[0]
  sid = re.sub(r'-[0-9]+$', '', sid)
  sid = re.sub(r'-repl$', '-replicate', sid)
  #print(sid, RR)

  if sid not in idx_old:
    idx_old[sid] = {'R1': [], 'R2': []}
  
  idx_old[sid][RR].append(_.replace('data/gsam','/home/r361003/mnt/neuro-genomic-1-ro2/gsam'))



idx_new = {}
for _ in new:
  if _.find("_R1_") != -1:
    RR = "R1"
  else:
    RR = "R2"

  sid = _.split("/")[-1].split("_NGS20")[0]
  sid = re.sub(r'^[0-9]+-', '', sid)
  sid = re.sub(r'-repl$', '-replicate', sid)
  print(sid, RR, _)

  if sid not in idx_new:
    idx_new[sid] = {'R1': [], 'R2': []}
  
  idx_new[sid][RR].append(_.replace('data/gsam','/home/r361003/mnt/neuro-genomic-1-ro2/gsam'))


  # for _ in data['R1']:
  #   print("R1: " + _)
  # for _ in data['R2']:
  #   print("R2: " + _)
  # print("")


with open("/tmp/fq-merge.sh","w") as ffh:
  for sid in sids:
    #print(sid)
    data = None
    # old sample; no '-new' found
    if sid.find("-new") == -1: #
      data = idx_old[sid]
    else:
      data = idx_new[sid.replace("-new","")]
    
    #print(data)
    assert(len(data['R1']) == len(data['R2']))
    
    if len(data['R1']) > 1:
      out1 = "cat '"+"' '".join(data['R1'])+ "' > "+"'/mnt/data-ngen0001/G-SAM.EGA/arriba/arriba_v2.1.0/fastq/"+sid+"_R1.fastq.gz'"
      out2 = "cat '"+"' '".join(data['R2'])+ "' > "+"'/mnt/data-ngen0001/G-SAM.EGA/arriba/arriba_v2.1.0/fastq/"+sid+"_R2.fastq.gz'"
      #print("\n")
      ffh.write(out1 + "\n")
      ffh.write(out2 + "\n")
      ffh.write("\n")
      
    else:
      out1 = "ln -s '"+data['R1'][0]+"' '/mnt/data-ngen0001/G-SAM.EGA/arriba/arriba_v2.1.0/fastq/"+sid+"_R1.fastq.gz'"
      out2 = "ln -s '"+data['R2'][0]+"' '/mnt/data-ngen0001/G-SAM.EGA/arriba/arriba_v2.1.0/fastq/"+sid+"_R2.fastq.gz'"
      
      ffh.write(out1 + "\n")
      ffh.write(out2 + "\n")
      ffh.write("\n")




with open("/tmp/arrib.sh","w") as fh:
  for sid in sids:
    cmd1 = [
    'nice',
    './run_arriba.sh',
    'STAR_index_hs37d5viral_GENCODE19/',
    'GENCODE19.gtf',
    'hs37d5viral.fa',
    'database/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz',
    'database/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz',
    'database/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3', '8',
    "'/mnt/data-ngen0001/G-SAM.EGA/arriba/arriba_v2.1.0/fastq/"+sid+"_R1.fastq.gz'",
    "'/mnt/data-ngen0001/G-SAM.EGA/arriba/arriba_v2.1.0/fastq/"+sid+"_R2.fastq.gz'"
    ]
    
    cmd2 = ['mv', "'fusions.tsv'", "'"+ sid + "_fusions.tsv'"]
    cmd3 = ['mv', "'fusions.discarded.tsv'",  "'" + sid + "_fusions.discarded.tsv'"]
    cmd4 = ['gzip', "'" + sid + "_fusions.discarded.tsv'"]
    
    #print(" ".join(cmd1) + " ; " +  " ".join(cmd2) + " ; " + " ".join(cmd3) + " ; " + " ".join(cmd4))
    fh.write(" ".join(cmd1) + " ; " +  " ".join(cmd2) + " ; " + " ".join(cmd3) + " ; " + " ".join(cmd4) + "\n")




