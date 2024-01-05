#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a list of representative transcripts per gene
"""
import pandas as pd
import sys
import subprocess
#braker_filename="/home/harry/Documents/manchester/global/sandbox/C480.braker.gtf"
#exonerate_filename="/home/harry/Documents/manchester/global/sandbox/C480.AF293.exonerate.gtf"
#collapsed_filename="/home/harry/Documents/manchester/global/sandbox/C480.braker-collapsed.gtf"
#clean_filename="/home/harry/Documents/manchester/global/sandbox/C480.braker-clean.gtf"

braker_filename=sys.argv[1]
exonerate_filename=sys.argv[2]
collapsed_filename=sys.argv[3]
clean_filename=sys.argv[4]


braker=pd.read_csv(braker_filename, sep='\t', header=None)
gene_id=braker[braker.iloc[:,2]=="gene"].iloc[:,8]
# Collapse and store braker results
to_keep=[]
for gid in gene_id:
    subset=braker[(braker.iloc[:,8].str.contains(gid + ".t*"))
                  &(braker.iloc[:,2]=="transcript")]
    best_id=subset.index[subset.iloc[:,5]==subset[5].max()].tolist()[0]
    to_keep.append(best_id)
    
collapsed=braker.iloc[to_keep]
collapsed.to_csv(collapsed_filename, sep='\t', header=None, index=False)

# Remove interesecting results using BEDtools
bedtools_cline="intersectBed -v -a " + collapsed_filename + " -b " + exonerate_filename
bedproc=subprocess.run(bedtools_cline, capture_output=True, shell=True)
bedres=bedproc.stdout.decode('utf-8')


output = open(clean_filename, "w")
output.write(bedres)
output.close()


