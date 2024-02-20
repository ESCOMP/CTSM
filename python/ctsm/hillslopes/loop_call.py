#!/usr/bin/env python
# coding: utf-8

# Import modules

import sys 
import string 
import subprocess
import copy
import argparse

# ------------------------------------------------------- #
#
#  Loop over script calls when chunking data
# 
# ------------------------------------------------------- #
parser = argparse.ArgumentParser(description='Loop script call')
parser.add_argument("snum", help="script flag", nargs='?',type=int,default=0)
parser.add_argument("-o","--overwrite", help="overwrite", action="store_true",default=False)
args = parser.parse_args()

snum = args.snum
if snum not in [1]:
    print('only 1 script currently supported')
    print('snum must be 1')
    stop

# Combine individual gridcell files into chunks
if snum == 1:
    script_to_call = './combine_gridcell_files.py'


# totalChunks must match value used in scripts
totalChunks = 36
    
istart,iend = 1,totalChunks+1
for i in range(istart,iend):
    print('===============================================')
    print('Beginning chunk {:d}'.format(i),flush=True)
    
    if args.overwrite:
        command=[script_to_call, str(i),'-o']
    else:
        command=[script_to_call, str(i)]

    try:
        x=subprocess.call(command,stderr=subprocess.PIPE)
    except:
        continue
