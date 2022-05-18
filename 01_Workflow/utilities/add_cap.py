#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 14:27:18 2021

@author: huy
"""

import numpy as np
import os
import glob

#%% Load full protein - we grab the C of previous residue for ACE
# Find the full protein pdb
full_pdb = glob.glob('rall*.pdb')[0]

with open(full_pdb,'r') as f:
    full_protein = f.readlines()

C_line = {}
O_line = {}
CA_line = {}

for idx, line in enumerate(full_protein):
    if ' C  ' in line:
        C_line[int(line[22:26])] = idx
    if ' CA ' in line:
        CA_line[int(line[22:26])] = idx
    if ' O  ' in line:
        O_line[int(line[22:26])] = idx
        
#print(C_line)
#print(CA_line)
#print(O_line)
    
#%% List all residue segments
listing = glob.glob('r*-*-*-b.pdb') # r00-54[t]-58[t]-b.pdb
listing.sort()
print(listing)

#%%
for l in listing:
                 # r00-54[t]-58[t]-    c.pdb
    new_f_name = l.split('.')[0][:-1]+'c.pdb'
    print(new_f_name)
    first = False
    last = False
    # parse file name
    caps = l.split('-')
    if 't' in caps[1]:
        first = True
    if 't' in caps[2]:
        last = True

    with open(l,'r') as f:
        seg_protein = f.readlines()
    # Replace OXT XXX as N   NME
    for idx, line in enumerate(seg_protein):
        if 'OXT' in line:
            seg_protein[idx] = line[:13] + 'N   NME' + line[20:]
            last_res = line[22:26]
        if not first:
            if 'ATOM      1' in line:
                seg_protein[idx] = 'ATOM      4' + line[11:]

    for i in [3,2,1]: # Remove hydrogens:
        if (' H ' in seg_protein[i]) or (' H1' in seg_protein[i]) or (' H2' in seg_protein[i]) or (' H3' in seg_protein[i]): 
            seg_protein.remove(seg_protein[i])
    # for line in seg_protein:
    #     print(line.strip('\n'))
            
    if not first:
        # Add atom 0 
        first_res = int(l.split('-')[1].replace('t',''))-1
        print(first_res)
        print(full_protein[C_line[first_res]].strip('\n'))
        first_line =  'ATOM      1  CH3 ACE     1' + full_protein[CA_line[first_res]][26:]
        second_line = 'ATOM      2  C   ACE     1' + full_protein[C_line[first_res]][26:]
        third_line =  'ATOM      3  O   ACE     1' + full_protein[O_line[first_res]][26:]
        print(first_line)
        seg_protein.insert(0, third_line)
        seg_protein.insert(0, second_line)
        seg_protein.insert(0, first_line)
    if not last:
        # Take the CA coordinate from next res
        next_res = int(l.split('-')[2].replace('t',''))+1
        next_line =  'ATOM      0  CH3 NME  ' + last_res + full_protein[CA_line[next_res]][26:]
        seg_protein.insert((len(seg_protein)-2), next_line)
    with open(new_f_name,'w') as f:
        f.writelines(seg_protein)
#for line in seg_protein:
#    print(line.strip('\n'))
    # with open()
        
