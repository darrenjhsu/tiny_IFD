
import mdtraj
import os, sys, glob
import numpy as np

case_name = sys.argv[1]
ref_file = sys.argv[2]

ref = mdtraj.load(f'../{case_name}/{ref_file}')
print(f'{case_name:16s}', end=' ')
num_docked = len(glob.glob(f'../{case_name}/docked_pdb/*.pdb'))
for i in range(1, num_docked+1):
    try:
        mol = mdtraj.load(f'../{case_name}/docked_pdb/{case_name}_{i}.pdb')
        rmsd = np.sqrt(np.mean(((ref.xyz[0]*10 - mol.xyz[0]*10)**2).sum(1)))
        print(f'{rmsd:5.2f}', end=' ')
    except:
        print(f'99999', end=' ')
print() 
