
# Broken protein and ligand for mdgx
import os
import parmed
r_list = os.listdir()
r_inpcrd = [x for x in r_list if 'c.inpcrd' in x]
r_prmtop = [x for x in r_list if 'c.prmtop' in x]

r_inpcrd.sort()
r_prmtop.sort()

rs = [x.split('.')[0] for x in r_inpcrd]
#print(r_list)
#print(r_inpcrd)
#print(r_prmtop)
print(f'Combining these segments: {rs}')

ps = [] # ps is protein segment
for ii in range(len(rs)):
    if ii == 0:
        ps = parmed.load_file(f'{rs[ii]}.prmtop',f'{rs[ii]}.inpcrd')
    else:
        ps += parmed.load_file(f'{rs[ii]}.prmtop',f'{rs[ii]}.inpcrd')

print(ps)
ps.save(f'apo_continue_nodS.inpcrd', overwrite=True)
ps.save(f'apo_continue_nodS.prmtop', overwrite=True)
ps.save(f'apo_continue_nodS.pdb', overwrite=True)

# Check if disulfide bonds are needed
with open('apo_continue_nodS.pdb', 'r') as f:
    cont = f.readlines()

# Generate fifth_antechamber.in

import numpy as np

SG_resid = []
SG_positions = []
for line in cont:
    if 'SG  CYS' in line:
        SG_resid.append(int(line[22:26]))
        SG_positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
SG_positions = np.array(SG_positions)
try:
    dS = np.where(np.sqrt(np.sum((SG_positions[:,:,None] - SG_positions[:,:,None].T)**2, axis=1))<3)
except:
    dS = ([0], [0])

resid_CYS2CYX = []

with open('fifth_antechamber.in', 'w') as f:
    f.write('## Autogenerated by assemble.py ##\n')
    f.write('source leaprc.protein.ff14SB #Source leaprc file for ff14SB protein force field\n')
    f.write('mol = loadpdb apo_continue_halfdS.pdb\n')
    for Si, Sj in zip(dS[0], dS[1]):
        if Si < Sj: # Only in this situation we make disulfide bonds
            # print(SG[Si], SG[Sj])
            resid_CYS2CYX.append(SG_resid[Si])
            resid_CYS2CYX.append(SG_resid[Sj])
            f.write(f'bond mol.{SG_resid[Si]}.SG mol.{SG_resid[Sj]}.SG\n')
    f.write('savepdb mol apo_continue.pdb\n')
    f.write('saveamberparm mol apo_continue.prmtop apo_continue.inpcrd\n')


# Change apo_continue_nodS.pdb
new_cont = []
for line in cont:
    if 'END' in line:
        new_cont.append(line)
    elif int(line[22:26]) in resid_CYS2CYX and 'HG' in line:
        pass # Do nothing = delete this atom
    elif int(line[22:26]) in resid_CYS2CYX:
        new_cont.append(line[:17] + 'CYX' + line[20:])  # Change from CYS to CYX
    else:
        new_cont.append(line) # Just copy paste old content
# for line in new_cont:    
    # print(line, end='')
with open(f'apo_continue_halfdS.pdb', 'w') as f:
    f.writelines(new_cont)
