
import glob
import numpy as np
import os, sys
from shutil import copyfile

class PDB:
    def __init__(self, target, center, lig_num_atoms=60):
        
        # center is a N*3-element coord array like this: [[x, y, z], [x, y, z] ...]

        self.target = target

        with open(f'{self.target}_h.pdb', 'r') as f:
            self.target_pdb = f.readlines()

        self.lig = "the ligand center"
        self.lig_xyz = center
        self.lig_num_atoms = lig_num_atoms

        self.target_pdb = parse_PDB(self.target_pdb)
        self.target_idx = np.array([int(x[1]) for x in self.target_pdb if x[0] == 'ATOM'])
        self.target_xyz = [x[6:9] for x in self.target_pdb if x[0] == 'ATOM']
        self.target_xyz = np.array([list(map(float, x)) for x in self.target_xyz])
        self.target_resid = np.array([int(x[12]) for x in self.target_pdb if x[0] == 'ATOM'])

        self.target_resid_chain_map = {}
        for idx, atom in zip(self.target_resid,self.target_pdb):
            self.target_resid_chain_map[idx] = atom[4]
        #print(self.target_resid_chain_map)


        self.dist_matrix = np.sqrt(((self.lig_xyz[:,:,None] - self.target_xyz[:,:,None].T)**2).sum(1))

    def determine_close_residues(self, thres=6.7, verbose=True):
        self.thres = thres
        self.close_atoms = np.where(np.sum(self.dist_matrix < self.thres, axis=0) > 0)
        self.close_res = np.unique(self.target_resid[self.close_atoms])
        if verbose:
            print(f'There are {len(self.close_res)} residues in {self.target} closer than {self.thres} A to {self.lig}')
            print(f'There are {len([x for x in self.target_resid if x in self.close_res])} atoms (including hydrogen) in these residues')
            print('')

    def determine_threshold(self, target_atoms=850, verbose=False):
        self.included_atoms = 0
        self.thres = 4.9
        while (self.included_atoms < target_atoms):
            self.thres += 0.1
            self.determine_close_residues(thres=self.thres, verbose=False)
            self.included_atoms = self.join_segments(verbose=False)
        # In case included atoms is larger than 928, revert one step
        if self.included_atoms > 925:
            self.thres -= 0.1
            self.determine_close_residues(thres=self.thres, verbose=False)
            self.included_atoms = self.join_segments(verbose=False)
        #print(self.segments)
        # print(f'There are {len(self.close_res)} residues in {self.target} closer than {self.thres} A to {self.lig}')
        # print(f'There are {len([x for x in self.target_resid if x in self.close_res])} atoms (including hydrogen) in these residues')
        print(f'Determined threshold to be {self.thres:.1f} A')
        self.print_joined_segments()

    def join_segments(self, verbose=True):
        self.cont_res = list(self.close_res)
        #print(self.cont_res)
        self.segments = []
        self.segment = []
        for ii in range(len(self.close_res)):
            #print(f'This resid (no. {ii}) is {self.close_res[ii]}')
            if ii == 0:
                self.segment.append(self.close_res[ii])
                current_chain = self.target_resid_chain_map[self.close_res[ii]]
            elif self.close_res[ii] - self.close_res[ii-1] == 1:
                if self.target_resid_chain_map[self.close_res[ii]] != current_chain:
                    # If the next res is in a different chain, break the chain
                    self.segments.append(self.segment)
                    self.segment = [self.close_res[ii]]
                    current_chain = self.target_resid_chain_map[self.close_res[ii]]
                else:
                    self.segment.append(self.close_res[ii])
            elif self.close_res[ii] - self.close_res[ii-1] == 2:
                if self.target_resid_chain_map[self.close_res[ii]] != current_chain:
                    # Then there's no need to include anything in between
                    # End the last segment and start a new one
                    self.segments.append(self.segment)
                    self.segment = [self.close_res[ii]]
                    #self.cont_res.append(self.close_res[ii]-1)
                    current_chain = self.target_resid_chain_map[self.close_res[ii]]
                else:
                    # Otherwise this forms a continuous segment
                    self.segment.append(self.close_res[ii]-1)
                    self.cont_res.append(self.close_res[ii]-1)
                    self.segment.append(self.close_res[ii])
                    self.cont_res.append(self.close_res[ii])


            elif self.close_res[ii] - self.close_res[ii-1] == 3:
                if self.target_resid_chain_map[self.close_res[ii]] != current_chain:
                    # Then there's no need to include anything in between
                    # End the last segment and start a new one
                    self.segments.append(self.segment)
                    self.segment = [self.close_res[ii]]
                    #self.cont_res.append(self.close_res[ii]-2)
                    current_chain = self.target_resid_chain_map[self.close_res[ii]]
                else:
                    # Otherwise this forms a continuous segment
                    self.segment.append(self.close_res[ii]-2)
                    self.cont_res.append(self.close_res[ii]-2)
                    self.segment.append(self.close_res[ii]-1)
                    self.cont_res.append(self.close_res[ii]-1)
                    self.segment.append(self.close_res[ii])
                    self.cont_res.append(self.close_res[ii])



            else:
                self.segments.append(self.segment)
                self.segment = [self.close_res[ii]]
                current_chain = self.target_resid_chain_map[self.close_res[ii]]

        if len(self.segment) > 0:
            self.segments.append(self.segment)
            self.segment = []
        self.cont_res.sort()

        self.est_receptor_atoms = len([x for x in self.target_resid if x in self.cont_res])
        self.est_cap_atoms = 12 * len(self.segments)
        self.est_lig_atoms = self.lig_num_atoms
        self.est_included_atoms = self.est_receptor_atoms + self.est_cap_atoms + self.est_lig_atoms
        if verbose:
            self.print_joined_segments()

        return self.est_included_atoms

    def print_joined_segments(self):

        print(f'Determined segments: ',end='')
        for seg in self.segments:
            if len(seg) > 1:
                print(f'{seg[0]} to {seg[-1]}', end=', ')
            else:
                print(f'{seg[0]}', end=', ')
        print('')

        # print(self.cont_res)
        print(f'There are {len(self.cont_res)} residues that will be included in the core')
        print(f'Rough estimate of how many atoms there will be: {self.est_included_atoms}')
        print(f'  with {self.est_receptor_atoms} receptor atoms, {self.est_cap_atoms} cap atoms, and {self.est_lig_atoms} ligand atoms')


    def write_all_PDBs(self, dry_run=True, path='./', write_script=False):
        try:
            os.makedirs(path)
        except:
            pass
        self.file_names = []
        for ii in range(len(self.segments)):
            first_flag = False
            last_flag = False
            if ii < len(self.segments) - 1:
                if self.segments[ii][-1] == self.segments[ii+1][0] - 1:
                    last_flag = True
            if ii > 0:
                if self.segments[ii][0] == self.segments[ii-1][-1] + 1:
                    first_flag = True
            if ii == 0:
                if self.segments[ii][0] == 1:
                    first_flag = True
                elif self.target_resid_chain_map[self.segments[ii][0]] != self.target_resid_chain_map[self.segments[ii][0]-1]:
                    first_flag = True # Example: this residue is chain B resid 36, and last residue is chain A resid 35
            if ii == len(self.segments)-1:
                if self.segments[ii][-1] == np.max(self.target_resid):
                    last_flag = True
                elif self.target_resid_chain_map[self.segments[ii][-1]] != self.target_resid_chain_map[self.segments[ii][-1]+1]:
                    last_flag = True # Example: this residue is chain B resid 300 and next residue is chain C resid 301
            if first_flag and last_flag:
                print(f'Segment {self.segments[ii][0]}t to {self.segments[ii][-1]}t')
            elif first_flag:
                print(f'Segment {self.segments[ii][0]}t to {self.segments[ii][-1]}')
            elif last_flag:
                print(f'Segment {self.segments[ii][0]} to {self.segments[ii][-1]}t')
            else:
                print(f'Segment {self.segments[ii][0]} to {self.segments[ii][-1]}')
            self.write_PDB(ii, dry_run = dry_run, path = path, first=first_flag, last=last_flag)
        self.write_PDB('all', dry_run = dry_run, path = path)
        # pass
        if write_script:
            self.write_all_scripts(path)

    def write_PDB(self, seg, dry_run=True, path='./', first=False, last=False):
        outputPDB = []
        counter = 1
        # if type(seg) == list:
        #     seg_list = seg
        #     seg_text = f'{seg:02d}'
        if type(seg) == int:
            seg_list = self.segments[seg]
            seg_text = f'{seg:02d}'
        elif seg == 'all':
            seg_list = np.unique(self.target_resid)
            seg_text = f'{seg}'
        elif type(seg) == list:
            raise NotImplementedError()

        if dry_run == True:
            for res, atom in zip(self.target_resid, self.target_pdb):
                if res in seg_list:
                    if atom[1] != '0':
                        print(atom_line(res, atom, counter),end='')
                        counter += 1
        else: # Write to file

            if first and last:
                file_name = f'r{seg_text}-{seg_list[0]}t-{seg_list[-1]}t.pdb'
            elif first:
                file_name = f'r{seg_text}-{seg_list[0]}t-{seg_list[-1]}.pdb'
            elif last:
                file_name = f'r{seg_text}-{seg_list[0]}-{seg_list[-1]}t.pdb'
            else:
                file_name = f'r{seg_text}-{seg_list[0]}-{seg_list[-1]}.pdb'
            self.file_names.append(file_name)
            with open(f'{path}/{file_name}','w') as f:
                for res, atom in zip(self.target_resid, self.target_pdb):
                    if res in seg_list:
                        if atom[1] != '0':
                            f.write(atom_line(res, atom, counter))
                            # print(atom_line(res, atom, counter),end='')
                            counter += 1
        # pass # Use the re-indexed resids as the new resid. Store the old ID in either occupancy or beta column

    def write_init_antechamber_script(self, path='./'):
        assert len(self.file_names) > 0, "First output pdb files [write_all_PDBs(path)] then do this!"
        with open(f'{path}/init_antechamber.in','w') as f:
            f.write('### Autogenerated by prepContRes ###\n')
            f.write('source leaprc.protein.ff14SB #Source leaprc file for ff14SB protein force field\n')
            for file_name in self.file_names:
                if 'all' not in file_name:
                    fbase = file_name.split('.')[0]
                    # print(fbase)
                    f.write(f'mol = loadpdb {fbase}.pdb\n')
                    f.write(f'savepdb mol {fbase}-b.pdb\n')


    def write_second_python_cap_script(self, path='./'):
        copyfile('../../01_Workflow/utilities/add_cap.py', f'{path}/add_cap.py')

    def write_third_antechamber_script(self, path='./'):
        assert len(self.file_names) > 0, "First output pdb files [write_all_PDBs(path)] then do this!"
        with open(f'{path}/third_antechamber.in','w') as f:
            f.write('### Autogenerated by prepContRes ###\n')
            f.write('source leaprc.protein.ff14SB #Source leaprc file for ff14SB protein force field\n')
            for file_name in self.file_names:
                if 'all' not in file_name:
                    fbase = file_name.split('.')[0]
                    f.write(f'mol = loadpdb {fbase}-c.pdb\n')
                    f.write(f'savepdb mol {fbase}-d.pdb\n')
                    f.write(f'saveamberparm mol {fbase}-c.prmtop {fbase}-c.inpcrd\n')


    def write_assemble_script(self, path='./'):
        copyfile('../../01_Workflow/utilities/assemble.py', f'{path}/assemble.py')

    def write_fix_ca_script(self, path='./'):
        copyfile('../../01_Workflow/utilities/fix_ca.py', f'{path}/fix_ca.py')
        copyfile('../../01_Workflow/utilities/fix_ca_openmm.py', f'{path}/fix_ca_openmm.py')
    
    def write_all_scripts(self, path='./'):
        self.write_init_antechamber_script(path)
        self.write_second_python_cap_script(path)
        self.write_third_antechamber_script(path)
        self.write_assemble_script(path)
        self.write_fix_ca_script(path)
        with open(f'{path}/overall_script.sh','w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('tleap -f init_antechamber.in > init_antechamber.out\n')
            f.write('python add_cap.py\n')
            f.write('tleap -f third_antechamber.in > third_antechamber.out\n')
            f.write('python assemble.py\n')
            f.write('tleap -f fifth_antechamber.in > fifth_antechamber.out\n')
            f.write('python fix_ca.py\n')
            f.write('python fix_ca_openmm.py\n')

def atom_line(res, atom, counter): # Atom is a line in self.target_pdb
    O = (int(atom[1]) - int(atom[1]) % 1000)/1000 + float((int(atom[5]) - int(atom[5])%100)/ 100 / 100)
    B = int(atom[1]) % 1000 + float(int(atom[5])%100 / 100)
    if len(atom[2]) < 4:
        a2 = ' ' + atom[2]
    else:
        a2 = atom[2]
    # print(type(a2))
    if len(atom[3]) < 4:
        a3 = atom[3] + ' '
    else:
        a3 = atom[3]

        #         'ATOM'       index        atom name    resname      chain        resid          ASSS      X          Y             Z        O,B-> ori. resid  rest
    return f'{atom[0]:<6s}{str(counter):>5s} {a2:<4s} {a3:>4s}{atom[4]:1s}{str(res):>4s}    {atom[6]:>8s}{atom[7]:>8s}{atom[8]:>8s}{O:6.2f}{B:6.2f}{atom[11]}'

def parse_PDB(PDB_content):
    splitPDB = []
    for line in PDB_content:
        if 'ATOM' in line:
            splitPDB.append([line[:6].strip(),    # 0. 'ATOM'
                             line[6:11].strip(),  # 1. index
                             line[12:16].strip(), # 2. atom name
                             line[16:20].strip(), # 3. residue name
                             line[21],            # 4. chain name
                             line[22:27].strip(), # 5. residue number
                             line[30:38].strip(), # 6. X
                             line[38:46].strip(), # 7. Y
                             line[46:54].strip(), # 8. Z
                             line[54:60].strip(), # 9. O
                             line[60:66].strip(), # 10. B
                             line[66:]            # 11. Rest of the line
                            ]) # B
            #The rest of the line we don't concern right now

    # Re-index resids using N to deal with resids like 60A, 60B etc.
    counter = 0
    for atom in splitPDB:
        if atom[2] == 'N':
            counter += 1
        atom.append(counter)                      # 12. reindexed resid

    return splitPDB




