

# Predict water positions
import time
import sys # To use sys.argv


import mdtraj
import numpy as np
import math


def BWHB(x1, x2, print_dist=False): # Bridge Water Hydrogen Bond
    d = np.sqrt(((x1*10-x2*10)**2).sum())
    if print_dist:
        print(f'{d:.3f}', end=' \t')
    return (1 / (1 + (d/2.6)**6)) / 0.58

def rotation3D(v, axis, degrees):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. Using the Euler-Rodrigues formula:
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    theta = degrees * math.pi / 180
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rot_mat = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    return np.dot(rot_mat, v)

def calcAngle(A, B, C):
    vec1 = np.array(A) - np.array(B)
    vec2 = np.array(C) - np.array(B)
    return np.arccos(np.dot(vec1, vec2) / np.linalg.norm(vec1) / np.linalg.norm(vec2))*180/ np.pi, vec1, vec2

def prepVec(v, length):
    v = v / np.linalg.norm(v)
    v.flatten()
    return v * length
    
def outputMockPDB(fname, p, e='X'):
    with open(fname,'w') as f:
#     f.write(f'{len(DOL)}\n\n')
        for ii in range(len(p)):
    #         ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
            f.write(f'ATOM    {ii+1:3d}  {e}   {e*3} A   1    {p[ii][0]*10:8.3f}{p[ii][1]*10:8.3f}{p[ii][2]*10:8.3f}  0.00  0.00            \n')
        f.write('END\n')
        
def outputMolecularMockPDB(fname, p, e=None, resname='HOH'):
    assert len(e) == len(p[0]), "Input: p as a N_mol * N_atom * 3 matrix, e as a list of elements with length N_atom"
    resname += ' '*3  # Fill in resnames
    e2 = [x + ' '*3 for x in e]
    with open(fname,'w') as f:
#     f.write(f'{len(DOL)}\n\n')
        for jj in range(len(p)):
            for ii in range(len(p[jj])):
    #         ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
                f.write(f'ATOM    {jj*len(p[jj])+ii+1:3d}  {e2[ii][:3]} {resname[:3]} A{jj:4d}    {p[jj][ii][0]*10:8.3f}{p[jj][ii][1]*10:8.3f}{p[jj][ii][2]*10:8.3f}  0.00  0.00            \n')
        f.write('END\n')



def calcWater(inItem, outFile=None, outResname='HOH', printTiming=False, printResult=False, processFrame=0, returnNumWater=False, returnHBscore=False):
    
    # processFrame takes a number, a list of numbers, or -1 for all frames


    if printTiming:
        t0 = time.time()

    if type(inItem) == str:
        comp = mdtraj.load(inItem)

    if type(inItem) == mdtraj.Trajectory:
        comp = inItem

    if printTiming:
        t1 = time.time()
        print(f'Load pdb: {t1-t0:.3f} s')
        t0 = time.time()
    
    protein_residue_list = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
        'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
        'SER', 'THR', 'TYR', 'VAL', 'TRP', 
        'ACE', 'NME', 'HID', 'HIE', 'HIP', 'WAT', 'HOH', 'TIP3'] # For additional caps and protonation states of HIS
    
    atom_vdW_radii = { # in nanometers
        'H': 0.107, 
        'C': 0.170, 
        'N': 0.155, 
        'O': 0.152, 
        'F': 0.147,
        'S': 0.180,
        'P': 0.180,
        'Cl': 0.175,
        'Br': 0.185,
        'I': 0.198
    }
    
    
    atom_list = {}
    bond_to_symbol = {}
    bond_to_index = {}
    atom_type = {} 
    atom_in_residue = {}
    NOS_index = []
    H_index = []
    ALL_vdW_radii = []
    
    for atom in comp.top.atoms:
    #     print(atom.index, atom)
        atom_list[atom.index] = atom.element.symbol
        atom_in_residue[atom.index] = atom.residue.name
        bond_to_symbol[atom.index] = []
        bond_to_index[atom.index] = []
        atom_type[atom.index] = []
        ALL_vdW_radii.append(atom_vdW_radii[atom.element.symbol])
        if atom.element.symbol in ['N','O','S']:
            NOS_index.append(atom.index)
        elif atom.element.symbol == 'H':
            H_index.append(atom.index)
            
    NOS_index = np.array(NOS_index)
    H_index = np.array(H_index)
    ALL_vdW_radii = np.array(ALL_vdW_radii)
    
    for bond in comp.top.bonds:
        bond_to_symbol[bond.atom1.index].append(bond.atom2.element.symbol)
        bond_to_symbol[bond.atom2.index].append(bond.atom1.element.symbol)
        bond_to_index[bond.atom1.index].append(bond.atom2.index)
        bond_to_index[bond.atom2.index].append(bond.atom1.index)
        
    
    
    for ii in range(len(atom_list)):
        bond_to_symbol[ii] = np.array(bond_to_symbol[ii])
        bond_to_index[ii] = np.array(bond_to_index[ii])
        sort_order = np.argsort(bond_to_symbol[ii])
        if len(sort_order) > 1:
    #         print(sort_order)
            bond_to_symbol[ii] = bond_to_symbol[ii][sort_order]
            bond_to_index[ii] = bond_to_index[ii][sort_order]
    #     print(ii, atom_list[ii], bond_to_symbol[ii], bond_to_index[ii])
        bond_to_symbol[ii] = list(bond_to_symbol[ii])
        bond_to_index[ii] = list(bond_to_index[ii])
        
    # Determine its atom type from bonds
    for ii in range(len(atom_list)):
        if atom_in_residue[ii] in protein_residue_list: # this is protein atom, treat with protein rules
            if atom_list[ii] == 'H':
                if bond_to_symbol[ii] == ['O']:        # This is a hydroxyl H, generate donor
                    atom_type[ii] = 'Hydroxyl H'
                    atom_type[ii] = 'OH'
                elif bond_to_symbol[ii] == ['N']:      # This is amino H, generate donor
                    atom_type[ii] = 'Amino H'
                    atom_type[ii] = 'NH'
                elif bond_to_symbol[ii] == ['C']:      # This is aliphatic H
                    atom_type[ii] = 'Aliphatic H'
                    atom_type[ii] = 'CH'
                elif bond_to_symbol[ii] == ['S']:      # This is thiol H, generate donor
                    atom_type[ii] = 'Thiol H'
                    atom_type[ii] = 'SH'
            elif atom_list[ii] == 'O':
                if bond_to_symbol[ii] == ['C']:        # This is carbonyl, generate two acceptor
                    atom_type[ii] = 'Carbonyl O'
                    atom_type[ii] = 'CO'
                    if bond_to_symbol[bond_to_index[ii][0]] == ['C','O','O']: # Carboxylate
                        atom_type[ii] = 'Carboxylate O'
                        atom_type[ii] = 'cCO'
                elif bond_to_symbol[ii] == ['C','H']:    # This is hydroxyl O, generate two acceptor, one donor
                    atom_type[ii] = 'Hydroxyl O'
                    atom_type[ii] = 'CHO'
                elif bond_to_symbol[ii] == ['H','H']:    # This is water O, generate two acceptor, two donor
                    atom_type[ii] = 'Water O'
                    atom_type[ii] = 'HHO'
            elif atom_list[ii] == 'N':
                if bond_to_symbol[ii] == ['C','C','H']:    # This is either the amide N (one donor), NME N, HIS N
                    atom_type[ii] = 'Amino N 1 donor'
                    atom_type[ii] = 'CCHN'
                elif bond_to_symbol[ii] == ['C','H','H']:    # This is ASN N, GLN N, or ARG N. Generate 2 donors each
                    atom_type[ii] = 'Amino N 2'
                    atom_type[ii] = 'CCHN'
                elif bond_to_symbol[ii] == ['C','C']:    # This is imine N (in Histidine and heterocycle cpds), one acceptor
                    atom_type[ii] = 'Imine N'
                    atom_type[ii] = 'CCN'
                elif bond_to_symbol[ii] == ['C','C','C']: # This is proline N. It has no H-bond
                    atom_type[ii] = 'Proline N'
                    atom_type[ii] = 'CCCN'
            elif atom_list[ii] == 'S':
                if bond_to_symbol[ii] == ['C','H']:    # This is thiol S, maybe generate acceptors??
                    atom_type[ii] = 'Thiol S'       
                    atom_type[ii] = 'CHS'
                elif bond_to_symbol[ii] == ['C','C']:    # This is thioester S, maybe generate acceptors??
                    atom_type[ii] = 'Thioester S'
                    atom_type[ii] = 'CCS'
    #         if atom_list[ii] != 'C' and atom_type[ii] == []:    
    #             print(ii, atom_list[ii], atom_type[ii], bond_to_symbol[ii])
        else: # Ligand residues
            if atom_list[ii] == 'H':
                if bond_to_symbol[ii] == ['O']:        # This is a hydroxyl H, generate donor
                    atom_type[ii] = 'Hydroxyl H'
                    atom_type[ii] = 'OH'
                elif bond_to_symbol[ii] == ['N']:      # This is amino H, generate donor
                    atom_type[ii] = 'Amino H'
                    atom_type[ii] = 'NH'
                elif bond_to_symbol[ii] == ['C']:      # This is aliphatic H
                    atom_type[ii] = 'Aliphatic H'
                    atom_type[ii] = 'CH'
                elif bond_to_symbol[ii] == ['S']:      # This is thiol H, generate donor
                    atom_type[ii] = 'Thiol H'
                    atom_type[ii] = 'SH'
                elif bond_to_symbol[ii] == ['P']:      # This is phosphine H, generate donor
                    atom_type[ii] = 'Phosphine H'
                    atom_type[ii] = 'PH'
    
            elif atom_list[ii] == 'O':
                if bond_to_symbol[ii] == ['C']:        # This is carbonyl, generate two acceptor
                    atom_type[ii] = 'Carbonyl O'
                    atom_type[ii] = 'CO'
                    if bond_to_symbol[bond_to_index[ii][0]] == ['C','O','O']: # Carboxylate group
                        atom_type[ii] = 'Carboxylate O'
                        atom_type[ii] = 'cCO'
                if bond_to_symbol[ii] == ['N']:        # This is nitro, generate one acceptor
                    atom_type[ii] = 'Nitro O'
                    atom_type[ii] = 'NO'
                elif bond_to_symbol[ii] == ['S']:        # This is carbonyl, generate two acceptor
                    atom_type[ii] = 'Sulfonyl / Sulfoxide O'
                    atom_type[ii] = 'SO'
                elif bond_to_symbol[ii] == ['P']:        # This is carbonyl, generate two acceptor
                    atom_type[ii] = 'Phosphine oxide / Phosphone O'
                    atom_type[ii] = 'PO'
                elif bond_to_symbol[ii] == ['C','H']:    # This is hydroxyl O, generate two acceptor, one donor
                    atom_type[ii] = 'Hydroxyl O'
                    atom_type[ii] = 'CHO'
                elif bond_to_symbol[ii] == ['C','C']:    # This is ether, generate one acceptor
                    atom_type[ii] = 'Ether O'
                    atom_type[ii] = 'CCO'
                elif bond_to_symbol[ii] == ['H','H']:    # This is water O, generate two acceptor, two donor
                    atom_type[ii] = 'Water O'
                    atom_type[ii] = 'HHO'
                elif bond_to_symbol[ii] == ['H','O']:    # This is water O, generate two acceptor, two donor
                    atom_type[ii] = 'Peroxide O'
                    atom_type[ii] = 'HOO'
                elif bond_to_symbol[ii] == ['C','O']:    # This is water O, generate two acceptor, two donor
                    atom_type[ii] = 'Peroxide O'
                    atom_type[ii] = 'COO'
            elif atom_list[ii] == 'N':
                if bond_to_symbol[ii] == ['C']:          # Cyanate / nitrile. Nitrile is an acceptor but cyanate is not
                    atom_type[ii] = 'Nitrile N'          # Default is nitrile unless the connecting C bonds to an O
                    atom_type[ii] = 'nCN'
                    for temp in bond_to_index[ii]:
                        if bond_to_symbol[temp] == 'O':
                            atom_type[ii] = 'Cyanate N'  # in which case this is an cyanate
                            atom_type[ii] = 'cCN'
    
                if bond_to_symbol[ii] == ['N']:          # Azide
                    atom_type[ii] = 'Azide N, acceptor'
                    atom_type[ii] = 'NN'
                if bond_to_symbol[ii] == ['C','C']:      # Isonitrile / sp2 nitrogen
                    atom_type[ii] = 'sp2 N acceptor'              # Default is sp2 nitroten
                    atom_type[ii] = 'aCCN'
                    for temp in bond_to_index[ii]:
                        if len(bond_to_symbol[temp]) == 1:
                            atom_type[ii] = 'Isonitrile N'  # If there's an sp carbon, than this is an isonitrile N
                            atom_type[ii] = 'iCCN'
    
                if bond_to_symbol[ii] == ['C','H']:      # Primary ketimine / aldimine
                    atom_type[ii] = 'Ketimine / aldimine N'
                    atom_type[ii] = 'CHN'
                if bond_to_symbol[ii] == ['C','N']:      # Azo / Azide
                    atom_type[ii] = 'Azo / Azide starting N'
                    atom_type[ii] = 'CNN'
                if bond_to_symbol[ii] == ['C','O']:      # Nitroso
                    atom_type[ii] = 'Nitroso / Oxime N'
                    atom_type[ii] = 'CON'
                if bond_to_symbol[ii] == ['N','N']:      # Azide
                    atom_type[ii] = 'Azide N, base of acceptor'
                    atom_type[ii] = 'NNN'
                if bond_to_symbol[ii] == ['O','O']:      # Nitrite
                    atom_type[ii] = 'Nitrite N'
                    atom_type[ii] = 'OON'
                if bond_to_symbol[ii] == ['C','C','C']:      # 
                    atom_type[ii] = 'sp3 N'
                    atom_type[ii] = 'aCCCN'
                    for temp in bond_to_index[ii]:
                        if len(bond_to_symbol[temp]) == 3:
                            atom_type[ii] = 'sp2 N'
                            atom_type[ii] = 'iCCCN'
    
                if bond_to_symbol[ii] == ['C','C','H']:      # 
                    atom_type[ii] = 'sp3 N'
                    atom_type[ii] = 'adCCHN'
                    for temp in bond_to_index[ii]:
                        if len(bond_to_symbol[temp]) == 3:
                            atom_type[ii] = 'sp2 N donor'
                            atom_type[ii] = 'dCCHN'
                if bond_to_symbol[ii] == ['C','H','H']:      # 
                    atom_type[ii] = 'sp3 N'
                    atom_type[ii] = 'adCHHN'
                    for temp in bond_to_index[ii]:
                        if len(bond_to_symbol[temp]) == 3:
                            atom_type[ii] = 'sp2 N donors'  
                            atom_type[ii] = 'dCHHN'
                                               
                if bond_to_symbol[ii] == ['H','H','H']:      # 
                    atom_type[ii] = 'Ammonia N, donors + acceptors'
                    atom_type[ii] = 'HHHN'
                if bond_to_symbol[ii] == ['C','O','O']:      # 
                    atom_type[ii] = 'Nitro N'
                    atom_type[ii] = 'COON'
                if bond_to_symbol[ii] == ['O','O','O']:      # 
                    atom_type[ii] = 'Nitrate / Nitric acid'
                    atom_type[ii] = 'OOON'
                if bond_to_symbol[ii] == ['C','C','C','C']:      # 
                    atom_type[ii] = '4 deg ammonium'
                    atom_type[ii] = 'CCCCN'
                if bond_to_symbol[ii] == ['H','H','H','H']:      # 
                    atom_type[ii] = 'Ammonium'
                    atom_type[ii] = 'HHHHN'
    
            elif atom_list[ii] == 'S':
                if bond_to_symbol[ii] == ['C']:        # This is thioketone, generate two acceptor
                    atom_type[ii] = 'Thioketone S / Isothiocyanite S / Thial S '
                    atom_type[ii] = 'CS'
                elif bond_to_symbol[ii] == ['C','H']:    # This is thiol S, maybe generate acceptors??
                    atom_type[ii] = 'Thiol S'       
                    atom_type[ii] = 'CHS'
                elif bond_to_symbol[ii] == ['C','C']:    # This is thioester S, maybe generate acceptors??
                    atom_type[ii] = 'Thioether S / Thiocyanate S'
                    atom_type[ii] = 'CCS'
                elif bond_to_symbol[ii] == ['C','S']:    # Disulfide
                    atom_type[ii] = 'Disulfide S'
                    atom_type[ii] = 'CSS'
                elif bond_to_symbol[ii] == ['C','C','O']: # Sulfoxide
                    atom_type[ii] = 'Sulfoxide S'
                    atom_type[ii] = 'CCOS'
                elif bond_to_symbol[ii] == ['C','C','O','O']: # Sulfone
                    atom_type[ii] = 'Sulfone S'
                    atom_type[ii] = 'CCOOS'
                elif bond_to_symbol[ii] == ['C','O','O']: # Sulfinic acid
                    atom_type[ii] = 'Sulfinic acid S'
                    atom_type[ii] = 'COOS'
                elif bond_to_symbol[ii] == ['C','O','O','O']: # Sulfonic acid / Ester
                    atom_type[ii] = 'Sulfonic acid / ester S'
                    atom_type[ii] = 'COOOS'
            elif atom_list[ii] == 'P':
                if bond_to_symbol[ii] == ['C','C','C']:
                    atom_type[ii] = 'CCCP'
                if bond_to_symbol[ii] == ['C','C','H']:
                    atom_type[ii] = 'CCCP'
                if bond_to_symbol[ii] == ['C','H','H']:
                    atom_type[ii] = 'CCCP'
                if bond_to_symbol[ii] == ['H','O','O','O']:
                    atom_type[ii] = 'HOOOP'
                if bond_to_symbol[ii] == ['O','O','O','O']:
                    atom_type[ii] = 'OOOOP'

    if printTiming: 
        t1 = time.time()
        print(f'Parse atom types: {t1-t0:.3f} s')
        t0 = time.time()


    if processFrame == -1:
        frames = range(len(comp))
    elif type(processFrame) == int:
        frames = [processFrame]
    elif type(processFrame) == list:
        frames = processFrame

    for frame in frames:
        NOS_coordinate = []
        H_coordinate = []
        ALL_coordinate = comp.xyz[frame]
        NOS_coordinate = comp.xyz[frame][NOS_index]
    #    print(NOS_coordinate.shape)
        H_coordinate = comp.xyz[frame][H_index]
    #    NOS_coordinate = np.array(NOS_coordinate)
    #    H_coordinate = np.array(H_coordinate)
        
        # Donor generation
        DOL = [] # Donor Oxygen Location
        DOL_index = []
        for ii in range(len(atom_list)):
            if atom_type[ii] in ['OH', 'NH', 'SH']: 
                vec = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]] 
                vec /= np.linalg.norm(vec)
                vec = vec.flatten()
                if atom_type[ii] == 'OH':
                    vec *= 0.18
                elif atom_type[ii] == 'NH':
                    vec *= 0.2
                elif atom_type[ii] == 'SH':
                    vec *= 0.245  # Temporary
                DOL.append(comp.xyz[frame][ii] + vec)
                DOL_index.append(ii)
        # print(np.array(DOL))
        DOL = np.array(DOL)
        DOL_index = np.array(DOL_index)
        # # Write out a test xyz
        # with open('test_donor.pdb','w') as f:
        # #     f.write(f'{len(DOL)}\n\n')
        #     for ii in range(len(DOL)):
        # #         ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
        #         f.write(f'ATOM    {ii+1:3d}  X   XXX A{ii:4d}    {DOL[ii][0]*10:8.3f}{DOL[ii][1]*10:8.3f}{DOL[ii][2]*10:8.3f}  0.00  0.00            \n')
        #     f.write('END\n')
        #     
        
        
        # Acceptor generation
        counter = 0
        AHL = [] # Acceptor hydrogen Location
        AHL_index = []
        for ii in range(len(atom_list)):
            if atom_type[ii] in ['CO', 'cCO']: # Carbonyl
        #         counter += 1
                vec1 = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]] 
                # print(ii)
                # print(bond_to_index[bond_to_index[ii][0]])
                for temp in bond_to_index[bond_to_index[ii][0]]:
                    if temp != ii:
                        vec2 = comp.xyz[frame][temp] - comp.xyz[frame][bond_to_index[ii][0]] 
                        rot_axis = np.cross(vec1.flatten(), vec2.flatten())
                        # print(rot_axis)
                        break
                vec = rotation3D(vec1, rot_axis, 60)
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.194))
                AHL_index.append(ii)
                vec = rotation3D(vec1, rot_axis, -60)
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.194))
                AHL_index.append(ii)
            if atom_type[ii] in ['NO']: # Nitro, very similar to CO
        #         counter += 1
                vec1 = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]] 
                # print(ii)
                # print(bond_to_index[bond_to_index[ii][0]])
                for temp in bond_to_index[bond_to_index[ii][0]]:
                    if temp != ii:
                        vec2 = comp.xyz[frame][temp] - comp.xyz[frame][bond_to_index[ii][0]] 
                        rot_axis = np.cross(vec1.flatten(), vec2.flatten())
                        # print(rot_axis)
                        break
                vec = rotation3D(vec1, rot_axis, 60)
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.233))
                AHL_index.append(ii)
                vec = rotation3D(vec1, rot_axis, -60)
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.233))
                AHL_index.append(ii)
        #     if atom_type[ii] in ['cCO']: # Carboxylate
        #         vec = comp.xyz[frame][bond_to_index[ii][0]] - comp.xyz[frame][bond_to_index[bond_to_index[ii][0]][0]]
        #         AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.192))
        #         AHL_index.append(ii)
            if atom_type[ii] in ['CHO']: # Hydroxyl
                vec1 = comp.xyz[frame][bond_to_index[ii][1]] - comp.xyz[frame][ii] # OH
                vec2 = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]] # CH
                vec = rotation3D(vec1, vec2, 180)
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.203))
                AHL_index.append(ii)
                
            if (atom_type[ii] in ['aCCN','CCN']): # mostly histidine and nitrogens in heterocycles
                counter += 1
        #         vec = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii]] 
                # print(ii)
                # print(bond_to_index[bond_to_index[ii][0]])
                avgCC = np.mean(comp.xyz[frame][bond_to_index[ii]],axis=0)
                vec = comp.xyz[frame][ii] - avgCC
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.21))
                AHL_index.append(ii)
            if atom_type[ii] in ['CCO']: # Ether
                counter += 1
        #         vec = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii]] 
                # print(ii)
                # print(bond_to_index[bond_to_index[ii][0]])
                avgCC = np.mean(comp.xyz[frame][bond_to_index[ii]],axis=0)
                vec = comp.xyz[frame][ii] - avgCC
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.21))
                AHL_index.append(ii)
            if atom_type[ii] in ['SO']: # Sulfoxide
                vec = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]]
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.196))
                AHL_index.append(ii)
            if atom_type[ii] in ['PO']: # Phosphine
                vec = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]]
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.188))
                AHL_index.append(ii)
            if atom_type[ii] in ['nCN']: # Nitrile
                vec = comp.xyz[frame][ii] - comp.xyz[frame][bond_to_index[ii][0]]
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.21))
                AHL_index.append(ii)
            if atom_type[ii] in ['aCCCN', 'adCCHN', 'adCHHN', 'HHHN']: # Various acceptors for amines
                avgCC = np.mean(comp.xyz[frame][bond_to_index[ii]],axis=0)
                vec = comp.xyz[frame][ii] - avgCC
                AHL.append(comp.xyz[frame][ii] + prepVec(vec, 0.213))
                AHL_index.append(ii)
        
        
        AHL = np.array(AHL)
        AHL_index = np.array(AHL_index)
        
        # # Write out a test xyz
        # with open('test_acceptor.pdb','w') as f:
        # #     f.write(f'{len(DOL)}\n\n')
        #     for ii in range(len(AHL)):
        # #         ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
        #         f.write(f'ATOM    {ii+1:3d}  Z   ZZZ A{ii:4d}    {AHL[ii][0]*10:8.3f}{AHL[ii][1]*10:8.3f}{AHL[ii][2]*10:8.3f}  0.00  0.00            \n')
        #     f.write('END\n')
        #     
    
        if printTiming: 
            t1 = time.time()
            print(f'Generate donors / acceptors: {t1-t0:.3f} s')
            t0 = time.time()
        
        # Calculate DOL -> N, O, S distance
        DOL_NOS_dist = np.min(np.sqrt(((DOL[:,:,None] - NOS_coordinate[:,:,None].T)**2).sum(1)),axis=1)
        # print(DOL_NOS_dist.shape)
        AHL_H_dist = np.min(np.sqrt(((AHL[:,:,None] - H_coordinate[:,:,None].T)**2).sum(1)),axis=1)
        # print(AHL_H_dist.shape)
        # import matplotlib.pyplot as plt
        # plt.hist(DOL_NOS_dist, bins=30)
        # plt.hist(AHL_H_dist, bins=30)
        # np.nonzero(AHL_H_dist > 0.1)[0]
        AHL_inter = AHL[np.nonzero(AHL_H_dist > 0.12)[0]]
        AHL_inter_index = AHL_index[np.nonzero(AHL_H_dist > 0.12)[0]] # Atom indices of acceptor atoms
        DOL_inter = DOL[np.nonzero(DOL_NOS_dist > 0.12)[0]]
        DOL_inter_index = DOL_index[np.nonzero(DOL_NOS_dist > 0.12)[0]]
        # outputMockPDB('test_AHL_inter.pdb',AHL_inter)
        # outputMockPDB('test_DOL_inter.pdb',DOL_inter)
        
        # Calculate inter-donor/acceptor distances
        DOL_DOL = np.sqrt(((DOL_inter[:,:,None] - DOL_inter[:,:,None].T)**2).sum(1))
        np.fill_diagonal(DOL_DOL, 1)
        AHL_DOL = np.sqrt(((AHL_inter[:,:,None] - DOL_inter[:,:,None].T)**2).sum(1))
        # np.fill_diagonal(AHL_DOL, 1)
        AHL_AHL = np.sqrt(((AHL_inter[:,:,None] - AHL_inter[:,:,None].T)**2).sum(1))
        np.fill_diagonal(AHL_AHL, 1)
        
        # Qualify donor/acceptors based on distances
        DOL_DOL_q = np.where(DOL_DOL < 0.14)
        AHL_DOL_q = np.where(AHL_DOL < 0.28)
        AHL_AHL_q = np.where(AHL_AHL < 0.28)
        
        if printTiming: 
            t1 = time.time()
            print(f'Qualify donors / acceptors: {t1-t0:.3f} s')
            t0 = time.time()
        
        
        
        # Place water molecules. Note this is ordered - AHL+AHL takes priority. 
        # This will happen in redundant water elimination and hydrogen placements
        
        temp_water = []
        temp_water_index = []
        temp_water_type = []
        for ii in range(len(AHL_AHL_q[0])):
            temp_water.append((AHL_inter[AHL_AHL_q[0][ii]]+AHL_inter[AHL_AHL_q[1][ii]])/2)
            temp_water_index.append([AHL_inter_index[AHL_AHL_q[0][ii]], AHL_inter_index[AHL_AHL_q[1][ii]]])
            temp_water_type.append(0)
        for ii in range(len(AHL_DOL_q[0])):
            temp_water.append((AHL_inter[AHL_DOL_q[0][ii]]+DOL_inter[AHL_DOL_q[1][ii]])/2)
            temp_water_index.append([AHL_inter_index[AHL_DOL_q[0][ii]], DOL_inter_index[AHL_DOL_q[1][ii]]])
            temp_water_type.append(1)
        for ii in range(len(DOL_DOL_q[0])):
            temp_water.append((DOL_inter[DOL_DOL_q[0][ii]]+DOL_inter[DOL_DOL_q[1][ii]])/2)
            temp_water_index.append([DOL_inter_index[DOL_DOL_q[0][ii]], DOL_inter_index[DOL_DOL_q[1][ii]]])
            temp_water_type.append(2)
        temp_water = np.array(temp_water)
        temp_water_index = np.array(temp_water_index)
        temp_water_type = np.array(temp_water_type)
        
        # outputMockPDB('temp_water.pdb',temp_water,e='O')
        
        if printTiming: 
            t1 = time.time()
            print(f'Place template oxygens: {t1-t0:.3f} s')
            t0 = time.time()
        
        # Test for waters that are too close
        ALL_WAT = np.sqrt(((ALL_coordinate[:,:,None] - temp_water[:,:,None].T)**2).sum(1)) - ALL_vdW_radii[:,None] - 0.06
        
        qualified_water = temp_water[np.where(np.all(ALL_WAT > 0, axis=0))]
        qualified_water_index = temp_water_index[np.where(np.all(ALL_WAT > 0, axis=0))]
        qualified_water_type = temp_water_type[np.where(np.all(ALL_WAT > 0, axis=0))]
        
        # outputMockPDB('qualified_water.pdb',qualified_water,e='O')
        
        if printTiming:
            t1 = time.time()
            print(f'Place qualified oxygens: {t1-t0:.3f} s')
            t0 = time.time()
        
        # Orient each qualifying water
        full_water = []
        full_water_index = []
        for ii in range(len(qualified_water)):
            this_water = np.array([qualified_water[ii], [0,0,0], [0,0,0]])
            if qualified_water_type[ii] == 0: # A two donor water - H's of water point to acceptors
                this_water[1] = this_water[0] + prepVec(comp.xyz[frame][qualified_water_index[ii][0]] - this_water[0], 0.09572) # Acceptor 1
                this_water[2] = this_water[0] + prepVec(comp.xyz[frame][qualified_water_index[ii][1]] - this_water[0], 0.09572) # Acceptor 2
                # Adjust angles to 104.5 deg
                thisAngle, vec1, vec2 = calcAngle(comp.xyz[frame][qualified_water_index[ii][0]], this_water[0], comp.xyz[frame][qualified_water_index[ii][1]])
                thisAxis = np.cross(vec1, vec2)
                this_water[1] = this_water[0] + rotation3D(this_water[1] - this_water[0], thisAxis, (thisAngle - 104.5)/2)
                this_water[2] = this_water[0] + rotation3D(this_water[2] - this_water[0], thisAxis, -(thisAngle - 104.5)/2)
            if qualified_water_type[ii] == 1: # A acceptor / donor water
                this_water[1] = this_water[0] + prepVec(comp.xyz[frame][qualified_water_index[ii][0]] - this_water[0], 0.09572) # Acceptor 1
                # Create second hydrogen and rotate that away (so the lone pair faces the donor)
                thisAngle, vec1, vec2 = calcAngle(comp.xyz[frame][qualified_water_index[ii][0]], this_water[0], comp.xyz[frame][qualified_water_index[ii][1]])
                thisAxis = np.cross(vec1, vec2)
                this_water[2] = this_water[0] + rotation3D(this_water[1] - this_water[0], thisAxis, 104.5)
                this_water[2] = this_water[0] + rotation3D(this_water[2] - this_water[0], this_water[1] - this_water[0], 120)
                #parallel_NH = comp.xyz[frame][qualified_water_index[ii][1]] - comp.xyz[frame][bond_to_index[qualified_water_index[ii][1]]]
                #this_water[2] = this_water[0] + prepVec(parallel_NH, 0.09572) # Acceptor 2
            full_water.append(this_water)
            full_water_index.append(qualified_water_index[ii])
        # print(full_water_index)
        
        # outputMolecularMockPDB('Full_water.pdb', full_water, e=['O','H','H'])
       
        if printTiming: 
            t1 = time.time()
            print(f'Place full waters: {t1-t0:.3f} s')
            t0 = time.time()
        
        
        # Geometry check
        geo_water = []
        geo_water_O = []
        geo_water_index = []
        lowerLength = {'O':0.161, 'N': 0.171, 'S':0.229}
        upperLength = {'O':0.230, 'N': 0.240, 'S':0.312}
        for ii in range(len(full_water)):
            # print(f'WATER {ii}')
            if qualified_water_type[ii] == 0: # Test two distances and two angles
                testLength = np.linalg.norm(full_water[ii][1] - comp.xyz[frame][qualified_water_index[ii][0]])
                # print(testLength)
        #         if atom_list[qualified_water_index[ii][0]] in ['O','S','N']:
                if (testLength < lowerLength[atom_list[qualified_water_index[ii][0]]]) or \
                   (testLength > upperLength[atom_list[qualified_water_index[ii][0]]]):
                    # print(f'Water {ii} failed due to distance between H1 and acceptor')
                    continue
                testLength = np.linalg.norm(full_water[ii][2] - comp.xyz[frame][qualified_water_index[ii][1]])
                # print(testLength)
        #         if atom_list[qualified_water_index[ii][0]] in ['O','S','N']:
                if (testLength < lowerLength[atom_list[qualified_water_index[ii][1]]]) or \
                   (testLength > upperLength[atom_list[qualified_water_index[ii][1]]]):
                    # print(f'Water {ii} failed due to distance between H2 and acceptor')
                    continue
                testAngle, _, _ = calcAngle(full_water[ii][0], full_water[ii][1], comp.xyz[frame][qualified_water_index[ii][0]])
                # print(testAngle)
                if testAngle < 120:
                    # print(f'Water {ii} failed due to angle between O, H1 and acceptor')
                    continue
                testAngle, _, _ = calcAngle(full_water[ii][0], full_water[ii][2], comp.xyz[frame][qualified_water_index[ii][1]])
                # print(testAngle)
                if testAngle < 120:
                    # print(f'Water {ii} failed due to angle between O, H2 and acceptor')
                    continue
            if qualified_water_type[ii] == 1: 
                # Test 1. distance between OH and acceptor, 
                #      2. angle between (donor, donor H, O), and 
                #      3. angle between (accpetor expected H, acceptor, H)
                testLength = np.linalg.norm(full_water[ii][1] - comp.xyz[frame][qualified_water_index[ii][0]])
                # print(testLength)
        
                if (testLength < lowerLength[atom_list[qualified_water_index[ii][0]]]) or \
                   (testLength > upperLength[atom_list[qualified_water_index[ii][0]]]):
                    # print(f'Water {ii} failed due to distance between H1 and acceptor')
                    continue
                # print(f'Qualified H, {qualified_water_index[ii][1]} (H) bonds to {bond_to_index[qualified_water_index[ii][1]][0]}')
                testAngle, _, _ = calcAngle(comp.xyz[frame][bond_to_index[qualified_water_index[ii][1]][0]], 
                                            comp.xyz[frame][qualified_water_index[ii][1]],
                                            full_water[ii][0])
                if testAngle < 120:
                    # print(f'Water {ii} failed due to angle between donor, donor H and water O')
                    continue
                for jj in np.where(AHL_index == qualified_water_index[ii][0])[0]:
                    testAngle, _, _ = calcAngle(AHL[jj],
                                                comp.xyz[frame][qualified_water_index[ii][0]],
                                                full_water[ii][1])
                    if testAngle < 45:
        #                 print(f'passed angle test')
                        break
        #             print(f'jjjj {testAngle}')
                else:
                    continue
            geo_water.append(full_water[ii])
            geo_water_O.append(full_water[ii][0])
            geo_water_index.append(full_water_index[ii])
        geo_water = np.array(geo_water)    
        geo_water_O = np.array(geo_water_O)
        geo_water_index = np.array(geo_water_index)
        
        # print(geo_water_index)
    
        if printTiming:
            t1 = time.time()
            print(f'Water geometry check: {t1-t0:.3f} s')
            t0 = time.time()
        
        # outputMolecularMockPDB('Geo_water.pdb', geo_water, e=['O','H','H'])
        
        # Eliminate duplicate O
        accept_water = []
        GEO_GEO = np.sqrt(((geo_water_O[:,:,None] - geo_water_O[:,:,None].T)**2).sum(1))
        np.fill_diagonal(GEO_GEO, 1)
        # print(GEO_GEO)
        for ii in range(len(GEO_GEO)):
            if np.any(GEO_GEO[ii][:ii] < 0.249):
                GEO_GEO[ii] = 1
                GEO_GEO[:,ii] = 1
            else:
                accept_water.append(ii)
        # print(GEO_GEO)
        # print(accept_water)
        final_water = geo_water[accept_water]
        final_water_O = geo_water_O[accept_water]
        final_water_index = geo_water_index[accept_water]
        
        
        if returnHBscore:
            HB_score = 0
            for ii in range(len(final_water)):

                atom1 = comp.top.atom(final_water_index[ii][0])
                atom2 = comp.top.atom(final_water_index[ii][1])
                
                if atom1.residue.name == 'LIG' or atom2.residue.name == 'LIG':
                    if atom1.residue.name == 'LIG' and atom2.residue.name == 'LIG':
                        pass
                    else:
                        if atom1.element.symbol == 'H':
                            # print(atom1)
                            # print(comp.top.atom(bond_to_index[final_water_index[ii][0]][0]), end=' \t')
                            btom1_index = bond_to_index[final_water_index[ii][0]][0]
                        else:
                            # print(atom1, end=' \t')
                            btom1_index = final_water_index[ii][0]
                        HB_score += BWHB(ALL_coordinate[btom1_index], final_water_O[ii], print_dist=False)
                        # print(f'{HB(ALL_coordinate[btom1_index], final_water_O[ii], print_dist=False):.3f}')
                        if atom2.element.symbol == 'H':
                            # print(comp.top.atom(bond_to_index[final_water_index[ii][1]][0]), end=' \t')
                            btom2_index = bond_to_index[final_water_index[ii][1]][0]
                        else:
                            # print(atom2, end=' \t')
                            btom2_index = final_water_index[ii][1]
                        HB_score += BWHB(ALL_coordinate[btom2_index], final_water_O[ii], print_dist=False)
                        # print(f'{HB(ALL_coordinate[btom2_index], final_water_O[ii], print_dist=True):.3f}')
            try:
                HB_score_list.append(HB_score)
            except:
                HB_score_list = [HB_score]
        # print()
        
        # print(final_water_index)
       
        if printTiming: 
            t1 = time.time()
            print(f'Final water determination: {t1-t0:.3f} s')
            t0 = time.time()
        
        if outFile is not None:
            outputMolecularMockPDB(outFile, final_water, e=['O','H1','H2'], resname=outResname)
        
        if printTiming:
            t1 = time.time()
            print(f'Output final file: {t1-t0:.3f} s')
            t0 = time.time()
    
        if printResult:
            print(f'For this input we can place {len(final_water)} waters')

        if returnNumWater:
            try:
                numWater.append(len(final_water))
            except:
                numWater = [len(final_water)]

    if returnNumWater:
        return numWater
    
    if returnHBscore:
        return HB_score_list
        
