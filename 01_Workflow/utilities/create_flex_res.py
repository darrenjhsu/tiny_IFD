
import numpy as np
import sys


def atom_dict(res):
    if res == 'VAL':
        return {'CA','CB','CG1','CG2'}
    if res == 'LEU':
        return {'CA','CB','CG','CD1','CD2'}
    if res == 'ILE':
        return {'CA','CB','CG2','CG1','CD1'} 
    if res == 'MET':
        return {'CA','CB','CG','SD','CE'} 
    if res == 'PHE':
        return {'CA','CB','CG','CD1','CE1','CZ','CE2','CD2'} 
    if res == 'TRP':
        return {'CA','CB','CG','CD1','NE1','HE1','CE2','CZ2','CH2','CZ3','CE3','CD2'} 
    if res == 'SER':
        return {'CA','CB','OG','HG'}
    if res == 'THR':
        return {'CA','CB','CG2','OG1','HG1'}
    if res == 'CYS':
        return {'CA','CB','SG','HG'}
    if res == 'TYR':
        return {'CA','CB','CG','CD1','CE1','CZ','CE2','CD2','OH','HH'} 
    if res == 'ASN':
        return {'CA','CB','CG','OD1','ND2','HD21','HD22'}
    if res == 'GLN':
        return {'CA','CB','CG','CD','OE1','NE2','HE21','HE22'}
    if res == 'ASP':
        return {'CA','CB','CG','OD1','OD2'}
    if res == 'GLU':
        return {'CA','CB','CG','CD','OE1','OE2'}
    if res == 'LYS':
        return {'CA','CB','CG','CD','CE','NZ','HZ1','HZ2','HZ3'}
    if res == 'ARG':
        return {'CA','CB','CG','CD','NE','HE','CZ','NH1','HH11','HH12','NH2','HH21','HH22'}
    if res == 'HIS' or res == 'HIE':
        return {'CA','CB','CG','ND1','CE1','NE2','HE2','CD2'} 

def generate_flex_res(resname, resid, coor, chain='A'):
    # coor is a dict like {'CA': [0.000, 1.000, 2.000]}
    if resid is None:
        resid = 9999
    if len(chain) > 1:
        chain = chain[0]

    # Based on the residue, determine the flex docking part
    if 0:
        pass
    elif resname == 'VAL': # Valine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG1 = f"{coor['CG1'][0]:>8s}{coor['CG1'][1]:>8s}{coor['CG1'][2]:>8s}"
        CG2 = f"{coor['CG2'][0]:>8s}{coor['CG2'][1]:>8s}{coor['CG2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
ATOM      3  CG1 {resname:3s} {chain}{resid:4d}    {CG1}  0.00  0.00    +0.000 C 
ATOM      4  CG2 {resname:3s} {chain}{resid:4d}    {CG2}  0.00  0.00    +0.000 C 
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''
    
    elif resname == 'LEU': # Leucine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD1 = f"{coor['CD1'][0]:>8s}{coor['CD1'][1]:>8s}{coor['CD1'][2]:>8s}"
        CD2 = f"{coor['CD2'][0]:>8s}{coor['CD2'][1]:>8s}{coor['CD2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C 
ATOM      4  CD1 {resname:3s} {chain}{resid:4d}    {CD1}  0.00  0.00    +0.000 C 
ATOM      5  CD2 {resname:3s} {chain}{resid:4d}    {CD2}  0.00  0.00    +0.000 C 
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'ILE': # Isoleucine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG2 = f"{coor['CG2'][0]:>8s}{coor['CG2'][1]:>8s}{coor['CG2'][2]:>8s}"
        CG1 = f"{coor['CG1'][0]:>8s}{coor['CG1'][1]:>8s}{coor['CG1'][2]:>8s}"
        CD1 = f"{coor['CD1'][0]:>8s}{coor['CD1'][1]:>8s}{coor['CD1'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
ATOM      3  CG2 {resname:3s} {chain}{resid:4d}    {CG2}  0.00  0.00    +0.000 C 
BRANCH   2   4
ATOM      4  CG1 {resname:3s} {chain}{resid:4d}    {CG1}  0.00  0.00    +0.000 C 
ATOM      5  CD1 {resname:3s} {chain}{resid:4d}    {CD1}  0.00  0.00    +0.000 C 
ENDBRANCH   2   4
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'MET': # Methionine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        SD = f"{coor['SD'][0]:>8s}{coor['SD'][1]:>8s}{coor['SD'][2]:>8s}"
        CE = f"{coor['CE'][0]:>8s}{coor['CE'][1]:>8s}{coor['CE'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C 
BRANCH   3   4
ATOM      4  SD  {resname:3s} {chain}{resid:4d}    {SD}  0.00  0.00    +0.000 S 
ATOM      5  CE  {resname:3s} {chain}{resid:4d}    {CE}  0.00  0.00    +0.000 C 
ENDBRANCH   3   4
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'PHE': # Phenylalanine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD1 = f"{coor['CD1'][0]:>8s}{coor['CD1'][1]:>8s}{coor['CD1'][2]:>8s}"
        CE1 = f"{coor['CE1'][0]:>8s}{coor['CE1'][1]:>8s}{coor['CE1'][2]:>8s}"
        CZ = f"{coor['CZ'][0]:>8s}{coor['CZ'][1]:>8s}{coor['CZ'][2]:>8s}"
        CE2 = f"{coor['CE2'][0]:>8s}{coor['CE2'][1]:>8s}{coor['CE2'][2]:>8s}"
        CD2 = f"{coor['CD2'][0]:>8s}{coor['CD2'][1]:>8s}{coor['CD2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 A
ATOM      4  CD1 {resname:3s} {chain}{resid:4d}    {CD1}  0.00  0.00    +0.000 A
ATOM      5  CE1 {resname:3s} {chain}{resid:4d}    {CE1}  0.00  0.00    +0.000 A 
ATOM      6  CZ  {resname:3s} {chain}{resid:4d}    {CZ}  0.00  0.00    +0.000 A 
ATOM      7  CE2 {resname:3s} {chain}{resid:4d}    {CE2}  0.00  0.00    +0.000 A 
ATOM      8  CD2 {resname:3s} {chain}{resid:4d}    {CD2}  0.00  0.00    +0.000 A 
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'TRP': # Tryptophan
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD1 = f"{coor['CD1'][0]:>8s}{coor['CD1'][1]:>8s}{coor['CD1'][2]:>8s}"
        NE1 = f"{coor['NE1'][0]:>8s}{coor['NE1'][1]:>8s}{coor['NE1'][2]:>8s}"
        HE1 = f"{coor['HE1'][0]:>8s}{coor['HE1'][1]:>8s}{coor['HE1'][2]:>8s}"
        CE2 = f"{coor['CE2'][0]:>8s}{coor['CE2'][1]:>8s}{coor['CE2'][2]:>8s}"
        CZ2 = f"{coor['CZ2'][0]:>8s}{coor['CZ2'][1]:>8s}{coor['CZ2'][2]:>8s}"
        CH2 = f"{coor['CH2'][0]:>8s}{coor['CH2'][1]:>8s}{coor['CH2'][2]:>8s}"
        CZ3 = f"{coor['CZ3'][0]:>8s}{coor['CZ3'][1]:>8s}{coor['CZ3'][2]:>8s}"
        CE3 = f"{coor['CE3'][0]:>8s}{coor['CE3'][1]:>8s}{coor['CE3'][2]:>8s}"
        CD2 = f"{coor['CD2'][0]:>8s}{coor['CD2'][1]:>8s}{coor['CD2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 A
ATOM      4  CD1 {resname:3s} {chain}{resid:4d}    {CD1}  0.00  0.00    +0.000 A
ATOM      5  NE1 {resname:3s} {chain}{resid:4d}    {NE1}  0.00  0.00    +0.000 N
ATOM      6  HE1 {resname:3s} {chain}{resid:4d}    {HE1}  0.00  0.00    +0.000 HD 
ATOM      7  CE2 {resname:3s} {chain}{resid:4d}    {CE2}  0.00  0.00    +0.000 A 
ATOM      8  CZ2 {resname:3s} {chain}{resid:4d}    {CZ2}  0.00  0.00    +0.000 A 
ATOM      9  CH2 {resname:3s} {chain}{resid:4d}    {CH2}  0.00  0.00    +0.000 A 
ATOM     10  CZ3 {resname:3s} {chain}{resid:4d}    {CZ3}  0.00  0.00    +0.000 A 
ATOM     11  CE3 {resname:3s} {chain}{resid:4d}    {CE3}  0.00  0.00    +0.000 A 
ATOM     12  CD2 {resname:3s} {chain}{resid:4d}    {CD2}  0.00  0.00    +0.000 A 
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'SER': # Serine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        OG = f"{coor['OG'][0]:>8s}{coor['OG'][1]:>8s}{coor['OG'][2]:>8s}"
        HG = f"{coor['HG'][0]:>8s}{coor['HG'][1]:>8s}{coor['HG'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  OG  {resname:3s} {chain}{resid:4d}    {OG}  0.00  0.00    +0.000 OA
ATOM      4  HG  {resname:3s} {chain}{resid:4d}    {HG}  0.00  0.00    +0.000 HD
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'THR': # Threonine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG2 = f"{coor['CG2'][0]:>8s}{coor['CG2'][1]:>8s}{coor['CG2'][2]:>8s}"
        OG1 = f"{coor['OG1'][0]:>8s}{coor['OG1'][1]:>8s}{coor['OG1'][2]:>8s}"
        HG1 = f"{coor['HG1'][0]:>8s}{coor['HG1'][1]:>8s}{coor['HG1'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
ATOM      3  CG2 {resname:3s} {chain}{resid:4d}    {CG2}  0.00  0.00    +0.000 C 
BRANCH   2   4
ATOM      4  OG1 {resname:3s} {chain}{resid:4d}    {OG1}  0.00  0.00    +0.000 OA
ATOM      5  HG1 {resname:3s} {chain}{resid:4d}    {HG1}  0.00  0.00    +0.000 HD
ENDBRANCH   2   4
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'CYS': # Cystine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        SG = f"{coor['SG'][0]:>8s}{coor['SG'][1]:>8s}{coor['SG'][2]:>8s}"
        HG = f"{coor['HG'][0]:>8s}{coor['HG'][1]:>8s}{coor['HG'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  SG  {resname:3s} {chain}{resid:4d}    {SG}  0.00  0.00    +0.000 S
ATOM      4  HG  {resname:3s} {chain}{resid:4d}    {HG}  0.00  0.00    +0.000 HD
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''


    elif resname == 'TYR': # Tyrosine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD1 = f"{coor['CD1'][0]:>8s}{coor['CD1'][1]:>8s}{coor['CD1'][2]:>8s}"
        CE1 = f"{coor['CE1'][0]:>8s}{coor['CE1'][1]:>8s}{coor['CE1'][2]:>8s}"
        CZ = f"{coor['CZ'][0]:>8s}{coor['CZ'][1]:>8s}{coor['CZ'][2]:>8s}"
        CE2 = f"{coor['CE2'][0]:>8s}{coor['CE2'][1]:>8s}{coor['CE2'][2]:>8s}"
        CD2 = f"{coor['CD2'][0]:>8s}{coor['CD2'][1]:>8s}{coor['CD2'][2]:>8s}"
        OH = f"{coor['OH'][0]:>8s}{coor['OH'][1]:>8s}{coor['OH'][2]:>8s}"
        HH = f"{coor['HH'][0]:>8s}{coor['HH'][1]:>8s}{coor['HH'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 A
ATOM      4  CD1 {resname:3s} {chain}{resid:4d}    {CD1}  0.00  0.00    +0.000 A
ATOM      5  CE1 {resname:3s} {chain}{resid:4d}    {CE1}  0.00  0.00    +0.000 A 
ATOM      6  CZ  {resname:3s} {chain}{resid:4d}    {CZ}  0.00  0.00    +0.000 A 
ATOM      7  CE2 {resname:3s} {chain}{resid:4d}    {CE2}  0.00  0.00    +0.000 A 
ATOM      8  CD2 {resname:3s} {chain}{resid:4d}    {CD2}  0.00  0.00    +0.000 A
BRANCH   6   9 
ATOM      9  OH  {resname:3s} {chain}{resid:4d}    {OH}  0.00  0.00    +0.000 OA
ATOM     10  HH  {resname:3s} {chain}{resid:4d}    {HH}  0.00  0.00    +0.000 HD
ENDBRANCH   6   9
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    
    elif resname == 'ASN': # Asparginine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        OD1 = f"{coor['OD1'][0]:>8s}{coor['OD1'][1]:>8s}{coor['OD1'][2]:>8s}"
        ND2 = f"{coor['ND2'][0]:>8s}{coor['ND2'][1]:>8s}{coor['ND2'][2]:>8s}"
        HD21 = f"{coor['HD21'][0]:>8s}{coor['HD21'][1]:>8s}{coor['HD21'][2]:>8s}"
        HD22 = f"{coor['HD22'][0]:>8s}{coor['HD22'][1]:>8s}{coor['HD22'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C 
ATOM      4  OD1 {resname:3s} {chain}{resid:4d}    {OD1}  0.00  0.00    +0.000 OA
BRANCH   3   5
ATOM      5  ND2 {resname:3s} {chain}{resid:4d}    {ND2}  0.00  0.00    +0.000 N
ATOM      6  HD21{resname:3s} {chain}{resid:4d}    {HD21}  0.00  0.00    +0.000 HD
ATOM      7  HD22{resname:3s} {chain}{resid:4d}    {HD22}  0.00  0.00    +0.000 HD
ENDBRANCH   3   5
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'GLN': # Glutamine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD = f"{coor['CD'][0]:>8s}{coor['CD'][1]:>8s}{coor['CD'][2]:>8s}"
        OE1 = f"{coor['OE1'][0]:>8s}{coor['OE1'][1]:>8s}{coor['OE1'][2]:>8s}"
        NE2 = f"{coor['NE2'][0]:>8s}{coor['NE2'][1]:>8s}{coor['NE2'][2]:>8s}"
        HE21 = f"{coor['HE21'][0]:>8s}{coor['HE21'][1]:>8s}{coor['HE21'][2]:>8s}"
        HE22 = f"{coor['HE22'][0]:>8s}{coor['HE22'][1]:>8s}{coor['HE22'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C
BRANCH   3   4 
ATOM      4  CD  {resname:3s} {chain}{resid:4d}    {CD}  0.00  0.00    +0.000 C 
ATOM      5  OE1 {resname:3s} {chain}{resid:4d}    {OE1}  0.00  0.00    +0.000 OA
BRANCH   4   6
ATOM      6  NE2 {resname:3s} {chain}{resid:4d}    {NE2}  0.00  0.00    +0.000 N
ATOM      7  HE21{resname:3s} {chain}{resid:4d}    {HE21}  0.00  0.00    +0.000 HD
ATOM      8  HE22{resname:3s} {chain}{resid:4d}    {HE22}  0.00  0.00    +0.000 HD
ENDBRANCH   4   6
ENDBRANCH   3   4
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'ASP': # Aspartate
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        OD1 = f"{coor['OD1'][0]:>8s}{coor['OD1'][1]:>8s}{coor['OD1'][2]:>8s}"
        OD2 = f"{coor['OD2'][0]:>8s}{coor['OD2'][1]:>8s}{coor['OD2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C 
ATOM      4  OD1 {resname:3s} {chain}{resid:4d}    {OD1}  0.00  0.00    +0.000 OA
ATOM      5  OD2 {resname:3s} {chain}{resid:4d}    {OD2}  0.00  0.00    +0.000 OA
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'GLU': # Glutamate
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD = f"{coor['CD'][0]:>8s}{coor['CD'][1]:>8s}{coor['CD'][2]:>8s}"
        OE1 = f"{coor['OE1'][0]:>8s}{coor['OE1'][1]:>8s}{coor['OE1'][2]:>8s}"
        OE2 = f"{coor['OE2'][0]:>8s}{coor['OE2'][1]:>8s}{coor['OE2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C
BRANCH   3   4 
ATOM      4  CD  {resname:3s} {chain}{resid:4d}    {CD}  0.00  0.00    +0.000 C 
ATOM      5  OE1 {resname:3s} {chain}{resid:4d}    {OE1}  0.00  0.00    +0.000 OA
ATOM      6  OE2 {resname:3s} {chain}{resid:4d}    {OE2}  0.00  0.00    +0.000 OA
ENDBRANCH   3   4
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'LYS': # Lysine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD = f"{coor['CD'][0]:>8s}{coor['CD'][1]:>8s}{coor['CD'][2]:>8s}"
        CE = f"{coor['CE'][0]:>8s}{coor['CE'][1]:>8s}{coor['CE'][2]:>8s}"
        NZ = f"{coor['NZ'][0]:>8s}{coor['NZ'][1]:>8s}{coor['NZ'][2]:>8s}"
        HZ1 = f"{coor['HZ1'][0]:>8s}{coor['HZ1'][1]:>8s}{coor['HZ1'][2]:>8s}"
        HZ2 = f"{coor['HZ2'][0]:>8s}{coor['HZ2'][1]:>8s}{coor['HZ2'][2]:>8s}"
        HZ3 = f"{coor['HZ3'][0]:>8s}{coor['HZ3'][1]:>8s}{coor['HZ3'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C 
BRANCH   3   4
ATOM      4  CD  {resname:3s} {chain}{resid:4d}    {CD}  0.00  0.00    +0.000 C 
BRANCH   4   5
ATOM      5  CE  {resname:3s} {chain}{resid:4d}    {CE}  0.00  0.00    +0.000 C 
BRANCH   5   6
ATOM      6  NZ  {resname:3s} {chain}{resid:4d}    {NZ}  0.00  0.00    +0.000 N
ATOM      7  HZ1 {resname:3s} {chain}{resid:4d}    {HZ1}  0.00  0.00    +0.000 HD
ATOM      8  HZ2 {resname:3s} {chain}{resid:4d}    {HZ2}  0.00  0.00    +0.000 HD
ATOM      9  HZ3 {resname:3s} {chain}{resid:4d}    {HZ3}  0.00  0.00    +0.000 HD
ENDBRANCH   5   6
ENDBRANCH   4   5
ENDBRANCH   3   4
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'ARG': # Arginine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        CD = f"{coor['CD'][0]:>8s}{coor['CD'][1]:>8s}{coor['CD'][2]:>8s}"
        NE = f"{coor['NE'][0]:>8s}{coor['NE'][1]:>8s}{coor['NE'][2]:>8s}"
        HE = f"{coor['HE'][0]:>8s}{coor['HE'][1]:>8s}{coor['HE'][2]:>8s}"
        CZ = f"{coor['CZ'][0]:>8s}{coor['CZ'][1]:>8s}{coor['CZ'][2]:>8s}"
        NH1 = f"{coor['NH1'][0]:>8s}{coor['NH1'][1]:>8s}{coor['NH1'][2]:>8s}"
        HH11 = f"{coor['HH11'][0]:>8s}{coor['HH11'][1]:>8s}{coor['HH11'][2]:>8s}"
        HH12 = f"{coor['HH12'][0]:>8s}{coor['HH12'][1]:>8s}{coor['HH12'][2]:>8s}"
        NH2 = f"{coor['NH2'][0]:>8s}{coor['NH2'][1]:>8s}{coor['NH2'][2]:>8s}"
        HH21 = f"{coor['HH21'][0]:>8s}{coor['HH21'][1]:>8s}{coor['HH21'][2]:>8s}"
        HH22 = f"{coor['HH22'][0]:>8s}{coor['HH22'][1]:>8s}{coor['HH22'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 C 
BRANCH   3   4
ATOM      4  CD  {resname:3s} {chain}{resid:4d}    {CD}  0.00  0.00    +0.000 C 
BRANCH   4   5
ATOM      5  NE  {resname:3s} {chain}{resid:4d}    {NE}  0.00  0.00    +0.000 N 
ATOM      6  HE  {resname:3s} {chain}{resid:4d}    {HE}  0.00  0.00    +0.000 HD
ATOM      7  CZ  {resname:3s} {chain}{resid:4d}    {CZ}  0.00  0.00    +0.000 C
ATOM      8  NH1 {resname:3s} {chain}{resid:4d}    {NH1}  0.00  0.00    +0.000 N
ATOM      9  HH11{resname:3s} {chain}{resid:4d}    {HH11}  0.00  0.00    +0.000 HD
ATOM     10  HH12{resname:3s} {chain}{resid:4d}    {HH12}  0.00  0.00    +0.000 HD
ATOM     11  NH2 {resname:3s} {chain}{resid:4d}    {NH2}  0.00  0.00    +0.000 N
ATOM     12  HH21{resname:3s} {chain}{resid:4d}    {HH21}  0.00  0.00    +0.000 HD
ATOM     13  HH22{resname:3s} {chain}{resid:4d}    {HH22}  0.00  0.00    +0.000 HD
ENDBRANCH   4   5
ENDBRANCH   3   4
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''

    elif resname == 'HIS' or resname == 'HIE': # Histidine
        assert set(coor.keys()) == atom_dict(resname)
        CA = f"{coor['CA'][0]:>8s}{coor['CA'][1]:>8s}{coor['CA'][2]:>8s}"
        CB = f"{coor['CB'][0]:>8s}{coor['CB'][1]:>8s}{coor['CB'][2]:>8s}"
        CG = f"{coor['CG'][0]:>8s}{coor['CG'][1]:>8s}{coor['CG'][2]:>8s}"
        ND1 = f"{coor['ND1'][0]:>8s}{coor['ND1'][1]:>8s}{coor['ND1'][2]:>8s}"
        CE1 = f"{coor['CE1'][0]:>8s}{coor['CE1'][1]:>8s}{coor['CE1'][2]:>8s}"
        NE2 = f"{coor['NE2'][0]:>8s}{coor['NE2'][1]:>8s}{coor['NE2'][2]:>8s}"
        HE2 = f"{coor['HE2'][0]:>8s}{coor['HE2'][1]:>8s}{coor['HE2'][2]:>8s}"
        CD2 = f"{coor['CD2'][0]:>8s}{coor['CD2'][1]:>8s}{coor['CD2'][2]:>8s}"
        res_output = f'''BEGIN_RES {resname:3s} {chain}{resid:4d}
ROOT
ATOM      1  CA  {resname:3s} {chain}{resid:4d}    {CA}  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  CB  {resname:3s} {chain}{resid:4d}    {CB}  0.00  0.00    +0.000 C 
BRANCH   2   3
ATOM      3  CG  {resname:3s} {chain}{resid:4d}    {CG}  0.00  0.00    +0.000 A
ATOM      4  ND1 {resname:3s} {chain}{resid:4d}    {ND1}  0.00  0.00    +0.000 NA
ATOM      5  CE1 {resname:3s} {chain}{resid:4d}    {CE1}  0.00  0.00    +0.000 A 
ATOM      6  NE2 {resname:3s} {chain}{resid:4d}    {NE2}  0.00  0.00    +0.000 A 
ATOM      7  HE2 {resname:3s} {chain}{resid:4d}    {HE2}  0.00  0.00    +0.000 HD 
ATOM      8  CD2 {resname:3s} {chain}{resid:4d}    {CD2}  0.00  0.00    +0.000 A
ENDBRANCH   2   3
ENDBRANCH   1   2
END_RES {resname:3s} {chain}{resid:4d}
'''


    return res_output


class PDBQT:
    def __init__(self, pdbqt_file):
        self.pdbqt_file = pdbqt_file
        self.rigid_file = pdbqt_file.replace('.pdbqt','_rigid.pdbqt')
        self.flex_file = pdbqt_file.replace('.pdbqt','_flex.pdbqt')
        with open(self.pdbqt_file, 'r') as f:
            self.pdbqt = f.readlines()
        self.pdbqt = parse_PDBQT(self.pdbqt)

        self.xyz = [x[6:9] for x in self.pdbqt if x[0] == 'ATOM']
        self.xyz = np.array([list(map(float, x)) for x in self.xyz])
        self.resid = np.array([int(x[12]) for x in self.pdbqt if x[0] == 'ATOM'])
        self.non_flex_resid = np.unique([int(x[12]) for x in self.pdbqt if x[3] in ['GLY','ALA','PRO']])
        print(f'These residues are not considered flexible: {self.non_flex_resid}')


    def report_flex_res(self, center, n_flex=3, n_consider=6):
        # Pick a center, and report n nearest residues to be flexible, considering a larger number of residues
        if n_flex > n_consider:
            n_consider = n_flex
        center = np.array(center) 
        dist_matrix = np.sqrt(((self.xyz - center)**2).sum(1))
        print(dist_matrix.shape)
        close_res = []

        thres = 0
        while len(close_res) < n_consider:
            thres += 0.02
            close_atoms = np.where(dist_matrix < thres)
            close_dist = {}
            for atom in close_atoms[0]:
                try:
                    if dist_matrix[atom] < close_dist[self.resid[atom]]:
                        close_dist[self.resid[atom]] = dist_matrix[atom]
                except:
                    close_dist[self.resid[atom]] = dist_matrix[atom]
            close_res = np.unique(self.resid[close_atoms])
            close_res = [x for x in close_res if x not in self.non_flex_resid]
            #print(f'{thres:.3f}, {close_res}, {close_dist}')

        print('We will consider these residues:')
        res_list = []
        for res in close_res:
            for atom in self.pdbqt:
                if atom[2] == 'CA' and int(atom[12]) == res: 
                    res_list.append(atom)
                    print(atom)

        # First pick charged residues, then polar, then nonpolar. Within each, sort by distance
        priority = np.zeros(len(res_list))

        for idx, res in enumerate(res_list):
            if res[3] in ['ASP','GLU','LYS','ARG','HIS','HIE']:
                priority[idx] = 8 - close_dist[int(res[12])]
            elif res[3] in ['SER','THR','CYS','TYR','ASN','GLN']:
                priority[idx] = 5 - close_dist[int(res[12])]
            elif res[3] in ['VAL','LEU','ILE','MET','PHE','TRP']:
                priority[idx] = -close_dist[int(res[12])]
            else:
                priority[idx] = -np.inf

        rank = np.argsort(priority)[::-1]
        print('We will use these residues:')
        picked_res = []
        for ii in rank[:n_flex]:
            print(f'{priority[ii]:.3f}, {res_list[ii]}')
            picked_res.append(res_list[ii][12])

        return picked_res
 
    def generate_rigid_flex_files(self, center, rigid_file=None, flex_file=None, n_consider=6, n_flex=3):
        if rigid_file is not None:
            self.rigid_file = rigid_file
        if flex_file is not None:
            self.flex_file = flex_file

        picked_res = self.report_flex_res(center=center, n_consider=n_consider, n_flex=n_flex)

        rigid = []
        flex = ''  # File content
        flex_coor = {}
        flex_resdata = {}
        for res in picked_res:
            flex_coor[res] = {}
            flex_resdata[res] = {}
        

        for atom in self.pdbqt:
            if (int(atom[12]) not in picked_res) or (atom[2] in ['N','H','C','O','H1','H2','H3','OXT']):
                rigid.append(atom)
            else:
                flex_coor[atom[12]][atom[2]] = [atom[6], atom[7], atom[8]]
                print(f'flex_coor[{atom[12]}][{atom[2]}] = [{atom[6]}, {atom[7]}, {atom[8]}]')
                flex_resdata[atom[12]]['resname'] = atom[3]
                flex_resdata[atom[12]]['chain'] = atom[4]

        for res in np.sort(picked_res):
            print(flex_coor[res])
            print(atom_dict(flex_resdata[res]['resname']))
            flex += generate_flex_res(flex_resdata[res]['resname'], res, flex_coor[res], flex_resdata[res]['chain'])
        counter = 1
        for atom in rigid: 
            atom[1] = counter
            atom[5] = atom[12]
            counter += 1
         
        with open(self.rigid_file,'w') as f:
            for l in rigid:
                f.write(atom_line(l[5], l, l[1]))
        with open(self.flex_file,'w') as f:
            f.write(flex)
#        print(rigid)
#        print(flex)

def atom_line(res, atom, counter): # Atom is a line in self.target_pdb
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
    return f'{atom[0]:<6s}{str(counter):>5s} {a2:<4s} {a3:>4s}{atom[4]:1s}{str(res):>4s}    {atom[6]:>8s}{atom[7]:>8s}{atom[8]:>8s}{atom[9]:>6s}{atom[10]:>6s}{atom[11]}'


def parse_PDBQT(PDB_content):
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
        atom.append(counter)                      # 12. reindexed resid (1-based)

    return splitPDB



pdbqt = PDBQT(sys.argv[1])
print(pdbqt)
#picked_res = pdbqt.report_flex_res([10.967,-32.321,-10.801], n_consider=6, n_flex=3)
#print(picked_res)
pdbqt.generate_rigid_flex_files([float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])], n_consider=6, n_flex=3)


