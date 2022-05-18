import os
import sys
import glob
import json
import copy
import mdtraj
import numpy as np
import time
import pandas as pd
import pickle
import mdtraj as md
import multiprocessing as mp
import re
from calcWater_func import *
from scipy.spatial import distance_matrix
import networkx as nx
from spyrmsd import rmsd as srmsd

try:
    import cupy as cp
    cudaExists = True
    import kernel
except ImportError as e:
    cudaExists = False
    print("Can't load CuPy, contact fingerprint will not run")
    
# sys.path.insert(1, os.path.join(sys.path[0], '../test_contacts/contacts/contacts/'))
import importlib

ligand_map = {'A': 0, 'C': 0, 'N': 1, 'NA': 1, 'O': 2, 'OA': 2, 'F': 3, 'P': 4, 'S': 5, 'SA': 5, 'CL': 6,
              'BR': 7, 'I': 8, 'H': 9}
protein_map = {'A': 0, 'C': 0, 'N': 1, 'NA': 1, 'O': 2, 'OA': 2, 'S': 3, 'SA': 3, 'H': 4}

vdW_map = {'C': 1.9, 'O': 1.7, 'N': 1.8, 'P': 2.1, 'S': 2.0, 'F': 1.5, 'Cl': 1.8, 'Br': 2.0, 'I': 2.2, 'H': 0.0}
metals = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ha', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'PB', 'Bi', 'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Lv'] # These have vdW radii = 1.2 A

# import MDAnalysis as mda
# import MDAnalysis.analysis.rms


def ALIGN_A_RMSD_B(P, Q, A, B, lgc=None, return_both_kinds=False): #lgc is ligandGraphConfig

    # P is the one to be aligned (N * 3)
    # Q is the ref (N * 3)
    # A is the list of index to be considered for alignment (protein) (N * 1)
    # B is the list of index to calculate RMSD (ligand) (N * 1)
    # Returns rmsd between subset P[B] and Q[B]
    PU = P[A] # Get subset
    QU = Q[A] # Get subset
    PC = PU - PU.mean(axis=0) # Center points
    QC = QU - QU.mean(axis=0) # Center points
    # Kabsch method
    C = np.dot(np.transpose(PC), QC)
    V, S, W = np.linalg.svd(C,full_matrices=False)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    P = P - PU.mean(axis=0) # Move all points
    Q = Q - QU.mean(axis=0) # Move all points
    P = np.dot(P, U) # Rotate P

    # TODO: Use spyrmsd to compute symmetry-corrected RMSD!
    if return_both_kinds and (lgc is not None):
        lgc_rmsd = srmsd.symmrmsd(P[B], Q[B], lgc['an'], lgc['an'], lgc['G'], lgc['G'])
        diff = P[B] - Q[B]
        N = len(P[B])
        raw_rmsd = np.sqrt((diff * diff).sum() / N)
        return lgc_rmsd, P + QU.mean(axis=0), raw_rmsd
    elif lgc is not None:
        lgc_rmsd = srmsd.symmrmsd(P[B], Q[B], lgc['an'], lgc['an'], lgc['G'], lgc['G'])
        return lgc_rmsd, P + QU.mean(axis=0)
    else:
        # Return raw rmsd
        diff = P[B] - Q[B]
        N = len(P[B])
        raw_rmsd = np.sqrt((diff * diff).sum() / N)
        return raw_rmsd, P + QU.mean(axis=0) 

def pureRMSD(P, Q):
    # Assume P and Q are aligned first
    diff = P - Q
    N = len(P)
    return np.sqrt((diff * diff).sum() / N)


def ReadLog(fname):
    with open(fname,'r') as f:
        cont = f.readlines()

    Step = []
    Time = []
    Temperature = []
    Etot = []
    EKtot = []
    EPtot = []
    Bond = []
    Angle = []
    Dihedral = []
    Elec = []
    vdW = []
    Solvent = []

    start_print = False
    for line in cont:
        if '(6.) RESULTS' in line:
            start_print = True
        if 'Averages over' in line:
            start_print = False
        if start_print == True:
            if 'Step' in line:
                ele = line.split()
                Step.append(int(ele[1]))
                Time.append(float(ele[3]))
            if 'Temperature' in line:
                Temperature.append(float(line.split(':')[1].strip()))
            if 'Etot' in line:
                ele = line.split()
                Etot.append(float(ele[1]))
                EKtot.append(float(ele[3]))
                EPtot.append(float(ele[5]))
            if 'Bond' in line:
                ele = line.split()
                Bond.append(float(ele[1]))
                Angle.append(float(ele[3]))
                Dihedral.append(float(ele[5]))
            if 'Elec' in line:
                ele = line.split()
                Elec.append(float(ele[1]))
                vdW.append(float(ele[3]))
                Solvent.append(float(ele[5]))
    #         print(line.strip('\n'))
    # output = np.array([Step, Time, Temperature, Etot, EKtot, EPtot, Bond, Angle, Dihedral, Elec, vdW, Solvent]).T
    output = np.array([Step, Time, Temperature, Etot, EKtot, EPtot, Bond, Angle, Dihedral, Elec, vdW, Solvent]).T[1:]
    output = pd.DataFrame(output, columns=['Step', 'Time', 'Temperature', 'Etot', 'EKtot', 'EPtot', 'Bond', 'Angle', 'Dihedral', 'Elec', 'vdW', 'Solvent'])
    return output


def Distribute_Lig_pose_RMSD(arg):
    obj, crystalComp, templateLigandGrid, ligandAtomConfig, ligandGraphConfig = arg
    obj.calculateRMSD(crystalComp, templateLigandGrid, ligandAtomConfig, ligandGraphConfig)
    return obj

def readMDCRD(fname, natom):
    with open(fname, 'r') as f:
        f.readline()
        cont = f.read()
    xyz = list(map(float,cont.split()))
    return np.array(xyz).reshape(-1, natom, 3)

def minmax(xyz):
    return np.min(xyz, axis=0) - 3, np.max(xyz, axis=0) + 3
def grid(xyz, width=0.2):
    mm = minmax(xyz)
    ming = mm[0] // width * width
    maxg = mm[1] // width * width + width * 2
    gx_ = np.arange(ming[0], maxg[0], width)
    gy_ = np.arange(ming[1], maxg[1], width)
    gz_ = np.arange(ming[2], maxg[2], width)
    
    gx, gy, gz = np.meshgrid(gx_, gy_, gz_)
    gxyz = np.array([gx.flatten(), gy.flatten(), gz.flatten()]).T
    return gxyz

def create_ref_grid(xyz, vdW=1.8, width=0.2):
    assert type(vdW) == float or len(vdW) == len(xyz), "vdW must be a single value or a list of the shape xyz"
    crude_grid = grid(xyz, width)
    if type(vdW) == float:
        dist = distance_matrix(crude_grid, xyz) - vdW
    elif len(vdW) == len(xyz):
        dist = distance_matrix(crude_grid, xyz) - np.array(vdW)
    fine_grid = crude_grid[np.any(dist < 0, axis=1)]
    return fine_grid
def overlap_mol(test_xyz, ref_xyz=None, test_vdW=1.8, ref_vdW=1.8, ref_grid=None, ref_width=0.2):
    assert ref_xyz is not None or ref_grid is not None
    assert type(test_vdW) == float or len(test_vdW) == len(test_xyz)
    assert type(ref_vdW) == float or len(ref_vdW) == len(ref_xyz)
    if ref_grid is None:
        ref_grid = create_ref_grid(ref_xyz, ref_vdW, ref_width)
    if type(test_vdW) == float:
        dist = distance_matrix(ref_grid, test_xyz) - test_vdW
    elif len(test_vdW) == len(test_xyz):
        print(test_vdW)
        dist = distance_matrix(ref_grid, test_xyz) - np.array(test_vdW)
    ratio = np.sum(np.any(dist < 0, axis=1)) / len(dist)
    return ratio, ref_grid


class mdgxTrajectory:
    def __init__(self, trjFile, outFile, rstFile, ioutfm):
        self.trjFile = trjFile
        self.outFile = outFile
        self.rstFile = rstFile
        self.ioutfm = ioutfm
        self.success = False
        self.hasRMSD = False
        self.hasProteinRMSD = False
        self.hasProteinCoreRMSD = False
        self.hasTemporalRMSD = False
        self.hasTemplateOverlap = False
        self.hasTrjFile = False
        self.hasOutFile = False
        self.output = None
        self.RMSD = []
        self.proteinRMSD = []
        self.proteinCoreRMSD = []
        self.temporalRMSD = []
        self.hasSASA = False
        self.hasHBHC = False
        self.hasBridgeWaterHB = False
        self.contactScore = []
        self.hasContactScore = False
        self.ligandTrajectory = None # We store this information for contact and clustering calculation
        self.ligandTrajectoryH = None # We store this information for contact and clustering calculation
        self.hasLigandTrajectory = False
#         self.trjLength = 0
#         self.rmsd = []
    
    def calculateRMSD(self, prmtop, crystalComp, lig_len, ligand_res, dockCenter, templateLigandGrid=None, 
                      ligandAtomConfig=None, ligandGraphConfig=None): 
        # This function also logs the ligand trajectory even if RMSD calculation fails
        try:
            if (not self.hasRMSD):
                Profile = np.random.random() < 0.00
                if Profile:
                    t0 = time.time()

                if self.ioutfm == 1: # Binary netcdf trajectory
                    # comp = mda.Universe(prmtop, self.trjFile)
                    comp = mdtraj.load_netcdf(self.trjFile, top=prmtop)[:-1] # Adjust for mismatch between output and traj
                elif self.ioutfm == 0: # ASCII MDCRD trajectory, loading is slower but partial trajs can be loaded
                    systemLen = crystalComp.n_atoms
                    # comp = readMDCRD(self.trjFile, systemLen)
                    #comp = readMDCRD(self.trjFile, systemLen)[:-1] 
                    comp = mdtraj.load_mdcrd(self.trjFile, top=prmtop)[:-1] # Adjust for mismatch between output and traj
                if Profile:
                    t1 = time.time()
                self.trjLength = len(comp)
    
                if Profile:
                    t2 = time.time()


                R = []   # RMSD
                rR = []  # raw RMSD
                PR = []  # Protein RMSD
                PCR = [] # Protein core RMSD
                LT = []
                systemLen = crystalComp.n_atoms
                referenceXYZ = crystalComp.xyz[0]
       
                ligand_heavy_atoms = crystalComp.top.select(f"residue {ligand_res} and not symbol H")
                protein_heavy_atoms = crystalComp.top.select(f"not residue {ligand_res} and not symbol H")
                core_res = np.where(np.sqrt(((crystalComp.xyz[0][crystalComp.top.select('name CA and protein')]*10 - np.array(dockCenter))**2).sum(1)) < 9)[0]
                if len(core_res) > 0:
                    core_idx = crystalComp.top.select(f'residue {" ".join([str(x) for x in core_res])} and not symbol H')
                else:
                    print("Core residue assignment failed!")
                    core_res = None
                    core_idx = None
                #print(f'ligand heavy atoms index: {ligand_heavy_atoms}')
                for xyz in comp.xyz:
                    R_this, LT_this, rR_this = ALIGN_A_RMSD_B(xyz*10, referenceXYZ*10,
                                            range(0, systemLen-lig_len), ligand_heavy_atoms, ligandGraphConfig, return_both_kinds=True)
                                            #range(0, systemLen-lig_len), range((systemLen-lig_len), systemLen))
                    R.append(R_this)
                    rR.append(rR_this)
                    LT.append(LT_this)
                    PR_this, _ = ALIGN_A_RMSD_B(xyz*10, referenceXYZ*10,
                                                range(0, systemLen-lig_len), protein_heavy_atoms)
                    PR.append(PR_this)
                    if core_idx is not None:
                        PCR_this, _ = ALIGN_A_RMSD_B(xyz*10, referenceXYZ*10,
                                                    range(0, systemLen-lig_len), core_idx)
                        PCR.append(PCR_this)
                        #print("PCR actually worked!")
                self.RMSD = np.array(R)
                self.rawRMSD = np.array(rR)
                self.proteinRMSD = np.array(PR)
                self.proteinCoreRMSD = np.array(PCR)
                LT = np.array(LT)
                
                

                if Profile:
                    t3 = time.time()
                self.output = ReadLog(self.outFile)
                self.getSASA(comp, ligand_res, ligandAtomConfig) 
                self.getHBHC(comp, ligandAtomConfig)
                self.output['proteinRMSD'] = self.proteinRMSD
                self.output['proteinCoreRMSD'] = self.proteinCoreRMSD
                if Profile:
                    t4 = time.time()
                if Profile:
                    print(f'    Profiling: Loading comp {(t1-t0)*1000:.3f} ms | Superpose {(t2-t1)*1000:.3f} ms | RMSD {(t3-t2)*1000:.3f} ms | Read output {(t4-t3)*1000:.3f} ms ')
                self.hasRMSD = True
                self.hasProteinRMSD = True
                self.hasProteinCoreRMSD = True
                self.hasTrjFile = True
                self.hasOutFile = True

                 
                if not self.hasLigandTrajectory: # Also store ligand trajectories for contact analysis

                    lig = crystalComp.top.select(f"residue {ligand_res} and not symbol H")
                    ligH = crystalComp.top.select(f"residue {ligand_res}")
                    pro = crystalComp.top.select(f"not residue {ligand_res} and not symbol H")
                    ligHeavyElementsvdW = np.array([vdW_map[r.element.symbol] if r.element.symbol in vdW_map.keys() else 1.2 if r.element.symbol in metals else 1.8 for r in list(crystalComp.top.atoms)])[lig]
                    self.ligandTrajectory = LT[:,lig]
                    self.ligandTrajectoryH = LT[:, ligH]
                    self.hasLigandTrajectory = True
                    if templateLigandGrid is not None:
                        self.getTemplateOverlap(self.ligandTrajectory, templateLigandGrid)#, ligHeavyElementsvdW)
                    self.getTemporalRMSD(self.ligandTrajectory)
                
                if np.all([self.hasRMSD, self.hasProteinRMSD, self.hasProteinCoreRMSD, self.hasSASA, self.hasHBHC, self.hasTemporalRMSD]):
                    self.success = True

        except:
            pass
    
    def readOutput(self):
        try:
            self.output = ReadLog(self.outFile)
            self.hasOutFile = True
        except:
            pass

    def getSASA(self, mol, ligand_res, lac): #lac = ligandAtomConfig
                # Change of ligand SASA is ligand fraction of SASA in complex - free SASA, not free ligand SASA + apo SASA - holo SASA
                # Also separate them to CSP SASA and Polar SASA, for both protein, ligand and complex
        try:
            lig = mol.atom_slice(mol.top.select(f'resid {ligand_res}'))
            pro = mol.atom_slice(mol.top.select(f'not resid {ligand_res}'))
            
            ligandAtomSASA  = mdtraj.shrake_rupley(lig, mode='atom', n_sphere_points=200)
            proteinAtomSASA = mdtraj.shrake_rupley(pro, mode='atom', n_sphere_points=200)
            complexAtomSASA = mdtraj.shrake_rupley(mol, mode='atom', n_sphere_points=200)
            self.output['apoProteinSASA'] = proteinAtomSASA.sum(1)
            self.output['freeLigandSASA'] = ligandAtomSASA.sum(1) 
            self.output['complexSASA'] = complexAtomSASA.sum(1)
            self.output['allChangeSASA'] = self.output['complexSASA'] - self.output['apoProteinSASA'] - self.output['freeLigandSASA']
            #print('Done with regular sasa')
            #print(complexAtomSASA.shape)
            self.output['ligandInComplexSASA'] = complexAtomSASA[:, lac['lig']].sum(1)
            self.output['proteinInComplexSASA'] = complexAtomSASA[:, lac['pro']].sum(1)
            self.output['ligandChangeSASA'] = self.output['ligandInComplexSASA'] - self.output['freeLigandSASA']
            #print('Done with ligand sasa')
            self.output['proteinChangeSASA'] = self.output['proteinInComplexSASA'] - self.output['apoProteinSASA']
            self.output['ligandCSPSASA'] = complexAtomSASA[:, lac['lig_CSP']].sum(1)
            self.output['proteinCSPSASA'] = complexAtomSASA[:, lac['pro_CSP']].sum(1) 
            self.output['complexCSPSASA'] = self.output['ligandCSPSASA'] + self.output['proteinCSPSASA']
            #print('Done with CSP sasa')
            self.output['ligandCHSASA'] = complexAtomSASA[:, lac['lig_CH']].sum(1)
            self.output['proteinCHSASA'] = complexAtomSASA[:, lac['pro_CH']].sum(1)
            self.output['complexCHSASA'] = self.output['ligandCHSASA'] + self.output['proteinCHSASA']
            #print('Done with CH sasa')
            self.output['ligandPolarSASA'] = complexAtomSASA[:, lac['lig_polar']].sum(1)
            self.output['proteinPolarSASA'] = complexAtomSASA[:, lac['pro_polar']].sum(1)
            self.output['complexPolarSASA'] = self.output['ligandPolarSASA'] + self.output['proteinPolarSASA']
            #print('Done with polar sasa')
            #if coreRes is not None:
                #print(f'core residues are {coreRes}')
            #    self.output['proteinCoreSASA'] = proteinResSASA[:, coreRes].sum(1)
            #    self.output['complexCoreSASA'] = complexResSASA[:, coreRes].sum(1)
            #    self.output['changeCoreSASA'] = self.output['complexCoreSASA'] - self.output['proteinCoreSASA'] - self.output['ligandSASA']
            self.hasSASA = True
        except:
            print("SASA failed!") 


    def HC(self, xyz, C_lig, C_pro, d0=3.8, detail=False):
    
        # xyz is mol.xyz[N]
        dist = distance_matrix(xyz[C_lig], xyz[C_pro])
        # print(dist)
        score = (1/1.5 * (d0 +2.0 - dist))
        score[dist < (d0 + 0.5)] = 1
        score[dist > (d0 + 2)] = 0
        if detail:
            print(C_lig[np.where(score > 0)[0]], C_pro[np.where(score > 0)[1]])
        return np.sum(score)/2
    
    def HB(self, x1, x2, xh=None, use_angle=False):
        if xh is None:
            use_angle = False
        dist_term = (1 / (1 + (np.sqrt(((x1-x2)**2).sum())/2.6)**6) / 0.58)
        if use_angle:
            angle = calcAngle(x1, xh, x2)[0]
            #print(angle)
            angle = min(angle, 180-angle)
            angle_term = max(0, min(1, 1 - (angle - 30) / 50))
            #print(angle, dist_term, angle_term)
            return dist_term * angle_term
        else:
            return dist_term


    def getHBHC(self, mol, lac): # lac = ligandAtomConfig
                # TODO: Use the new definition for CSP contacts
        try:
            AllHB = []
            AllHC = []
            #C_lig = mol.top.select(f'residue {ligand_res} and element C')
            #C_pro = mol.top.select(f'not residue {ligand_res} and element C')
            lig_atoms = lac['lig']
            #lig_atoms = mol.top.select(f'resid {ligand_res}')

            for idx, frame in enumerate(mol):
                score = 0
                hbonds = mdtraj.baker_hubbard(frame, freq=0.1, distance_cutoff=0.25)
                for hbond in hbonds:
                    if (hbond[0] in lig_atoms) and (hbond[2] in lig_atoms): # Intraligand HB doesn't count
                        pass
                    elif (hbond[0] not in lig_atoms) and (hbond[2] not in lig_atoms): # Intraprotein HB doesn't count
                        pass
                    else:
                        #print(mol.xyz[idx][hbond[0]]*10, mol.xyz[idx][hbond[2]]*10, xh=mol.xyz[idx][hbond[1]]*10, use_angle=True)
                        score += self.HB(mol.xyz[idx][hbond[0]]*10, mol.xyz[idx][hbond[2]]*10, xh=mol.xyz[idx][hbond[1]]*10, use_angle=False)
                AllHB.append(score)
                AllHC.append(self.HC(mol.xyz[idx]*10, lac['lig_allCSP'], lac['pro_allCSP']))
            self.output['HBond'] = np.array(AllHB)
            self.output['Contact'] = np.array(AllHC)
            self.hasHBHC = True
        except:
            print('getHBHC failed!')

    def getBridgeWaterHB(self, mol):
        try:
            self.output['BridgeWaterHB'] = calcWater(mol, printResult=False, processFrame=-1, returnHBscore=True)
            self.hasBridgeWaterHB = True
        except:
            print('getBridgeWaterHB failed!')

    def getTemplateOverlap(self, LT, tlg, lhev=1.8): # tlg is templateLigandGrid, lhev is ligHeavyElementsvdW
        try:
            templateOverlap = []
            for xyz in LT:
                templateOverlap.append(overlap_mol(xyz, ref_grid=tlg, test_vdW=lhev)[0])
            #Profile = np.random.random() < 0.01
            #if Profile:
            #    print(f'Profiling of templateOverlap:')
            #    print(templateOverlap)
            self.output['templateOverlap'] = np.array(templateOverlap)
            self.hasTemplateOverlap = True
            #print('getTemplateOverlap succeeded!')
            
        except:
            print('getTemplateOverlap failed!')


    def getTemporalRMSD(self, LT):
        try:
            dist = distance_matrix(LT.reshape(len(LT), -1), LT.reshape(len(LT), -1))
            tRMSD = [np.std(dist[i, max(0, i-5):min(len(dist), i+5)]) for i in range(len(dist))]
            self.output['temporalRMSD'] = np.array(tRMSD)
            self.hasTemporalRMSD = True
        except:
            print('getTemporalRMSD failed!')
 
    def readLigTraj(self, prmtop, initialPose, ligand_res): 
        # This function just logs the ligand trajectory
        # initialPose is for getting the systemLen
        if not self.hasLigandTrajectory:
            if self.ioutfm == 1: # Binary netcdf trajectory
                initialComp = mdtraj.load_netcdf(self.trjFile, top=prmtop)
            elif self.ioutfm == 0: # ASCII MDCRD trajectory, loading is slower
                initialComp = mdtraj.load(self.initialPose, top=prmtop)
                systemLen = initialComp.n_atoms
                comp = readMDCRD(self.trjFile, systemLen)

            lig = initialComp.top.select(f"residue {ligand_res} and not symbol H")
            ligH = initialComp.top.select(f"residue {ligand_res}")
            pro = initialComp.top.select(f"not residue {ligand_res} and not symbol H")
            self.ligandTrajectory = LT[:,lig]
            self.ligandTrajectoryH = LT[:, ligH]
            self.hasLigandTrajectory = True
    
    
        
    
    def inheritMdgxTrajectory(self, TRJ):
        self.trjFile = TRJ.trjFile
        self.outFile = TRJ.outFile
        self.rstFile = TRJ.rstFile
        self.success = TRJ.success
        self.hasRMSD = TRJ.hasRMSD
        self.hasProteinRMSD = TRJ.hasProteinRMSD
        self.hasProteinCoreRMSD = TRJ.hasProteinCoreRMSD
        self.hasTemplateOverlap = TRJ.hasTemplateOverlap
        self.hasTemporalRMSD = TRJ.hasTemporalRMSD
#         if self.hasRMSD:
        self.RMSD = TRJ.RMSD
        try:
            self.rawRMSD = TRJ.rawRMSD
        except:
            pass
        self.proteinRMSD = TRJ.proteinRMSD
        self.proteinCoreRMSD = TRJ.proteinCoreRMSD
        self.temporalRMSD = TRJ.temporalRMSD
        self.hasSASA = TRJ.hasSASA
        self.hasHBHC = TRJ.hasHBHC
        self.hasBridgeWaterHB = TRJ.hasBridgeWaterHB 
        self.ioutfm = TRJ.ioutfm
        self.hasTrjFile = TRJ.hasTrjFile
        self.hasOutFile = TRJ.hasOutFile
        self.output = TRJ.output
        try:
            self.hasContactScore = TRJ.hasContactScore
            self.contactScore = TRJ.contactScore
        except:
            self.hasContactScore = False
            self.contactScore = None
#         try:
        self.ligandTrajectory = TRJ.ligandTrajectory
        self.ligandTrajectoryH = TRJ.ligandTrajectoryH
        self.hasLigandTrajectory = TRJ.hasLigandTrajectory
#         self.proteinCheck = TRJ.proteinCheck
#         except:
#             self.ligandTrajectory = None
#             self.hasLigandTrajectory = False
#         self.RMSD = TRJ.RMSD
        
        
class Pose:
    def __init__(self, name, rank, successQR=False, ligand_res=None, settings=None, 
                 folderMetadata=None, simulationPrefixes=['EM','QR','MD'],
                 dockCenter=None): 
        self.poseName = name
        self.ligandName = name.split('_')[0]
        self.nrep = {}
        self.length = {}
        self.writeCrdInterval = {}
        self.ioutfm = {}
        self.timeStep = {}
        self.successQR = successQR # If QR succeeds all subsequent rounds succeed.
        self.ligand_res = ligand_res
        self.simulationPrefixes = simulationPrefixes
        self.dockCenter = dockCenter
        self.outSuffix = {}
        self.rstSuffix = {}
        self.crdSuffix = {}
        for simPrefix in self.simulationPrefixes:
            self.nrep[simPrefix] = int(settings[simPrefix]['N-rep'])
            self.length[simPrefix] = int(settings[simPrefix]['cntrl']['nstlim'])
            self.writeCrdInterval[simPrefix] = int(settings[simPrefix]['cntrl']['ntwx'])
            self.ioutfm[simPrefix] = int(settings[simPrefix]['cntrl']['ioutfm'])
            self.timeStep[simPrefix] = float(settings[simPrefix]['cntrl']['dt'])
            if 'EX' in simPrefix:
                self.nrep[simPrefix] = self.nrep['MD']
            self.outSuffix[simPrefix] = settings[simPrefix]['files']['-osf']
            self.rstSuffix[simPrefix] = settings[simPrefix]['files']['-rsf']
            self.crdSuffix[simPrefix] = settings[simPrefix]['files']['-xsf']
        self.rootFolder = folderMetadata['rootFolder']
        self.referenceFolder = folderMetadata['referenceFolder']
        self.structureFolder = folderMetadata['structureFolder']
        self.simulationFolder = folderMetadata['simulationFolder']
        self.inpcrdFolder = folderMetadata['inpcrdFolder']
        self.prmtopFolder = folderMetadata['prmtopFolder']

        self.traj = {}
        for simPrefix in self.simulationPrefixes:
            self.traj[simPrefix] = []
            if simPrefix == 'EM' or simPrefix == 'QR':
                pass
            elif not self.successQR:
                continue
            for ii in range(self.nrep[simPrefix]):

                trjFile = f'{self.simulationFolder}/{self.poseName}/{simPrefix}_R{ii+1}{self.crdSuffix[simPrefix]}'
                outFile = f'{self.simulationFolder}/{self.poseName}/{simPrefix}_R{ii+1}{self.outSuffix[simPrefix]}'
                rstFile = f'{self.simulationFolder}/{self.poseName}/{simPrefix}_R{ii+1}{self.rstSuffix[simPrefix]}'
                self.traj[simPrefix].append(mdgxTrajectory(trjFile, outFile, rstFile, self.ioutfm[simPrefix]))
        # Also include initial pose
        self.initialPose = f'{self.inpcrdFolder}/{name}.inpcrd'
        self.prmtop = f'{self.prmtopFolder}/{self.ligandName}.prmtop'

        if not os.path.exists(self.prmtop):
            self.prmtop = f'{folderMetadata["prmtopFolder"]}/complex.prmtop'

        self.initialHasRMSD = False
        self.lig_len = 0
#         self.initialRMSD = []
        
    def calculateRMSD(self, crystalComp, templateLigandGrid=None, 
                      ligandAtomConfig=None, ligandGraphConfig=None): 
        # calculate RMSD for the initial pose and then for each traj in EM, QR and MD
        try:
            t0 = time.time()
            if not self.initialHasRMSD:
                
                comp = mdtraj.load(self.initialPose, top=self.prmtop)
                self.lig_len = len(comp.top.select(f"residue {self.ligand_res}"))
                self.sys_len = comp.n_atoms
                ligand_heavy_atoms = crystalComp.top.select(f"residue {self.ligand_res} and not symbol H")
                comp.superpose(crystalComp, frame=0, atom_indices=range(0,self.sys_len - self.lig_len - 10))
                self.initialRMSD, _ = ALIGN_A_RMSD_B(comp.xyz[0]*10, crystalComp.xyz[0]*10, 
                                                     range(0, self.sys_len - self.lig_len), ligand_heavy_atoms, ligandGraphConfig)
                                                    #range(0, systemLen-lig_len), range((systemLen-lig_len), systemLen))
                #self.initialRMSD = mdtraj.rmsd(comp, crystalComp, frame=0, atom_indices=range(crystalComp.n_atoms-self.lig_len, crystalComp.n_atoms))*10
                 
                self.initialHasRMSD = True
                print(f'  Pose {self.poseName} has an initial RMSD of {self.initialRMSD}', end = ' ')
                if self.initialRMSD < 2.5:
                    print('            VVV')
                else:
                    print(' ')
            for simPrefix in self.simulationPrefixes:
                for ii in self.traj[simPrefix]:
                    ii.calculateRMSD(self.prmtop, crystalComp, self.lig_len, self.ligand_res, self.dockCenter, templateLigandGrid, ligandAtomConfig, ligandGraphConfig)
            t1 = time.time()
            displayString = f'  Done with {self.poseName:18s}: '
            for simPrefix in self.simulationPrefixes:
                displayString += f'{len(self.traj[simPrefix]):3d} {simPrefix:3s},'
            displayString += f'  {t1-t0:.3f} s'
            print(displayString)
        except:
            print(f'{self.poseName} failed!!')

    def readOutput(self):
        for simPrefix in self.simulationPrefixes:
            for ii in self.traj[simPrefix]:
                ii.readOutput()

    def readLigTraj(self):
        t0 = time.time()
        if self.lig_len == 0:
            comp = mdtraj.load(self.initialPose, top=self.prmtop)
            self.lig_len = len(comp.top.select(f"residue {self.ligand_res}"))
        for simPrefix in self.simulationPrefixes:
            for ii in self.traj[simPrefix]:
                ii.readLigTraj(self.prmtop, self.initialPose, self.ligand_res)
        t1 = time.time()
        print(f'Reading this pose took {t1-t0:.3f} s')
        
    def gatherMagnitude(self):
        # returns EM frames / EM length / QR frames / QR length / MD frames / MD length (in ns)
        frame = {}
        length = {}
        for simPrefix in self.simulationPrefixes:
            frame[simPrefix] =  len(self.traj[simPrefix]) * self.length[simPrefix] / self.writeCrdInterval[simPrefix]
            length[simPrefix] = len(self.traj[simPrefix]) * self.length[simPrefix] * self.timeStep[simPrefix] / 1000
        return frame, length
    

    def inheritPose(self, POSE):
        self.poseName = POSE.poseName
        self.ligandName = POSE.ligandName
        self.successQR = POSE.successQR
        self.ligand_res = POSE.ligand_res
        self.simulationPrefixes = POSE.simulationPrefixes
        self.dockCenter = POSE.dockCenter
        self.outSuffix = POSE.outSuffix
        self.rstSuffix = POSE.rstSuffix
        self.crdSuffix = POSE.crdSuffix
        self.nrep = {}
        for simPrefix in self.simulationPrefixes:
            self.nrep[simPrefix] = POSE.nrep[simPrefix]
            self.length[simPrefix] = POSE.length[simPrefix]
            self.writeCrdInterval[simPrefix] = POSE.writeCrdInterval[simPrefix]
            self.ioutfm[simPrefix] = POSE.ioutfm[simPrefix]
            self.timeStep[simPrefix] = POSE.timeStep[simPrefix]
            if 'EX' in simPrefix:
                self.nrep[simPrefix] = self.nrep['MD']
        self.rootFolder = POSE.rootFolder
        self.referenceFolder = POSE.referenceFolder
        self.structureFolder = POSE.structureFolder
        self.simulationFolder = POSE.simulationFolder
        self.inpcrdFolder = POSE.inpcrdFolder
        self.prmtopFolder = POSE.prmtopFolder
        self.traj = {}
        for simPrefix in self.simulationPrefixes:
            self.traj[simPrefix] = []
            if simPrefix == 'EM' or simPrefix == 'QR':
                pass
            elif not self.successQR:
                continue
            for ii in range(self.nrep[simPrefix]):
                trjFile = f'{self.simulationFolder}/{self.poseName}/{simPrefix}_R{ii+1}{self.crdSuffix}'
                outFile = f'{self.simulationFolder}/{self.poseName}/{simPrefix}_R{ii+1}{self.outSuffix}'
                rstFile = f'{self.simulationFolder}/{self.poseName}/{simPrefix}_R{ii+1}{self.rstSuffix}'
                self.traj[simPrefix].append(mdgxTrajectory(trjFile, outFile, rstFile, self.ioutfm[simPrefix]))
                self.traj[simPrefix][ii].inheritMdgxTrajectory(POSE.traj[simPrefix][ii])

        self.initialPose = POSE.initialPose
        self.prmtop = POSE.prmtop
        self.initialHasRMSD = POSE.initialHasRMSD
        if self.initialHasRMSD:
            self.lig_len = POSE.lig_len
            self.initialRMSD = POSE.initialRMSD


        

class Ligand:
    def __init__(self, ligandName, poseNames, ligand_res, settings, success, folderMetadata, simulationPrefixes, dockCenter=None, templateLigandFile=None):

        self.simulationPrefixes = simulationPrefixes
        self.dockCenter = dockCenter
        self.templateLigandFile = templateLigandFile
        self.ligandName = ligandName

        self.numPoses = len(poseNames)
        self.poseNames = poseNames
        #print(f'Here are the poseNames: {self.poseNames}')
        #print(f'{success["QR"]}')
        self.ligand_res = ligand_res # This may be '-1'. See below for error handling
        try:
            self.poseRanks = [x.split('_')[1] for x in poseNames]
        except:
            self.poseRanks = [re.sub('[a-zA-Z]', '', x) for x in poseNames]
        #print(f'Here are the poseRanks: {self.poseRanks}')
        self.settings = settings
        self.folderMetadata = folderMetadata
        self.nrep = {}
        self.qualifiedList = {}
        self.qualifiedTruth = {}
        for simPrefix in self.simulationPrefixes:
            self.nrep[simPrefix] = int(settings[simPrefix]['N-rep'])
            if 'EX' in simPrefix:
                self.nrep[simPrefix] = int(settings['MD']['N-rep'])
            self.qualifiedList[simPrefix] = [x for x in self.poseNames if x in success['QR']]
#         print(f'{len(self.qualifiedList)} of the {self.numPoses} for ligand {self.ligandName} qualified')
            self.qualifiedTruth[simPrefix] = [(x in self.qualifiedList[simPrefix]) for x in self.poseNames]
#             print(success[simPrefix])
        print(self.qualifiedList)
        print(self.qualifiedTruth)
        self.prmtop = f'{folderMetadata["prmtopFolder"]}/{self.ligandName}.prmtop'

        if not os.path.exists(self.prmtop):
            self.prmtop = f'{folderMetadata["prmtopFolder"]}/complex.prmtop'
        

        self.hasProcessedRMSD = False
        self.hasProcessedContact = False
        
            
        # Also include crystal structure
        self.crystalPose = f'{folderMetadata["referenceFolder"]}/{self.ligandName}_0A.inpcrd'
        self.initialPose = f'{folderMetadata["inpcrdFolder"]}/{self.ligandName}_0.inpcrd'

        # Error handling
        if not os.path.exists(self.crystalPose):
            self.crystalPose = f'{folderMetadata["referenceFolder"]}/complex.inpcrd'
        if not os.path.exists(self.initialPose):
            self.initialPose = glob.glob(folderMetadata["inpcrdFolder"] + '*')[0]

        initialComp = mdtraj.load(self.initialPose, top=self.prmtop)
        if self.ligand_res == '-1':
            self.ligand_res = str(initialComp.n_residues-1)

        self.ligandAtomConfig = self.getLigandAtomConfig(initialComp)

        self.ligandGraphConfig = self.getLigandGraphConfig(initialComp)

        print(f'We have determined for this ligand its residue number is {self.ligand_res}')

        print(self.templateLigandFile)

        templateLigandComp = mdtraj.load(self.templateLigandFile)

        templateLigandElements = [vdW_map[r.element.symbol] if r.element.symbol in vdW_map.keys() else 1.2 if r.element.symbol in metals else 1.8 for r in list(templateLigandComp.top.atoms)]
        
        self.templateLigandGrid = create_ref_grid(templateLigandComp.xyz[0] * 10)#, templateLigandElements)
        print(self.templateLigandGrid.shape)

        # Create poses
        
        self.Poses = {} 

        for ii in range(self.numPoses):
            try:
                self.Poses[self.poseRanks[ii]] = Pose(self.poseNames[ii], self.poseRanks[ii], self.qualifiedTruth['QR'][ii], 
                                                      self.ligand_res, settings, folderMetadata, 
                                                      simulationPrefixes=self.simulationPrefixes,
                                                      dockCenter = self.dockCenter) 
            except:
                self.Poses[self.poseRanks[ii]] = Pose(self.poseNames[ii], self.poseRanks[ii], self.qualifiedTruth['MD'][ii], 
                                                      self.ligand_res, settings, folderMetadata, 
                                                      simulationPrefixes=self.simulationPrefixes,
                                                      dockCenter = self.dockCenter)
            
        # For stats
        self.frame = {}
        self.length = {}

        
        self.LTA = LigandTrajectoryAggregate()

    def getLigandAtomConfig(self, mol):
        protein_residue_list = [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 
            'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
            'SER', 'THR', 'TYR', 'VAL', 'TRP', 
            'ACE', 'NME', 'HID', 'HIE', 'HIP', 
            'WAT', 'HOH', 'TIP3'] # For additional caps and protonation states of HIS
        pro = []
        pro_allCSP = []
        pro_polar = [] # Includes N, O, and H connected to N and O
        pro_CSP = []   # Includes C, S and P that are not connected to N, or O
        pro_CH = []    # Includes C and H where C is not connected to N or O
        lig = []
        lig_allCSP = []
        lig_polar = [] # Includes F, N, O, and H connected to N and O
        lig_CSP = []   # Includes C, S and P that are not connected to F, N, or O
        lig_CH = []
        lig_halo = []  # Includes F, Cl, and Br
        
        atom_list = {}
        bond_to_symbol = {}
        bond_to_index = {}
        
        for atom in mol.top.atoms:
            atom_list[atom.index] = atom.element.symbol
            bond_to_symbol[atom.index] = []
            bond_to_index[atom.index] = []
        
        for bond in mol.top.bonds:
            bond_to_symbol[bond.atom1.index].append(bond.atom2.element.symbol)
            bond_to_symbol[bond.atom2.index].append(bond.atom1.element.symbol)
            bond_to_index[bond.atom1.index].append(bond.atom2.index)
            bond_to_index[bond.atom2.index].append(bond.atom1.index)
        
        for ii in range(len(atom_list)):
            bond_to_symbol[ii] = np.array(bond_to_symbol[ii])
            bond_to_index[ii] = np.array(bond_to_index[ii])
            sort_order = np.argsort(bond_to_symbol[ii])
            if len(sort_order) > 1:
                bond_to_symbol[ii] = bond_to_symbol[ii][sort_order]
                bond_to_index[ii] = bond_to_index[ii][sort_order]
            bond_to_symbol[ii] = list(bond_to_symbol[ii])
            bond_to_index[ii] = list(bond_to_index[ii])
        
            
        for atom in mol.top.atoms:
            if atom.residue.name in protein_residue_list: # part of a protein
                pro.append(atom.index)
                if atom.element.symbol in ['N','O']:
                    pro_polar.append(atom.index)
                elif atom.element.symbol == 'H':
                    if bond_to_symbol[atom.index][0] in ['N', 'O']:
                        pro_polar.append(atom.index)
                    elif bond_to_symbol[atom.index][0] == 'C':
                        if not np.any([x in ['N','O'] for x in bond_to_symbol[bond_to_index[atom.index][0]]]):
                            pro_CH.append(atom.index)
                elif atom.element.symbol in ['C','S','P']:
                    pro_allCSP.append(atom.index)
                    if not np.any([x in ['N','O'] for x in bond_to_symbol[atom.index]]):
                        pro_CH.append(atom.index)
                        pro_CSP.append(atom.index)
            else: # part of a ligand
                lig.append(atom.index)
                if atom.element.symbol in ['N','O','F']:
                    lig_polar.append(atom.index)
                elif atom.element.symbol == 'H':
                    if bond_to_symbol[atom.index][0] in ['N','O','F']:
                        lig_polar.append(atom.index)
                    elif bond_to_symbol[atom.index][0] == 'C':
                        if not np.any([x in ['N','O','F'] for x in bond_to_symbol[bond_to_index[atom.index][0]]]):
                            lig_CH.append(atom.index)
                elif atom.element.symbol in ['C','S','P']:
                    lig_allCSP.append(atom.index)
                    if not np.any([x in ['N','O','F'] for x in bond_to_symbol[atom.index]]):
                        lig_CH.append(atom.index)
                        lig_CSP.append(atom.index)
                elif atom.element.symbol in ['Cl','Br','I']:
                    lig_halo.append(atom.index)
                
        molConfig = {'pro': np.array(pro), 'pro_polar': np.array(pro_polar), 
                     'pro_allCSP': np.array(pro_allCSP), 'pro_CSP': np.array(pro_CSP), 'pro_CH': np.array(pro_CH),
                     'lig': np.array(lig), 'lig_polar': np.array(lig_polar), 
                     'lig_allCSP': np.array(lig_allCSP), 'lig_CSP': np.array(lig_CSP), 
                     'lig_CH': np.array(lig_CH), 'lig_halo': np.array(lig_halo)}

        #print('molConfig dump:')
        #print(molConfig)

        return molConfig
 
    def getLigandGraphConfig(self, mol):
        if self.ligand_res == '-1':
            self.ligand_res = str(crystalComp.n_residues-1)
        lig_index = mol.top.select(f'resid {self.ligand_res} and not element H')
        lig_index_h = mol.top.select(f'resid {self.ligand_res}')
        mol_lig = mol.atom_slice(lig_index) 
        mol_lig_h = mol.atom_slice(lig_index_h)
        an = np.array([atom.element.atomic_number for atom in mol_lig.top.atoms])
        an_h = np.array([atom.element.atomic_number for atom in mol_lig_h.top.atoms])
        G = nx.Graph()
        G.add_nodes_from(range(len(mol_lig.xyz[0])))
        bond_list = []
        for bond in mol_lig.top.bonds:
            bond_list.append([bond.atom1.index, bond.atom2.index])
        G.add_edges_from(bond_list)
        G_h = nx.Graph()
        G_h.add_nodes_from(range(len(mol_lig_h.xyz[0])))
        bond_list = []
        for bond in mol_lig_h.top.bonds:
            bond_list.append([bond.atom1.index, bond.atom2.index])
        G_h.add_edges_from(bond_list)

        graphConfig = {'lig_index': lig_index, 'lig_index_h': lig_index_h,
                       'an': an, 'an_h': an_h, 'G': G, 'G_h': G_h} 
        return graphConfig
        
        # lig_index is indices of ligand atom without H, lig_index_h is with
        # an is "atomic number", and G is the "graph"
        # These are for the spyrmsd calculations



    def calculateRMSD(self): # Invoke calculateRMSD in each pose
        t0 = time.time()
        print(f'Calculating RMSD for ligand {self.ligandName} ...')
#// mdtraj
        crystalComp = mdtraj.load(self.crystalPose, top=self.prmtop)
        if self.ligand_res == '-1':
            self.ligand_res = str(crystalComp.n_residues-1)
#// MDAnalysis
#         crystalComp = mda.Universe(self.prmtop, self.crystalPose)
        for ii in range(self.numPoses):
            self.Poses[self.poseRanks[ii]].calculateRMSD(crystalComp)
        t1 = time.time()
        return (t1-t0)
    
    
    def calculateRMSD_parallel(self):
        t0 = time.time()
        if self.hasProcessedRMSD:
            print(f'{self.ligandName} has processed RMSD, skip.')
            t1 = time.time()

        else:

            print(f'Calculating RMSD for ligand {self.ligandName} ...')
            print(f'self.crystalPose is at {self.crystalPose}')
            try:
                crystalComp = mdtraj.load(self.crystalPose, top=self.prmtop)
                print(f'Successfully loaded crystal pose at {self.crystalPose}')
            except:
                crystalComp = mdtraj.load(self.initialPose, top=self.prmtop)
            if self.ligand_res == '-1':
                self.ligand_res = str(crystalComp.n_residues-1)

            pool = mp.Pool()
            lr = pool.map(Distribute_Lig_pose_RMSD, ((self.Poses[self.poseRanks[ii]], crystalComp, self.templateLigandGrid, self.ligandAtomConfig, self.ligandGraphConfig) for ii in range(self.numPoses)))
            pool.close()
            pool.join()
            for ii in range(self.numPoses):
                self.Poses[self.poseRanks[ii]].inheritPose(lr[ii])
            self.hasProcessedRMSD = True
            t1 = time.time()
        return (t1-t0)


    def calculateRMSD_parallel2(self, pool):
        t0 = time.time()
        if self.hasProcessedRMSD:
            print(f'{self.ligandName} has processed RMSD, skip.')
            t1 = time.time()

        else:

            print(f'Calculating RMSD for ligand {self.ligandName} ...')
    #// mdtraj
            crystalComp = mdtraj.load(self.crystalPose, top=self.prmtop)
    #// MDAnalysis
    #         crystalComp = mda.Universe(self.prmtop, self.crystalPose)
            lr = pool.map(Distribute_Lig_pose_RMSD, ((self.Poses[self.poseRanks[ii]], crystalComp) for ii in range(self.numPoses)))
    #         print(lr[0].trajMD[0].RMSD)
            for ii in range(self.numPoses):
                self.Poses[self.poseRanks[ii]].inheritPose(lr[ii])
            self.hasProcessedRMSD = True
            t1 = time.time()
        return (t1-t0)
    
    def readOutput(self):
        for ii in range(self.numPoses):
            self.Poses[self.poseRanks[ii]].readOutput()

    def readLigTraj(self):
        for ii in range(self.numPoses):
            self.Poses[self.poseRanks[ii]].readLigTraj()
        
    def gatherMagnitude(self):
        print(f'In ligand {self.ligandName}, we have {self.numPoses} poses')
        self.frame = {}
        self.length = {}
        for ii in range(self.numPoses):
            print(f'In pose {ii}')
            f, l = self.Poses[self.poseRanks[ii]].gatherMagnitude()
            if ii == 0:
                for simPrefix in self.simulationPrefixes:
                    self.frame[simPrefix] = f[simPrefix]
                    self.length[simPrefix] = l[simPrefix]
            else:
                for simPrefix in self.simulationPrefixes:
                    self.frame[simPrefix] += f[simPrefix]
                    self.length[simPrefix] += l[simPrefix]
                
        self.numberOfSimulations = {}
        for simPrefix in self.simulationPrefixes:
            if simPrefix == 'EM' or simPrefix == 'QR':
                self.numberOfSimulations[simPrefix] = self.nrep[simPrefix] * self.numPoses
            else:
                self.numberOfSimulations[simPrefix] = self.nrep[simPrefix] * len(self.qualifiedList[simPrefix])
#         print(self.numEM, self.frameEM, self.lengthEM, self.numQR, self.frameQR, self.lengthQR, self.numMD, self.frameMD, self.lengthMD)
        return self.frame, self.length, self.numberOfSimulations
#         return np.array([self.numEM, self.frameEM, self.lengthEM, self.numQR, self.frameQR, self.lengthQR, self.numMD, self.frameMD, self.lengthMD],dtype=float)
    
    def gatherRMSD(self, mode='MD', trajStartCut=5, trajEndCut=0): # Discard trajStartCut frames from the beginning and trajEndCut frames from the end for every trajectory
        RMSD = []
        output = []
        reference = []
        for ii in range(self.numPoses):
            for traj in self.Poses[self.poseRanks[ii]].traj[mode]:
                if traj.success:
                    frames = range(trajStartCut, len(traj.RMSD)-trajEndCut)
                    RMSD.append(traj.RMSD[trajStartCut:(len(traj.RMSD)-trajEndCut)])
                    output.append(traj.output[trajStartCut:(len(traj.RMSD)-trajEndCut)])
                    reference.append([(traj.trjFile, idx) for idx in frames])
        try:
            RMSD = np.array(RMSD).flatten()
        except:
            pass
        if RMSD.dtype != float:
            print('We have ragged RMSD arrays, use lists to fix this')
            RMSD = np.array(sum([list(x) for x in RMSD], []))
        try:
            output = pd.concat(output)
        except:
            pass
        try:
            reference = [y for x in reference for y in x]
        except:
            pass

                


        return RMSD, output, reference
    
    def prepareLTA(self):
        # First prepare some invariant things
        
        crystalComp = mdtraj.load(self.crystalPose, top=self.prmtop)
        
        # Protein selection AFTER removing H
        self.proteinAtomSelection = np.array([ 45, 106, 107, 167, 168, 170, 175, 176, 177, 178, 179, 180, 181, 182, 232, 259, 381, 386, 387, 388])
        pro = crystalComp.top.select(f"not residue {self.ligand_res} and not symbol H")
        self.LTA.proteinCoordinate = crystalComp.xyz[0][pro][self.proteinAtomSelection]*10
#         self.LTA.proteinType = np.array([protein_map[x.element.symbol.upper()] for x in np.array(list(crystalComp.top.atoms))[pro][self.proteinAtomSelection]])
        self.LTA.proteinType = np.arange(len(self.proteinAtomSelection))
        
        lig = crystalComp.top.select(f"residue {self.ligand_res} and not symbol H")
        self.LTA.baseLigandType = np.array([ligand_map[x.element.symbol.upper()] for x in np.array(list(crystalComp.top.atoms))[lig]])
            
        for ii in self.Poses:
            if self.Poses[ii].successQR:
                for jj in self.Poses[ii].trajEM:
                    if jj.hasLigandTrajectory:
                        self.LTA.ligandTrajectory.append(jj.ligandTrajectory)
                        self.LTA.trajectoryAddress.append(jj)
                        self.LTA.trajectoryLength.append(len(jj.ligandTrajectory))
                        self.LTA.Nsim += 1
                for jj in self.Poses[ii].trajQR:
                    if jj.hasLigandTrajectory:
                        self.LTA.ligandTrajectory.append(jj.ligandTrajectory)
                        self.LTA.trajectoryAddress.append(jj)
                        self.LTA.trajectoryLength.append(len(jj.ligandTrajectory))
                        self.LTA.Nsim += 1
                for jj in self.Poses[ii].trajMD:
                    if jj.hasLigandTrajectory:
                        self.LTA.ligandTrajectory.append(jj.ligandTrajectory)
                        self.LTA.trajectoryAddress.append(jj)
                        self.LTA.trajectoryLength.append(len(jj.ligandTrajectory))
                        self.LTA.Nsim += 1

        self.LTA.trajectoryLength = np.array(self.LTA.trajectoryLength)
    
    def distributeLTA(self):
        if self.LTA.hasFeatures:
            cumLength = 0
            for traj, length in zip(self.LTA.trajectoryAddress, self.LTA.trajectoryLength):
                traj.contactScore = (self.LTA.features[cumLength:cumLength+length] > 0).sum(1)
                traj.features = self.LTA.features[cumLength:cumLength+length]
                traj.hasContactScore = True
                cumLength += length
            self.hasProcessedContact = True
        else:
            raise('First calculate the contact then distribute')
        self.LTA.trajectoryAddress = []
        self.LTA.trajectoryLength = []
        self.LTA.ligandTrajectory = []
        self.LTA.features = []
    
    def inheritLigand(self, LIG):
        try:
            self.simulationPrefixes = LIG.simulationPrefixes
        except:
            self.simulationPrefixes = ['EM','QR','MD']
        self.dockCenter = LIG.dockCenter
        self.templateLigandFile = LIG.templateLigandFile
        self.templateLigandGrid = LIG.templateLigandGrid
        self.ligandName = LIG.ligandName
        self.numPoses = LIG.numPoses
        self.poseNames = LIG.poseNames
        self.poseRanks = LIG.poseRanks
        self.ligandAtomConfig = LIG.ligandAtomConfig
        self.ligandGraphConfig = LIG.ligandGraphConfig
        self.ligand_res = LIG.ligand_res
        self.nrep = {}
        for simPrefix in simulationPrefixes:
            self.nrep[simPrefix] = LIG.nrep[simPrefix]
            if 'EX' in simPrefix:
                self.nrep[simPrefix] = LIG.nrep['MD']
        self.qualifiedList = LIG.qualifiedList
        self.qualifiedTruth = LIG.qualifiedTruth

        self.prmtop = LIG.prmtop
        self.crystalPose = LIG.crystalPose
        self.initialPose = LIG.initialPose
        self.hasProcessedRMSD = LIG.hasProcessedRMSD
        self.hasProcessedContact = LIG.hasProcessedContact
        
        self.Poses = {}
        for ii in range(self.numPoses):
            self.Poses[self.poseRanks[ii]] = Pose(self.poseNames[ii], self.poseRanks[ii], self.qualifiedTruth[ii], self.settings, self.folderMetadata)
            self.Poses[self.poseRanks[ii]].inheritPose(LIG.Poses[self.poseRanks[ii]])

        try:
            self.LTA = LIG.LTA
        except:
            self.LTA = LigandTrajectoryAggregate()

            
            
class MolecularDynamicsRefinement:
    def __init__(self, root='./',simulationPrefixes=['EM','QR','MD'], ligand_res=-1, **kwargs):

        # Ligand residue usually is the last residue of the system
        
        self.simulationPrefixes = simulationPrefixes
        
        self.rootFolder = root
        self.referenceFolder = self.rootFolder + '/reference_structure/'
        self.structureFolder = self.rootFolder + '/Structure/'
        self.simulationFolder = self.rootFolder + '/Simulation/'
        self.planningFolder = self.structureFolder + '/planning/'
        self.configFolder = self.structureFolder + '/planning/'
        self.inpcrdFolder = self.structureFolder + '/inpcrd/'
        self.prmtopFolder = self.structureFolder + '/prmtop/'
        self.saveFolder = self.rootFolder + '/pickle/'
        self.overrideSuccess = False # This is for newly devleoped method where minimization is done with openmm
        self.inputFolder = None
 
        for k, v in kwargs.items(): # Override these folder path assignments
            if k == 'rootFolder':
                self.rootFolder = v
            elif k == 'referenceFolder':
                self.referenceFolder = v
            elif k == 'structureFolder':
                self.structureFolder = v
            elif k == 'simulationFolder':
                self.simulationFolder = v
            elif k == 'planningFolder':
                self.planningFolder = v
            elif k == 'configFolder':
                self.configFolder = v
            elif k == 'inpcrdFolder':
                self.inpcrdFolder = v
            elif k == 'prmtopFolder':
                self.prmtopFolder = v
            elif k == 'saveFolder':
                self.saveFolder = v
            elif k == 'overrideSuccess':
                self.overrideSuccess = v
            elif k == 'inputFolder':
                self.inputFolder = v


        
        self.folderMetadata = {
            'rootFolder': self.rootFolder,
            'referenceFolder': self.referenceFolder,
            'structureFolder': self.structureFolder,
            'simulationFolder': self.simulationFolder,
            'planningFolder': self.planningFolder,
            'configFolder': self.configFolder,
            'inpcrdFolder': self.inpcrdFolder,
            'prmtopFolder': self.prmtopFolder,
            'saveFolder': self.saveFolder
        }        
       
        #print(f'folderMetadata is {self.folderMetadata}')


        temp = os.listdir(self.inpcrdFolder)
        temp.sort()
        self.inpcrdFiles = [x for x in temp if '.inpcrd' in x]
        self.numPoses = len(self.inpcrdFiles)
        self.poseNames = [x.split('.')[0] for x in self.inpcrdFiles]
        
        temp = os.listdir(self.prmtopFolder)
        temp.sort()
        self.prmtopFiles = [x for x in temp if ('.prmtop' in x) and ('openmm' not in x)]
        self.numSystems = len(self.prmtopFiles)
        self.sysNames = [x.split('.')[0] for x in self.prmtopFiles]
        self.customSysNames = False
        for k, v in kwargs.items(): # Override these folder path assignments
            if k == 'sysNames':
                self.sysNames = [v]
                self.customSysNames = True


        self.ligand_res = str(ligand_res)

        with open(self.configFolder + 'simulation_config.json', 'r') as f:
            self.settings = json.load(f)
            
        # Load all successes
        if not self.overrideSuccess:
            self.success = {}
            for simPrefix in self.simulationPrefixes:
                try:
                    with open(self.planningFolder + simPrefix + '_success.txt','r') as f:
                        if self.customSysNames:
                            self.success[simPrefix] = [x.split('/')[-1] for x in f.readlines()[0].strip('\n').split(',') if self.sysNames[0] in x]
                        else:
                            self.success[simPrefix] = f.readlines()[0].strip('\n').split(',')
                        print(f'Success runs: {self.success[simPrefix]}')
                except:
                    self.success[simPrefix] = []
        else:
            print('Overriding success ... for openmm')
            self.success = {}
            for simPrefix in self.simulationPrefixes:
                self.success['QR'] = self.poseNames 
        
        self.Ligands = {}
        

        
        
    def createLigands(self):
        
        # Read the job description file
        if self.inputFolder is not None:
            try:
                Jobs = pd.read_csv(self.inputFolder + 'job_description.csv') 
            except:
                Jobs = None

        for sysName in self.sysNames:
            try:
                dockCenter = [Jobs.loc[Jobs['Job_name'] == sysName]['dockX'].values[0],
                              Jobs.loc[Jobs['Job_name'] == sysName]['dockY'].values[0],
                              Jobs.loc[Jobs['Job_name'] == sysName]['dockZ'].values[0]]
                print(f'For {sysName} the dock center is {dockCenter}')
            except:
                dockCenter = None
            try:
                templateLigandFile = self.inputFolder + Jobs.loc[Jobs['Job_name'] == sysName]['template_file'].values[0]
                print(f'Template ligand file is {templateLigandFile}')
            except:
                templateLigandFile = None 
            if len(self.sysNames) == 1: # Only one system, all poseNames belong to the sysNames
                poseName = [x for x in self.poseNames]
            else:
                poseName = [x for x in self.poseNames if sysName in x]
            self.Ligands[sysName] = Ligand(sysName, poseName, self.ligand_res, self.settings, self.success, self.folderMetadata,
                                           simulationPrefixes=self.simulationPrefixes, dockCenter=dockCenter, 
                                           templateLigandFile=templateLigandFile)
            
    def calculateRMSD(self): # invoke calculateRMSD in each ligand
        completedRMSD = 0
        for sysName in self.sysNames:
            print(f'MDR: Ligand {sysName}')
            t = self.Ligands[sysName].calculateRMSD()
#             completedRMSD += self.Ligands[sysName].gatherMagnitude()
            completedRMSD += 1
            if self.sumMD > 0:
                print(f'Done with {sysName} in {t:.3f} s ({(completedRMSD / self.numSystems * 100):.2f} %) | Estimated remaining: {(self.numSystems - completedRMSD) * t:.0f} s')

    def calculateRMSD_parallel(self): # invoke calculateRMSD_parallel in each ligand
        # First check if there are ligand pickles with processed RMSD - saves time

        try:
            self.loadLigands()
        except:
            pass
        completedRMSD = 0
        self.save()
        for sysName in self.sysNames:
            print(f'MDR: Ligand {sysName}')
            t = self.Ligands[sysName].calculateRMSD_parallel()
            self.saveOneLigand(sysName)
#             completedRMSD += self.Ligands[sysName].gatherMagnitude()
            completedRMSD += 1
            try:
                print(f'Done with {sysName} in {t:.3f} s ({(completedRMSD / self.numSystems * 100):.2f} %) | Estimated remaining: {(self.numSystems - completedRMSD) * t:.0f} s')
            except:
                pass
            
    def calculateRMSD_parallel2(self): # invoke calculateRMSD_parallel2 in each ligand
        self.save()
        try:
            self.loadLigands()
        except:
            pass
        completedRMSD = 0
        pool = mp.Pool(processes=16)
        for sysName in self.sysNames:
            print(f'MDR: Ligand {sysName}')
            t = self.Ligands[sysName].calculateRMSD_parallel2(pool)
            self.saveOneLigand(sysName)
#             completedRMSD += self.Ligands[sysName].gatherMagnitude()
            completedRMSD += 1
            if self.sumMD > 0:
                print(f'Done with {sysName} in {t:.3f} s ({(completedRMSD / self.numSystems * 100):.2f} %) | Estimated remaining: {(self.numSystems - completedRMSD) * t:.0f} s')
        pool.close()
        pool.join()
                
                
    def readOutput(self):
        for sysName in self.sysNames:
            self.Ligands[sysName].readOutput()
            print(f'Done with force re-reading output for ligand {sysName}')

    def readLigTraj(self):
        for sysName in self.sysNames:
            self.Ligands[sysName].readLigTraj()
            print(f'Done with force re-reading ligand trajectory for ligand {sysName}')
    
    def gatherMagnitude(self, mute=False):
        if not mute:
            print(f'In this dataset, there are:')
            print(f'{self.numSystems:5d} different molecules')
            print(f'{self.numPoses:5d} different poses')
        
#         accList = ([0, 0, 0, 0, 0, 0, 0, 0, 0])
#         for sysName in self.sysNames:
#             PD = self.Ligands[sysName].gatherMagnitude()
#             accList += PD
        self.numberOfSimulations = {}
        self.length = {}
        self.frame = {}
        for idx, sysName in enumerate(self.sysNames):
            print(f'In number {idx} and sysName {sysName}')
            f, l, n = self.Ligands[sysName].gatherMagnitude()
            if idx == 0:
                for simPrefix in self.simulationPrefixes:
                    self.numberOfSimulations[simPrefix] = n[simPrefix]
                    self.length[simPrefix] = l[simPrefix]
                    self.frame[simPrefix] = f[simPrefix]
            else:
                for simPrefix in self.simulationPrefixes:
                    self.numberOfSimulations[simPrefix] += n[simPrefix]
                    self.length[simPrefix] += l[simPrefix]
                    self.frame[simPrefix] += f[simPrefix]
                    
#         self.sumEM = accList[0]
#         self.frameEM = accList[1]
#         self.lengthEM = accList[2]
#         self.sumQR = accList[3]
#         self.frameQR = accList[4]
#         self.lengthQR = accList[5]
#         self.sumMD = accList[6]
#         self.frameMD = accList[7]
#         self.lengthMD = accList[8]
        if not mute:
            for simPrefix in self.simulationPrefixes:
                if simPrefix == 'QR':
                    print(f'{self.numberOfSimulations[simPrefix]:6.0f} QR   | {self.frame[simPrefix]:8.0f} frames | {self.length[simPrefix]:10.2f} nanoseconds | {len(self.success[simPrefix]):5d} qualified poses')
                else:
                    print(f'{self.numberOfSimulations[simPrefix]:6.0f} {simPrefix:3s}  | {self.frame[simPrefix]:8.0f} frames | {self.length[simPrefix]:10.2f} nanoseconds')
#             print(f'{self.sumMD:6.0f} molecular dynamics   | {self.frameMD:8.0f} frames | {self.lengthMD:10.2f} nanoseconds')

    def processLTA(self):
        for sysName in self.sysNames:
            if self.Ligands[sysName].hasProcessedContact:
                print(f'{sysName} has processed contacts, skip.')
            else:
                self.Ligands[sysName].prepareLTA()
                self.Ligands[sysName].LTA.calculateContact()
                self.Ligands[sysName].distributeLTA()
    
    def save(self, fname='MDR_dump.pkl'):
        try:
            os.mkdir(self.saveFolder)
        except:
            pass # Folder exists
        # Only saves MDR-level variables
        varlist = [self.rootFolder, self.referenceFolder, self.structureFolder, self.simulationFolder,
                   self.planningFolder, self.inpcrdFolder, self.prmtopFolder, self.saveFolder, self.folderMetadata,
                   self.inpcrdFiles, self.numPoses, self.poseNames, self.prmtopFiles, self.numSystems,
                   self.sysNames, self.ligand_res, self.settings, self.success, self.overrideSuccess, self.simulationPrefixes,
                   self.inputFolder]
        with open(f'{self.saveFolder}/{fname}', 'wb') as f:
            pickle.dump(varlist, f)
            
    def load(self, fname):
        with open(fname, 'rb') as f:
            self.rootFolder, self.referenceFolder, self.structureFolder, self.simulationFolder, \
               self.planningFolder, self.inpcrdFolder, self.prmtopFolder, self.saveFolder, self.folderMetadata, \
               self.inpcrdFiles, self.numPoses, self.poseNames, self.prmtopFiles, self.numSystems, \
               self.sysNames, self.ligand_res, self.settings, self.success, self.overrideSuccess, self.simulationPrefixes, \
               self.inputFolder = \
                    pickle.load(f)

    def saveOneLigand(self, sysName):
        with open(f'{self.saveFolder}/{sysName}.pkl', 'wb') as f:
            pickle.dump(self.Ligands[sysName], f)
            
    def saveLigands(self):
        for sysName in self.sysNames:
            with open(f'{self.saveFolder}/{sysName}.pkl', 'wb') as f:
                pickle.dump(self.Ligands[sysName], f)
                
    def loadLigands(self):
        for sysName in self.sysNames:
            with open(f'{self.saveFolder}/{sysName}.pkl', 'rb') as f:
                self.Ligands[sysName] = pickle.load(f)
                print(f'Loading {sysName} successful')
    
    def inheritMDR(self, MDR): # This clears current MDR!
        
        try:
            self.simulationPrefixes = MDR.simulationPrefixes
        except:
            self.simulationPrefixes = ['EM','QR','MD']
        
        self.rootFolder = MDR.rootFolder
        self.referenceFolder = MDR.referenceFolder
        self.structureFolder = MDR.structureFolder
        self.simulationFolder = MDR.simulationFolder
        self.planningFolder = MDR.planningFolder
        self.inpcrdFolder = MDR.inpcrdFolder
        self.prmtopFolder = MDR.prmtopFolder
        self.saveFolder = MDR.saveFolder
        self.overrideSuccess = MDR.overrideSuccess 
        self.folderMetadata = MDR.folderMetadata
        
        self.inpcrdFiles = MDR.inpcrdFiles
        self.numPoses = MDR.numPoses
        self.poseNames = MDR.poseNames
        
        self.prmtopFiles = MDR.prmtopFiles
        self.numSystems = MDR.numSystems
        self.sysNames = MDR.sysNames

        self.settings = MDR.settings
            
        self.QR_success = MDR.QR_success
            
        self.Ligands = {}
        
        for sysName in self.sysNames:
            poseName = [x for x in self.poseNames if sysName in x]
            self.Ligands[sysName] = Ligand(sysName, poseName, self.settings, self.QR_success, self.folderMetadata, simulationPrefixes=self.simulationPrefixes)
            self.Ligands[sysName].inheritLigand(MDR.Ligands[sysName])
       


 
            
class LigandTrajectoryAggregate: #
    def __init__(self):
        
        # These are to be filled at the Ligand level in Ligand.prepareLTA()
        
        # Nsim includes EM, QR and MD
        self.Nsim = 0
        
        # aggregated ligand trajectory, [(Nsim * Nframe) * N_lig_atoms * 3]
        # Stored as a list of numpy arrays, converted to a huge cp array at runtime
        self.ligandTrajectory = []
        
        # Ligand types, [(Nsim * Nframe) * N_lig_atoms]
        self.baseLigandType = []
#         self.ligandType = [] # This will simply be a repeat tile of baseLigandType
        
        # This will be something like np.linspace(0, len(self.ligandType), self.Nsim+1)
#         self.ligandSeparation = []
        
        # Just N_prot_atoms * 3 or * 1 arrays
        self.proteinCoordinate = None
        self.proteinType = None
        
        # References to the mdgxTrajectory objects (to put things back into place)
        self.trajectoryAddress = []
        
        # Used in distribution back to mdgxTrajectory
        self.trajectoryLength = []
        
        # Return from kernel
        self.features = None
        self.hasFeatures = None
        
    def calculateContact(self):
        t0 = time.time()
        LTraj = np.concatenate(self.ligandTrajectory).reshape((-1, 3))
        LType = np.tile(self.baseLigandType, np.sum(self.trajectoryLength))
        LSep = np.linspace(0, len(LType), np.sum(self.trajectoryLength)+1, dtype=int)
        
        t1 = time.time()        
        x_ligand = cp.array(LTraj[:,0],dtype=cp.float32)
        y_ligand = cp.array(LTraj[:,1],dtype=cp.float32)
        z_ligand = cp.array(LTraj[:,2],dtype=cp.float32)
        x_protein = cp.array(self.proteinCoordinate[:,0])
        y_protein = cp.array(self.proteinCoordinate[:,1])
        z_protein = cp.array(self.proteinCoordinate[:,2])
        t_protein = cp.array(self.proteinType, dtype=cp.int32)
        types_ligand = cp.array(LType, dtype=cp.int32)
        offset = cp.array(LSep, dtype=cp.int32)
        nbins = 1
        binsize = 3.4

        t2 = time.time()
        feat = kernel.compute(x_ligand, y_ligand, z_ligand, types_ligand, offset,
                   x_protein, y_protein, z_protein, t_protein,
                   cutoff=nbins*binsize,binsize=binsize,nbins=nbins,
                   n_receptor_types=len(t_protein),max_ligand_atoms=(len(self.baseLigandType)+31)//32*32)
        
#         print(self.feat.shape)
#         print(type(self.feat))
        t3 = time.time()
#         self.features = feat
        self.features = cp.asnumpy(feat)
        self.hasFeatures = True
        t4 = time.time()
        print(f'Timing: Load {(t1-t0)*1000:.2f} ms, set {(t2-t1)*1000:.2f} ms, calc {(t3-t2)*1000:.2f} ms, convert {(t4-t3)*1000:.2f} ms')

        
        
