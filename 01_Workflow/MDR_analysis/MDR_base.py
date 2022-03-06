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

# import MDAnalysis as mda
# import MDAnalysis.analysis.rms


def ALIGN_A_RMSD_B(P, Q, A, B):
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
    diff = P[B] - Q[B]
    N = len(P[B])
    return np.sqrt((diff * diff).sum() / N), P + QU.mean(axis=0) 

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
    obj, crystalComp = arg
    obj.calculateRMSD(crystalComp)
    return obj

def readMDCRD(fname, natom):
    with open(fname, 'r') as f:
        f.readline()
        cont = f.read()
    xyz = list(map(float,cont.split()))
    return np.array(xyz).reshape(-1, natom, 3)


class mdgxTrajectory:
    def __init__(self, trjFile, outFile, rstFile, ioutfm):
        self.trjFile = trjFile
        self.outFile = outFile
        self.rstFile = rstFile
        self.ioutfm = ioutfm
        self.success = False
        self.hasRMSD = False
        self.hasTrjFile = False
        self.hasOutFile = False
        self.output = None
        self.RMSD = []
        self.contactScore = []
        self.hasContactScore = False
        self.ligandTrajectory = None # We store this information for contact and clustering calculation
        self.ligandTrajectoryH = None # We store this information for contact and clustering calculation
        self.hasLigandTrajectory = False
#         self.trjLength = 0
#         self.rmsd = []
    
    def calculateRMSD(self, prmtop, crystalComp, lig_len, ligand_res): 
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
                    comp = readMDCRD(self.trjFile, systemLen)[:-1] 
                if Profile:
                    t1 = time.time()
                self.trjLength = len(comp)
    
                if Profile:
                    t2 = time.time()


                R = []
                LT = []
                systemLen = crystalComp.n_atoms
                referenceXYZ = crystalComp.xyz[0]
       
                ligand_heavy_atoms = crystalComp.top.select(f"residue {ligand_res} and not symbol H")
                #print(f'ligand heavy atoms index: {ligand_heavy_atoms}')
                if self.ioutfm == 1:
                    for xyz in comp.xyz:
                        R_this, LT_this = ALIGN_A_RMSD_B(xyz*10, referenceXYZ*10,
                                                range(0, systemLen-lig_len), ligand_heavy_atoms)
                                                #range(0, systemLen-lig_len), range((systemLen-lig_len), systemLen))
                        R.append(R_this)
                        LT.append(LT_this)
                elif self.ioutfm == 0:
                    for xyz in comp:
                        R_this, LT_this = ALIGN_A_RMSD_B(xyz, referenceXYZ*10,
                                                range(0, systemLen-lig_len), ligand_heavy_atoms)
                                                #range(0, systemLen-lig_len), range((systemLen-lig_len), systemLen))
                        R.append(R_this)
                        LT.append(LT_this)
                self.RMSD = np.array(R)
                LT = np.array(LT)
                
                

                if Profile:
                    t3 = time.time()
                self.output = ReadLog(self.outFile)
                if Profile:
                    t4 = time.time()
                if Profile:
                    print(f'    Profiling: Loading comp {(t1-t0)*1000:.3f} ms | Superpose {(t2-t1)*1000:.3f} ms | RMSD {(t3-t2)*1000:.3f} ms | Read output {(t4-t3)*1000:.3f} ms ')
                self.hasRMSD = True
                self.hasTrjFile = True
                self.hasOutFile = True
                self.success = True
                 
                if not self.hasLigandTrajectory: # Also store ligand trajectories for contact analysis

                    lig = crystalComp.top.select(f"residue {ligand_res} and not symbol H")
                    ligH = crystalComp.top.select(f"residue {ligand_res}")
                    pro = crystalComp.top.select(f"not residue {ligand_res} and not symbol H")
                    self.ligandTrajectory = LT[:,lig]
                    self.ligandTrajectoryH = LT[:, ligH]
                    self.hasLigandTrajectory = True

        except:
            pass
    
    def readOutput(self):
        try:
            self.output = ReadLog(self.outFile)
            self.hasOutFile = True
        except:
            pass
    
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
    
    
    def getContactScore(self, prmtop, ligand_res):
        t0 = time.time()
        if self.ioutfm == 1:
            comp = md.load(self.trjFile, top=prmtop)
        elif self.ioutfm == 0:
            comp = md.load_mdcrd(self.trjFile, top=prmtop)
#         print(comp)
        t1 = time.time()
        pro = comp.top.select(f"not residue {ligand_res} and not symbol H")
        lig = comp.top.select(f"residue {ligand_res} and not symbol H")
        close_atoms = np.array([ 45, 106, 107, 167, 168, 170, 175, 176, 177, 178, 179, 180, 181, 182, 232, 259, 381, 386, 387, 388])
        
#         self.proteinCoordinates = comp.xyz[0][pro][close_atoms]*10
#         self.ligandTrajectory2 = comp.xyz[:,lig]*10
        
        x_ligand = cp.array(comp.xyz[:,lig,0].flatten()*10)
        y_ligand = cp.array(comp.xyz[:,lig,1].flatten()*10)
        z_ligand = cp.array(comp.xyz[:,lig,2].flatten()*10)
        x_protein = cp.array(comp.xyz[0][pro][close_atoms][:,0]*10)
        y_protein = cp.array(comp.xyz[0][pro][close_atoms][:,1]*10)
        z_protein = cp.array(comp.xyz[0][pro][close_atoms][:,2]*10)

        t_protein = cp.array(cp.arange(len(close_atoms)), dtype=cp.int32)
#         types_protein = cp.array([protein_map[x.element.symbol.upper()] for x in np.array(list(comp.top.atoms))[pro]], dtype=cp.int32)
        types_ligand = cp.tile(cp.array([ligand_map[x.element.symbol.upper()] for x in np.array(list(comp.top.atoms))[lig]], dtype=cp.int32), comp.n_frames)
        offset = cp.linspace(0,len(types_ligand),comp.n_frames+1,dtype=cp.int32)
        nbins = 1
        binsize = 3.4

        t2 = time.time()
        feat = kernel.compute(x_ligand, y_ligand, z_ligand, types_ligand, offset,
                   x_protein, y_protein, z_protein, t_protein,
                   cutoff=nbins*binsize,binsize=binsize,nbins=nbins,
                   n_receptor_types=len(t_protein), max_ligand_atoms=(len(lig)+31)//32*32)
        print(feat.shape)
        t3 = time.time()
        self.features = cp.asnumpy(feat)
#         print(feat2.shape)
        t35 = time.time()
#         self.contactScore = feat2.reshape(len(comp), -1, len(pro))
#         self.contactScore = contactScore.get()
#         self.hasContactScore = True
        t4 = time.time()
        print(f'Timing: Load {(t1-t0)*1000:.2f} ms, set {(t2-t1)*1000:.2f} ms, calc {(t3-t2)*1000:.2f} ms, convert {(t4-t3)*1000:.2f} ms')
        print(f'Timing: feat get {(t35-t3)*1000:.2f} ms, reshape {(t4-t35)*1000:.2f} ms')
        
    
    def inheritMdgxTrajectory(self, TRJ):
        self.trjFile = TRJ.trjFile
        self.outFile = TRJ.outFile
        self.rstFile = TRJ.rstFile
        self.success = TRJ.success
        self.hasRMSD = TRJ.hasRMSD
#         if self.hasRMSD:
        self.RMSD = TRJ.RMSD

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
    def __init__(self, name, rank, successQR=False, ligand_res=None, settings=None, folderMetadata=None, simulationPrefixes=['EM','QR','MD']):
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
        
    def calculateRMSD(self, crystalComp): # calculate RMSD for the initial pose and then for each traj in EM, QR and MD
        t0 = time.time()
        if not self.initialHasRMSD:
            
#// mdtraj
            comp = mdtraj.load(self.initialPose, top=self.prmtop)
            self.lig_len = len(comp.top.select(f"residue {self.ligand_res}"))
            self.sys_len = comp.n_atoms
            ligand_heavy_atoms = crystalComp.top.select(f"residue {self.ligand_res} and not symbol H")
            comp.superpose(crystalComp, frame=0, atom_indices=range(0,self.sys_len - self.lig_len - 10))
            self.initialRMSD, _ = ALIGN_A_RMSD_B(comp.xyz[0]*10, crystalComp.xyz[0]*10, 
                                                 range(0, self.sys_len - self.lig_len), ligand_heavy_atoms)
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
                ii.calculateRMSD(self.prmtop, crystalComp, self.lig_len, self.ligand_res)
        t1 = time.time()
        displayString = f'  Done with {self.poseName:18s}: '
        for simPrefix in self.simulationPrefixes:
            displayString += f'{len(self.traj[simPrefix]):3d} {simPrefix:3s},'
        displayString += f'  {t1-t0:.3f} s'
        print(displayString)
        
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
#         frameQR =  len(self.trajQR) * self.lengthQR / self.writeCrdIntervalQR
#         lengthQR = len(self.trajQR) * self.lengthQR * self.timeStepQR / 1000
#         frameMD =  len(self.trajMD) * self.lengthMD / self.writeCrdIntervalMD
#         lengthMD = len(self.trajMD) * self.lengthMD * self.timeStepMD / 1000
#         print(frameEM, lengthEM, frameQR, lengthQR, frameMD, lengthMD)
        return frame, length
    

    def inheritPose(self, POSE):
        self.poseName = POSE.poseName
        self.ligandName = POSE.ligandName
        self.successQR = POSE.successQR
        self.ligand_res = POSE.ligand_res
        self.simulationPrefixes = POSE.simulationPrefixes
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
#         self.nrepEM = POSE.nrepEM
#         self.nrepQR = POSE.nrepQR
#         self.nrepMD = POSE.nrepMD
#         self.lengthEM = POSE.lengthEM
#         self.lengthQR = POSE.lengthQR
#         self.lengthMD = POSE.lengthMD
#         self.writeCrdIntervalEM = POSE.writeCrdIntervalEM
#         self.writeCrdIntervalQR = POSE.writeCrdIntervalQR
#         self.writeCrdIntervalMD = POSE.writeCrdIntervalMD
#         self.ioutfmEM = POSE.ioutfmEM
#         self.ioutfmQR = POSE.ioutfmQR
#         self.ioutfmMD = POSE.ioutfmMD
#         self.timeStepEM = POSE.timeStepEM
#         self.timeStepQR = POSE.timeStepQR
#         self.timeStepMD = POSE.timeStepMD
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

#         self.trajEM = []
#         for ii in range(self.nrepEM):
#             trjFile = f'{self.simulationFolder}/{self.poseName}/EM_R{ii+1}{self.crdSuffix}'
#             outFile = f'{self.simulationFolder}/{self.poseName}/EM_R{ii+1}{self.outSuffix}'
#             rstFile = f'{self.simulationFolder}/{self.poseName}/EM_R{ii+1}{self.rstSuffix}'
#             self.trajEM.append(mdgxTrajectory(trjFile, outFile, rstFile, self.ioutfmEM))
#         self.trajQR = []
#         for ii in range(self.nrepQR):
#             trjFile = f'{self.simulationFolder}/{self.poseName}/QR_R{ii+1}{self.crdSuffix}'
#             outFile = f'{self.simulationFolder}/{self.poseName}/QR_R{ii+1}{self.outSuffix}'
#             rstFile = f'{self.simulationFolder}/{self.poseName}/QR_R{ii+1}{self.rstSuffix}'
#             self.trajQR.append(mdgxTrajectory(trjFile, outFile, rstFile, self.ioutfmEM))
#         self.trajMD = []
#         if self.successQR:    
#             for ii in range(self.nrepMD):
#                 trjFile = f'{self.simulationFolder}/{self.poseName}/MD_R{ii+1}{self.crdSuffix}'
#                 outFile = f'{self.simulationFolder}/{self.poseName}/MD_R{ii+1}{self.outSuffix}'
#                 rstFile = f'{self.simulationFolder}/{self.poseName}/MD_R{ii+1}{self.rstSuffix}'
#                 self.trajMD.append(mdgxTrajectory(trjFile, outFile, rstFile, self.ioutfmEM))
        self.initialPose = POSE.initialPose
        self.prmtop = POSE.prmtop
        self.initialHasRMSD = POSE.initialHasRMSD
        if self.initialHasRMSD:
            self.lig_len = POSE.lig_len
            self.initialRMSD = POSE.initialRMSD


        

class Ligand:
    def __init__(self, ligandName, poseNames, ligand_res, settings, success, folderMetadata, simulationPrefixes):

        self.simulationPrefixes = simulationPrefixes
        self.ligandName = ligandName

        self.numPoses = len(poseNames)
        self.poseNames = poseNames
        print(f'Here are the poseNames: {self.poseNames}')
        print(f'{success["QR"]}')
        self.ligand_res = ligand_res # This may be '-1'. See below for error handling
        try:
            self.poseRanks = [x.split('_')[1] for x in poseNames]
        except:
            self.poseRanks = [re.sub('[a-zA-Z]', '', x) for x in poseNames]
        print(f'Here are the poseRanks: {self.poseRanks}')
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

        print(f'We have determined for this ligand its residue number is {self.ligand_res}')

        # Create poses
        
        self.Poses = {} 

        for ii in range(self.numPoses):
            try:
                self.Poses[self.poseRanks[ii]] = Pose(self.poseNames[ii], self.poseRanks[ii], self.qualifiedTruth['QR'][ii], self.ligand_res, 
                                                      settings, folderMetadata, simulationPrefixes=self.simulationPrefixes)
            except:
                self.Poses[self.poseRanks[ii]] = Pose(self.poseNames[ii], self.poseRanks[ii], self.qualifiedTruth['MD'][ii], self.ligand_res, 
                                                      settings, folderMetadata, simulationPrefixes=self.simulationPrefixes)

            
        # For stats
        self.frame = {}
        self.length = {}

        
        self.LTA = LigandTrajectoryAggregate()

        
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

            pool = mp.Pool(processes=20)
            lr = pool.map(Distribute_Lig_pose_RMSD, ((self.Poses[self.poseRanks[ii]], crystalComp) for ii in range(self.numPoses)))
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
                

            
            
#         self.frameEM = accList[0] 
#         self.lengthEM = accList[1]
#         self.frameQR = accList[2]
#         self.lengthQR = accList[3]
#         self.frameMD = accList[4]
#         self.lengthMD = accList[5]
        self.numberOfSimulations = {}
        for simPrefix in self.simulationPrefixes:
            if simPrefix == 'EM' or simPrefix == 'QR':
                self.numberOfSimulations[simPrefix] = self.nrep[simPrefix] * self.numPoses
            else:
                self.numberOfSimulations[simPrefix] = self.nrep[simPrefix] * len(self.qualifiedList[simPrefix])
#         print(self.numEM, self.frameEM, self.lengthEM, self.numQR, self.frameQR, self.lengthQR, self.numMD, self.frameMD, self.lengthMD)
        return self.frame, self.length, self.numberOfSimulations
#         return np.array([self.numEM, self.frameEM, self.lengthEM, self.numQR, self.frameQR, self.lengthQR, self.numMD, self.frameMD, self.lengthMD],dtype=float)
    
    def gatherRMSD(self, mode='MD'):
        RMSD = []
        for ii in range(self.numPoses):
            RMSD.append(self.Poses[self.poseRanks[ii]].gatherRMSD())
        return RMSD
    
    
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
        self.ligandName = LIG.ligandName
        self.numPoses = LIG.numPoses
        self.poseNames = LIG.poseNames
        self.poseRanks = LIG.poseRanks
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
       
        print(f'folderMetadata is {self.folderMetadata}')


        temp = os.listdir(self.inpcrdFolder)
        temp.sort()
        self.inpcrdFiles = [x for x in temp if '.inpcrd' in x]
        self.numPoses = len(self.inpcrdFiles)
        self.poseNames = [x.split('.')[0] for x in self.inpcrdFiles]
        
        temp = os.listdir(self.prmtopFolder)
        temp.sort()
        self.prmtopFiles = [x for x in temp if '.prmtop' in x]
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
        if len(self.sysNames) == 1: # Only one system, all poseNames belong to the sysNames
            for sysName in self.sysNames:
                poseName = [x for x in self.poseNames]
                self.Ligands[sysName] = Ligand(sysName, poseName, self.ligand_res, self.settings, self.success, self.folderMetadata,
                                               simulationPrefixes=self.simulationPrefixes)
        else:
            for sysName in self.sysNames:
                poseName = [x for x in self.poseNames if sysName in x]
                self.Ligands[sysName] = Ligand(sysName, poseName, self.ligand_res, self.settings, self.success, self.folderMetadata,
                                               simulationPrefixes=self.simulationPrefixes)
            
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
                   self.sysNames, self.ligand_res, self.settings, self.success, self.overrideSuccess, self.simulationPrefixes]
        with open(f'{self.saveFolder}/{fname}', 'wb') as f:
            pickle.dump(varlist, f)
            
    def load(self, fname):
        with open(fname, 'rb') as f:
            self.rootFolder, self.referenceFolder, self.structureFolder, self.simulationFolder, \
               self.planningFolder, self.inpcrdFolder, self.prmtopFolder, self.saveFolder, self.folderMetadata, \
               self.inpcrdFiles, self.numPoses, self.poseNames, self.prmtopFiles, self.numSystems, \
               self.sysNames, self.ligand_res, self.settings, self.success, self.overrideSuccess, self.simulationPrefixes = \
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

        
        
