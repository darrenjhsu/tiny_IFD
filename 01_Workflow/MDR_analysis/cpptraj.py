
import os, glob
import numpy as np
import pandas as pd


class cpptrajEnergy:
    def __init__(self, rootFolder='.',simType=['EM','QR','MD']):
        self.rootFolder = rootFolder + '/'
        self.resultFolder = self.rootFolder + 'out/'
        
        self.result = {}
        
        self.poseList = [x.split('/')[-1] for x in glob.glob(self.resultFolder + '*')]
        self.poseList.sort()
        self.ligandList = np.unique([x.split('_')[0] for x in self.poseList])
        self.ligandList.sort()
       
        print(f'In cpptraj: self.poseList is {self.poseList} and self.ligandList is {self.ligandList}')

        for ligand in self.ligandList:
            self.result[ligand] = {}
        
        this_ligand = ''
        
        for pose in self.poseList:
            if pose.split('_')[0] != this_ligand:
                print('')
                this_ligand = pose.split('_')[0]
                print(pose.split('_')[0], end=' ')
            print(pose.split('_')[-1], end = ' ')
            ligand = pose.split('_')[0]
            self.result[ligand][pose.split('_')[1]] = {}
            traj_list = glob.glob(self.resultFolder + pose + '/*')
            for t in simType:
                self.result[ligand][pose.split('_')[1]][t] = {}
                self.result[ligand][pose.split('_')[1]][t]['apo'] = {}
                self.result[ligand][pose.split('_')[1]][t]['lig'] = {}
                self.result[ligand][pose.split('_')[1]][t]['holo'] = {}
            for traj in traj_list:
                fileName = traj.split('/')[-1]
                if fileName.split('_')[0] in simType:
                    # print(traj)

                    self.result[ligand][pose.split('_')[1]][fileName.split('_')[0]][fileName.split('_')[-1].split('.')[0]][int(fileName.split('_')[1].replace('R',''))] = pd.read_table(traj, delim_whitespace=True)
                    
                

            # print(len(traj_list))
            # print(traj_list[0])
            
def append_cpptraj_to_MDR(MDR, CE):                
    # Append energy
    fullLength = 0
    for ligand in MDR.Ligands.keys():
        print(ligand)
        for pose in MDR.Ligands[ligand].Poses.keys():
            for idx, traj in enumerate(MDR.Ligands[ligand].Poses[pose].traj['MD']):
                # print(traj.outFile)
                trajEnergy = CE.result[ligand][pose]['MD']
                trajNum = idx+1
                try:
                    traj.output['apovdW14'] = trajEnergy['apo'][trajNum]['ENE_00001[vdw14]'][:-1]
                    traj.output['apovdW'] = trajEnergy['apo'][trajNum]['ENE_00001[vdw]'][:-1]
                    traj.output['ligvdW14'] = trajEnergy['lig'][trajNum]['ENE_00010[vdw14]'][:-1]
                    traj.output['ligvdW'] = trajEnergy['lig'][trajNum]['ENE_00010[vdw]'][:-1]
                    traj.output['holovdW14'] = trajEnergy['holo'][trajNum]['ENE_00019[vdw14]'][:-1]
                    traj.output['holovdW'] = trajEnergy['holo'][trajNum]['ENE_00019[vdw]'][:-1]

                    traj.output['extravdW14'] = traj.output['holovdW14'] - traj.output['apovdW14'] - traj.output['ligvdW14']
                    traj.output['extravdW'] = traj.output['holovdW'] - traj.output['apovdW'] - traj.output['ligvdW']

                    traj.output['apoelec14'] = trajEnergy['apo'][trajNum]['ENE_00001[elec14]'][:-1]
                    traj.output['apoelec'] = trajEnergy['apo'][trajNum]['ENE_00001[elec]'][:-1]
                    traj.output['ligelec14'] = trajEnergy['lig'][trajNum]['ENE_00010[elec14]'][:-1]
                    traj.output['ligelec'] = trajEnergy['lig'][trajNum]['ENE_00010[elec]'][:-1]
                    traj.output['holoelec14'] = trajEnergy['holo'][trajNum]['ENE_00019[elec14]'][:-1]
                    traj.output['holoelec'] = trajEnergy['holo'][trajNum]['ENE_00019[elec]'][:-1]
                    traj.output['extraelec14'] = traj.output['holoelec14'] - traj.output['apoelec14'] - traj.output['ligelec14']
                    traj.output['extraelec'] = traj.output['holoelec'] - traj.output['apoelec'] - traj.output['ligelec']

                    traj.output['ligBond'] = trajEnergy['lig'][trajNum]['ENE_00010[bond]'][:-1]
                    traj.output['ligAngle'] = trajEnergy['lig'][trajNum]['ENE_00010[angle]'][:-1]
                    traj.output['ligDihe'] = trajEnergy['lig'][trajNum]['ENE_00010[dih]'][:-1]
                    traj.output['apoBond'] = trajEnergy['apo'][trajNum]['ENE_00001[bond]'][:-1]
                    traj.output['apoAngle'] = trajEnergy['apo'][trajNum]['ENE_00001[angle]'][:-1]
                    traj.output['apoDihe'] = trajEnergy['apo'][trajNum]['ENE_00001[dih]'][:-1]
                except:
                    pass
                fullLength += len(traj.output['apovdW'])
    print(fullLength)
