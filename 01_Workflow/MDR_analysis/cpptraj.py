
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
            
            
