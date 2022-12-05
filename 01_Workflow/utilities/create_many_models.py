
import pickle
import numpy as np
import glob, os
import pandas as pd
import os, sys, gc, glob
  
sys.path.insert(1, os.path.join(sys.path[0], '../../01_Workflow/MDR_analysis/'))
#print(sys.path)
from MDR_base import *
import MDR_plot
import MDR_cluster
from MDR_cluster import gimme_best_pose
import cpptraj
from cpptraj import cpptrajEnergy
from xgboost_util import *
from xgboost_config import create_config_list
from xgboost import XGBClassifier

job = sys.argv[1]
config_csv_fname = sys.argv[2]

config = create_config_list(config_csv_fname)


pathConfig = {
            'root': f'../../05_Refinement/{job}/',
            'planningFolder': '../../05_Refinement/script/',
            'configFolder': '../../01_Workflow/utilities/',
            'saveFolder': f'../{job}/pickle/',
            'sysNames': job,
            'overrideSuccess': True,
            'inputFolder': '../../02_Input/'
}


MDR = MolecularDynamicsRefinement(simulationPrefixes=['MD'],
        **dict(pathConfig))
MDR.createLigands()
MDR.gatherMagnitude()
MDR.loadLigands()
#MDR.calculateRMSD_parallel()


RMSD = []
OUTPUT = []
for ligand in MDR.Ligands.keys():
    # print(ligand)
    for pose in MDR.Ligands[ligand].Poses.keys():
        try:
            for idx, traj in enumerate(MDR.Ligands[ligand].Poses[pose].traj['MD']):
                if traj.hasRMSD:
                    RMSD.append(traj.RMSD[10:])
                    OUTPUT.append(traj.output[10:])
        except:
            pass

try:
    RMSD = np.array(RMSD).flatten()
except:
    pass
if RMSD.dtype != float:
    print("We will fix this!")
    RMSD = np.array(sum([list(x) for x in RMSD],[]))
    print(RMSD)
try:
    OUTPUT = pd.concat(OUTPUT)
except:
    pass


norm_cols = config[0]['norm_cols']
glob_cols = config[0]['glob_cols']
raw_cols = config[0]['raw_cols']
normal_data = pd.concat([normalize(OUTPUT[norm_cols]), OUTPUT[glob_cols], OUTPUT[raw_cols]], axis=1)
    
conc_data, conc_rmsd = concentrate(normal_data, RMSD, thres=2.5, ratio=0.01)
X, yv, y, _ = clean_train_data(conc_data, conc_rmsd, dataset='train')

os.makedirs('models',exist_ok=True) 

for conf in config:
    n_estimators = conf['n_estimators']
    max_depth = conf['max_depth']
    gamma = conf['gamma']
    colsample_bytree = conf['colsample_bytree']
    subsample = conf['subsample']
    suffix = conf['suffix']
    
        
    model = XGBClassifier(n_estimators=n_estimators, max_depth=max_depth, gamma=gamma, 
                          colsample_bytree=colsample_bytree, subsample=subsample, verbosity=0)
    
    model.fit(X, y)
    print(model)
    model.save_model(f"models/xgboost_{suffix}.json") 
