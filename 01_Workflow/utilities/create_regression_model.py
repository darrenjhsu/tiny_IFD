
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
from xgboost_config import xgboost_config
from xgboost import XGBRegressor

job = sys.argv[1]
conf = xgboost_config()
n_estimators = conf['n_estimators']
max_depth = conf['max_depth']
gamma = conf['gamma']
norm_cols = conf['norm_cols']
glob_cols = conf['glob_cols']
raw_cols = conf['raw_cols']
suffix = conf['suffix']

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

MDR.calculateRMSD_parallel()


RMSD = []
OUTPUT = []
for ligand in MDR.Ligands.keys():
    # print(ligand)
    for pose in MDR.Ligands[ligand].Poses.keys():
        for idx, traj in enumerate(MDR.Ligands[ligand].Poses[pose].traj['MD']):
            if traj.hasRMSD:
                RMSD.append(traj.RMSD[10:])
                OUTPUT.append(traj.output[10:])

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

normal_data = pd.concat([normalize(OUTPUT[norm_cols]), OUTPUT[glob_cols], OUTPUT[raw_cols]], axis=1)

#conc_data, conc_rmsd = concentrate(normal_data, RMSD, thres=2.5, ratio=0.01)
conc_data, conc_rmsd = normal_data, RMSD
X, yv, y, _ = clean_train_data(conc_data, conc_rmsd, dataset='train')



model = XGBRegressor(n_estimators=n_estimators, max_depth=max_depth, gamma=gamma, verbosity=2)

model.fit(X, yv)
print(model)

model.save_model(f"xgboost_regression_{suffix}.json") 
