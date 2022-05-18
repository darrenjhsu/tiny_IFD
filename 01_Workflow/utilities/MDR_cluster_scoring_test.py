import os, sys, gc, glob

sys.path.insert(1, os.path.join(sys.path[0], '../MDR_analysis/'))

#print(sys.path)
from MDR_base import *
import MDR_plot
import MDR_cluster
from MDR_cluster import gimme_best_pose
import cpptraj
from cpptraj import cpptrajEnergy, append_cpptraj_to_MDR


pathConfig = {
        'root': sys.argv[1],
        'planningFolder': '../../05_Refinement/script/',
        'configFolder': '../../01_Workflow/utilities/',
        'saveFolder': './pickle/',
        'sysNames': sys.argv[1].split('/')[-1],
        'overrideSuccess': True,
        'inputFolder': '../../02_Input/'
}


MDR = MolecularDynamicsRefinement(simulationPrefixes=['MD'],
        **dict(pathConfig))
MDR.createLigands()
#print(MDR.Ligands)
MDR.gatherMagnitude()

MDR.loadLigands()
exit()

