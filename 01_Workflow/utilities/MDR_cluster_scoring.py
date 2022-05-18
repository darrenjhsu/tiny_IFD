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

MDR.calculateRMSD_parallel()

#print(MDR.Ligands)

#print(MDR.Ligands[sys.argv[1].split('/')[-1]].Poses)

min_RMSD = np.inf
for p in MDR.Ligands[sys.argv[1].split('/')[-1]].Poses:
    for t in MDR.Ligands[sys.argv[1].split('/')[-1]].Poses[p].traj['MD']:
        if len(t.RMSD) == 0:
            pass
        elif np.min(t.RMSD) < min_RMSD:
            min_RMSD = np.min(t.RMSD)
            print(min_RMSD)


CE = cpptrajEnergy(rootFolder=f'{sys.argv[1]}/cpptraj')
try:
    append_cpptraj_to_MDR(MDR, CE)
except:
    pass

MDR.saveLigands()


#try:
#    t = gimme_best_pose(MDR, ligand=sys.argv[1].split('/')[-1], top_select=5, plot=False,
#                    #metric="traj.output['extravdW'] + traj.output['extraCoul'] + traj.output['Solvent']",
#                    metric="traj.output['extravdW']",
#                    min_size_multiplier=1, cluster_min_samples=5, eps=0.5, speed=20000,
#                    filter_dist=False, filter_dist_thres=10, show_pose=False, useLigandHTraj=False,
#                    rank=1,outputPDB=True, ligand_res=-1, timing=True)
#
#    print(t)
#except:
#    pass
