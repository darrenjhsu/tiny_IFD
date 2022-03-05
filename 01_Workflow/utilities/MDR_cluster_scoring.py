import os, sys, gc, glob

sys.path.insert(1, os.path.join(sys.path[0], '../MDR_analysis/'))

#print(sys.path)
from MDR_base import *
import MDR_plot
import MDR_cluster
from MDR_cluster import gimme_best_pose
import cpptraj
from cpptraj import cpptrajEnergy


pathConfig = {
        'root': sys.argv[1],
        'planningFolder': '../../05_Refinement/script/',
        'configFolder': '../../01_Workflow/utilities/',
        'saveFolder': './pickle/',
        'sysNames': sys.argv[1].split('/')[-1]
}


MDR = MolecularDynamicsRefinement(simulationPrefixes=['EM','QR','MD'],
        **dict(pathConfig))
MDR.createLigands()
print(MDR.Ligands)
MDR.gatherMagnitude()

MDR.calculateRMSD_parallel()

print(MDR.Ligands)

print(MDR.Ligands[sys.argv[1].split('/')[-1]].Poses)

min_RMSD = np.inf
for p in MDR.Ligands[sys.argv[1].split('/')[-1]].Poses:
    for t in MDR.Ligands[sys.argv[1].split('/')[-1]].Poses[p].traj['MD']:
        if len(t.RMSD) == 0:
            pass
        elif np.min(t.RMSD) < min_RMSD:
            min_RMSD = np.min(t.RMSD)
            print(min_RMSD)


#CE = cpptrajEnergy(rootFolder=f'{sys.argv[1]}/cpptraj')
#
##print(CE.result)
#
## Append energy
#fullLength = 0
#for ligand in MDR.Ligands.keys():
#    print(ligand)
#    for pose in MDR.Ligands[ligand].Poses.keys():
#        for idx, traj in enumerate(MDR.Ligands[ligand].Poses[pose].traj['MD']):
#            # print(traj.outFile)
#            trajEnergy = CE.result[ligand][pose]['MD']
#            trajNum = idx+1
#            try:
#                extravdW = \
#                trajEnergy['holo'][trajNum]['ENE_00019[vdw14]'] + trajEnergy['holo'][trajNum]['ENE_00019[vdw]'] - \
#                trajEnergy['apo'][trajNum]['ENE_00001[vdw14]'] - trajEnergy['apo'][trajNum]['ENE_00001[vdw]'] - \
#                trajEnergy['lig'][trajNum]['ENE_00010[vdw14]'] - trajEnergy['lig'][trajNum]['ENE_00010[vdw]']
#                ligvdW = trajEnergy['lig'][trajNum]['ENE_00010[vdw14]'] + trajEnergy['lig'][trajNum]['ENE_00010[vdw]']
#                traj.output['extravdW'] = extravdW[:-1]
#                traj.output['ligvdW'] = ligvdW[:-1]
#                traj.output['ligAngle'] = trajEnergy['lig'][trajNum]['ENE_00010[angle]'][:-1]
#                traj.output['ligBond'] = trajEnergy['lig'][trajNum]['ENE_00010[bond]'][:-1]
#                traj.output['ligDihe'] = trajEnergy['lig'][trajNum]['ENE_00010[dih]'][:-1]
#                extraCoul = \
#                trajEnergy['holo'][trajNum]['ENE_00019[elec14]'] + trajEnergy['holo'][trajNum]['ENE_00019[elec]'] - \
#                trajEnergy['apo'][trajNum]['ENE_00001[elec14]'] - trajEnergy['apo'][trajNum]['ENE_00001[elec]'] - \
#                trajEnergy['lig'][trajNum]['ENE_00010[elec14]'] - trajEnergy['lig'][trajNum]['ENE_00010[elec]']
#                traj.output['extraCoul'] = extraCoul[:-1]
#            except:
#                pass
#            fullLength += len(extravdW)
#print(fullLength)

t = gimme_best_pose(MDR, ligand=sys.argv[1].split('/')[-1], top_select=5, plot=False,
                    #metric="traj.output['extravdW'] + traj.output['extraCoul'] + traj.output['Solvent']",
                    metric="traj.output['vdW']",
                    min_size_multiplier=1, cluster_min_samples=5, eps=0.7, speed=5000,
                    filter_dist=False, filter_dist_thres=10, show_pose=False, useLigandHTraj=False,
                    rank=1,outputPDB=True, ligand_res=-1, timing=True)

print(t)
