import os
import sys
import json
import copy
import mdtraj
import numpy as np
import time
import pandas as pd
import pickle
import mdtraj as md
import multiprocessing as mp
try:
    import cupy as cp
    cudaExists = True
except ImportError as e:
    cudaExists = False
    print("Can't load CuPy, fall back to numba")
from numba import jit, prange
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from scipy.stats import rankdata



@jit(nopython=True)
def pureRMSD(P, Q):
    # Assume P and Q are aligned first
    diff = P - Q
    N = len(P)
    return np.sqrt((diff * diff).sum() / N)

@jit(nopython=True,parallel=True)
def pureRMSDrow(P, Q):
    # Assume P and Q are aligned first, P is a ref, Q is an array of P-like matrices
    rmsd = np.zeros(len(Q))
    diff = np.zeros_like(Q)
    for ii in prange(len(Q)):
#         diff[ii] = P - Q[ii] 
        N = len(P)
#         rmsd[ii] = np.sqrt((diff[ii] * diff[ii]).sum() / N)
        rmsd[ii] = np.sqrt(((P-Q[ii]) * (P-Q[ii])).sum() / N)
    return rmsd

@jit(nopython=True,parallel=True)
def pureRMSDself(P):
    # Assume P and Q are aligned first, P is a ref, Q is an array of P-like matrices
    rmsd = np.zeros((len(P),len(P)))
    lenP = len(P[0])
    for ii in prange(len(P)):
#         if ii % 50 == 0:
#             print(ii)
        for jj in range(len(P)):
            Q = P[ii] - P[jj]
            rmsd[ii,jj] = np.sqrt(np.sum(Q * Q) / lenP)
    return rmsd

def pureRMSDcupy(P):
    Pcp = cp.array(P)
#     print(P.shape)
    rmsd = cp.zeros((len(P),len(P)))
#     print(rmsd.shape)
    lenP = len(P[0])
    for ii in range(len(P)):
#         if ii % 400 == 0:
#             print(ii)
        Q = Pcp - Pcp[ii]
#         print(Q.shape)
        rmsd[ii] = cp.sqrt(cp.sum(cp.sum(Q * Q, axis=2), axis=1) / lenP)
    rmsd_local = rmsd.get()
    del rmsd
    return rmsd_local
    
# traj = MDR.Ligands['Mpro-x10959'].Poses['3'].traj['MD'][0]
# len_traj = len(traj.ligandTrajectory)
# rmsd_cluster = np.zeros((len_traj,len_traj))
# for ii in range(len_traj):
#     for jj in range(len_traj):
#         rmsd_cluster[ii,jj] = pureRMSD(traj.ligandTrajectory[ii],traj.ligandTrajectory[jj])

def gimme_best_pose(MDR, ligand='Mpro-x10959', metric='vdW', ligand_res=56, filter_dist=False, filter_dist_thres=2.5, 
                    top_select=5, plot=True, cluster_min_samples=3, eps=0.4, min_size_multiplier=1, speed=15000, show_pose=False,
                    simRound='MD',rank=None, outputPDB=False, returnComp=True, returnFullActualComp=False, 
                    useLigandHTraj=True, timing=False, sizeEffect=False, warmUpFrames=10, pre_classifier_model=None,
                    use_consensus_ranking=False, discard_pose=['0']):
        
        
    simRound = simRound
    ligand_res = str(ligand_res)
    warmup = warmUpFrames

    if ligand_res == '-1':
        initialComp = mdtraj.load(MDR.Ligands[ligand].initialPose, top=MDR.Ligands[ligand].prmtop)
        ligand_res = str(initialComp.n_residues-1)
    print(f'We have determined for this ligand its residue number is {ligand_res}')

    #if not filter_dist:
    #    filter_dist = True
    #    filter_dist_thres = 500
    
    if timing:
        t0 = time.time()
    
    # perform full ligand sample clustering to pick best poses based on the metric
    trajAgg = []
    rmsdAgg = []
    metricAgg = []

    for pp in MDR.Ligands[ligand].Poses.keys():
        if pp in discard_pose:
            continue
        poseNum = str(pp)

        if MDR.Ligands[ligand].Poses[poseNum].successQR:

            for traj in MDR.Ligands[ligand].Poses[poseNum].traj[simRound]:
                if traj.hasTrjFile and traj.hasRMSD:
                    # Find out the minimum length recorded
                    minLength = np.min([len(traj.ligandTrajectoryH), len(eval(metric).values)])
                    if useLigandHTraj:
                        trajAgg.append(traj.ligandTrajectoryH[warmup:minLength])
                    else:
                        trajAgg.append(traj.ligandTrajectory[warmup:minLength])
                    rmsdAgg.append(traj.RMSD[warmup:minLength])
                    # metricAgg.append(traj.output[metric].values[warmup:minLength])
                    metricAgg.append(eval(metric).values[warmup:minLength]) # Let's try this
    trajAgg = np.array(np.concatenate(trajAgg))
    rmsdAgg = np.array(np.concatenate(rmsdAgg)).flatten()
    metricAgg = np.array(np.concatenate(metricAgg)).flatten()

    len_traj = len(trajAgg)

    stride = max(len_traj // speed,1)
    print(f'Stride factor is {stride} (number of frames: {len_traj//stride})')
    
    if timing:
        t1 = time.time()
        print(f'Aggregating traj / rmsd took {t1-t0:.3f} s')
        t0 = time.time()
    
    # if cudaExists:
    #     rmsd_cluster = pureRMSDcupy(trajAgg[::stride])
    # else:
    #     rmsd_cluster = pureRMSDself(trajAgg[::stride])
    print(f'RMSD calculation done on {len(trajAgg[::stride])} frames.')
    
    if timing:
        t1 = time.time()
        print(f'Calculating cluster rmsd took {t1-t0:.3f} s')
        t0 = time.time()
    
    
    cluster_min_samples = cluster_min_samples
    # clustering = DBSCAN(eps=eps, min_samples=cluster_min_samples, metric='precomputed').fit(rmsd_cluster)
    
    clustering = DBSCAN(eps=eps * np.sqrt(trajAgg.shape[1]), min_samples=cluster_min_samples).fit(trajAgg[::stride].reshape(len(trajAgg[::stride]), -1))

    print("Clustering done")

    if timing:
        t1 = time.time()
        print(f'DBSCAN took {t1-t0:.3f} s')
        t0 = time.time()
    
    simLigTraj = []
    simTraj = []
    simOutput = []
    simOutputExtravdW = []
    simOutputChangeSASA = []
    simLocation = []
    simCumLength = [0] # This line needs check
    simCumLength = []
    
    if pre_classifier_model is not None:
        X_test = None

    for pp in MDR.Ligands[ligand].Poses.keys():
        if pp in discard_pose:
            continue
        poseNum = str(pp)

        if MDR.Ligands[ligand].Poses[poseNum].successQR:

            for jj in range(0, len(MDR.Ligands[ligand].Poses[poseNum].traj[simRound])):
                for simType in [simRound]:

                    traj = MDR.Ligands[ligand].Poses[poseNum].traj[simType][jj]
                    if traj.hasRMSD and traj.hasTrjFile:
                        minLength = np.min([len(traj.ligandTrajectoryH), len(eval(metric).values)])
                        simLigTraj.append(traj.ligandTrajectoryH[warmup:minLength])
                        simTraj.append(traj.RMSD[warmup:minLength])
                        simOutput.append(eval(metric).values[warmup:minLength])
                        simOutputExtravdW.append(traj.output['extravdW'].values[warmup:minLength])
                        simOutputChangeSASA.append(traj.output['changeSASA'].values[warmup:minLength])
                        if pre_classifier_model is not None:
                            if X_test is None:
                                X_test = traj.output.iloc[warmup:minLength]
                            else:
                                X_test = pd.concat([X_test, traj.output.iloc[warmup:minLength]])
                        simLocation.append(traj)
                        if len(simCumLength) > 0:
                            simCumLength.append(simCumLength[-1]+len(traj.RMSD[warmup:minLength]))
                        else:
                            simCumLength = [len(traj.RMSD[warmup:minLength])]
    simLigTraj = np.concatenate(simLigTraj)
    simTraj = np.array(np.concatenate(simTraj)).flatten()
    simOutput = np.array(np.concatenate(simOutput)).flatten()
    simOutputExtravdW = np.array(np.concatenate(simOutputExtravdW).flatten())
    simOutputChangeSASA = np.array(np.concatenate(simOutputChangeSASA).flatten())
    simCumLength = np.array(simCumLength)

    
    
    if timing:
        t1 = time.time()
        print(f'Aggregating traj / rmsd, second round, took {t1-t0:.3f} s')
        t0 = time.time()
    
    
    
    if plot:
        plt.figure(figsize=(12,8))
        plt.scatter(simTraj[::stride],clustering.labels_,s=5)
        plt.yticks(np.sort(np.unique(clustering.labels_)))
        plt.ylabel('Cluster index',fontsize=15)
        plt.xlabel('RMSD (Å)',fontsize=15)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.xlim([1, 15])
#         plt.ylim([-1, 250])


    # For additional classificaiton
    if pre_classifier_model is not None:
        
        from xgboost import XGBClassifier
        logreg = XGBClassifier()
        logreg.load_model(pre_classifier_model)
        
        X_test = (X_test - np.mean(X_test, axis=0)) / np.std(X_test, axis=0)
        X_mask = ~np.any(np.abs(X_test) > 5, axis=1)
        # print(X_test.shape, X_mask.shape)
        X_test = X_test[X_mask]
        
        print(len(X_test))
        # feature_cols = ['Etot', 'EKtot', 'EPtot', 'Bond', 'Angle', 'Dihedral', 'Elec', 'vdW', 
        #         'Solvent', 'extravdW', 'ligvdW', 'ligAngle', 'ligBond', 'ligDihe', 
        #         'extraCoul', 'proteinRMSD', 'changeSASA']
        # feature_cols = ['Dihedral', 'Elec', 'vdW', 
        #         'Solvent', 'extravdW', 'ligvdW', 'ligAngle', 'ligBond', 'ligDihe', 
        #         'extraCoul', 'proteinRMSD', 'changeSASA']
        feature_cols = ['Dihedral', 'Elec', 'vdW', 
                'Solvent', 'extravdW', 'ligvdW', 'ligAngle', 'ligBond', 'ligDihe', 
                'extraCoul', 'proteinRMSD', 'changeSASA', 'HBond', 'Contact']
        X_test = X_test[feature_cols]
        y_proba = logreg.predict_proba(X_test)[::,1]
        y_pred = y_proba > 0.5
        all_pick = simTraj[X_mask][np.where(y_pred == True)[0]]
        good_pick = np.sum(all_pick < 2.5)

        print(f'ML model predicted {len(all_pick)} picks!')
        print(f'Of which, there are {good_pick} ({good_pick / len(all_pick) * 100:.1f} %) picks with <2.5 A RMSD.')
    
        plt.figure(dpi=150)
        plt.hist(simTraj[X_mask][y_pred], bins=40)
        plt.show()
        
        trajML = trajAgg[X_mask][y_pred]
        
        if cudaExists:
            trajML_cluster = pureRMSDcupy(trajML)
        else:
            trajML_cluster = pureRMSDself(trajML)
        print(f'ML RMSD calculation done on {len(trajML)} frames.')
        MLclustering = DBSCAN(eps=2, min_samples=len(trajML)//100, metric='precomputed').fit(trajML_cluster)
        # MLclustering = DBSCAN(eps=2.5, min_samples=2).fit(X_test[y_pred])
        print(f'There are {len(np.unique(MLclustering.labels_))} clusters in ML picks.')
        
        MLrepresentative = {}
        for idx, ii in enumerate(np.unique(MLclustering.labels_)):
            loc = np.where(MLclustering.labels_ == ii)[0]
            print(f'Cluster {idx}: lowest output = {np.min(simOutput[X_mask][y_pred][loc]):.3f}, median output = {np.median(simOutput[X_mask][y_pred][loc]):.3f}, median proba = {np.median(y_proba[y_pred][loc]):.3f}, median RMSD = {np.median(all_pick[loc]):.3f}, RMSD of lowest output {all_pick[loc][np.argmin(simOutput[X_mask][y_pred][loc])]:.3f}')
            # MLrepresentative[ii] = np.argmax(y_proba[y_pred][loc])
            MLrepresentative[ii] = np.argmax(simOutput[X_mask][y_pred][loc])
            # MLrepresentative[ii] = np.argsort(y_proba[y_pred][loc])[len(y_proba[y_pred][loc])//2]
        
        
    results = []
    representative = {}
    
    # Create distance filter
    
    usingCrystalComp = False
    if filter_dist:
        try:
            crystalComp = mdtraj.load(MDR.Ligands[ligand].crystalPose, top=MDR.Ligands[ligand].prmtop)
            usingCrystalComp = True
        except:
            if rank is not None:
                crystalComp = mdtraj.load(f'{MDR.inpcrdFolder}/{ligand}_{rank}.inpcrd', top=MDR.Ligands[ligand].prmtop)
                print(f'Using ranked reference')
                print(f'No crystal comp found - using pose {rank} instead. RMSD will be meaningless')
            else:
                crystalComp = mdtraj.load(f'{MDR.inpcrdFolder}/rank{rank}.inpcrd', top=MDR.Ligands[ligand].prmtop)
                print(f'Using rank {rank} as reference')
                print(f'No crystal comp found - using pose {rank} instead. RMSD will be meaningless')

        pro = crystalComp.top.select(f'not residue {ligand_res} and not symbol H')
        close_atoms = np.array([ 45, 106, 107, 167, 168, 170, 175, 176, 177, 178, 179, 180, 181, 182, 232, 259, 381, 386, 387, 388])
        COM_active_site = crystalComp.xyz[0][pro][close_atoms].mean(0)
    else:
        try:
            crystalComp = mdtraj.load(MDR.Ligands[ligand].crystalPose, top=MDR.Ligands[ligand].prmtop)
            usingCrystalComp = True
        except:
            pass


    if timing:
        t1 = time.time()
        print(f'Creating filter took {t1-t0:.3f} s')
        t0 = time.time()
    
    
    # Filter results
    if use_consensus_ranking:
        extravdW = []
        changeSASA = []        
        for idx, ii in enumerate(np.unique(clustering.labels_)):

            if np.sum(clustering.labels_ == ii) >= min_size_multiplier*cluster_min_samples and ii >= 0:
                selected = simLigTraj[::stride][clustering.labels_ == ii]
                centroid = np.mean(selected,axis=0)
                distance_to_centroid = pureRMSDrow(centroid, selected)
                representative[ii] = np.argmin(distance_to_centroid)
                extravdW.append(np.mean(simOutputExtravdW[::stride][clustering.labels_ == ii]))
                changeSASA.append(np.mean(simOutputChangeSASA[::stride][clustering.labels_ == ii]))
            else:
                extravdW.append(1e8)
                changeSASA.append(1e8)
        # print(extravdW)
        extravdW_rank = rankdata(extravdW, method='min').astype(int)
        changeSASA_rank = rankdata(changeSASA, method='min').astype(int)
        consensusRank = extravdW_rank + changeSASA_rank - 2
        # for idx, ii in enumerate(np.unique(clustering.labels_)):
        #     print(ii, extravdW_rank[idx], changeSASA_rank[idx])
        for idx, ii in enumerate(np.unique(clustering.labels_)):
            if np.sum(clustering.labels_ == ii) >= min_size_multiplier*cluster_min_samples and ii >= 0:
                results.append([ii, 
                                np.sum(clustering.labels_ == ii), 
                                np.mean(simTraj[::stride][clustering.labels_ == ii]), 
                                # np.mean(simOutput[::stride][clustering.labels_ == ii]),
                                consensusRank[idx],
                                np.inf])
                # print([ii, 
                #                 np.sum(clustering.labels_ == ii), 
                #                 np.mean(simTraj[::stride][clustering.labels_ == ii]), 
                #                 # np.mean(simOutput[::stride][clustering.labels_ == ii]),
                #                 consensusRank[idx],
                #                 np.inf])
                
    else:
        for ii in np.unique(clustering.labels_):

            if np.sum(clustering.labels_ == ii) >= min_size_multiplier*cluster_min_samples and ii >= 0:
                selected = simLigTraj[::stride][clustering.labels_ == ii]
                centroid = np.mean(selected,axis=0)


                distance_to_centroid = pureRMSDrow(centroid, selected)
                representative[ii] = np.argmin(distance_to_centroid)
                if filter_dist:
                    ligandActiveSiteDistance = np.sqrt(((simLigTraj[::stride][clustering.labels_ == ii].mean(1) - \
                                                         COM_active_site*10)**2).sum(1)).mean(0)

                    if ligandActiveSiteDistance < filter_dist_thres:
                        try:
                            results.append([ii, 
                                            np.sum(clustering.labels_ == ii), 
                                            np.mean(simTraj[::stride][clustering.labels_ == ii]), 
                                            np.mean(simOutput[::stride][clustering.labels_ == ii]), 
                                            ligandActiveSiteDistance])
                        except:
                            print(f'Error at cluster {ii}, simTraj has the shape {simTraj.shape}, simOutput has the shape {simOutput.shape}')
                            raise

                else:
                    results.append([ii, 
                                    np.sum(clustering.labels_ == ii), 
                                    np.mean(simTraj[::stride][clustering.labels_ == ii]), 
                                    np.mean(simOutput[::stride][clustering.labels_ == ii]),
                                    np.inf])

                
#             if plot:
#                 plt.scatter(np.mean(simTraj[::stride][clustering.labels_ == ii]), ii, c='r',s=8)
            
    results = np.array(results)
    #print(results)

    if timing:
        t1 = time.time()
        print(f'Filtering results took {t1-t0:.3f} s')
        t0 = time.time()



    if plot:
        plt.figure(figsize=(8,6))
        minOutput = np.inf
        maxOutput = -np.inf
        for row in results:
            ii = row[0]

            plt.scatter(simTraj[::stride][clustering.labels_ == ii],
                        simOutput[::stride][clustering.labels_ == ii],s=2)
            
            plt.scatter(np.mean(simTraj[::stride][clustering.labels_ == ii]), 
                        np.mean(simOutput[::stride][clustering.labels_ == ii]), 
                        s = np.sum(clustering.labels_ == ii), 
                        color = [0,0,1,0.05], edgecolor='r', linewidths=2.4)

            if np.mean(simOutput[::stride][clustering.labels_ == ii]) < minOutput:
                minOutput = np.mean(simOutput[::stride][clustering.labels_ == ii])
            if np.mean(simOutput[::stride][clustering.labels_ == ii]) > maxOutput:
                maxOutput = np.mean(simOutput[::stride][clustering.labels_ == ii])
            
            plt.xlabel('RMSD (Å)', fontsize=15)
            plt.ylabel('Raw vdW energy (kcal/mol)', fontsize=15)
            plt.tick_params(axis='both', which='major', labelsize=15)

        outputRange = (maxOutput - minOutput) * 0.2 + minOutput
        
        plt.plot([0, 15], [outputRange, outputRange], 'r', alpha=0.5, linewidth = 2)
        plt.xlim([0, 15])
            
    if plot:

        plt.figure(figsize=(11,5))
        plt.grid(alpha=0.3)
        plt.scatter(results[:,2], results[:,3], s = results[:,1], 
                    color = [0,0,1,0.05], edgecolor='r', linewidths=2.4)
        plt.xlabel('RMSD (Å)', fontsize=15)
        plt.ylabel('Raw vdW energy (kcal/mol)', fontsize=15)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.plot([0, 15], [outputRange, outputRange], 'r', alpha=0.5, linewidth = 2)
        plt.xlim([0, 15])
        
#         plt.figure(figsize=(12,8))
#         plt.scatter(results[:,4], results[:,3], s = results[:,1], color = [0,0,1,0.02], edgecolor='r', linewidths=1.5)
#         plt.xlabel('Dist (Å)', fontsize=15)
#         plt.ylabel('Raw vdW energy (kcal/mol)', fontsize=15)
#         plt.tick_params(axis='both', which='major', labelsize=15)

#         plt.figure(figsize=(12,8))
#         plt.scatter(results[:,4], results[:,2], s = results[:,1], color = [0,0,1,0.02], edgecolor='r', linewidths=1.5)
#         plt.xlabel('Dist (Å)', fontsize=15)
#         plt.ylabel('RMSD (Å)', fontsize=15)
#         plt.tick_params(axis='both', which='major', labelsize=15)

    
    if top_select == 'all':
        top_select = len(results)

        
    if sizeEffect:
        cluster_selection = results[:,0][np.argsort(results[:,3]-np.log(results[:,1])/np.log(np.e))[:top_select]].astype(int)
    else:
        cluster_selection = results[:,0][np.argsort(results[:,3])[:top_select]].astype(int)
        print(np.sort(results[:,3])[:top_select])
    print(f'Selected clusters are: {cluster_selection}')
    
    rmsd_selection = [simTraj[::stride][clustering.labels_ == ii][representative[ii]] for ii in cluster_selection]
    lowest_in_cluster = [np.min(simTraj[::stride][clustering.labels_ == ii]) for ii in cluster_selection]
    if not usingCrystalComp:
        if rank is not None:
            print(f'RMSD of the clusters : {rmsd_selection}' + 
                  f' NOTE: these values are meaningless - they are compared to predicted pose {rank}')
        else:
            print(f'RMSD of the clusters : {rmsd_selection}' +
                  f' NOTE: these values are meaningless - they are compared to predicted pose 1')
    elif min(rmsd_selection) < 2.5:
        print(f'RMSD of the clusters : {rmsd_selection} PASS')
    else:
        print(f'RMSD of the clusters : {rmsd_selection} ')
    print(f'Lowest RMSD in each cluster: {lowest_in_cluster}')
    if usingCrystalComp:
        print(f'Lowest possible median RMSD of a cluster : {np.min(results[:,2])}')
    
#     print(simCumLength)


    if timing:
        t1 = time.time()
        print(f'Plotting results took {t1-t0:.3f} s')
        t0 = time.time()

        
    # Output ML picked conformers
    if pre_classifier_model is not None:
        comp = {}
        actualComp = {}
        comp_selection = [simLigTraj[X_mask][y_pred][MLclustering.labels_ == ii][MLrepresentative[ii]] for ii in np.unique(MLclustering.labels_)]   
        loc_selection = [np.arange(len(simLigTraj))[X_mask][y_pred][MLclustering.labels_ == ii][MLrepresentative[ii]] for ii in np.unique(MLclustering.labels_)]
    #     print(loc_selection)
        for idx, ii in enumerate(np.unique(MLclustering.labels_)):
            if usingCrystalComp:
                ref_pose = mdtraj.load(f'{MDR.referenceFolder}/complex.inpcrd',top=f'{MDR.referenceFolder}/complex.prmtop')
            else:
                if rank is not None:
                    ref_pose = mdtraj.load(f'{MDR.inpcrdFolder}/rank{rank}.inpcrd', top=MDR.Ligands[ligand].prmtop)
                else:
                    ref_pose = mdtraj.load(f'{MDR.inpcrdFolder}/rank1.inpcrd', top=MDR.Ligands[ligand].prmtop)
            if usingCrystalComp:
                comp[ii] = mdtraj.load(f'{MDR.referenceFolder}/complex.inpcrd',top=f'{MDR.referenceFolder}/complex.prmtop')
            else:
                if rank is not None:
                    comp[ii] = mdtraj.load(f'{MDR.inpcrdFolder}/rank{rank}.inpcrd', top=MDR.Ligands[ligand].prmtop)
                else:
                    comp[ii] = mdtraj.load(f'{MDR.inpcrdFolder}/rank1.inpcrd', top=MDR.Ligands[ligand].prmtop)
            comp[ii].xyz[0][-len(comp_selection[0]):] = simLigTraj[X_mask][y_pred][MLclustering.labels_ == ii][MLrepresentative[ii]]/10
            comp[ii].superpose(ref_pose, atom_indices=range(0,700))
            actualLocation = np.where(simCumLength > loc_selection[idx])[0][0]
            actualFile = simLocation[actualLocation]
            actualComp[ii] = mdtraj.load_mdcrd(f'{actualFile.trjFile}',top=MDR.Ligands[ligand].prmtop)
            actualComp[ii].superpose(ref_pose, atom_indices=range(0,700))
            if actualLocation == 0:
                actualFrame = loc_selection[idx] + 10
            else:
                actualFrame = loc_selection[idx] - simCumLength[actualLocation-1] + 10
            actualComp[ii] = actualComp[ii][actualFrame]


        # Output PDB
        if outputPDB:
            try:
                os.mkdir('Refined_Poses')
            except:
                pass
            try:
                os.mkdir(f'Refined_Poses/{ligand}')
            except:
                pass
            for idx, ii in enumerate(np.unique(MLclustering.labels_)):
                # There's some bug with actualComp, so we output comp for now
                print(actualComp[ii][0])
                actualComp[ii][0].save(f'Refined_Poses/{ligand}/MLPose{idx}.pdb')
                # Save ligand only pdb
                lig_select = actualComp[ii][0].top.select(f'residue {ligand_res}')
                pro_select = actualComp[ii][0].top.select(f'not residue {ligand_res}')
                pro_select2 = actualComp[ii][0].top.select(f'not residue {ligand_res} and not resname NME ACE')
                actualComp[ii][0].atom_slice(lig_select).save(f'Refined_Poses/{ligand}/MLPose{idx}_lig.pdb') 
                actualComp[ii][0].atom_slice(pro_select).save(f'Refined_Poses/{ligand}/MLPose{idx}_pro.pdb')
                actualComp[ii][0].atom_slice(pro_select2).save(f'Refined_Poses/{ligand}/MLPose{idx}_pro2.pdb')
    #             comp[ii][0].save(f'Refined_Poses/{ligand}/Pose{idx}.pdb')
            # ref_pose[0].save(f'Refined_Poses/ref.pdb')
        
        
        

    # Set the comps for each selected conformers
    if returnComp:
        comp = {}
        actualComp = {}
        comp_selection = [simLigTraj[::stride][clustering.labels_ == ii][representative[ii]] for ii in cluster_selection]   
        loc_selection = [np.arange(len(simLigTraj))[::stride][clustering.labels_ == ii][representative[ii]] for ii in cluster_selection]
    #     print(loc_selection)
        for idx, ii in enumerate(cluster_selection):
            if usingCrystalComp:
                ref_pose = mdtraj.load(f'{MDR.referenceFolder}/complex.inpcrd',top=f'{MDR.referenceFolder}/complex.prmtop')
            else:
                if rank is not None:
                    ref_pose = mdtraj.load(f'{MDR.inpcrdFolder}/rank{rank}.inpcrd', top=MDR.Ligands[ligand].prmtop)
                else:
                    ref_pose = mdtraj.load(f'{MDR.inpcrdFolder}/rank1.inpcrd', top=MDR.Ligands[ligand].prmtop)
            if usingCrystalComp:
                comp[ii] = mdtraj.load(f'{MDR.referenceFolder}/complex.inpcrd',top=f'{MDR.referenceFolder}/complex.prmtop')
            else:
                if rank is not None:
                    comp[ii] = mdtraj.load(f'{MDR.inpcrdFolder}/rank{rank}.inpcrd', top=MDR.Ligands[ligand].prmtop)
                else:
                    comp[ii] = mdtraj.load(f'{MDR.inpcrdFolder}/rank1.inpcrd', top=MDR.Ligands[ligand].prmtop)
            comp[ii].xyz[0][-len(comp_selection[0]):] = simLigTraj[::stride][clustering.labels_ == ii][representative[ii]]/10
            comp[ii].superpose(ref_pose, atom_indices=range(0,700))
            actualLocation = np.where(simCumLength > loc_selection[idx])[0][0]
            actualFile = simLocation[actualLocation]
            actualComp[ii] = mdtraj.load_mdcrd(f'{actualFile.trjFile}',top=MDR.Ligands[ligand].prmtop)
            actualComp[ii].superpose(ref_pose, atom_indices=range(0,700))
            if actualLocation == 0:
                actualFrame = loc_selection[idx] + 10
            else:
                actualFrame = loc_selection[idx] - simCumLength[actualLocation-1] + 10
            actualComp[ii] = actualComp[ii][actualFrame]

        if timing:
            t1 = time.time()
            print(f'Fetching comps took {t1-t0:.3f} s')
            t0 = time.time()

        # Output PDB
        if outputPDB:
            try:
                os.mkdir('Refined_Poses')
            except:
                pass
            try:
                os.mkdir(f'Refined_Poses/{ligand}')
            except:
                pass
            for idx, ii in enumerate(cluster_selection):
                # There's some bug with actualComp, so we output comp for now
                print(actualComp[ii][0])
                actualComp[ii][0].save(f'Refined_Poses/{ligand}/Pose{idx}.pdb')
                # Save ligand only pdb
                lig_select = actualComp[ii][0].top.select(f'residue {ligand_res}')
                pro_select = actualComp[ii][0].top.select(f'not residue {ligand_res}')
                pro_select2 = actualComp[ii][0].top.select(f'not residue {ligand_res} and not resname NME ACE')
                actualComp[ii][0].atom_slice(lig_select).save(f'Refined_Poses/{ligand}/Pose{idx}_lig.pdb') 
                actualComp[ii][0].atom_slice(pro_select).save(f'Refined_Poses/{ligand}/Pose{idx}_pro.pdb')
                actualComp[ii][0].atom_slice(pro_select2).save(f'Refined_Poses/{ligand}/Pose{idx}_pro2.pdb')
    #             comp[ii][0].save(f'Refined_Poses/{ligand}/Pose{idx}.pdb')
            ref_pose[0].save(f'Refined_Poses/ref.pdb')


        if timing:
            t1 = time.time()
            print(f'Outputting PDB took {t1-t0:.3f} s')
            t0 = time.time()



        if show_pose: # Show selected poses
            view = nglview.NGLWidget() 
            view.camera = 'orthographic'
    #         print(loc_selection)
            c = {}


    #         ac = {}
            colors = ['red','orange','yellow','green','blue','purple']
            for idx, ii in enumerate(cluster_selection):

    #             ac[ii] = view.add_trajectory(actualComp[ii][actualFrame])
    #             ac[ii].clear()
    #             ac[ii].add_cartoon(selection="protein")
    #             ac[ii].add_licorice(selection="protein",opacity=0.3,color='white')
    # For poster
    #             comp[ii].xyz[0][-len(comp_selection[0]):][:,0] -= idx // 3 * 1.5
    #             comp[ii].xyz[0][-len(comp_selection[0]):][:,1] -= idx % 3

                c[ii] = view.add_trajectory(comp[ii])
                c[ii].clear()
                if idx == 0:

                    c[ii].add_cartoon(selection="protein")
                    c[ii].add_licorice(selection='33')
                    c[ii].add_licorice(selection="protein",opacity=0.25)
    #                 c[ii].add_surface(selection="protein", opacity=0.3) # Not often used

                c[ii].add_licorice(selection=ligand_res,color=colors[idx%6],opacity=0.6)
    #             c[ii].add_licorice(selection=ligand_res,color=colors[idx%6],opacity=0.8)
    #             c[ii].add_ball_and_stick(selection=ligand_res,opacity=1)

        # Reference pose
            c2 = view.add_trajectory(ref_pose[-1])
            c2.clear()
            c2.add_ball_and_stick(selection=ligand_res,width=0.5)

            if timing:
                t1 = time.time()
                print(f'Processing nglview took {t1-t0:.3f} s')
                t0 = time.time()


            #del rmsd_cluster
            return cluster_selection, rmsd_selection, view, comp, actualComp
        else:
            if returnFullActualComp:
                # Aggregate full actual comps into one file
                for ii in cluster_selection:
                    selected = np.where(clustering.labels_ == ii)[0]
                    print(f'For cluster {ii} there are {len(selected)} snapshots')
                    print(f'{selected[:20]}')

            else:
                return cluster_selection, rmsd_selection, None, comp, actualComp
    else: # Not returnComp
        return cluster_selection, rmsd_selection, lowest_in_cluster, None, None


