    def calculateRMSD(self, prmtop, crystalComp, lig_len):
        try:
            if (not self.hasRMSD):
                Profile = np.random.random() < 0.00
                if Profile:
                    t0 = time.time()

                if self.ioutfm == 1: # Binary netcdf trajectory
#// MDAnalysis
#                     comp = mda.Universe(prmtop, self.trjFile, in_memory=True)
#// mdtraj
                    comp = mdtraj.load_netcdf(self.trjFile, top=prmtop)
                elif self.ioutfm == 0: # ASCII MDCRD trajectory, loading is slower
#// MDAnalysis
#                     comp = mda.Universe(prmtop, self.trjFile, format='mdcrd', in_memory=True)
#// mdtraj
#                     comp = mdtraj.load_mdcrd(self.trjFile, top=prmtop)
#// Custom read
                    systemLen = crystalComp.n_atoms
                    comp = readMDCRD(self.trjFile, systemLen)
                if Profile:
                    t1 = time.time()
#// mdtraj
                self.trjLength = len(comp)
#// MDAnalysis
#                 self.trjLength = len(comp.trajectory) # mdanalysis
    
                if Profile:
                    t2 = time.time()
#// MDAnalysis
#                 R = mda.analysis.rms.RMSD(comp, crystalComp)
#                 R.run()
#                 self.RMSD = np.sqrt((R.rmsd.T[2])**2/lig_len*len(crystalComp.atoms))

#// Custom RMSD with MDAnalysis
#                 R = []
#                 systemLen = len(crystalComp.atoms)
#                 referenceXYZ = crystalComp.atoms.positions
#                 for ts in comp.trajectory:
#                     R.append(ALIGN_A_RMSD_B(comp.atoms.positions, referenceXYZ,
#                                             range(systemLen-lig_len-50, systemLen-lig_len), range((systemLen-lig_len), systemLen)))
#                 self.RMSD = np.array(R)

#// Custom RMSD with mdtraj
                R = []
                LT = []
                systemLen = crystalComp.n_atoms
                referenceXYZ = crystalComp.xyz[0]
        
                if self.ioutfm == 1:
                    for xyz in comp.xyz:
                        R_this, LT_this = ALIGN_A_RMSD_B(xyz*10, referenceXYZ*10,
                                                range(0, systemLen-lig_len), range((systemLen-lig_len), systemLen))
#                                                 range(systemLen-lig_len-50, systemLen-lig_len), range((systemLen-lig_len), systemLen))
                        R.append(R_this)
                        LT.append(LT_this)
                elif self.ioutfm == 0:
                    for xyz in comp:
#// Custom read
                        R_this, LT_this = ALIGN_A_RMSD_B(xyz, referenceXYZ*10,
                                                range(0, systemLen-lig_len), range((systemLen-lig_len), systemLen))
#                                                 range(systemLen-lig_len-50, systemLen-lig_len), range((systemLen-lig_len), systemLen))
                        R.append(R_this)
                        LT.append(LT_this)
#// mdtraj
#                         R.append(ALIGN_A_RMSD_B(xyz*10, referenceXYZ*10,
#                                                 range(systemLen-lig_len-50, systemLen-lig_len), range((systemLen-lig_len), systemLen)))
                self.RMSD = np.array(R)
                LT = np.array(LT)
                
                
#// Load with MDAnalysis but use mdtraj
#                 R = []
#                 systemLen = len(crystalComp.atoms)
#                 referenceXYZ = crystalComp.atoms.positions
# #                 print(referenceXYZ.dtype)
#                 for ts in comp.trajectory:
#                     R.append(alignment.rmsd_qcp(comp.atoms.positions, referenceXYZ))
#                 self.RMSD = np.sqrt((np.array(R)**2)*systemLen/lig_len)*10
                
#// mdtraj                
#                 self.RMSD = (np.sqrt((mdtraj.rmsd(comp, crystalComp, frame=0, atom_indices=range(crystalComp.n_atoms-400, crystalComp.n_atoms)))**2/lig_len*400)*10)

                if Profile:
                    t3 = time.time()
                self.output = ReadLog(self.outFile)
                if Profile:
                    t4 = time.time()
                self.hasRMSD = True
                if Profile:
                    print(f'    Profiling: Loading comp {(t1-t0)*1000:.3f} ms | Superpose {(t2-t1)*1000:.3f} ms | RMSD {(t3-t2)*1000:.3f} ms | Read output {(t4-t3)*1000:.3f} ms ')
                self.hasTrjFile = True
                self.hasOutFile = True
                 
                if not self.hasLigandTrajectory: # Also store ligand trajectories for contact analysis
#                     print('Now try logging lig traj')
                    lig = crystalComp.top.select("residue 56 and not symbol H")
                    ligH = crystalComp.top.select("residue 56")
                    pro = crystalComp.top.select("not residue 56 and not symbol H")
#                     if self.ioutfm == 1:
#                         self.ligandTrajectory = comp.xyz[:,lig]*10
#                         self.hasLigandTrajectory = True
#                     elif self.ioutfm == 0:
#                         self.ligandTrajectory = comp[:,lig]
#                         self.hasLigandTrajectory = True
#                     # else: Nothing changes
                    self.ligandTrajectory = LT[:,lig]
                    self.ligandTrajectoryH = LT[:, ligH]
                    proteinAtomSelection = np.array([ 45, 106, 107, 167, 168, 170, 175, 176, 177, 178, 179, 180, 181, 182, 232, 259, 381, 386, 387, 388])
#                     self.proteinCheck = LT[0,pro][proteinAtomSelection]
#                     self.totalTrajectory = LT
                    self.hasLigandTrajectory = True

        except:
            pass