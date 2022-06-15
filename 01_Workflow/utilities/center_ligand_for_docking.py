
import mdtraj, sys
import numpy as np

mol = mdtraj.load(sys.argv[1])

dock_center = np.array([float(sys.argv[2])/10.0, float(sys.argv[3])/10.0, float(sys.argv[4])/10.0])

mol.xyz[0] -= mol.xyz[0].mean(0)
mol.xyz[0] += dock_center

mol.save(sys.argv[5])
