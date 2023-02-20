import sys
from vina import Vina
import mdtraj
import numpy as np

# Usage: vina_dock.py <receptor_file> <ligand_file> <dockX> <dockY> <dockZ> <num_output_poses> <energy_range>
# energy_range can be float, 'inf', 'infinity', 'np.inf'; otherwise, defaults to None (let Vina set it)

mol = mdtraj.load(sys.argv[2].replace('.pdbqt','_d_c.pdb'))
rg = mdtraj.compute_rg(mol)[0] * 10
rgbox = rg * 2.9 # Based on a paper on optimal box size

print(f'The rg of ligand {sys.argv[2]} is {rg} A and therefore the box size is first determined to be {rgbox} A cubed.')

center = np.array([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])])

# Expand box a bit
max_dist = np.max(np.abs(mol.xyz[0] * 10 - center)) # This is per axis

if max_dist > rgbox / 2.0:
    rgbox = max_dist * 2.1
    print(f'The box side length is increased to {rgbox} A')

if rgbox > 40:
    print('Warning: rgbox is unexpectedly large, check your ligand')

#if rgbox < 12:
#    rgbox = 12
#    print(f'The box is too small - side length is increased to {rgbox} A')

v = Vina(sf_name='vina')

v.set_receptor(sys.argv[1])
v.compute_vina_maps(center=[float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])], box_size=[24, 24, 24])
#v.compute_vina_maps(center=center, box_size=[rgbox, rgbox, rgbox])
    
v.set_ligand_from_file(sys.argv[2])
    
    
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=2*int(sys.argv[6]))
if sys.argv[7] in ['inf', 'np.inf', 'infinity']:
    v.write_poses('docking.dlg', n_poses=int(sys.argv[6]), overwrite=True, energy_range=np.inf)
else:
    try:
        float(sys.argv[7])
        v.write_poses('docking.dlg', n_poses=int(sys.argv[6]), overwrite=True, energy_range=float(sys.argv[7]))
    except:
        print("sys.argv[7] is not a float, use no energy range")
        v.write_poses('docking.dlg', n_poses=int(sys.argv[6]), overwrite=True)
