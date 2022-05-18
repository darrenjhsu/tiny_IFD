import sys
from vina import Vina

v = Vina(sf_name='vina')

v.set_receptor(sys.argv[1])
v.compute_vina_maps(center=[float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])], box_size=[24, 24, 24])
    
v.set_ligand_from_file(sys.argv[2])
    
    
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=40)
v.write_poses('docking.dlg', n_poses=20, overwrite=True)
