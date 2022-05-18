
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import os

# argv: 1. a file including docked pose names (e.g. 1cqp-803-to-3e2m_1)


with open(sys.argv[1],'r') as f:
    cases = f.readlines()

print(f'In this file there are {len(cases)} poses to be minimized ...')

for case in cases:
    print(f'Minimizing {case.strip()} ...')
    pose = case.strip()
    ligand_name = pose.split('_')[0]
    inpcrd_name = pose + '.inpcrd'
    prmtop_name = ligand_name + '.prmtop'
    pdb_name =    pose + '.pdb'

    try:
        os.mkdir(f'../{ligand_name}/Structure/inpcrd_em')
    except:
        pass
    try:
        os.mkdir(f'../{ligand_name}/Structure/inpcrd_pdb')
    except:
        pass
    
     
    
    # Load sim system
    prmtop = AmberPrmtopFile(f'../{ligand_name}/Structure/prmtop/' + prmtop_name)
    inpcrd = AmberInpcrdFile(f'../{ligand_name}/Structure/inpcrd/' + inpcrd_name)
    system = prmtop.createSystem(nonbondedMethod=NoCutoff)
    integrator = LangevinIntegrator(0*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    simulation.minimizeEnergy()
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'../{ligand_name}/Structure/inpcrd_pdb/' + pdb_name, 'w'))
    
    #simulation.reporters.append(PDBReporter('inpcrd_em/' + pdb_name, 1))
    #simulation.reporters.append(StateDataReporter(stdout, 1, step=True,
    #        potentialEnergy=True, temperature=True))
    #simulation.step(1)
    
    
    with open(f'../{ligand_name}/Structure/inpcrd_pdb/' + pdb_name,'r') as p:
        cont = p.readlines()
        
    xyz = []
    for line in cont:
        if ('ATOM  ' in line) or ('HETATM' in line):
            # print(line)
            xyz.append([float(x) for x in line[30:54].split()])
    xyz = np.array(xyz)
    
    with open(f'../{ligand_name}/Structure/inpcrd_em/' + inpcrd_name,'w') as f:
        f.write('\n')
        f.write(f'{len(xyz):5d}{0.0:15.7e}\n')
        for idx, coor in enumerate(xyz):
            f.write(f'{coor[0]:12.7f}{coor[1]:12.7f}{coor[2]:12.7f}')
            if idx % 2 == 1:
                f.write('\n')
        f.write('\n')
