
def write_gpf(fname, rname, X, Y, Z):
    # fname is file name with relative path, e.g. ../03_Gridmaps/1lyz_X_Y_Z/1lyz.gpf
    # rname is receptor name, e.g. 1lyz
    with open(fname,'w') as f:
        f.write(f'''
npts 65 65 65                   # num.grid points in xyz
gridfld {rname}.fld             # grid_data_file
spacing 0.375                      # spacing(A)
receptor_types A C OA N P S HD NA     # receptor atom types
ligand_types C N OA HD A NA F SA S Cl P I Br Fe   # ligand atom types
receptor {rname}.pdbqt               # macromolecule
gridcenter {X:.3f} {Y:.3f} {Z:.3f}         # xyz coordinates for the center of the box.ATOM    352  CA  PRO P  39, shifted by -8 6 8
smooth 0.5                           # store minimum energy w/in rad(A)
map {rname}.C.map                    # atom-specific affinity map
map {rname}.N.map                    # atom-specific affinity map
map {rname}.OA.map                   # atom-specific affinity map
map {rname}.HD.map                    # atom-specific affinity map
map {rname}.NA.map                    # atom-specific affinity map
map {rname}.A.map                   # atom-specific affinity map
map {rname}.F.map                   # atom-specific affinity map
map {rname}.Cl.map                   # atom-specific affinity map
map {rname}.SA.map                   # atom-specific affinity map
map {rname}.S.map                   # atom-specific affinity map
map {rname}.P.map                   # atom-specific affinity map
map {rname}.I.map                   # atom-specific affinity map
map {rname}.Br.map                   # atom-specific affinity map
map {rname}.Fe.map                   # atom-specific affinity map
elecmap {rname}.e.map                # electrostatic potential map
dsolvmap {rname}.d.map              # desolvation potential map
dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant
''')


def write_prepgrid(fh, config, dname, rname, rfname):
    fh.write(f'''
        cd ../{dname}
        {config["reduceH"]} -NOFLIP ../../02_Input/{rfname} > {rname}_h.pdb
        obabel -ipdb {rname}_h.pdb -opdbqt -O {rname}.pdbqt --partialcharge gasteiger -xr -xp -xc
        {config["autogrid"]} -p {rname}.gpf
        ''')

def write_dock(fh, config, jname, dname, rname, lname, lfname):
    fh.write(f'''
        cd ../{jname}
        obabel -ipdb ../../02_Input/{lfname} -opdbqt -O {lname}.pdbqt -p --partialcharge eem
        obabel -ipdbqt {lname}.pdbqt -opdb -O {lname}.pdb -d
        lig_atoms=`obabel -ipdb {lname}.pdb -opdb -h | grep "ATOM" | wc -l`
        echo "Ligand {lname} has $lig_atoms atoms including hydrogens"
        {config["autodock"]} -nrun 50 -autostop 1 -nev 3000000 -ffile ../../03_Gridmaps/{dname}/{rname}.fld -lfile {lname}.pdbqt
        python ../../01_Workflow/utilities/parseDock.py . {jname}
        mkdir -p docked_pdb
        #mkdir -p docked_pdb_h
        for i in `ls docked`;
        do
            obabel -ipdbqt docked/${{i}} -opdb -Odocked_pdb/${{i%.*}}.pdb -d
        done
        mkdir -p ../../05_Refinement/{jname}/Structure/docked_pdb_h
        cd ../../03_Gridmaps/{dname}
        python assembly.py $lig_atoms
        cd assembly
        sh overall_script.sh
        cd ..
        cd ../../04_Docking/{jname}
        ''')

def write_assembly(fname, rname, X, Y, Z):
    with open(fname, 'w') as f:
        f.write(f'''
import glob
import numpy as np
import os, sys
from shutil import copyfile
sys.path.insert(0,'../../01_Workflow/utilities')
from assembleCore import PDB, atom_line, parse_PDB 

try:
    lig_num_atoms = int(sys.argv[1])
except:
    lig_num_atoms = 60

proj = PDB('{rname}', np.array([[{X}, {Y}, {Z}]]), lig_num_atoms)

proj.determine_threshold(target_atoms=860)
proj.write_all_PDBs(dry_run = False, path=f'./assembly/', write_script=True)
''')


def write_antechamber(fh, jname, lname, dname):
    fh.write(f'''
        mkdir -p ../{jname}
        mkdir -p ../{jname}/Structure
        cd ../{jname}/Structure
        sh ../../../01_Workflow/utilities/pipeline_extra.sh ../../../04_Docking/{jname}/{lname.split('/')[-1]} > antechamber.log
        mv {lname.split('/')[-1].split('.')[0]} lig_param
        cd ..
        ''')

def write_prepareComplex(fh, jname, lname, dname):
    fh.write(f''' 
        # First convert all the docked_pdb files in 04_Docking
        cd ../{jname}/Structure
        cd ../../../04_Docking/{jname}/
        for i in `ls docked_pdb`;
        do
            obabel -ipdb docked_pdb/${{i}} -opdb -O ../../05_Refinement/{jname}/Structure/docked_pdb_h/${{i%.*}}.pdb -h
        done
        cd ../../05_Refinement/{jname}/Structure
        python ../../../01_Workflow/utilities/prepareComplex.py ../../../03_Gridmaps/{dname}/assembly/ lig_param/ docked_pdb_h/ {jname}
        cd ../
        ''')

def write_prepareComplex_openmm(fh, jname, lname, dname):
    fh.write(f''' 
        cd ../../05_Refinement/{jname}/Structure
        python ../../../01_Workflow/utilities/prepareComplex_openmm.py ../../../03_Gridmaps/{dname}/assembly/ lig_param/ docked_pdb_h/ {jname}
        cd ../
        ''')

def write_cpptraj(fname, jname, parallelcpptraj):
    with open(fname, 'w') as f:
        f.write(f'''
import glob, os, sys
inpcrd = [x.split('/')[-1].split('.')[0] for x in glob.glob('../../Structure/inpcrd/*')]
inpcrd.sort()
# Get num of res
with open(f'../../reference_structure/complex.prmtop') as f:
    cont = f.readlines()
n_res = 0
count_res_flag = 0
for line in cont:
    if 'RESIDUE_LABEL' in line:
        count_res_flag = 1
        continue
    if 'RESIDUE_POINTER' in line:
        break
    if 'FORMAT' not in line and count_res_flag:
        n_res += len(line.split())

# make dirs
for pose in inpcrd:
    os.makedirs(f'../in/{{pose}}', exist_ok=True)
    os.makedirs(f'../out/{{pose}}', exist_ok=True)

process_these_inpcrd = inpcrd
num_script_files = min({parallelcpptraj}, len(process_these_inpcrd))
for ii in range(0, num_script_files):
    try:
        os.remove(f'process_cpptraj_{{ii}}.sh')
    except:
        pass

this_ligand = ''
for idx, pose in enumerate(process_these_inpcrd):
    with open(f'process_cpptraj_{{idx % num_script_files}}.sh','a') as g:
        if pose.split('_')[0] != this_ligand:
            print('')
            this_ligand = pose.split('_')[0]
            print(pose.split('_')[0], end=' ')
        print(pose.split('_')[-1], end = ' ')
        lig = pose.split('_')[0]
        traj_list = [x.split('/')[-1].split('.')[0] for x in glob.glob(f'../../Simulation/{{pose}}/*.nc')]
        for traj in traj_list:
            g.write(f'~/Tools/amber_rhel8_2/bin/cpptraj -i ../in/{{pose}}/{{traj}}.in\\n')
            with open(f'../in/{{pose}}/{{traj}}.in','w') as f:
                f.write(f'parm ../../Structure/prmtop/{{this_ligand}}.prmtop\\n')
                f.write(f'trajin ../../Simulation/{{pose}}/{{traj}}.nc\\n')
                f.write(f'energy :1-{{n_res-1}} out ../out/{{pose}}/{{traj}}_apo.dat\\n')
                f.write(f'energy :{{n_res}} out ../out/{{pose}}/{{traj}}_lig.dat\\n')
                f.write(f'energy out ../out/{{pose}}/{{traj}}_holo.dat\\n')

with open('process_all_cpptraj.sh','w') as f:
    #f.write(f'for i in {{{{0..{{num_script_files-1}}}}}}; do jsrun -n 1 sh process_cpptraj_${{{{i}}}}.sh > output_${{{{i}}}}.log & done\\n')
    f.write(f'for i in {{{{0..{{num_script_files-1}}}}}}; do sh process_cpptraj_${{{{i}}}}.sh > output_${{{{i}}}}.log & done\\n')
    f.write('wait\\n')
''')

def write_em(fh, jname):
    fh.write(f'''
        cd ../{jname}/Structure
        for i in `seq 1 20`
        do
            echo "Minimizing energy for {jname}_${{i}} ..."
            python ../../../01_Workflow/utilities/em.py {jname}_${{i}}
        done
        cd ../
        ''')

def write_em_list(fh, jname):
    for i in range(1, 21):
        fh.write(f'{jname}_{i}\n')
        
def write_MDR_analysis(fname, jname):
    with open(fname, 'w') as f:
        f.write(f'python ../../01_Workflow/utilities/MDR_cluster_scoring.py ../../05_Refinement/{jname}')
