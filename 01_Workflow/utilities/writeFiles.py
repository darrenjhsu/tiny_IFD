

def write_preppdbqt(fh, config, dname, rname, rfname, dockX, dockY, dockZ, rigid=True):
    fh.write(f'''
        cd ../{dname}
        obabel -ipdb ../../02_Input/{rfname} -opdb -O {rname}_d.pdb -d
        tleap -f fix_protein.in > fix_protein.out
        python ../../01_Workflow/utilities/restore_residues.py {rname}
        obabel -ipdb {rname}_restored_h.pdb -opdb -O {rname}_restored_e.pdb # Just to get elements
        obabel -ipdb {rname}_restored_e.pdb -opdbqt -O {rname}.pdbqt --partialcharge gasteiger -xr -xp -xc''')
    if not rigid:
        fh.write(f'''
        python ../../01_Workflow/utilities/create_flex_res.py {rname}.pdbqt {dockX} {dockY} {dockZ}
        ''')


def write_vina_dock(fh, config, jname, dname, rname, lname, lfname, dockX, dockY, dockZ):
    fh.write(f'''
        cd ../{jname}
        obabel -ipdb ../../02_Input/{lfname} -opdb -O {lname}_d.pdb -d # Remove hydrogens that may come with the ligand XRD
        python ../../01_Workflow/utilities/center_ligand_for_docking.py {lname}_d.pdb {dockX} {dockY} {dockZ} {lname}_d_c.pdb # Center the ligand file to be used in docking
        obabel -ipdb {lname}_d.pdb -opdbqt -O {lname}_d.pdbqt -p --partialcharge eem # Get an extra pdbqt (atom indices may change) for calculating RMSD in the subsequent analyses
        obabel -ipdb {lname}_d_c.pdb -opdbqt -O {lname}.pdbqt -p --partialcharge eem # Get the pdbqt file for docking
        obabel -ipdbqt {lname}_d.pdbqt -opdb -O {lname}.pdb -d # Get the extra pdb for calculating RMSD in the subseuqent analyses
        lig_atoms=`obabel -ipdb {lname}.pdb -opdb -h | grep "ATOM" | wc -l`
        echo "Ligand {lname} has $lig_atoms atoms including hydrogens"
        python ../../01_Workflow/utilities/vina_dock.py ../../03_PrepProtein/{dname}/{rname}.pdbqt {lname}.pdbqt {dockX} {dockY} {dockZ}
        python ../../01_Workflow/utilities/parseVinaDock.py . {jname}
        mkdir -p docked_pdb
        #mkdir -p docked_pdb_h
        for i in `ls docked`;
        do
            obabel -ipdbqt docked/${{i}} -opdb -Odocked_pdb/${{i%.*}}.pdb -d
        done
        mkdir -p ../../05_Refinement/{jname}/Structure/docked_pdb_h
        cd ../../03_PrepProtein/{dname}
        mkdir -p assembly
        python assembly.py $lig_atoms > assembly.log
        cd assembly
        sh overall_script.sh
        cd ..
        cd ../../04_Docking/{jname}
        ''')

def write_vina_flex_dock(fh, config, jname, dname, rname, lname, lfname, dockX, dockY, dockZ):
    fh.write(f'''
        cd ../{jname}
        obabel -ipdb ../../02_Input/{lfname} -opdb -O {lname}_d.pdb -d
        obabel -ipdb {lname}_d.pdb -opdbqt -O {lname}.pdbqt -p --partialcharge eem
        obabel -ipdbqt {lname}.pdbqt -opdb -O {lname}.pdb -d
        lig_atoms=`obabel -ipdb {lname}.pdb -opdb -h | grep "ATOM" | wc -l`
        echo "Ligand {lname} has $lig_atoms atoms including hydrogens"
        python ../../01_Workflow/utilities/vina_flex_dock.py ../../03_PrepProtein/{dname}/{rname}.pdbqt {lname}.pdbqt {dockX} {dockY} {dockZ}
        python ../../01_Workflow/utilities/parseVinaFlexDock.py . {jname} ../../03_PrepProtein/{dname}
        mkdir -p docked_pdb
        #mkdir -p docked_pdb_h
        for i in `ls docked`;
        do
            obabel -ipdbqt docked/${{i}} -opdb -Odocked_pdb/${{i%.*}}.pdb -d
        done
        mkdir -p ../../05_Refinement/{jname}/Structure/docked_pdb_h
        cd ../../03_PrepProtein/{dname}
        mkdir -p assembly
        for i in `seq 1 20`;
        do
            mkdir -p flex_assembly_${{i}}
        done
        python flex_assembly.py $lig_atoms > flex_assembly.log
        cd assembly
        sh overall_script.sh
        cd ..
        for i in `seq 1 20`;
        do
            cd flex_assembly_${{i}}
            sh overall_script.sh
            cd ..
        done

        cd ../../04_Docking/{jname}
        ''')


def write_fix_protein(fname, rname):
    with open(fname, 'w') as f:
        f.write(f'''
source leaprc.protein.ff14SB #Source leaprc file for ff14SB protein force field
mol = loadpdb {rname}_d.pdb
savepdb mol {rname}_full_h.pdb
''')


def write_assembly(fname, rname, X, Y, Z, aMode='atoms', aThres=920):
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

proj = PDB('{rname}_restored_h.pdb', np.array([[{X}, {Y}, {Z}]]), lig_num_atoms)
proj.determine_threshold(target_{aMode}={aThres})
proj.write_all_PDBs(dry_run = False, path=f'./assembly/', write_script=True)
''')

def write_flex_assembly(fname, rname, jname, X, Y, Z):
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

proj = PDB('{rname}_restored_h.pdb', np.array([[{X}, {Y}, {Z}]]), lig_num_atoms)

proj.determine_threshold(target_atoms=920)
proj.write_all_PDBs(dry_run = False, path=f'./assembly/', write_script=True)
for ii in range(1, 21):
    try:
        proj.write_all_PDBs(dry_run = False, path=f'./flex_assembly_{{ii}}/', write_script=True, overwrite_coor_file=f'flex/{jname}_{{ii}}.pdbqt')
    except:
        pass
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
        python ../../../01_Workflow/utilities/prepareComplex.py ../../../03_PrepProtein/{dname}/assembly/ lig_param/ docked_pdb_h/ {jname}
        cd ../reference_structure
        obabel -ipdb ligand.pdb -opdb -O ligand_d.pdb -d
        obabel -ipdb ligand_d.pdb -omol -O ligand.mol
        cd ../
        ''')

def write_prepareFlexComplex(fh, jname, lname, dname):
    fh.write(f''' 
        # First convert all the docked_pdb files in 04_Docking
        cd ../{jname}/Structure
        cd ../../../04_Docking/{jname}/
        for i in `ls docked_pdb`;
        do
            obabel -ipdb docked_pdb/${{i}} -opdb -O ../../05_Refinement/{jname}/Structure/docked_pdb_h/${{i%.*}}.pdb -h
        done
        cd ../../05_Refinement/{jname}/Structure
        python ../../../01_Workflow/utilities/prepareFlexComplex.py ../../../03_PrepProtein/{dname}/ lig_param/ docked_pdb_h/ {jname}
        cd ../
        ''')

def write_prepareComplex_openmm(fh, jname, lname, dname):
    fh.write(f''' 
        cd ../../05_Refinement/{jname}/Structure
        python ../../../01_Workflow/utilities/prepareComplex_openmm.py ../../../03_PrepProtein/{dname}/assembly/ lig_param/ docked_pdb_h/ inpcrd/ {jname}
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
            g.write(f'~/Tools/amber_andes/bin/cpptraj -i ../in/{{pose}}/{{traj}}.in\\n')
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
        for i in `seq 0 20`
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

def write_create_xgboost_model(fname, jname):
    with open(fname, 'w') as f:
        f.write(f'python ../../01_Workflow/utilities/create_model.py {jname}')

def write_create_many_xgboost_models(fname, jname):
    with open(fname, 'w') as f:
        f.write(f'python ../../01_Workflow/utilities/create_many_models.py {jname} ../../01_Workflow/utilities/xgboost_config_list.csv')

def write_create_xgboost_regression_model(fname, jname):
    with open(fname, 'w') as f:
        f.write(f'python ../../01_Workflow/utilities/create_regression_model.py {jname}')
