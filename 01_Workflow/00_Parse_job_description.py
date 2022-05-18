
import numpy as np
import pandas as pd
import os
from os.path import exists
import argparse
import shutil

from utilities.writeFiles import *

# Parse arguments

parser = argparse.ArgumentParser(description="Generate related work shell scripts")
parser.add_argument('--parallelGrid', nargs='?', type=int, default=1, help="Split grid prep to N chunks for parallel execution")
parser.add_argument('--parallelDock', nargs='?', type=int, default=1, help="Split docking to N chunks for parallel execution")
parser.add_argument('--parallelAntechamber', nargs='?', type=int, default=1, help="Split ligand parametrization to N chunks for parallel execution")
parser.add_argument('--parallelPrepareComplex', nargs='?', type=int, default=1, help="Split complex preparation to N chunks for parallel execution")
parser.add_argument('--parallelEM', nargs='?', type=int, default=1, help="Split energy minimization to N chunks for parallel execution")
parser.add_argument('--parallelcpptraj', nargs='?', type=int, default=1, help="Parallelize CPU cpptraj calculation")
parser.add_argument('--parallelMDR', nargs='?', type=int, default=1, help="Split post processing of MD simulations to N chunks")
parser.add_argument('--parallelXGBoost', nargs='?', type=int, default=1, help="Split XGBoost training of features to N chunks")
parser.add_argument('--dockingMode', type=str, default='rigid', help="Whether to use flexible docking (changes some residue orientations) (rigid|flex)")
parser.add_argument('--config', nargs='?', type=str, default=None, help="File logging locations of executables")

args = parser.parse_args()
print(args)

config = {}
config['parallelGrid'] = args.parallelGrid
config['parallelDock'] = args.parallelDock
config['parallelAntechamber'] = args.parallelAntechamber
config['parallelPrepareComplex'] = args.parallelPrepareComplex
config['parallelEM'] = args.parallelEM
config['parallelcpptraj'] = args.parallelcpptraj
config['parallelMDR'] = args.parallelMDR
config['parallelXGBoost'] = args.parallelXGBoost
config['dockingMode'] = args.dockingMode


os.makedirs('../03_Gridmaps', exist_ok=True)
os.makedirs('../04_Docking', exist_ok=True)
os.makedirs('../05_Refinement', exist_ok=True)
os.makedirs('../06_Analysis', exist_ok=True)


if args.config is not None:
    with open(args.config,'r') as f:
        cont = f.readlines()
    for line in cont:
        try:
            print(line)
            if line.split('=')[0] in ['parallelGrid', 'parallelDock', 'parallelAntechamber', 'parallelPrepareComplex', 'parallelEM', 'parallelcpptraj', 'parallelMDR', 'parallelXGBoost']:
                try:
                    config[line.split('=')[0]] = int(line.split('=')[1].strip())
                except:
                    pass
            else:
                config[line.split('=')[0]] = line.split('=')[1].strip()
        except:
            pass

# Short hand for flexibleDocking:
if config['dockingMode'] == 'rigid':
    rigid = True
elif config['dockingMode'] == 'flex':
    rigid = False
else:
    raise ValueError('dockingMode must be either rigid or flex')

print(config)

# Load job description

Jobs = pd.read_csv('../02_Input/job_description.csv')#, column=['Receptor_file_name','Ligand_file_name','dockX','dockY','dockZ'])
Jobs['Receptor_name'] = Jobs['Receptor_file_name'].str.split('/').str[-1].str.split('.').str[0] # What kind of awful code is this??
Jobs['Ligand_name'] = Jobs['Ligand_file_name'].str.split('/').str[-1].str.split('.').str[0] 

print(Jobs)

# Check all files are in folders
Rec = Jobs['Receptor_file_name'].unique()
Lig = Jobs['Ligand_file_name'].unique()
#Ref = Jobs['Ref_receptor_file_name'].unique()

print(Rec, Lig)

for R in Rec:
    if not exists(f'../02_Input/{R}'):
        print(f'Receptor {R} does not exist!')
        raise FileNotFoundError
for L in Lig:
    if not exists(f'../02_Input/{L}'):
        print(f'Ligand {L} does not exist!')
        raise FileNotFoundError
#for R in Ref:
#    if not exists(f'../02_Input/{R}'):
#        print(f'Reference receptor {R} does not exist!')
#        raise FileNotFoundError

print("All files in job description exist")


# Get unique receptor / center combos so we only need to process each receptor once
# Processing includes writing gpf files, prepgrid script (parallelizable), and assembly script

Receptors = Jobs.drop(columns=['Ligand_file_name', 'Job_name', 'Ligand_name'])
Receptors = Receptors.drop_duplicates()
print(Receptors)


try:
    os.mkdir('../03_Gridmaps/script')
except:
    pass
fh = []  # File handle array

g = open('../03_Gridmaps/script/andes_prepDock.sh', 'w')
g.write('#!/bin/bash\n')
g.write(f'#SBATCH -N 1\n')
g.write('#SBATCH -t 1:00:00\n')
g.write('#SBATCH -A STF006\n\n')
g.write('echo "`date`: Job starts"\n')
g.write('source ~/.andesrc\n\n\n')

for ii in range(config['parallelGrid']):
    fh.append(open(f'../03_Gridmaps/script/prepDock{ii}.sh','w'))
    g.write(f'srun -n1 -N1 -c1 -s sh prepDock{ii}.sh & \n')

g.write('wait\n\n')
g.write('echo "`date`: All done!\n\n"')
g.close()

for idx, row in Receptors.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    os.makedirs(f'../03_Gridmaps/{dname}', exist_ok=True)
    write_preppdbqt(fh[idx % config['parallelGrid']], config, dname, row["Receptor_name"], row["Receptor_file_name"], row["dockX"], row["dockY"], row["dockZ"], rigid)
    write_fix_protein(f'../03_Gridmaps/{dname}/fix_protein.in', row["Receptor_name"])
    write_assembly(f'../03_Gridmaps/{dname}/assembly.py',row["Receptor_name"], row["dockX"], row["dockY"], row["dockZ"])
for ii in range(config['parallelGrid']):
    fh[ii].close()


# Write docking script. To parallelize this part it requires an active allocation.

try:
    os.mkdir('../04_Docking/script')
except:
    pass

fh = []  # File handle array

g = open('../04_Docking/script/andes_dock.sh', 'w')
g.write('#!/bin/bash\n')
g.write(f'#SBATCH -N {config["parallelDock"]}\n')
g.write('#SBATCH -t 1:00:00\n')
g.write('#SBATCH -A STF006\n\n')
g.write('echo "`date`: Job starts"\n')
g.write('source ~/.andesrc\n\n\n')

for ii in range(config['parallelDock']):
    fh.append(open(f'../04_Docking/script/dock{ii}.sh','w'))
    g.write(f'srun -n1 -N1 -c32 sh dock{ii}.sh > dock{ii}.log & \n')

g.write('wait\n\n')
g.write('echo "`date`: All done!\n\n"')
g.close()


with open(f'../04_Docking/script/check_docking.sh','w') as f:
    for idx, row in Jobs.iterrows():
        dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
        lname = row["Ligand_name"]
        jname = row["Job_name"]
        os.makedirs(f'../04_Docking/{jname}', exist_ok=True)
        if rigid:
            write_vina_dock(fh[idx % config['parallelDock']], config, jname, dname, row["Receptor_name"], row["Ligand_name"], row["Ligand_file_name"], row['dockX'], row['dockY'], row['dockZ'])
        else:
            write_vina_flex_dock(fh[idx % config['parallelDock']], config, jname, dname, row["Receptor_name"], row["Ligand_name"], row["Ligand_file_name"], row['dockX'], row['dockY'], row['dockZ'])
            write_flex_assembly(f'../03_Gridmaps/{dname}/flex_assembly.py',row["Receptor_name"], jname, row["dockX"], row["dockY"], row["dockZ"])
        f.write(f'python ../../01_Workflow/utilities/check_dock.py {jname} {lname}.pdb\n')

for ii in range(config['parallelDock']):
    fh[ii].close()



# The heavy-weight jobs need multiple andes nodes


# Write antechamber and complex preparation script

try:
    os.makedirs('../05_Refinement/script')
except:
    pass

# Write antechamber script

fh = []
g = open('../05_Refinement/script/andes_antechamber.sh', 'w')
g.write('#!/bin/bash\n')
g.write(f'#SBATCH -N {config["parallelAntechamber"]}\n')
g.write('#SBATCH -t 2:00:00\n')
g.write('#SBATCH -A STF006\n\n')
g.write('echo "`date`: Job starts"\n')
g.write('source ~/.andesrc\n\n\n')

for ii in range(config['parallelAntechamber']):
    fh.append(open(f'../05_Refinement/script/antechamber{ii}.sh','w'))
    g.write(f'srun -n1 -N1 -c32 sh antechamber{ii}.sh & \n')

g.write('wait\n\n')
g.write('echo "`date`: All done!\n\n"')
g.close()

for idx, row in Jobs.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    lname = row["Ligand_file_name"]
    jname = row["Job_name"]
    try:
        os.makedirs(f'../05_Refinement/{jname}/Structure/')
    except:
        pass
    write_antechamber(fh[idx % config['parallelAntechamber']], jname, lname, dname)

for ii in range(config['parallelDock']):
    fh[ii].close()

# Write complex assembly script

fh = []
gh = []
if rigid:
    g = open('../05_Refinement/script/andes_prepareComplex.sh', 'w')
    g.write('#!/bin/bash\n')
    g.write(f'#SBATCH -N 1\n')
    g.write('#SBATCH -t 1:00:00\n')
    g.write('#SBATCH -A STF006\n\n')
    g.write('echo "`date`: Job starts"\n')
    g.write('source ~/.andesrc\n\n\n')
else:
    ih = open('../05_Refinement/script/andes_prepareFlexComplex.sh', 'w')
    ih.write('#!/bin/bash\n')
    ih.write(f'#SBATCH -N 1\n')
    ih.write('#SBATCH -t 1:00:00\n')
    ih.write('#SBATCH -A STF006\n\n')
    ih.write('echo "`date`: Job starts"\n')
    ih.write('source ~/.andesrc\n\n\n')

for ii in range(config['parallelPrepareComplex']):
    if rigid:
        fh.append(open(f'../05_Refinement/script/prepareComplex{ii}.sh','w'))
        g.write(f'srun -n1 -N1 -c1 -s sh prepareComplex{ii}.sh & \n')
    else:
        gh.append(open(f'../05_Refinement/script/prepareFlexComplex{ii}.sh','w'))
        ih.write(f'srun -n1 -N1 -c1 -s sh prepareFlexComplex{ii}.sh & \n')

if rigid:
    g.write('wait\n\n\n')
else:
    ih.write('wait\n\n\n')

for idx, row in Jobs.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    lname = row["Ligand_file_name"]
    jname = row["Job_name"]
    os.makedirs(f'../05_Refinement/{jname}/Structure/', exist_ok=True)
    os.makedirs(f'../05_Refinement/{jname}/Simulation/', exist_ok=True)
    if rigid:
        write_prepareComplex(fh[idx % config['parallelPrepareComplex']], jname, lname, dname)
    else:
        write_prepareFlexComplex(gh[idx % config['parallelPrepareComplex']], jname, lname, dname)

for ii in range(config['parallelPrepareComplex']):
    if rigid:
        fh[ii].close()
    else:
        gh[ii].close()

# Write complex assembly prmtop for openmm
fh = []

for ii in range(config['parallelPrepareComplex']):
    fh.append(open(f'../05_Refinement/script/prepareComplex_openmm{ii}.sh','w'))
    if rigid:
        g.write(f'srun -n1 -N1 -c1 -s sh prepareComplex_openmm{ii}.sh & \n')
    else: 
        # This is actually wrong, we need a prepareFlexComplex_openmm thing
        ih.write(f'srun -n1 -N1 -c1 -s sh prepareComplex_openmm{ii}.sh & \n')

if rigid:
    g.write('wait\n\n')
    g.write('echo "`date`: All done!\n\n"')
    g.close()
else:
    ih.write('wait\n\n')
    ih.write('echo "`date`: All done!\n\n"')
    ih.close()

for idx, row in Jobs.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    lname = row["Ligand_file_name"]
    jname = row["Job_name"]
    write_prepareComplex_openmm(fh[idx % config['parallelPrepareComplex']], jname, lname, dname)

for ii in range(config['parallelPrepareComplex']):
    fh[ii].close()

fh = []
g = open('../05_Refinement/script/andes_EM.sh', 'w')
g.write('#!/bin/bash\n')
g.write(f'#SBATCH -N {config["parallelEM"]}\n')
g.write('#SBATCH -t 1:00:00\n')
g.write('#SBATCH -A STF006\n\n')
g.write('echo "`date`: Job starts"\n')
g.write('source ~/.andesrc\n\n\n')
for ii in range(config['parallelEM']):
    fh.append(open(f'../05_Refinement/script/energyMinimization{ii}.sh','w'))
    g.write(f'srun -n1 -N1 -c32 sh energyMinimization{ii}.sh & \n')

g.write('wait\n\n')
g.write('echo "`date`: All done!\n\n"')
g.close()

for idx, row in Jobs.iterrows():
    jname = row["Job_name"]
    write_em(fh[idx % config['parallelEM']], jname)

for ii in range(config['parallelEM']):
    fh[ii].close()

# Copy some files 
shutil.copyfile('utilities/planEM.sh', '../05_Refinement/script/planEM.sh')
shutil.copyfile('utilities/planMD_openmm.sh', '../05_Refinement/script/planMD_openmm.sh')
shutil.copyfile('utilities/check_atom_num.sh', '../05_Refinement/script/check_atom_num.sh')
shutil.copyfile('utilities/check_sims.py', '../05_Refinement/script/check_sims.py')


# Write cpptraj script

with open(f'../05_Refinement/script/gen_cpptraj_script.sh','w') as f,\
     open(f'../05_Refinement/script/run_cpptraj_script.sh','w') as g,\
     open(f'../05_Refinement/script/run_cpptraj_script_andes.sh','w') as h:
    h.write('#!/bin/bash\n')
    h.write(f'#SBATCH -N {config["parallelcpptraj"]}\n')
    h.write('#SBATCH -t 1:00:00\n')
    h.write('#SBATCH -A STF006\n\n')

    h.write('echo "`date`: Job starts"\n')
    h.write('source ~/.andesrc\n\n\n')
    batch_counter = 0
    for idx, row in Jobs.iterrows():
        jname = row["Job_name"]
        os.makedirs(f'../05_Refinement/{jname}/cpptraj', exist_ok=True)
        os.makedirs(f'../05_Refinement/{jname}/cpptraj/in', exist_ok=True)
        os.makedirs(f'../05_Refinement/{jname}/cpptraj/out', exist_ok=True)
        os.makedirs(f'../05_Refinement/{jname}/cpptraj/script', exist_ok=True)
        write_cpptraj(f'../05_Refinement/{jname}/cpptraj/script/generate_cpptraj_inputs.py', jname, config['parallelcpptraj'])
        f.write(f'cd ../{jname}/cpptraj/script\n')
        f.write(f'python generate_cpptraj_inputs.py\n')
        f.write(f'cd ../../\n')
        g.write(f'echo processing {jname}\n')
        g.write(f'cd ../{jname}/cpptraj/script\n')
        g.write(f'sh process_all_cpptraj.sh\n')
        g.write(f'wait\n')
        g.write(f'cd ../../\n')
        if idx % config["parallelcpptraj"] == 0:
            if idx > 0:
                h.write('wait\n\n')
            batch_counter += 1
            h.write(f'echo "Batch {batch_counter}"\n')
        h.write(f'cd ../{jname}/cpptraj/script\n')
        h.write(f'srun -n1 -N1 -c32 sh process_all_cpptraj.sh &\n')
        h.write(f'cd ../../\n')
    if idx % config["parallelcpptraj"] != 0:
        h.write('wait\n\n')
    h.write('echo "`date`: All jobs done!"\n\n')
# Write analysis script

try:
    os.makedirs('../06_Analysis/script')
except:
    pass


with open('../06_Analysis/script/MDR_analysis_all.sh', 'w') as f,\
     open('../06_Analysis/script/MDR_analysis_andes.sh', 'w') as g: # A script for andes batch job
    g.write('#!/bin/bash\n')
    g.write(f'#SBATCH -N {config["parallelMDR"]}\n')
    g.write('#SBATCH -t 1:00:00\n')
    g.write('#SBATCH -A STF006\n\n')

    g.write('echo "`date`: Job starts"\n')
    g.write('source ~/.andesrc\n\n\n')
    batch_counter = 0
    for idx, row in Jobs.iterrows():
        jname = row["Job_name"]
        try:
            os.makedirs(f'../06_Analysis/{jname}')
        except:
            pass
        write_MDR_analysis(f'../06_Analysis/{jname}/MDR_analysis.sh', jname)
        f.write(f'cd ../{jname}\n')
        f.write(f'echo {jname}\n')
        f.write(f'sh MDR_analysis.sh > MDR.log\n')
        if idx % config["parallelMDR"] == 0:
            if idx > 0:
                g.write('wait\n\n')
            batch_counter += 1
            g.write('date\n')
            g.write(f'echo "Batch {batch_counter}"\n')
        g.write(f'cd ../{jname}\n')
        g.write('srun -n1 -N1 -c32 sh MDR_analysis.sh > MDR.log &\n')
    if idx % config["parallelMDR"] != 0:
        g.write('wait\n\n')
    g.write('echo "`date`: All jobs done!"\n\n')


# Individual xgboost classification models
with open('../06_Analysis/script/create_xgboost_model_all.sh', 'w') as f,\
     open('../06_Analysis/script/create_xgboost_model_andes.sh', 'w') as g: # A script for andes batch job
    g.write('#!/bin/bash\n')
    g.write(f'#SBATCH -N {config["parallelXGBoost"]}\n')
    g.write('#SBATCH -t 1:00:00\n')
    g.write('#SBATCH -A STF006\n\n')

    g.write('echo "`date`: Job starts"\n')
    g.write('source ~/.andesrc\n\n\n')
    batch_counter = 0
    for idx, row in Jobs.iterrows():
        jname = row["Job_name"]
        try:
            os.makedirs(f'../06_Analysis/{jname}')
        except:
            pass
        write_create_xgboost_model(f'../06_Analysis/{jname}/create_model.sh', jname)
        f.write(f'cd ../{jname}\n')
        f.write(f'echo {jname}\n')
        f.write(f'sh create_model.sh > create_model.log\n')
        if idx % config["parallelXGBoost"] == 0:
            if idx > 0:
                g.write('wait\n\n')
            batch_counter += 1
            g.write('date\n')
            g.write(f'echo "Batch {batch_counter}"\n')
        g.write(f'cd ../{jname}\n')
        g.write('srun -n1 -N1 -c32 sh create_model.sh > create_model.log &\n')
    if idx % config["parallelXGBoost"] != 0:
        g.write('wait\n\n')
    g.write('echo "`date`: All jobs done!"\n\n')

## Individual xgboost regression models
#with open('../06_Analysis/script/create_xgboost_regression_model_all.sh', 'w') as f,\
#     open('../06_Analysis/script/create_xgboost_regression_model_andes.sh', 'w') as g: # A script for andes batch job
#    g.write('#!/bin/bash\n')
#    g.write(f'#SBATCH -N {N_andes}\n')
#    g.write('#SBATCH -t 1:00:00\n')
#    g.write('#SBATCH -A STF006\n\n')
#
#    g.write('echo "`date`: Job starts"\n')
#    g.write('source ~/.andesrc\n\n\n')
#    batch_counter = 0
#    for idx, row in Jobs.iterrows():
#        jname = row["Job_name"]
#        try:
#            os.makedirs(f'../06_Analysis/{jname}')
#        except:
#            pass
#        write_create_xgboost_regression_model(f'../06_Analysis/{jname}/create_regression_model.sh', jname)
#        f.write(f'cd ../{jname}\n')
#        f.write(f'echo {jname}\n')
#        f.write(f'sh create_regression_model.sh > create_regression_model.log\n')
#        if idx % N_andes == 0:
#            if idx > 0:
#                g.write('wait\n\n')
#            batch_counter += 1
#            g.write('date\n')
#            g.write(f'echo "Batch {batch_counter}"\n')
#        g.write(f'cd ../{jname}\n')
#        g.write('srun -n1 -N1 -c32 sh create_regression_model.sh > create_regression_model.log &\n')
#    if idx % N_andes != 0:
#        g.write('wait\n\n')
#    g.write('echo "`date`: All jobs done!"\n\n')

