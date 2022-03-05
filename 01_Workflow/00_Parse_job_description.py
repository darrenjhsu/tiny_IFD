
import numpy as np
import pandas as pd
import os
from os.path import exists
import argparse

from utilities.writeFiles import *

# Parse arguments

parser = argparse.ArgumentParser(description="Generate related work shell scripts")
parser.add_argument('--parallelGrid', nargs='?', type=int, default=1, help="Split grid prep to N chunks for parallel execution")
parser.add_argument('--parallelDock', nargs='?', type=int, default=1, help="Split docking to N chunks for parallel execution")
parser.add_argument('--parallelAntechamber', nargs='?', type=int, default=1, help="Split ligand parametrization to N chunks for parallel execution")
parser.add_argument('--parallelEM', nargs='?', type=int, default=1, help="Split energy minimization to N chunks for parallel execution")
parser.add_argument('--reduceH', nargs='?', type=str, default='.', help="Location of reduce executable (a tool that adds H)")
parser.add_argument('--autogrid', nargs='?', type=str, default='.', help="Location of autogrid4 executable")
parser.add_argument('--autodock', nargs='?', type=str, default='.', help="Location of autodock executable")
parser.add_argument('--config', nargs='?', type=str, default=None, help="File logging locations of executables")

args = parser.parse_args()
print(args)

config = {}
config['autogrid'] = args.autogrid
config['autodock'] = args.autodock
config['reduceH'] = args.reduceH
config['parallelGrid'] = args.parallelGrid
config['parallelDock'] = args.parallelDock
config['parallelAntechamber'] = args.parallelAntechamber
config['parallelEM'] = args.parallelAntechamber

if args.config is not None:
    with open(args.config,'r') as f:
        cont = f.readlines()
    for line in cont:
        try:
            print(line)
            if line.split('=')[0] in ['parallelGrid', 'parallelDock', 'parallelAntechamber', 'parallelEM']:
                try:
                    config[line.split('=')[0]] = int(line.split('=')[1].strip())
                except:
                    pass
            else:
                config[line.split('=')[0]] = line.split('=')[1].strip()
        except:
            pass


for item in ['reduceH','parallelGrid','parallelDock']:
    assert config[item] != '.', f'autogrid executable cannot be " {config[item]} "'

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
for ii in range(config['parallelGrid']):
    fh.append(open(f'../03_Gridmaps/script/prepDock{ii}.sh','w'))


for idx, row in Receptors.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    try:
        #print(f'Trying to make dir ../03_Gridmap/{dname}')
        os.mkdir(f'../03_Gridmaps/{dname}')
    except:
        pass
        #print('mkdir failed!!!!!!!!!')
    write_gpf(f'../03_Gridmaps/{dname}/{row["Receptor_name"]}.gpf',
              row["Receptor_name"], row["dockX"], row["dockY"], row["dockZ"])
    write_prepgrid(fh[idx % config['parallelGrid']], config, dname, row["Receptor_name"], row["Receptor_file_name"])
    write_assembly(f'../03_Gridmaps/{dname}/assembly.py',row["Receptor_name"], row["dockX"], row["dockY"], row["dockZ"])
for ii in range(config['parallelGrid']):
    fh[ii].close()


# Write docking script. To parallelize this part it requires an active allocation.

try:
    os.mkdir('../04_Docking/script')
except:
    pass
fh = []  # File handle array
for ii in range(config['parallelDock']):
    fh.append(open(f'../04_Docking/script/dock{ii}.sh','w'))


for idx, row in Jobs.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    lname = row["Ligand_name"]
    jname = row["Job_name"]
    try:
        #print(f'Trying to make dir ../03_Gridmap/{dname}')
        os.mkdir(f'../04_Docking/{jname}')
    except:
        pass
        #print('mkdir failed!!!!!!!!!')
    write_dock(fh[idx % config['parallelDock']], config, jname, dname, row["Receptor_name"], row["Ligand_name"], row["Ligand_file_name"])

for ii in range(config['parallelDock']):
    fh[ii].close()


# Write antechamber and complex preparation script

try:
    os.makedirs('../05_Refinement/script')
except:
    pass

# Write antechamber script

fh = []
for ii in range(config['parallelAntechamber']):
    fh.append(open(f'../05_Refinement/script/antechamber{ii}.sh','w'))

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
for ii in range(config['parallelAntechamber']):
    fh.append(open(f'../05_Refinement/script/prepareComplex{ii}.sh','w'))

for idx, row in Jobs.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    lname = row["Ligand_file_name"]
    jname = row["Job_name"]
    try:
        os.makedirs(f'../05_Refinement/{jname}/Structure/')
    except:
        pass
    write_prepareComplex(fh[idx % config['parallelAntechamber']], jname, lname, dname)

for ii in range(config['parallelAntechamber']):
    fh[ii].close()

# Write complex assembly prmtop for openmm
fh = []
for ii in range(config['parallelAntechamber']):
    fh.append(open(f'../05_Refinement/script/prepareComplex_openmm{ii}.sh','w'))

for idx, row in Jobs.iterrows():
    dname = f'{row["Receptor_name"]}_{row["dockX"]}_{row["dockY"]}_{row["dockZ"]}'
    lname = row["Ligand_file_name"]
    jname = row["Job_name"]
    try:
        os.makedirs(f'../05_Refinement/{jname}/Structure/')
    except:
        pass
    write_prepareComplex_openmm(fh[idx % config['parallelAntechamber']], jname, lname, dname)

for ii in range(config['parallelAntechamber']):
    fh[ii].close()

fh = []
for ii in range(config['parallelEM']):
    fh.append(open(f'../05_Refinement/script/energyMinimization{ii}.sh','w'))

for idx, row in Jobs.iterrows():
    jname = row["Job_name"]
    write_em(fh[idx % config['parallelEM']], jname)

for ii in range(config['parallelEM']):
    fh[ii].close()

# Write analysis script

try:
    os.makedirs('../06_Analysis/script')
except:
    pass

with open('../06_Analysis/script/MDR_analysis_all.sh', 'w') as f:
    for idx, row in Jobs.iterrows():
        jname = row["Job_name"]
        try:
            os.makedirs(f'../06_Analysis/{jname}')
        except:
            pass
        write_MDR_analysis(f'../06_Analysis/{jname}/MDR_analysis.sh', jname)
        f.write(f'cd ../{jname}\n')
        f.write(f'sh MDR_analysis.sh\n')
