
import os, glob
import sys
import json
import numpy as np
import pandas as pd

if len(sys.argv) < 2: # print usage
    print("Usage: python Plan_simulation.py <config_file>")
    exit()


Jobs = pd.read_csv('../../02_Input/job_description.csv')#, column=['Receptor_file_name','Ligand_file_name','dockX','dockY','dockZ'])

inpcrd_list = []
prmtop_list = []

for idx, row in Jobs.iterrows():
    jname = row['Job_name']
    inpcrd_dir = f'../{jname}/Structure/inpcrd' #lig_base + '_complexes'
    prmtop_dir = f'../{jname}/Structure/prmtop'

    inpcrd_listing = glob.glob(f'{inpcrd_dir}/*')
    prmtop_listing = glob.glob(f'{prmtop_dir}/*')
    inpcrd_listing.sort()
    prmtop_listing.sort()

    #print(inpcrd_listing, prmtop_listing)

    # Create a corresponding prmtop_list
    prmtop_list += [prmtop_listing[0] for x in inpcrd_listing]
    inpcrd_list += inpcrd_listing

#print(list(zip(prmtop_list, inpcrd_list)))
num_sys = len(inpcrd_list)
sim_folder = ['../' + x.split('/')[1] + '/Simulation/' + x.split('/')[-1].split('.')[0] for x in inpcrd_list] # ../Job_name/Simulation/rank1 etc
#print(num_sys)

# load json config file
with open(sys.argv[1], 'r') as f:
    settings = json.load(f)

try:
    PERF = settings["Perf_assumption"]
except:
    PERF = 250000.0 # steps / minute for each simulation

# Plan simulation - we want to use all cards and minimize rounds of simulations

EM_Nrep = settings["EM"]["N-rep"]
EM_concurrency = 80 - 80 % EM_Nrep
EM_Nsim = num_sys * EM_Nrep
#EM_Nround = np.ceil(EM_Nsim / NGPU / 80)
EM_min_per_sim = settings["EM"]["cntrl"]["nstlim"] / PERF
EM_max_round = 120 / EM_min_per_sim
EM_Nsim_when_max_round = EM_Nsim / EM_max_round
EM_GPU_when_max_round = EM_Nsim_when_max_round / EM_concurrency

if EM_Nsim <= 80:
    EM_GPU_when_max_round = 1
    EM_max_round = 1
    EM_Nsim_when_max_round = EM_Nsim  
elif EM_GPU_when_max_round < 6:
    EM_GPU_when_max_round = 6
    EM_max_round = np.ceil(EM_Nsim / EM_GPU_when_max_round / EM_concurrency)
    EM_Nsim_when_max_round = np.ceil(EM_Nsim / EM_GPU_when_max_round / EM_concurrency) * EM_concurrency

    if EM_concurrency < 80:
        print(f'\n\n  Notice: Since N-rep is set to {EM_Nrep} and not a factor of 80, the concurrency for EM round is reduced to {EM_concurrency}')

print(f'''
*** EM - Energy minimization:
  There are {num_sys} unique systems whose energy we need to minimize.
  We are running {EM_Nrep} independent simulation(s) for each system, totalling {EM_Nsim} simulations.

  Minimal time to complete this round is to run on {np.ceil(EM_Nsim / EM_concurrency / 6) * 6} GPU cards or {np.ceil(EM_Nsim / EM_concurrency / 6)} nodes.
  Each GPU card takes {EM_concurrency} simulations.
  This would take {EM_min_per_sim:.2f} minutes and {np.ceil(EM_Nsim / EM_concurrency / 6) * EM_min_per_sim / 60:.2f} node hours.

  Alternatively, to use the minimal amount of nodes (and keep sim time < 2 hours), 
  you could run on {EM_GPU_when_max_round} GPU for {EM_max_round} rounds, having each GPU run {EM_Nsim_when_max_round} simulations.
  This would take {EM_max_round * EM_min_per_sim:.2f} minutes and {EM_GPU_when_max_round / 6 * EM_max_round * EM_min_per_sim / 60:.2f} node hours.
''')

NGPU = input("How many GPUs do you want to use to complete this task? ")
NGPU = int(NGPU)

#exit()
    # write EM part
systems_assignment = np.linspace(0,num_sys,NGPU+1,dtype=int)

script_mode = "EM"

for ii in range(NGPU):
    with open(f'mdgxGPU_EM_{ii}.in','w') as f:
        f.write(f'&files\n')
        for key in settings[script_mode]["files"]:
          f.write(f'  {key}  {settings[script_mode]["files"][key]}\n')
        f.write(f'&end\n\n')
    
        f.write(f'&cntrl\n')
        for key in settings[script_mode]["cntrl"]:
          f.write(f'  {key} = {settings[script_mode]["cntrl"][key]},\n')
        f.write(f'&end\n\n')
    
        f.write(f'&pptd\n')
        for key in settings[script_mode]["pptd"]:
          f.write(f'  {key} = {settings[script_mode]["pptd"][key]},\n')
        for jj in range(systems_assignment[ii], systems_assignment[ii+1]):
            if EM_Nrep > 1:
                f.write(f'  oligomer -p {prmtop_list[jj]} -c {inpcrd_list[jj]} -o {sim_folder[jj]}/EM -x {sim_folder[jj]}/EM -r {sim_folder[jj]}/EM N-rep {EM_Nrep}\n')
            else:
                f.write(f'  oligomer -p {prmtop_list[jj]} -c {inpcrd_list[jj]} -o {sim_folder[jj]}/EM_R1 -x {sim_folder[jj]}/EM_R1 -r {sim_folder[jj]}/EM_R1\n')
            # print(f'blah blah for system {jj}')
    # Loop through the oligmer script
        f.write(f'&end\n\n')
    

with open(f'EM_script.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    for ii in sim_folder:
        f.write(f'mkdir -p {ii}\n')

if NGPU == 1: # Also write local script for execution
    with open(f'EM_local.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('sh EM_script.sh\n')
        f.write('mdgx.cuda -O -i mdgxGPU_EM_0.in -Reckless &\n')
        f.write('wait')

with open(f'Ssubmit_EM.sh','w') as f:
    f.write(f'''#!/bin/bash
#BSUB -P STF006
#BSUB -W 2:00
#BSUB -nnodes {int(np.ceil(NGPU / 6))}
#BSUB -J mdgx_test
module load gcc/9.3.0 cuda/11.0.3 cmake readline zlib bzip2 boost netcdf-c netcdf-cxx netcdf-fortran parallel-netcdf  openblas netlib-lapack fftw

~/miniconda/bin/conda init bash
source ~/.bashrc
conda activate amber

sh EM_script.sh

for i in {{0..{NGPU-1}}};
do
  jsrun -n 1 -g 1 -a 1 -c 1 --smpiargs="off" mdgx.cuda -O -i mdgxGPU_EM_${{i}}.in -Reckless &
done

wait
''')
