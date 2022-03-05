
import os
import sys
import json
import numpy as np

if len(sys.argv) < 2: # print usage
    print("Usage: python Plan_simulation.py <config_file>")
    exit()


# print json
with open(sys.argv[1], 'r') as f:
    settings = json.load(f)

try:
    PERF = settings["Perf_assumption"]
except:
    PERF = 250000.0 # steps / minute for each simulation

#lig_base = 'Z1530724813'
#lig_prefix = lig_base + '_1_T1'
#lig_suffix = '_H_bcc_sim'
#param_dir = lig_base + '_parameter_tests'
inpcrd_dir = '../Structure/inpcrd' #lig_base + '_complexes'
prmtop_dir = '../Structure/prmtop'

inpcrd_listing = os.listdir('../'+inpcrd_dir)
prmtop_listing = os.listdir('../'+prmtop_dir)
inpcrd_listing.sort()
prmtop_listing.sort()


# Create a corresponding prmtop_list
inpcrd_list = [x.split('.')[0] for x in inpcrd_listing]


# Check EX3 round is successful
simulation_dir = '../../Simulation'



if os.path.isfile('EX3_success.txt'):
    print('Found existing EX3 tally, use as is!')
    with open('EX3_success.txt','r') as f:
        EX3_success = f.readlines()[0].split(',')
#    print(EX3_success)
else:
    EX3_success = []
    counter = 0
    for idx, ii in enumerate(inpcrd_list):
        EX3_success_this = True
        EX3_out_this = os.listdir(f'{simulation_dir}/{ii}')
        EX3_out_this = [x for x in EX3_out_this if 'EX3' in x]
        EX3_out_this = [x for x in EX3_out_this if 'out' in x]
        #print(EX3_out_this)
        EX3_out_this.sort()
        for jj in EX3_out_this:
            Temp = 0.0
            with open(f'{simulation_dir}/{ii}/{jj}','r') as f:
                cont = f.readlines()
            for line in cont:
                if 'Temperature' in line:
                    Temp = line.split(':')[-1]
            if np.abs(float(Temp)) > 1000:
                EX3_success_this = False
                break
            if EX3_success_this:
                EX3_success.append(f"{ii}/{jj.split('.')[0]}")
                print(f'{counter:5d}: Ligand {ii}\'s EX3 round {jj} is successfully tested')
            else:
                print(f'{counter:5d}: Ligand {ii} failed test')
            counter += 1
    with open('EX3_success.txt','w') as f:
        f.writelines(','.join(EX3_success))

inpcrd_list = EX3_success
inpcrd_ul = [x.split('_')[0] for x in inpcrd_list] # ul is unique list
prmtop_ul = [x.split('.')[0] for x in prmtop_listing if 'prmtop' in x]
prmtop_list = [x + '.prmtop' for x in inpcrd_ul]

# Check all inpcrd has corresponding prmtop in the folder
for ii in inpcrd_ul:
    #print(ii.split('_')[0])
    if ii.split('_')[0] not in prmtop_ul:
        print(f'There is no {ii}.prmtop in {prmtop_dir} folder!')
        exit()

print("File check passed. All inpcrd have corresponding prmtop.")

num_sys = len(inpcrd_ul)

print(f'Num sys: {num_sys}')
#exit()


# EX4 round - the big one

EX4_Nrep = settings["EX4"]["N-rep"]
if EX4_Nrep > 1:
    print("Currently can't do multiples in the EX rounds ... exiting")
    exit()
EX4_concurrency = 80 - 80 % EX4_Nrep
EX4_Nsim = num_sys * EX4_Nrep
#EX4_Nround = np.ceil(EX4_Nsim / NGPU / 80)
EX4_min_per_sim = settings["EX4"]["cntrl"]["nstlim"] / PERF
EX4_max_round = np.floor(120 / EX4_min_per_sim)
EX4_GPU_when_max_round = np.ceil(EX4_Nsim / EX4_concurrency / EX4_max_round / 6) * 6
EX4_Nsim_when_max_round = np.ceil(EX4_Nsim / EX4_GPU_when_max_round / EX4_Nrep) * EX4_Nrep

if EX4_GPU_when_max_round <= 6:
    EX4_GPU_when_max_round = 6
    EX4_max_round = np.ceil(EX4_Nsim / EX4_GPU_when_max_round / EX4_concurrency)
    EX4_Nsim_when_max_round = EX4_Nsim / EX4_GPU_when_max_round

try:
    EX4_max_wait = settings["EX4"]["Max-wait"]
    EX4_set_max_wait = True
    EX4_max_round_max_wait = np.floor(EX4_max_wait / EX4_min_per_sim)
    EX4_GPU_when_max_wait = np.ceil(EX4_Nsim / EX4_concurrency / EX4_max_round_max_wait / 6) * 6
    EX4_Nsim_when_max_wait = np.ceil(EX4_Nsim / EX4_GPU_when_max_wait / EX4_Nrep) * EX4_Nrep
    EX4_max_round_max_wait = np.ceil(EX4_Nsim_when_max_wait / EX4_concurrency)
except:
    EX4_set_max_wait = False

if EX4_concurrency < 80:
    print(f'\n\n  Notice: Since N-rep is set to {EX4_Nrep} and not a factor of 80, the concurrency for EX4 round is reduced to {EX4_concurrency}')
print(f'''
*** EX4 - Molecular Dynamics:  
  We will run dynamics of the {num_sys} unique systems.
  We are running {EX4_Nrep} independent simulation(s) for each system, totalling {EX4_Nsim} simulations.

  Minimal time to complete this round is to run on {np.ceil(EX4_Nsim / EX4_concurrency / 6) * 6} GPU cards or {np.ceil(EX4_Nsim / EX4_concurrency / 6)} nodes.
  Each GPU card takes {EX4_concurrency} simulations.
  This would take {EX4_min_per_sim:.2f} minutes and {np.ceil(EX4_Nsim / EX4_concurrency / 6) * EX4_min_per_sim / 60:.2f} node hours.

  Alternatively, to use the minimal amount of nodes (and keep sim time < 2 hours), 
  you could run on {EX4_GPU_when_max_round} GPU for {EX4_max_round} rounds, having each GPU run {EX4_Nsim_when_max_round} simulations.
  This would take {EX4_max_round * EX4_min_per_sim:.2f} minutes and {EX4_GPU_when_max_round / 6 * EX4_max_round * EX4_min_per_sim / 60:.2f} node hours.''')

if EX4_set_max_wait:
    print(f'''
  Finally, to keep wait time less than {EX4_max_wait} minutes,
  you could run on {EX4_GPU_when_max_wait} GPU for {EX4_max_round_max_wait} rounds, having each GPU run {EX4_Nsim_when_max_wait} simulations.
  This would take {EX4_max_round_max_wait * EX4_min_per_sim:.2f} minutes and {EX4_GPU_when_max_wait / 6 * EX4_max_round_max_wait * EX4_min_per_sim / 60:.2f} node hours.
''')

if settings["EX4"]["Vary-wait"]:
    perf_table = []
    print('  Here is a plan table for your reference:')
    print(f'    Wait, NGPU, Sim/GPU, Nodes, Node Hours')
    for EX4_max_wait in range(2,480,15):
        EX4_set_max_wait = True
        EX4_max_round_max_wait = np.floor(EX4_max_wait / EX4_min_per_sim)
        if EX4_max_round_max_wait == 0:
            continue
        EX4_GPU_when_max_wait = np.ceil(EX4_Nsim / EX4_concurrency / EX4_max_round_max_wait / 6) * 6
        EX4_Nsim_when_max_wait = np.ceil(EX4_Nsim / EX4_GPU_when_max_wait / EX4_Nrep) * EX4_Nrep
        EX4_max_round_max_wait = np.ceil(EX4_Nsim_when_max_wait / EX4_concurrency)
        EX4_node_hours = EX4_GPU_when_max_wait / 6 * EX4_max_round_max_wait * EX4_min_per_sim / 60
        #if EX4_max_wait == EX4_max_round_max_wait * EX4_min_per_sim:
        print(f'    {EX4_max_wait:4.0f}, {EX4_GPU_when_max_wait:4.0f}, {EX4_Nsim_when_max_wait:7.0f}, {EX4_GPU_when_max_wait/6:5.0f}, {EX4_node_hours:10.2f}')

NGPU = input("How many GPUs do you want to use to complete this task? ")
NGPU = int(NGPU)

# Write EX4 part
systems_assignment = np.linspace(0,num_sys,NGPU+1,dtype=int)

script_mode = "EX4"

for ii in range(NGPU):
    
    with open(f'mdgxGPU_EX4_{ii}.in','w') as f:
        f.write(f'&files\n')
        for key in settings[script_mode]["files"]:
          f.write(f'  {key}  {settings[script_mode]["files"][key]}\n')
        f.write(f'&end\n\n')
    
        f.write(f'&cntrl\n')
        for key in settings[script_mode]["cntrl"]:
          f.write(f'  {key} =  {settings[script_mode]["cntrl"][key]},\n')
        f.write(f'&end\n\n')
    
        f.write(f'&pptd\n')
        for key in settings[script_mode]["pptd"]:
          f.write(f'  {key} =  {settings[script_mode]["pptd"][key]},\n')
        for jj in range(systems_assignment[ii], systems_assignment[ii+1]):
            #out_base = inpcrd_list[jj].split('/')[0]
            #round_base = inpcrd_list[jj].split('_')[-1].split('.')[0]
            out_base = inpcrd_list[jj].replace("EX3","EX4")
            f.write(f'  oligomer -p {prmtop_dir}/{prmtop_list[jj]} -c {inpcrd_list[jj]}.rst -o {out_base} -x {out_base} -r {out_base}\n')
    # Loop through the oligmer script
        f.write(f'&end\n\n')



with open(f'Ssubmit_EX4.sh','w') as f:
    f.write(f'''#!/bin/bash
#BSUB -P STF006
#BSUB -W 2:00
#BSUB -nnodes {int(np.ceil(NGPU / 6))}
#BSUB -J mdgx_test
module load gcc/9.3.0 cuda/11.0.3 cmake readline zlib bzip2 boost netcdf netcdf-cxx4 netcdf-fortran parallel-netcdf openblas netlib-lapack fftw

~/miniconda/bin/conda init bash
source ~/.bashrc
conda activate amber


for i in {{0..{NGPU-1}}};
do
  jsrun -n 1 -g 1 -a 1 -c 1 --smpiargs="off" ~/Tools/amber_build_test/bin/mdgx.cuda -O -i mdgxGPU_EX4_${{i}}.in &
done
wait
''')



