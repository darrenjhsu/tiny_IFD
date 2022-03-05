
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


# Check EX2 round is successful
simulation_dir = '../../Simulation'



if os.path.isfile('EX2_success.txt'):
    print('Found existing EX2 tally, use as is!')
    with open('EX2_success.txt','r') as f:
        EX2_success = f.readlines()[0].split(',')
#    print(EX2_success)
else:
    EX2_success = []
    counter = 0
    for idx, ii in enumerate(inpcrd_list):
        EX2_success_this = True
        EX2_out_this = os.listdir(f'{simulation_dir}/{ii}')
        EX2_out_this = [x for x in EX2_out_this if 'EX2' in x]
        EX2_out_this = [x for x in EX2_out_this if 'out' in x]
        #print(EX2_out_this)
        EX2_out_this.sort()
        for jj in EX2_out_this:
            Temp = 0.0
            with open(f'{simulation_dir}/{ii}/{jj}','r') as f:
                cont = f.readlines()
            for line in cont:
                if 'Temperature' in line:
                    Temp = line.split(':')[-1]
            if np.abs(float(Temp)) > 1000:
                EX2_success_this = False
                break
            if EX2_success_this:
                EX2_success.append(f"{ii}/{jj.split('.')[0]}")
                print(f'{counter:5d}: Ligand {ii}\'s EX2 round {jj} is successfully tested')
            else:
                print(f'{counter:5d}: Ligand {ii} failed test')
            counter += 1
    with open('EX2_success.txt','w') as f:
        f.writelines(','.join(EX2_success))

inpcrd_list = EX2_success
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


# EX3 round - the big one

EX3_Nrep = settings["EX3"]["N-rep"]
if EX3_Nrep > 1:
    print("Currently can't do multiples in the EX rounds ... exiting")
    exit()
EX3_concurrency = 80 - 80 % EX3_Nrep
EX3_Nsim = num_sys * EX3_Nrep
#EX3_Nround = np.ceil(EX3_Nsim / NGPU / 80)
EX3_min_per_sim = settings["EX3"]["cntrl"]["nstlim"] / PERF
EX3_max_round = np.floor(120 / EX3_min_per_sim)
EX3_GPU_when_max_round = np.ceil(EX3_Nsim / EX3_concurrency / EX3_max_round / 6) * 6
EX3_Nsim_when_max_round = np.ceil(EX3_Nsim / EX3_GPU_when_max_round / EX3_Nrep) * EX3_Nrep

if EX3_GPU_when_max_round <= 6:
    EX3_GPU_when_max_round = 6
    EX3_max_round = np.ceil(EX3_Nsim / EX3_GPU_when_max_round / EX3_concurrency)
    EX3_Nsim_when_max_round = EX3_Nsim / EX3_GPU_when_max_round

try:
    EX3_max_wait = settings["EX3"]["Max-wait"]
    EX3_set_max_wait = True
    EX3_max_round_max_wait = np.floor(EX3_max_wait / EX3_min_per_sim)
    EX3_GPU_when_max_wait = np.ceil(EX3_Nsim / EX3_concurrency / EX3_max_round_max_wait / 6) * 6
    EX3_Nsim_when_max_wait = np.ceil(EX3_Nsim / EX3_GPU_when_max_wait / EX3_Nrep) * EX3_Nrep
    EX3_max_round_max_wait = np.ceil(EX3_Nsim_when_max_wait / EX3_concurrency)
except:
    EX3_set_max_wait = False

if EX3_concurrency < 80:
    print(f'\n\n  Notice: Since N-rep is set to {EX3_Nrep} and not a factor of 80, the concurrency for EX3 round is reduced to {EX3_concurrency}')
print(f'''
*** EX3 - Molecular Dynamics:  
  We will run dynamics of the {num_sys} unique systems.
  We are running {EX3_Nrep} independent simulation(s) for each system, totalling {EX3_Nsim} simulations.

  Minimal time to complete this round is to run on {np.ceil(EX3_Nsim / EX3_concurrency / 6) * 6} GPU cards or {np.ceil(EX3_Nsim / EX3_concurrency / 6)} nodes.
  Each GPU card takes {EX3_concurrency} simulations.
  This would take {EX3_min_per_sim:.2f} minutes and {np.ceil(EX3_Nsim / EX3_concurrency / 6) * EX3_min_per_sim / 60:.2f} node hours.

  Alternatively, to use the minimal amount of nodes (and keep sim time < 2 hours), 
  you could run on {EX3_GPU_when_max_round} GPU for {EX3_max_round} rounds, having each GPU run {EX3_Nsim_when_max_round} simulations.
  This would take {EX3_max_round * EX3_min_per_sim:.2f} minutes and {EX3_GPU_when_max_round / 6 * EX3_max_round * EX3_min_per_sim / 60:.2f} node hours.''')

if EX3_set_max_wait:
    print(f'''
  Finally, to keep wait time less than {EX3_max_wait} minutes,
  you could run on {EX3_GPU_when_max_wait} GPU for {EX3_max_round_max_wait} rounds, having each GPU run {EX3_Nsim_when_max_wait} simulations.
  This would take {EX3_max_round_max_wait * EX3_min_per_sim:.2f} minutes and {EX3_GPU_when_max_wait / 6 * EX3_max_round_max_wait * EX3_min_per_sim / 60:.2f} node hours.
''')

if settings["EX3"]["Vary-wait"]:
    perf_table = []
    print('  Here is a plan table for your reference:')
    print(f'    Wait, NGPU, Sim/GPU, Nodes, Node Hours')
    for EX3_max_wait in range(2,480,15):
        EX3_set_max_wait = True
        EX3_max_round_max_wait = np.floor(EX3_max_wait / EX3_min_per_sim)
        if EX3_max_round_max_wait == 0:
            continue
        EX3_GPU_when_max_wait = np.ceil(EX3_Nsim / EX3_concurrency / EX3_max_round_max_wait / 6) * 6
        EX3_Nsim_when_max_wait = np.ceil(EX3_Nsim / EX3_GPU_when_max_wait / EX3_Nrep) * EX3_Nrep
        EX3_max_round_max_wait = np.ceil(EX3_Nsim_when_max_wait / EX3_concurrency)
        EX3_node_hours = EX3_GPU_when_max_wait / 6 * EX3_max_round_max_wait * EX3_min_per_sim / 60
        #if EX3_max_wait == EX3_max_round_max_wait * EX3_min_per_sim:
        print(f'    {EX3_max_wait:4.0f}, {EX3_GPU_when_max_wait:4.0f}, {EX3_Nsim_when_max_wait:7.0f}, {EX3_GPU_when_max_wait/6:5.0f}, {EX3_node_hours:10.2f}')

NGPU = input("How many GPUs do you want to use to complete this task? ")
NGPU = int(NGPU)

# Write EX3 part
systems_assignment = np.linspace(0,num_sys,NGPU+1,dtype=int)

script_mode = "EX3"

for ii in range(NGPU):
    
    with open(f'mdgxGPU_EX3_{ii}.in','w') as f:
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
            out_base = inpcrd_list[jj].replace("EX2","EX3")
            f.write(f'  oligomer -p {prmtop_dir}/{prmtop_list[jj]} -c {inpcrd_list[jj]}.rst -o {out_base} -x {out_base} -r {out_base}\n')
    # Loop through the oligmer script
        f.write(f'&end\n\n')



with open(f'Ssubmit_EX3.sh','w') as f:
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
  jsrun -n 1 -g 1 -a 1 -c 1 --smpiargs="off" ~/Tools/amber_build_test/bin/mdgx.cuda -O -i mdgxGPU_EX3_${{i}}.in &
done
wait
''')



