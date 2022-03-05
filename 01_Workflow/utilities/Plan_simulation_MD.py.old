
import os
import sys
import json
import numpy as np

if len(sys.argv) < 2: # print usage
    print("Usage: python Plan_simulation.py <config_file>")
    exit()


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

# Check QR round is successful
simulation_dir = '../../Simulation'

if os.path.isfile('QR_success.txt'):
    print('Found existing QR tally, use as is!')
    with open('QR_success.txt','r') as f:
        QR_success = f.readlines()[0].split(',')
#    print(QR_success)
else:
    QR_success = []
    for idx, ii in enumerate(inpcrd_list):
        QR_success_this = True
        QR_out_this = os.listdir(f'{simulation_dir}/{ii}')
        QR_out_this = [x for x in QR_out_this if 'QR' in x]
        QR_out_this = [x for x in QR_out_this if 'out' in x]
        if len(QR_out_this) == 0:
            QR_success_this = False
            
        for jj in QR_out_this:
            Temp = 0.0
            with open(f'{simulation_dir}/{ii}/{jj}','r') as f:
                cont = f.readlines()
            for line in cont:
                if 'Temperature' in line:
                    Temp = line.split(':')[-1]
            if np.abs(float(Temp)) > 1000:
                QR_success_this = False
                break
        if QR_success_this:
            QR_success.append(ii)
            print(f'{idx:5d}: Ligand {ii} is successfully tested')
        else:
            print(f'{idx:5d}: Ligand {ii} failed test')
    with open('QR_success.txt','w') as f:
        f.writelines(','.join(QR_success))

inpcrd_list = QR_success
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

# print json
with open(sys.argv[1], 'r') as f:
    settings = json.load(f)

try:
    PERF = settings["Perf_assumption"]
except:
    PERF = 250000.0 # steps / minute for each simulation


# MD round - the big one

MD_Nrep = settings["MD"]["N-rep"]
MD_concurrency = 80 - 80 % MD_Nrep
MD_Nsim = num_sys * MD_Nrep
#MD_Nround = np.ceil(MD_Nsim / NGPU / 80)
MD_min_per_sim = settings["MD"]["cntrl"]["nstlim"] / PERF
MD_max_round = np.floor(120 / MD_min_per_sim)
MD_GPU_when_max_round = np.ceil(MD_Nsim / MD_concurrency / MD_max_round / 6) * 6
MD_Nsim_when_max_round = np.ceil(MD_Nsim / MD_GPU_when_max_round / MD_Nrep) * MD_Nrep

if MD_GPU_when_max_round <= 6:
    MD_GPU_when_max_round = 6
    MD_max_round = np.ceil(MD_Nsim / MD_GPU_when_max_round / MD_concurrency)
    MD_Nsim_when_max_round = MD_Nsim / MD_GPU_when_max_round

try:
    MD_max_wait = settings["MD"]["Max-wait"]
    MD_set_max_wait = True
    MD_max_round_max_wait = np.floor(MD_max_wait / MD_min_per_sim)
    MD_GPU_when_max_wait = np.ceil(MD_Nsim / MD_concurrency / MD_max_round_max_wait / 6) * 6
    MD_Nsim_when_max_wait = np.ceil(MD_Nsim / MD_GPU_when_max_wait / MD_Nrep) * MD_Nrep
    MD_max_round_max_wait = np.ceil(MD_Nsim_when_max_wait / MD_concurrency)
except:
    MD_set_max_wait = False

if MD_concurrency < 80:
    print(f'\n\n  Notice: Since N-rep is set to {MD_Nrep} and not a factor of 80, the concurrency for MD round is reduced to {MD_concurrency}')
print(f'''
*** MD - Molecular Dynamics:  
  We will run dynamics of the {num_sys} unique systems.
  We are running {MD_Nrep} independent simulation(s) for each system, totalling {MD_Nsim} simulations.

  Minimal time to complete this round is to run on {np.ceil(MD_Nsim / MD_concurrency / 6) * 6} GPU cards or {np.ceil(MD_Nsim / MD_concurrency / 6)} nodes.
  Each GPU card takes {MD_concurrency} simulations.
  This would take {MD_min_per_sim:.2f} minutes and {np.ceil(MD_Nsim / MD_concurrency / 6) * MD_min_per_sim / 60:.2f} node hours.

  Alternatively, to use the minimal amount of nodes (and keep sim time < 2 hours), 
  you could run on {MD_GPU_when_max_round} GPU for {MD_max_round} rounds, having each GPU run {MD_Nsim_when_max_round} simulations.
  This would take {MD_max_round * MD_min_per_sim:.2f} minutes and {MD_GPU_when_max_round / 6 * MD_max_round * MD_min_per_sim / 60:.2f} node hours.''')

if MD_set_max_wait:
    print(f'''
  Finally, to keep wait time less than {MD_max_wait} minutes,
  you could run on {MD_GPU_when_max_wait} GPU for {MD_max_round_max_wait} rounds, having each GPU run {MD_Nsim_when_max_wait} simulations.
  This would take {MD_max_round_max_wait * MD_min_per_sim:.2f} minutes and {MD_GPU_when_max_wait / 6 * MD_max_round_max_wait * MD_min_per_sim / 60:.2f} node hours.
''')

if settings["MD"]["Vary-wait"]:
    perf_table = []
    print('  Here is a plan table for your reference:')
    print(f'    Wait, NGPU, Sim/GPU, Nodes, Node Hours')
    for MD_max_wait in range(2,480,15):
        MD_set_max_wait = True
        MD_max_round_max_wait = np.floor(MD_max_wait / MD_min_per_sim)
        if MD_max_round_max_wait == 0:
            continue
        MD_GPU_when_max_wait = np.ceil(MD_Nsim / MD_concurrency / MD_max_round_max_wait / 6) * 6
        MD_Nsim_when_max_wait = np.ceil(MD_Nsim / MD_GPU_when_max_wait / MD_Nrep) * MD_Nrep
        MD_max_round_max_wait = np.ceil(MD_Nsim_when_max_wait / MD_concurrency)
        MD_node_hours = MD_GPU_when_max_wait / 6 * MD_max_round_max_wait * MD_min_per_sim / 60
        #if MD_max_wait == MD_max_round_max_wait * MD_min_per_sim:
        print(f'    {MD_max_wait:4.0f}, {MD_GPU_when_max_wait:4.0f}, {MD_Nsim_when_max_wait:7.0f}, {MD_GPU_when_max_wait/6:5.0f}, {MD_node_hours:10.2f}')

NGPU = input("How many GPUs do you want to use to complete this task? ")
NGPU = int(NGPU)

# Write MD part
systems_assignment = np.linspace(0,num_sys,NGPU+1,dtype=int)

script_mode = "MD"

for ii in range(NGPU):
    
    with open(f'mdgxGPU_MD_{ii}.in','w') as f:
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
            if MD_Nrep > 1:
                f.write(f'  oligomer -p {prmtop_dir}/{prmtop_list[jj]} -c {inpcrd_list[jj]}/EM_R1.rst -o {inpcrd_list[jj]}/MD -x {inpcrd_list[jj]}/MD -r {inpcrd_list[jj]}/MD N-rep {MD_Nrep}\n')
            else:
                f.write(f'  oligomer -p {prmtop_dir}/{prmtop_list[jj]} -c {inpcrd_list[jj]}/EM_R1.rst -o {inpcrd_list[jj]}/MD_R1 -x {inpcrd_list[jj]}/MD_R1 -r {inpcrd_list[jj]}/MD_R1\n')
    # Loop through the oligmer script
        f.write(f'&end\n\n')



with open(f'Ssubmit_MD.sh','w') as f:
    f.write(f'''#!/bin/bash
#BSUB -P STF006
#BSUB -W 2:00
#BSUB -nnodes {int(np.ceil(NGPU / 6))}
#BSUB -J mdgx_test
module load gcc/9.3.0 cuda/11.0.3 cmake readline zlib bzip2 boost netcdf-c netcdf-cxx netcdf-fortran parallel-netcdf openblas netlib-lapack fftw

~/miniconda/bin/conda init bash
source ~/.bashrc
conda activate amber


for i in {{0..{NGPU-1}}};
do
  jsrun -n 1 -g 1 -a 1 -c 1 --smpiargs="off" ~/Tools/amber_build_rhel8/bin/mdgx.cuda -O -i mdgxGPU_MD_${{i}}.in -Reckless &
done
wait
''')



