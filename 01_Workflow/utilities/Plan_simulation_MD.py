
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


num_sys = len(inpcrd_list)
sim_folder = ['../' + x.split('/')[1] + '/Simulation/' + x.split('/')[-1].split('.')[0] for x in inpcrd_list] # ../Job_name/Simulation/rank1 etc


if os.path.isfile('QR_success.txt'):
    print('Found existing QR tally, use as is!')
    with open('QR_success.txt','r') as f:
        QR_success = f.readlines()[0].split(',')
#    print(QR_success)
else:
    QR_success = []
    QR_success_count = 0
    print('Succ. / Total')
    for idx, ii in enumerate(sim_folder):
        QR_success_this = True
        QR_out_this = os.listdir(f'{ii}/')
        QR_out_this = [x for x in QR_out_this if 'QR' in x]
        QR_out_this = [x for x in QR_out_this if 'out' in x]
        if len(QR_out_this) == 0:
            QR_success_this = False

        for jj in QR_out_this:
            Temp = 0.0
            with open(f'{ii}/{jj}','r') as f:
                cont = f.readlines()
            for line in cont:
                if 'Temperature' in line:
                    Temp = line.split(':')[-1]
            if np.abs(float(Temp)) > 1000:
                QR_success_this = False
                break
        if QR_success_this:
            QR_success_count += 1
            QR_success.append(ii)
            print(f'{QR_success_count:5d} / {idx:5d}: Ligand {ii} is successfully tested')
        else:
            print(f'{QR_success_count:5d} / {idx:5d}: Ligand {ii} failed test')
    with open('QR_success.txt','w') as f:
        f.writelines(','.join(QR_success))


inpcrd_list = QR_success
#print(inpcrd_list)
sim_folder = ['../' + x.split('/')[1] + '/Simulation/' + x.split('/')[-1].split('.')[0] for x in inpcrd_list] # ../Job_name/Simulation/rank1 etc
prmtop_list =  [x.split('_')[0].replace('Simulation','Structure/prmtop') + '.prmtop' for x in inpcrd_list]
#print(sim_folder)
#print(prmtop_list)

num_sys = len(inpcrd_list)


# print json
with open(sys.argv[1], 'r') as f:
    settings = json.load(f)

try:
    PERF = settings["Perf_assumption"]
except:
    PERF = 250000.0 # steps / minute for each simulation

# MD round
MD_Nrep = settings["MD"]["N-rep"]
MD_concurrency = 80 - 80 % MD_Nrep
MD_Nsim = num_sys * MD_Nrep
#MD_Nround = np.ceil(MD_Nsim / NGPU / 80)
MD_min_per_sim = settings["MD"]["cntrl"]["nstlim"] / PERF
MD_max_round = np.floor(120 / MD_min_per_sim)
MD_GPU_when_max_round = np.ceil(MD_Nsim / MD_concurrency / MD_max_round / 6) * 6
MD_Nsim_when_max_round = np.ceil(MD_Nsim / MD_GPU_when_max_round / MD_Nrep) * MD_Nrep


if MD_Nsim <= 80:
    MD_GPU_when_max_round = 1
    MD_max_round = 1
    MD_Nsim_when_max_round = MD_Nsim  
if MD_GPU_when_max_round <= 6:
    MD_GPU_when_max_round = 6
    MD_max_round = np.ceil(MD_Nsim / MD_GPU_when_max_round / MD_concurrency)
    MD_Nsim_when_max_round = np.ceil(MD_Nsim / MD_GPU_when_max_round / 6) * 6

try:
    MD_max_wait = settings["MD"]["Max-wait"]
    MD_set_max_wait = True
    MD_max_round_max_wait = np.floor(MD_max_wait / MD_min_per_sim)
    MD_GPU_when_max_wait = np.ceil(MD_Nsim / MD_concurrency / MD_max_round_max_wait / 6) * 6
    MD_Nsim_when_max_wait = np.ceil(MD_Nsim / MD_GPU_when_max_wait)
    MD_max_round_max_wait = np.ceil(MD_Nsim_when_max_wait / MD_concurrency)
except:
    MD_set_max_wait = False

if MD_concurrency < 80:
    print(f'\n\n  Notice: Since N-rep is set to {MD_Nrep} and not a factor of 80, the concurrency for MD round is reduced to {MD_concurrency}')
print(f'''
*** MD - Molecular Dynamics:  
  We will test the stability of the {num_sys} unique systems.
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
    current_GPU_when_max_wait = None
    for MD_max_wait in range(5, 480, 1):
        MD_set_max_wait = True
        MD_max_round_max_wait = np.floor(MD_max_wait / MD_min_per_sim)
        if MD_max_round_max_wait == 0:
            continue
        MD_GPU_when_max_wait = np.ceil(MD_Nsim / MD_concurrency / MD_max_round_max_wait / 6) * 6
        MD_Nsim_when_max_wait = np.ceil(MD_Nsim / MD_GPU_when_max_wait)
        MD_max_round_max_wait = np.ceil(MD_Nsim_when_max_wait / MD_concurrency)
        MD_node_hours = MD_GPU_when_max_wait / 6 * MD_max_round_max_wait * MD_min_per_sim / 60
        if (np.abs(MD_max_wait - MD_max_round_max_wait * MD_min_per_sim) < 1) and MD_GPU_when_max_wait != current_GPU_when_max_wait:
            current_GPU_when_max_wait = MD_GPU_when_max_wait
            print(f'    {MD_max_wait:4.0f}, {MD_GPU_when_max_wait:4.0f}, {MD_Nsim_when_max_wait:7.0f}, {MD_GPU_when_max_wait/6:5.0f}, {MD_node_hours:10.2f}')



NGPU = input("How many GPUs do you want to use to complete this task? ")
NGPU = int(NGPU)
# write MD part
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
                f.write(f'  oligomer -p {prmtop_list[jj]} -c {sim_folder[jj]}/QR_R1.rst -o {sim_folder[jj]}/MD -x {sim_folder[jj]}/MD -r {sim_folder[jj]}/MD N-rep {MD_Nrep}\n')
            else:
                f.write(f'  oligomer -p {prmtop_list[jj]} -c {sim_folder[jj]}/QR_R1.rst -o {sim_folder[jj]}/MD_R1 -x {sim_folder[jj]}/MD_R1 -r {sim_folder[jj]}/MD_R1\n')
    # Loop through the oligmer script
        f.write(f'&end\n\n')


if NGPU == 1: # Also write local script for execution
    with open(f'MD_local.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('mdgx.cuda -O -i mdgxGPU_MD_0.in -Reckless &')
        f.write('wait')

with open(f'Ssubmit_MD.sh','w') as f:
    f.write(f'''#!/bin/bash
#BSUB -P STF006
#BSUB -W 2:00
#BSUB -nnodes {int(np.ceil(NGPU / 6))}
#BSUB -J mdgx_test
module load gcc/9.3.0 cuda/11.0.3 cmake readline zlib bzip2 boost netcdf-c netcdf-cxx netcdf-fortran parallel-netcdf  openblas netlib-lapack fftw

~/miniconda/bin/conda init bash
source ~/.bashrc
conda activate amber


for i in {{0..{NGPU-1}}};
do
  jsrun -n 1 -g 1 -a 1 -c 1 --smpiargs="off" mdgx.cuda -O -i mdgxGPU_MD_${{i}}.in -Reckless &
done
wait
''')


