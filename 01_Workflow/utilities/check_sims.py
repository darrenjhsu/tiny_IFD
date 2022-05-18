
import glob, os, sys

#sim_list = glob.glob('../*-*-*/Simulation/*/MD_R*.out')


lig_set = [x for x in os.listdir('../') if 'script' not in x]
print(f'\n\nThere are {len(lig_set)} IFD cases.\n\n')

#lig_set = list(set([x.split('/')[1] for x in sim_list]))

#print(lig_set)

total_sim = 0
total_failed_sim = 0
total_success_sim = 0

for lig in lig_set:
    num_sim = 0
    failed_sim = 0
    success_sim = 0
    finished_sim = 0
    lig_sim_list = glob.glob(f'../{lig}/Simulation/*/MD_R*.out') #[x for x in sim_list if lig in x]
    #print(lig_sim_list)
    for idx, sim in enumerate(lig_sim_list):
        #if (idx + 1) % 20 == 0:
        #    print(idx, end=' ')
        #print(sim)
        with open(sim, 'r') as f:
            cont = f.readlines()
        finished = False
        for line in cont[-36:]:
            if 'Root mean squared' in line:
                finished = True
                finished_sim += 1
            if 'Temperature' in line:
                if float(line.split(':')[-1]) > 500.0:
                    failed_sim += 1
                    #print(f'  {sim} failed')
                elif finished:
                    success_sim += 1
                else:
                    failed_sim += 1
                break
        else:
            failed_sim += 1
        num_sim += 1
    total_sim += num_sim
    total_failed_sim += failed_sim
    total_success_sim += success_sim
    #print()
    print(f'Ligand: {lig} | {num_sim} simulations | {finished_sim} finished | {success_sim} succeeded |', end=' ')
    if num_sim > 0:
        print(f'{failed_sim} failed ({failed_sim / num_sim * 100:.1f} %)')
    else:
        print(' ')
print(f'Total simulations: {total_sim}')
print(f'Total failed simulations: {total_failed_sim}')
print(f'Total success simulations: {total_success_sim}')
print(f'Success percentage: {total_success_sim / total_sim}')
