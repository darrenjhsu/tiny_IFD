
import sys, os

assert len(sys.argv) > 0, "Must supply the folder name containing the docking.dlg"
path = sys.argv[1]
jname = sys.argv[2]
if path[-1] == '/':
    path = path[:-1] # Remove the slash

with open(f'{path}/docking.dlg','r') as f:
    cont = f.readlines()
PDBQT = {}
Current_run = None
run_selection = {}
for line in cont:

    if 'Run:' in line:
        PDBQT[int(line.split()[1])] = []
        Current_run = int(line.split()[1])

    if 'DOCKED:' in line:
        PDBQT[Current_run].append(line.replace('DOCKED: ',''))

subrank = 1
counter = 1
while len(run_selection) < 20:
    current_len = len(run_selection)
    for line in cont:
        if 'RANKING' in line:
            if line.split()[1] == str(subrank):
                #run_selection[int(line.split()[0])] = int(line.split()[2])
                run_selection[counter] = int(line.split()[2])
                print(f'Selecting rank {line.split()[0]} and sub-rank {line.split()[1]} which is run {line.split()[2]}')
                counter += 1
    subrank += 1
    print(f"Right now there are {len(run_selection)} logged poses.")
    if len(run_selection) == current_len: # meaning the list did not grow, this will become infinite loop
        print(f"Let's break because list did not grow")
        break
    else:
        print(f"List grew, let's continue!")

try:
    os.mkdir(f'{path}/docked')
except:
    pass
try:
    for rank in range(1,21):
        with open(f'{path}/docked/{jname}_{rank}.pdbqt','w') as f:
            f.writelines(PDBQT[run_selection[rank]])
except:
    os.remove(f'{path}/docked/{jname}_{rank}.pdbqt')
