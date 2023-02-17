
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

    if 'MODEL' in line:
        PDBQT[int(line.split()[1])] = []
        Current_run = int(line.split()[1])

    PDBQT[Current_run].append(line)

try:
    os.mkdir(f'{path}/docked')
except:
    pass
try:
    for rank in range(1, len(PDBQT)+1):
        with open(f'{path}/docked/{jname}_{rank}.pdbqt','w') as f:
            f.writelines(PDBQT[rank])
except:
    os.remove(f'{path}/docked/{jname}_{rank}.pdbqt')
