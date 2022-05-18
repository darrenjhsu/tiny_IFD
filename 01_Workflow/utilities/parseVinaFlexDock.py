
import sys, os

assert len(sys.argv) > 0, "Must supply the folder name containing the docking.dlg"
path = sys.argv[1]
jname = sys.argv[2]
flex_path = sys.argv[3]
if path[-1] == '/':
    path = path[:-1] # Remove the slash

with open(f'{path}/docking.dlg','r') as f:
    cont = f.readlines()
PDBQT = {}
FLEX = {}
Current_run = None
run_selection = {}
for line in cont:

    if 'MODEL' in line:
        PDBQT[int(line.split()[1])] = []
        FLEX[int(line.split()[1])] = []
        Current_run = int(line.split()[1])
        record = 'model'
    if 'BEGIN_RES' in line:
        record = 'flex'
    if 'ENDMOL' in line:
        PDBQT[Current_run].append(line)
        record = 'model'
        continue 

    if record == 'model':
        PDBQT[Current_run].append(line)
    elif record == 'flex':
        FLEX[Current_run].append(line)

try:
    os.mkdir(f'{path}/docked')
except:
    pass
try:
    os.mkdir(f'{flex_path}/flex')
except:
    pass
try:
    for rank in range(1,21):
        with open(f'{path}/docked/{jname}_{rank}.pdbqt','w') as f:
            f.writelines(PDBQT[rank])
except:
    os.remove(f'{path}/docked/{jname}_{rank}.pdbqt')
try:
    for rank in range(1,21):
        with open(f'{flex_path}/flex/{jname}_{rank}.pdbqt','w') as f:
            f.writelines(FLEX[rank])
except:
    os.remove(f'{flex_path}/flex/{jname}_{rank}.pdbqt')
