
import parmed
import os, sys
import glob
import re
import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
#import multiprocessing as mp

ps_dir = sys.argv[1]
try:
    ps0 = parmed.load_file(f'{ps_dir}/apo_continue_fix_CA.prmtop',f'{ps_dir}/apo_continue.inpcrd')
except:
    try:
        ps0 = parmed.load_file(f'{ps_dir}/apo_continue.prmtop',f'{ps_dir}/apo_continue.inpcrd')
    except:
        print(f"{ps_dir}/apo_continue files not found! Exiting ...")
        exit()
lig_param_dir = sys.argv[2]
lig_dir = sys.argv[3]
if lig_dir[-1] == '/':
    lig_dir = lig_dir[:-1]
lig_param_suffix = '_H_bcc_sim' 
lig_base_name = sys.argv[4]
lig_names = os.listdir(lig_dir) # should be rank[1-20].pdb

ref_lig_name = glob.glob(f'{lig_param_dir}/*_H.pdb')
#print(ref_lig_name)

lig_prmtop_fname = [x for x in os.listdir(lig_param_dir) if f'{lig_param_suffix}.prmtop' in x][0]

output_inpcrd = 'inpcrd'
output_prmtop = 'prmtop'
output_reference = '../reference_structure'

try:
    os.mkdir(output_inpcrd)
except:
    pass
try:
    os.mkdir(output_prmtop)
except:
    pass
try:
    os.mkdir(output_reference)
except:
    pass

for ii in lig_names:
    try:
            print(f'Ligand {ii}')         # Test1_1.pdb
            lig_base = ii.split('_')[0]   # Test1
            #lig_base = 'complex'
            lig_id = ii.split('.')[0]     # Test1_1
            #lig_prmtop_fname = lig_base + lig_param_suffix + '.prmtop'
            #print(f'{lig_prmtop_fname}')
            ligand_structure = parmed.load_file(f'{lig_param_dir}/{lig_prmtop_fname}',
                                                f'{lig_dir}/{ii}')
            complex_structure = ps0 + ligand_structure
            print('Complex assembled')
            complex_structure.save(f'{output_inpcrd}/{lig_id}.inpcrd', overwrite=True)
            complex_structure.save(f'{output_prmtop}/{lig_base}.prmtop', overwrite=True)
    except:
        try:
            print(f"Ligand {ii} failed initial assembly, let's try something different ...")
            tautomer_file = int(re.sub(r'[^0-9]', '', glob.glob(f'{lig_param_dir}/tautomer*')[0]))
            try:
                os.mkdir('docked_pdb_rdkit')
            except:
                pass
            for ii in lig_names:
                lig_base = ii.split('_')[0]   # Test1
                lig_id = ii.split('.')[0]     # Test1_1
                print(ii)
                a = Chem.MolFromPDBFile(f'../../../04_Docking/{lig_base}/docked_pdb/{ii}')
                enumerator = rdMolStandardize.TautomerEnumerator()
                tauts = enumerator.Enumerate(a)
                Chem.MolToPDBFile(Chem.rdmolops.AddHs(tauts[tautomer_file], addCoords=True), f'docked_pdb_rdkit/{ii}')
                ligand_structure = parmed.load_file(f'{lig_param_dir}/{lig_prmtop_fname}',
                                                    f'docked_pdb_rdkit/{ii}')
                complex_structure = ps0 + ligand_structure
                print('Complex assembled')
                complex_structure.save(f'{output_inpcrd}/{lig_id}.inpcrd', overwrite=True)
                complex_structure.save(f'{output_prmtop}/{lig_base}.prmtop', overwrite=True)
        
            print(tautomer_file)
        except:
            print('Complex assembly failed')

# Finally prepare the reference
print(f'Reference_complex')         # rank1.pdb
#lig_base = ii.split('_')[0]   # 
lig_id = lig_base + '_0'     # This is for reference that we also simulate
#print(f'{lig_prmtop_fname}')
ligand_structure = parmed.load_file(f'{lig_param_dir}/{lig_prmtop_fname}',
                                    ref_lig_name[0])
complex_structure = ps0 + ligand_structure
print('Complex assembled')
ligand_structure.save(f'{output_reference}/ligand.inpcrd', overwrite=True)
ligand_structure.save(f'{output_reference}/ligand.prmtop', overwrite=True)
ps0.save(f'{output_reference}/protein.inpcrd', overwrite=True)
ps0.save(f'{output_reference}/protein.prmtop', overwrite=True)

complex_structure.save(f'{output_reference}/complex.inpcrd', overwrite=True)
complex_structure.save(f'{output_reference}/complex.prmtop', overwrite=True)
complex_structure.save(f'{output_inpcrd}/{lig_id}.inpcrd', overwrite=True)



