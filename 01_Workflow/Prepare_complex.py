
import parmed
import os
import multiprocessing as mp


ps0 = parmed.load_file('apo_continue_fix_CA.prmtop','apo_continue.inpcrd')

lig_dir = 'SDF_match_PDB'
lig_param_dir = 'lig_prmtop'
lig_param_suffix = '_H_bcc_sim' 

lig_names = os.listdir(lig_dir)

output_inpcrd = 'inpcrd'
output_prmtop = 'prmtop'


def assemble(ii):
    try:
        print(f'Ligand {ii}')         # Mpro-x10959_0.pdb
        lig_base = ii.split('_')[0]   # Mpro-x10959
        lig_id = ii.split('.')[0]     # Mpro-x10959_0
        lig_prmtop_fname = lig_base + lig_param_suffix + '.prmtop'
        print(f'{lig_prmtop_fname}')
        ligand_structure = parmed.load_file(f'{lig_param_dir}/{lig_prmtop_fname}',
                                            f'{lig_dir}/{ii}')
        complex_structure = ps0 + ligand_structure
        print('Complex assembled')
        complex_structure.save(f'{output_inpcrd}/{lig_id}.inpcrd', overwrite=True)
        complex_structure.save(f'{output_prmtop}/{lig_base}.prmtop', overwrite=True)
    except:
        print(f'Ligand {ii} failed')

        
if __name__ == '__main__':
    with mp.Pool(24) as p:
        p.map(assemble, lig_names)
