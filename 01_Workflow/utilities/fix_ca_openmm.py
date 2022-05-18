
# Parse prmtop file and generate a new prmtop file with C O N CA's masses set to 9.999E+03.

import numpy as np

with open('apo_continue.prmtop','r') as f:
    cont = f.readlines()


get_atom_name = False
get_atom_mass = False
get_residue_label = False
get_residue_pointer = False
atom_names = []
atom_mass = []
residue_label = []
residue_pointer = []
for line in cont:
    if 'ATOM_NAME' in line:
#        print(line)
        get_atom_name = True
        continue
    if 'CHARGE' in line:
        get_atom_name = False
    if 'MASS' in line:
        get_atom_mass = True
        continue
    if 'ATOM_TYPE_INDEX' in line:
        get_atom_mass = False
        continue
    if 'RESIDUE_LABEL' in line:
        get_residue_label = True
        continue
    if 'RESIDUE_POINTER' in line:
        get_residue_label = False
        get_residue_pointer = True
        continue
    if 'BOND_FORCE_CONSTANT' in line:
        get_residue_pointer = False
        continue

    if get_atom_name == True and not 'FORMAT' in line:
#        print(line)
#        print(len(line))
        for ii in range(0, len(line.strip('\n')),4):
#            print(line[ii:ii+4])
            atom_names.append(line[ii:ii+4])

    if get_atom_mass == True and not 'FORMAT' in line:
        for ii in range(0, len(line.strip('\n')),16):
            atom_mass.append(line[ii:ii+16])

    if get_residue_label == True and not 'FORMAT' in line:
        for ii in range(0, len(line.strip('\n')), 4):
#            print(line[ii:ii+4])
            residue_label.append(line[ii:ii+4])

    if get_residue_pointer == True and not 'FORMAT' in line:
        for ii in range(0, len(line.strip('\n')), 8):
#            print(line[ii:ii+8])
            residue_pointer.append(int(line[ii:ii+8]))
         
    
#print(atom_names)
#print(atom_mass)

print(len(residue_pointer))

residue_pointer = np.array(residue_pointer)
res_name = []
#res_name = [''] * len(atom_names)
for ii in range(len(residue_pointer)):
#    print(ii, residue_label[ii])
    if ii == len(residue_pointer) - 1:
        for jj in range(residue_pointer[ii]-1, len(atom_names)):
            res_name.append(residue_label[ii])
#            print(jj, residue_label[ii])
    else:
        for jj in range(residue_pointer[ii]-1, residue_pointer[ii+1]-1):
            res_name.append(residue_label[ii])
#            print(jj, residue_label[ii])
#print(len(res_name))

counter = 0
for idx, (ii, jj) in enumerate(zip(atom_names,atom_mass)):
#    print(idx, counter, ii, jj, res_name[idx])
    if ii == 'CA  ':
        if res_name[idx] != 'WAT ':
            atom_mass[counter] = '  0.00000000E+00'
#    print(counter, ii, atom_mass[counter])
    counter += 1


get_atom_name = False
get_atom_mass = False
counter = 0
for idx, line in enumerate(cont):
    if 'ATOM_NAME' in line:
#        print(line)
        get_atom_name = True
        continue
    if 'CHARGE' in line:
        get_atom_name = False
    if 'MASS' in line:
        get_atom_mass = True
        continue
    if 'ATOM_TYPE_INDEX' in line:
        get_atom_mass = False

    if get_atom_mass == True and not 'FORMAT' in line:
#        print(line)
        new_line = ''
        for ii in range(0, 5):
            try:
                new_line += atom_mass[counter]
#                print(counter)
                counter += 1
            except:
                pass
        new_line += '\n'
#        print(new_line)
        cont[idx] = new_line

with open('apo_continue_fix_CA_openmm.prmtop','w') as f:
    f.writelines(cont)
#print(cont)
