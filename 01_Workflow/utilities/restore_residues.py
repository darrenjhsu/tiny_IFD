import sys

fbase = sys.argv[1]

with open(f'{fbase}_d.pdb','r') as f, open(f'{fbase}_full_h.pdb','r') as g:
    refcont = f.readlines()
    tarcont = g.readlines()

ref_res = []
tar_res = []
for line in refcont:
    if 'ATOM' in line:
        if len(ref_res) == 0:
            ref_res.append(line[17:26])
        if line[17:26] != ref_res[-1]:
            ref_res.append(line[17:26])
for line in tarcont:
    if 'ATOM' in line:
        if len(tar_res) == 0:
            tar_res.append(line[17:26])
        if line[17:26] != tar_res[-1]:
            tar_res.append(line[17:26])
tar_ref = {}




for k, v in zip(tar_res, ref_res):
    print(k, v)
    tar_ref[k] = v


new_tarcont= []
for line in tarcont:
    if 'ATOM' in line:
        new_tarcont.append(line[:17] + tar_ref[line[17:26]] + line[26:])
    else:
        new_tarcont.append(line)

with open(f'{fbase}_restored_h.pdb', 'w') as f:
    f.writelines(new_tarcont)
