
import os, glob, sys
import numpy as np

results = glob.glob('docking_results/*')

result_table = []
for fname in results:
    with open(fname, 'r') as f:
        cont = f.readlines()
        result_table.append([float(x) for x in cont[0].split()[1:]])

result_table = np.array(result_table)

print(result_table)
print(np.sum(np.any(result_table[:,:2] < 2.5, axis=1)),
      np.sum(np.any(result_table[:,:5] < 2.5, axis=1)),
      np.sum(np.any(result_table < 2.5, axis=1)),
      len(result_table))
