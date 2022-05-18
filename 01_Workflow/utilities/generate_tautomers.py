import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import sys

a = Chem.MolFromPDBFile(sys.argv[1])

enumerator = rdMolStandardize.TautomerEnumerator()
tauts = enumerator.Enumerate(a)
for idx, mol in enumerate(tauts):
    Chem.MolToMolFile(Chem.rdmolops.AddHs(mol, addCoords=True), f'tautomer{idx:02d}.sdf')
