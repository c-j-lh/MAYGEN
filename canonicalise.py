from rdkit import Chem
from tqdm import *

with open("odir/CH5O2PS.smi", 'r') as inp:
	with open("odir/canon-CH5O2PS.smi", 'w') as out:
		for line in tqdm(inp):
			canon = Chem.MolToSmiles(Chem.MolFromSmiles(line), True)
			#print(canon)
			out.write(canon + '\n')