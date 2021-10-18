from rdkit import Chem
from rdkit.Chem import *
from rdkit.Chem import AllChem
from os import listdir
from re import *
import os



#print(listdir())
'''
def print(*args, **a):
	if type(args[0]) is Chem.rdChem.Atom:
		print(Chem.MolToSmiles(args[0])) 
		'''

initial = 'CH6O6P2'

assert 'N' not in initial
noO = int(search('O[0-9]*', initial)[0][1:])
noP = int(search('P[0-9]*', initial)[0][1:])
max = min(noO//3, noP)
ss = []
for i in range(2, max+1):
	init = initial
	init = sub('O[0-9]*', f'O{noO-3*i}', init).replace('O0','')
	#print(init)
	init = sub('P[0-9]*', f'P{noP-i}', init).replace('P0','')
	#print(init)
	init += f'N{i}'
	print(init)
	os.system(f'java -jar MAYGEN-1.8.jar -f "{init}" -v -t -smi -m 24 -o odir')
	ss += open(f'odir/{init}.smi', 'r').readlines()

#s = lambda i: f'[{Chem.MolToSmiles(i), [*map(s, i.GetAtoms())]}]  ' if isinstance(i,Chem.rdchem.Mol) else f'{i.GetIdx(), i.GetSymbol()} '

#a = open('odir/CH6N2.smi', 'r').readlines()
#print('read', a)
for s in ss:
	assert all(dig not in s for dig in '0123456789')
	#m = Chem.MolFromSmiles(s)
	
	### Approach SMILES
	#s = Chem.MolToSmiles(m)
	#iis = [i for i, atom in enumerate(s) if atom=='N']
	#i = iis[0]
	
	out = [[*string] for string in ss]
	#print(s)
	#for i in iis:

print('out', out)
for i in range(s.count('N')):
	temp = [t[:] for t in out]
	out = []
	for infix in ('OP(=O)O', 'P(=O)(O)O',  'OP(=O)(O)'):
		for l in temp:
			ll = [*l]
			print(l, ''.join(l))
			i = ''.join(l).index('N')
			ll[i:i+1] = infix # 'P(=O)(O)'  'OP(=O)'
			news = ''.join(ll)
			out.append(news)
	print(out,'\n==================\n')
print(i, s, *out)

s = set(Chem.MolToSmiles(Chem.MolFromSmiles(mol)) for mol in out)

for i, tt in enumerate(s):
	print(i,tt)
'''print(m, *m.GetAtoms())

qs = [i for i, atom in enumerate(m.GetAtoms()) if atom.GetAtomicNum()==7]
print(qs)
a = m.GetAtoms()[qs[1]] # 1,N
neis = [nei.GetEndAtomIdx() for nei in a.GetBonds()]

print(a.GetSymbol(), a.GetBonds(), a.GetNeighbors()[0].GetSymbol(), a.GetNeighbors()[-1].GetBonds())
n = Chem.ReplaceSubstructs(m, Chem.MolFromSmarts('N'), Chem.MolFromSmarts('O'))[0]
n = Chem.ReplaceSubstructs(n, Chem.MolFromSmarts('N'), Chem.MolFromSmarts('O'))[0]

template = Chem.MolFromSmiles('CNN')
AllChem.Compute2DCoords(template)
#AllChem.GenerateDepictionMatching2DStructure(n)

print(Chem.MolToSmiles(n))


from rdkit.Chem.Draw import rdMolDraw2D
#mol = Chem.MolFromSmiles('Cl[C@H](F)NC\C=C\C')
d = rdMolDraw2D.MolDraw2DCairo(250, 200) # or MolDraw2DSVG to get SVGs
#mol.GetAtomWithIdx(2).SetProp('atomNote', 'foo')
#mol.GetBondWithIdx(0).SetProp('bondNote', 'bar')
#d.drawOptions().addStereoAnnotation = True
#d.drawOptions().addAtomIndices = True

#n2 = Chem.AddHs(m)
#d.DrawMolecule(n2)
#d.FinishDrawing()
#d.WriteDrawingText('atom_annotation_1.png')  
'''

