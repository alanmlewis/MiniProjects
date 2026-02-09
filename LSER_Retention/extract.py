import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import Counter


data = pd.read_excel('LSER parameters Calculation.xlsx')
data = data[data['Literature']=='Abraham Absolv']

subset = data[["Name","SMILES","E","S","A","B","V"]].dropna()

subset.info()

i = 0
idx = []
ha_list = []
for string in subset["SMILES"].values:
	try:
		mol = Chem.MolFromSmiles(string)
		wt = Descriptors.ExactMolWt(mol)
	except:
		continue
	ha = np.max([atom.GetAtomicNum() for atom in mol.GetAtoms()])
	if wt > 80 and wt < 400 and ha < 26:
		idx.append(i)
		ha_list.append(ha)
	i += 1

print(len(idx))
subset = subset.iloc[idx]

print(Counter(ha_list))

subset.to_csv('data.csv',index=False)
