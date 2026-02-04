# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:55:50 2025

@author: RoisinMurphy
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import Descriptors

#smiles list
smiles_list = np.genfromtxt('SMILES.txt',dtype=str,comments='#######')
#converting smiles into rdkit mol objects
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
#creating a fingerprint generator
MorganFPGen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
#converts the molecules within the list into morgan fingerprint bit vectors
MorganFP = [MorganFPGen.GetFingerprint(mol) for mol in mols]
#converts the fingerprints into a 2D numpy array 
fp = np.array(MorganFP)
print(fp.shape)

features = [Descriptors.CalcMolDescriptors(mol) for mol in mols]
print(len(features),len(features[0]))

labels = list(features[0].keys())

print(labels)

values = [list(feature.values()) for feature in features]
values = np.array(values)
values[np.where(np.isnan(values))]=0

print(values)

np.savetxt('features.txt',values)
np.savetxt('labels.txt',labels,fmt='%s')




