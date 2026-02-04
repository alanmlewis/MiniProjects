#Libraries

import pandas as pd
import pubchempy as pcp
import time

# Reads the dataset

mols = pd.read_csv('SMRT_dataset.csv', delimiter=";")

# Ensure CIDs are valid numbers

valid_cids = mols['pubchem'].dropna().astype(str)
valid_cids = valid_cids[valid_cids.str.isdigit()].astype(int)

# Gets pubchem compounds without getting blocked (hopefully)

def get_compound(cid, max_retries=5):
    for attempt in range(max_retries):
        try:
            return pcp.Compound.from_cid(cid)
        except Exception as e:
            print(f"Attempt {attempt+1} failed for CID {cid}: {e}")
            time.sleep(2 ** attempt)  # Exponential backoff (1s, 2s, 4s, ...)
    print(f"Failed to fetch CID {cid} after {max_retries} attempts.")
    return None
compounds = [get_compound(cid) for cid in valid_cids]

# Removes None values

compounds = [c for c in compounds if c is not None]

# Generates SMILES

smiles = [comp.isomeric_smiles for comp in compounds]

# Exports SMILES 

df = pd.DataFrame(smiles)
print(df)
df.to_csv("SMILES_df.csv", header=False, index=False)