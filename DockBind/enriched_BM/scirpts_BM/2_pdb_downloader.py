import pandas as pd
import requests
from urllib.request import urlretrieve as ur
import os


df = pd.read_csv('../data/clean_BindingMOAD.csv')


TEMP = f'../pdb_files_not_aligned'
print('DOWNLOAD PDB FILES\n\n')
for i, uni, pdb in df[['UniProt_ID', 'Protein_ID']].itertuples():
    print(f'Staritng {i+1}/{len(df)}')
    if os.path.exists(f'{TEMP}/{uni}/{pdb}.pdb') == False:
        os.makedirs(f'{TEMP}/{uni}', exist_ok=True)
        ur(f'https://files.rcsb.org/download/{pdb}.pdb', filename=f'{TEMP}/{uni}/{pdb}.pdb')

