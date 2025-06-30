import os
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from tqdm import tqdm
from urllib.request import urlretrieve as ur

files = [i for i in os.listdir('.') if 'csv' in i]

df = pd.DataFrame()
for file in files:
    temp = pd.read_csv(file, skiprows=1)
    df = df.append(temp).reset_index(drop=True)

df = df.drop(columns=['Entry ID', 'Release Date', 'Unnamed: 7'])

cof = pd.read_csv('../known_cofactors.csv')
buf = pd.read_csv('../crystal_buffer_components.csv')
ama = pd.read_csv('../amino_acids.csv')

df['PDB ID'] = df['PDB ID'].fillna(method='ffill')
df['Refinement Resolution (Å)'] = df['Refinement Resolution (Å)'].fillna(method='ffill')
df.info()

df = df.dropna(subset=['Ligand SMILES', 'Ligand ID'])
df = df[~df['Ligand ID'].isin(buf.code.unique())]
df = df[~df['Ligand ID'].isin(cof.code.unique())]
df = df[~df['Ligand ID'].isin(ama.HET.unique())]
df = df[df['Ligand ID'].str.len() == 3].reset_index(drop=True)
df

uniprot = pd.read_csv(f'idmapping_2024_07_02.tsv', sep='\t')
uniprot = uniprot.drop_duplicates(subset='From').reset_index(drop=True)
pdb_uni = dict(zip(uniprot.From, uniprot.Entry))

# Uniprot mapping csv
for i, u in df[['PDB ID']].itertuples():
    print(f'{i+1}/{len(df)}')
    try:
        df.at[i, 'UniProt_ID'] = pdb_uni[u]
    except Exception as e:
        print(e)
        df.at[i, 'UniProt_ID'] = 'Not found'

df = df[df['UniProt_ID'] != 'Not found']

uniprots = df.UniProt_ID.unique()

bad_id = []
for i, ch_id, smi in df[['Ligand Name', 'Ligand SMILES']].itertuples():
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad_id.append(ch_id)
            
    except:
        continue

confirmed_pdb = df[~df['Ligand Name'].isin(bad_id)]

for i, pdb in enumerate(confirmed_pdb['PDB ID'].unique()):
    try:
        print(f'Staritng {i+1}/{len(df)}')
        uni = pdb_uni[pdb]
        os.makedirs(f'pdb_files_not_aligned/{uni}', exist_ok=True)
        ur(f'https://files.rcsb.org/download/{pdb}.pdb', filename=f'pdb_files_not_aligned/{uni}/{pdb}.pdb')
    except:
        print(f'{pdb} not found')
















