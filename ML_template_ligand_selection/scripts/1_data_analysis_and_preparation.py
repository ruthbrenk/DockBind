import os
import pandas as pd
from matplotlib import pyplot as plt

BASE = '..'
PATH_TO_DATA = f'{BASE}/data'
PATH_TO_OUTPUT = f'{BASE}/preprocessed_data'

df = pd.read_csv(f'../data/astex_selected_for_docking.csv', usecols=['UniProt_ID', 'PDB', 'Lig', 'Smiles_String'])
df['Template'] = df['PDB'] + '_' + df['Lig']

df = df.drop(columns=['PDB', 'Lig'])

grouped = df.groupby('UniProt_ID', sort=False)['Template'].apply(list).reset_index(name='Template')
smiles = df.groupby('UniProt_ID', sort=False)['Smiles_String'].apply(list).tolist()
grouped['Smiles_String'] = smiles

grouped['len'] = grouped['Template'].str.len()
grouped.sort_values(by='len', ascending=False).drop(columns='len').reset_index(drop=True)

uniprots = grouped['UniProt_ID'].tolist()
ligands = grouped['Template'].tolist()
smiles = grouped['Smiles_String'].tolist()

for i in range(len(uniprots)):
    if len(smiles[i]) > 1:
        os.makedirs(f'{PATH_TO_OUTPUT}', exist_ok=True)
        with open(f'{PATH_TO_OUTPUT}/{uniprots[i]}_Smiles_List_w_Names.txt', 'w') as output:
            for j in range(len(smiles[i])):
                output.write(smiles[i][j] + '  ' + ligands[i][j]+'\n')

