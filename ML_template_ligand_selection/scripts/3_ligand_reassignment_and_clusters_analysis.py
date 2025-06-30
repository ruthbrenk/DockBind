import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdFMCS
import os
import time

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def run_the_analysis(uniprot):
    smiles, ligand_name, clusters_keys, mcs_mol = load_required_files(uniprot)
    if len(clusters_keys) == 0:
        with open(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/Clusters_ID_with_matches.txt', 'a') as out:
            out.write('No relevant clusters have been found')
    else:
        find_matches(smiles, ligand_name, clusters_keys, mcs_mol)
        get_images(smiles, ligand_name, clusters_keys, mcs_mol)
        with open(f'{PATH_TO_CLUSTER_ANALYSIS}/Updated_UniProtID_list.txt', 'a') as up_uni:
            up_uni.write(uniprot+'\n')
            
    
def load_required_files(uniprot):
    '''This function generate the data required for further analysis. It loads:
            1. Smiles_List
            2. MCS_dict
        And outputs:
            1. Reference dictionary with Smiles as keys and Ligand_Name as value
            2. A list with Clusters_Keys and a list with the respective Mols object
                        '''
    with open(f'{PATH_TO_PREPROCESSED}/{uniprot}_Smiles_List_w_Names.txt','r') as f:
        smiles = []
        ligand_name = []
        for line in f:
            smiles.append(line.split('  ')[0])
            ligand_name.append(line.split('  ')[1][:-1])
    
    clusters = pd.read_csv(f'{PATH_TO_CLUSTERING_OUTPUT}/{uniprot}/{uniprot}_MCS_dict.csv', index_col=0).to_dict()
    clusters_keys = list(clusters.keys())
    mcs_mol = [Chem.MolFromSmarts(clusters[_][2]) for _ in clusters_keys]
        
    return smiles, ligand_name, clusters_keys, mcs_mol



    
def find_matches(smiles, ligand_name, clusters_keys, mcs_mol):
    '''This function reassign the SMILES to each cluster and saves two files:
        1. for each cluster it saves a txt file with the matched SMILES and the corresponding ligand name
        2. one file with'''

    for i in range(len(clusters_keys)):
        os.makedirs(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/{clusters_keys[i]}', exist_ok=True)
        counter = 0
        with open(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/{clusters_keys[i]}/Matches_Smiles_Ligand_Name.txt', 'w') as f:
            for lig in range(len(smiles)):
                if Chem.MolFromSmiles(smiles[lig]).HasSubstructMatch(mcs_mol[i]):
                    f.write(smiles[lig]+'  '+ligand_name[lig]+'\n')
                    counter = counter+1
        with open(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/Clusters_ID_with_matches.txt', 'a') as matched:
            matched.write('Cluster number '+clusters_keys[i]+' has matched molecules:  '+str(counter)+'\n')
                    
def get_images(smiles, ligand_name, clusters_keys, mcs_mol):

    mcss = Draw.MolsToGridImage(mcs_mol, subImgSize=(1000,1000))
    mcss.save(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/MCS_Structures.png')
    for i in range(len(clusters_keys)):
        with open(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/{clusters_keys[i]}/Matches_Smiles_Ligand_Name.txt', 'r') as f:
            mols = []
            for line in f:
                mols.append(Chem.MolFromSmiles(line.split('  ')[0]))
        highlights = [_.GetSubstructMatch(mcs_mol[i]) for _ in mols]
        try:
            img = Draw.MolsToGridImage(mols, highlightAtomLists=highlights, subImgSize=(1024,1024))
            img.save(f'{PATH_TO_CLUSTER_ANALYSIS}/{uniprot}/{clusters_keys[i]}/Ligands_with_MCS_highlighted.png')
        except Exception as e:
            print(e)
            print(f'Cannot save image for uniprot {uniprot} cluster number {clusters_keys[i]}')

            
BASE = '..'
PATH_TO_DATA = f'{BASE}/data'
PATH_TO_PREPROCESSED = f'{BASE}/preprocessed_data'
PATH_TO_CLUSTERING_OUTPUT = f'{BASE}/clustering_output'
PATH_TO_CLUSTER_ANALYSIS = f'{BASE}/clusters_analysis'

os.makedirs(PATH_TO_CLUSTERING_OUTPUT, exist_ok=True)

df = pd.read_csv(f'../data/astex_selected_for_docking.csv', usecols=['UniProt_ID', 'PDB', 'Lig', 'Smiles_String'])
df = df[df.UniProt_ID != 'P04818']
uniprots = df.UniProt_ID.unique()
len(uniprots)

for _ in uniprots:
    os.makedirs(f'{PATH_TO_CLUSTER_ANALYSIS}/{_}', exist_ok=True)

start_time = time.time()

with open(f'{PATH_TO_CLUSTER_ANALYSIS}/cluster_analysis_log.txt', 'a') as output_txt:
    for uniprot in uniprots:
        print(f'Starting with UniProt: {uniprot}')
        output_txt.write(f'Starting with uniprot: {uniprot}\n')
        run_the_analysis(uniprot)
        output_txt.write(f'{uniprot}: done successfully \n \n')
        print(f'UniProt: {uniprot} done')

    seconds = time.time() - start_time

    print(f'Time Taken: {time.strftime("%H:%M:%S",time.gmtime(seconds))}')
    output_txt.write(f'\n\n\nTIME TAKEN: {time.strftime("%H:%M:%S",time.gmtime(seconds))}')            

