import pandas as pd
import numpy as np
import time
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdFMCS, rdSubstructLibrary
from sklearn.cluster import AgglomerativeClustering
import os

BASE = '..'
PATH_TO_DATA = f'{BASE}/data'
PATH_TO_PREPROCESSED = f'{BASE}/preprocessed_data'
PATH_TO_CLUSTERING_OUTPUT = f'{BASE}/clustering_output'

def UPGMA_Clustering(uniprot, mol, dists):
    
    cluster = AgglomerativeClustering(n_clusters=2,
                                      compute_full_tree=True,
                                      metric='precomputed',
                                      linkage='average')
    cluster.fit(dists)

    NumMolList, MolDict = CalcSizeAndAssignment(cluster.children_,len(dists))

    MCSdict = DetermineRelevantMCS(len(dists), mol, cluster.children_, MolDict)

    return MCSdict, MolDict


def GetDistMat(uniprot, save=bool):
    '''This function will generate the distance matrix and a list of mol objects for each target.
    The DistMat will be saved in the respective target folder'''
    smiles_list = []
    ligand_names = []
    
    with open(f'{PATH_TO_PREPROCESSED}/{uniprot}_Smiles_List_w_Names.txt','r') as f:
        for line in f:
            smiles_list.append(line.split('  ')[0])
            ligand_names.append(line.split('  ')[1][:8])
    mol = [Chem.MolFromSmiles(smile) for smile in smiles_list]
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024) for x in mol]
    
    dists=np.zeros([len(smiles_list),len(smiles_list)])
    for i in range(len(smiles_list)):
        for j in range(i):
            dists[i,j]=1-DataStructs.TanimotoSimilarity(fps[i],fps[j])
            dists[j,i]=dists[i,j]
    if save==True:
        df1 = pd.DataFrame(dists, index=ligand_names, columns=ligand_names)
        df1.to_csv(f'{PATH_TO_CLUSTERING_OUTPUT}/{uniprot}/Distance_Matrix.csv', index=True)
    
    return mol, dists, ligand_names


def CalcSizeAndAssignment(children, Ndata):
    NumMolList = []
    MolDict = {}
    for i in range(len(children)):
        N = 0
        mols_assigned = []
        for j in range(len(children[i])):
            if children[i][j]<Ndata:
                N+=1
                mols_assigned.append(children[i][j])
            else:
                N+=NumMolList[children[i][j]-Ndata]
                mols_assigned+=MolDict[children[i][j]]
        NumMolList.append(N)
        MolDict[i+Ndata] = mols_assigned

    return NumMolList, MolDict


def DetermineRelevantMCS(Ndata,mol,children,MolDict):
    # filter out irrelevant clusters and calculate MCS on selected clusters
    currlayer=[Ndata*2-2]
    MCSdict={}
    while len(currlayer)>0:
        childlayer=[]
        for c in currlayer:
            if c>=Ndata:
                if len(MolDict[c])>=2:
                    print(c,' MCS iteration started')
                    moldata = []
                    for index in MolDict[c]:
                        moldata.append(mol[index])
                    print('Calculating FChembl score...')
                    fChembl,Smarts=MCSFromMollist(moldata, sslib)
                    print(f'Score calculated: {fChembl}')
                    moldata.clear()
                    if fChembl>=1e-3:
                        childlayer+=children[c-Ndata].tolist()
                    else:
                        MCSdict[c]=(fChembl,len(MolDict[c]),Smarts)
        currlayer=childlayer
        
    return MCSdict


def MCSFromMollist(mollist, sslib):
    MCSSmarts2=rdFMCS.FindMCS(mollist,
                              atomCompare=rdFMCS.AtomCompare.CompareAny,
                              bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                              ringMatchesRingOnly=False,
                              timeout=1).smartsString
    MCSSmarts=rdFMCS.FindMCS(mollist,
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=False,
                             timeout=1).smartsString
    if MCSSmarts2 == '':
        fChembl2 = 1
    else:
        qry = Chem.MolFromSmarts(MCSSmarts2)
        matches2 = sslib.GetMatches(qry, maxResults=5000)
        fChembl2 = round((len(matches2)+1)/(len(sslib)+2), 6)
    if MCSSmarts == '': 
        fChembl1 = 1
    else:
        qry = Chem.MolFromSmarts(MCSSmarts)
        matches1 = sslib.GetMatches(qry, maxResults=5000)
        fChembl1 = round((len(matches1)+1)/(len(sslib)+2), 6)
    if fChembl2 < fChembl1:
        fChembl1 = fChembl2
        MCSSmarts = MCSSmarts2
    return fChembl1, MCSSmarts

print('ChEMBL library loading...')
with open(f'{PATH_TO_DATA}/chembl34_ssslib.pkl','rb') as inf:
    sslib = pickle.load(inf)
print('Library loaded successfully')

os.makedirs(PATH_TO_CLUSTERING_OUTPUT, exist_ok=True)

df = pd.read_csv(f'../data/astex_selected_for_docking.csv', usecols=['UniProt_ID', 'PDB', 'Lig', 'Smiles_String'])
uniprots = df.UniProt_ID.unique()

for _ in uniprots:
    os.makedirs(f'{PATH_TO_CLUSTERING_OUTPUT}/{_}', exist_ok=True)

start_time = time.time()

for i, uniprot in enumerate(uniprots):
    
    t1 = time.time()
    print('Starting with UniProt: '+uniprot)
    with open(f'{PATH_TO_CLUSTERING_OUTPUT}/clustering_log.txt', 'a') as output_txt:
        output_txt.write('Starting with uniprot:'+ uniprot+'\n')

    mol, dists, ligand_names = GetDistMat(uniprot, True)

    MCSdict, MolDict = UPGMA_Clustering(uniprot, mol, dists)

    MCS = pd.DataFrame.from_dict(MCSdict)

    MCS.to_csv(f'{PATH_TO_CLUSTERING_OUTPUT}/{uniprot}/{uniprot}_MCS_dict.csv')

    temp_list = []
    keys = list(MolDict.keys())
    for key in MolDict.keys():
        temp_list.append(MolDict[key])
    #print(temp_list)
    with open(f'{PATH_TO_CLUSTERING_OUTPUT}/{uniprot}/{uniprot}_MolDict.txt', 'w') as f:
        for i in range(len(temp_list)):
            f.write(str(keys[i]))
            for j in range(len(temp_list[i])):
                f.write(f'  {ligand_names[temp_list[i][j]]}')
            f.write('\n')
    tf = time.time() - t1

    with open(f'{PATH_TO_CLUSTERING_OUTPUT}/clustering_log.txt', 'a') as output_txt:
        output_txt.write(f'{uniprot}: done successfully\nTIME: {time.strftime("%H:%M:%S",time.gmtime(tf))}\n\n')

    print(f'UniProt: {uniprot} done\n{time.strftime("%H:%M:%S",time.gmtime(tf))}')

seconds = time.time() - start_time

print(f'Total time: {time.strftime("%H:%M:%S",time.gmtime(seconds))}')
with open(f'{PATH_TO_CLUSTERING_OUTPUT}/clustering_log.txt', 'a') as output_txt:
    output_txt.write(f'\n\n\nTOTAL TIME: {time.strftime("%H:%M:%S",time.gmtime(seconds))}')
