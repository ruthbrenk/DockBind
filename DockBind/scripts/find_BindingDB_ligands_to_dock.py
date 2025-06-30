import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

print('Loading datasets...')
df = pd.read_csv('../data/selected_templates_updated_2024.csv')
ch = pd.read_csv('../data/BindingDB_activities_for_BindingMOAD_2024.csv')

uni_w_act = ch['UniProt (SwissProt) Primary ID of Target Chain'].drop_duplicates().to_list()
bm_ch = df[df.UniProt_ID.isin(uni_w_act)].reset_index(drop=True)

print('Starting SMILES check...')

bad_id = []
for i, ch_id, smi in ch[['BindingDB MonomerID', 'Ligand SMILES']].itertuples():
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad_id.append(ch_id)
            
    except:
        continue
print('Starting ligand analysis...')
to_append = []
for i, uni, template, smi in bm_ch[['UniProt_ID', 'Template', 'Smiles_String']].itertuples():
    if i > 5000 :
        continue
    try:
        print(f'{i+1}/{len(bm_ch)}')
        temp_df = ch[ch['UniProt (SwissProt) Primary ID of Target Chain'] == uni]
        temp_df = temp_df.drop_duplicates(subset=['BindingDB MonomerID'])

        ref_mol = Chem.MolFromSmiles(smi)
        ref_fps = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, 2048)

        for j, ch_id, ch_smi in temp_df[['BindingDB MonomerID', 'Ligand SMILES']].itertuples():

            if ch_id in bad_id:
                continue
            else:
                ch_mol = Chem.MolFromSmiles(ch_smi)

                ch_fps = AllChem.GetMorganFingerprintAsBitVect(ch_mol, 2, 2048)

                ts = round(DataStructs.TanimotoSimilarity(ref_fps, ch_fps), 2)

                to_append.append({'UniProt_ID': uni,
                             'BM_Template': template,
                             'BM_Smiles': smi,
                             'BindingDB_ID': ch_id,
                             'BindingDB_Smiles': ch_smi,
                             'Tanimoto_Similarity': ts
                            })
    except:
        continue

TS_df = pd.DataFrame(to_append)
TS_df.to_csv('../data/full_BindingMOAD_BindingDB_Tanimoto_similarity_1.csv', index=False)

ch_selected_for_docking = TS_df[TS_df.Tanimoto_Similarity >= 0.7].reset_index(drop=True)
ch_selected_for_docking.to_csv('../data/BindingDB_ligands_selected_for_docking_1.csv', index=False)

