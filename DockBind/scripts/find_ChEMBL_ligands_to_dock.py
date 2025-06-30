import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

print('Loading dataframes...')

df = pd.read_csv('../data/selected_templates_updated_2024.csv')
ch = pd.read_csv('../data/ChEMBL_activities_for_bindingMOAD_2024.csv')

uni_w_act = ch.accession.drop_duplicates().to_list()
bm_ch = df[df.UniProt_ID.isin(uni_w_act)].reset_index(drop=True)

print('Starting ligand analysis')

to_append = []
for i, uni, template, smi in bm_ch[['UniProt_ID', 'Template', 'Smiles_String']].itertuples():
    try:
        print(f'{i+1}/{len(bm_ch)}')
        temp_df = ch[ch.accession == uni]
        temp_df = temp_df.drop_duplicates(subset=['molecule_chembl_id'])

        ref_mol = Chem.MolFromSmiles(smi)
        ref_fps = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, 2048)

        for j, ch_id, ch_smi in temp_df[['molecule_chembl_id', 'canonical_smiles']].itertuples():
            ch_mol = Chem.MolFromSmiles(ch_smi)
            ch_fps = AllChem.GetMorganFingerprintAsBitVect(ch_mol, 2, 2048)

            ts = round(DataStructs.TanimotoSimilarity(ref_fps, ch_fps), 2)

            to_append.append({'UniProt_ID': uni,
                         'BM_Template': template,
                         'BM_Smiles': smi,
                         'ChEMBL_ID': ch_id,
                         'ChEMBL_Smiles': ch_smi,
                         'Tanimoto_Similarity': ts
                        })
    except Exception as e:
        print(e)
        
TS_df = pd.DataFrame(to_append)
TS_df.to_csv('../data/full_BindingMOAD_ChEMBL_Tanimoto_similarity.csv', index=False)

ch_selected_for_docking = TS_df[TS_df.Tanimoto_Similarity >= 0.7].reset_index(drop=True)
ch_selected_for_docking.to_csv('../data/ChEMBL_ligands_selected_for_docking.csv', index=False)
