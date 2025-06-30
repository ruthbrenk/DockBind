import pymysql
import pandas as pd
import numpy as np
import os
# from chembl_webresource_client.new_client import new_client
import chembl_webresource_client
from pymol import cmd
import requests
from dotenv import load_dotenv
import math
import statistics
from pymol import cmd

os.getcwd()

print('Starting BindingMOAD preprocessing', flush=True)
# BindingMOAD
bm = pd.read_csv(f'../data/every.csv', names=['Protein_Class', 'Protein_Family', 'Protein_ID',
                                           'Ligand_Name', 'Ligand_Validity', 'Affinity_Measure',
                                           'Relation', 'Affinity_Value', 'Affinity_Unit', 'Smiles_String'])

# Cofactors
cof = pd.read_csv(f'../data/known_cofactors.csv').dropna()

# Buffer components
cbc = pd.read_csv(f'../data/crystal_buffer_components.csv')

# PDB-UniProt mapping
uniprot = pd.read_csv(f'../data/pdb_chain_uniprot.csv',
                      skiprows=[0], usecols=['PDB', 'SP_PRIMARY']).drop_duplicates(subset=['PDB'], keep='first')
uniprot['PDB'] = uniprot['PDB'].str.upper()
pdb_uni = dict(zip(uniprot['PDB'], uniprot['SP_PRIMARY']))

# 1st fill PDB lines
bm = bm.iloc[1:]
bm['Protein_ID'].fillna(method='ffill', inplace=True)
bm = bm.dropna(subset=['Ligand_Name'])

# 2nd select onli valid ligands
valid_ligands = bm[bm['Ligand_Validity'] == 'valid']
valid_ligands[['Lig_Name', 'Chain', 'Resn']] = valid_ligands['Ligand_Name'].str.split(pat=':', expand=True)
valid_ligands = valid_ligands.drop(columns=['Protein_Class', 'Protein_Family', 'Ligand_Name', 'Ligand_Validity'])

# If PDB has ligands with binding affinity: remove all other entries
# Else: drop duplicates and keep one row of each ligand
valid_ligands = valid_ligands[valid_ligands['Affinity_Value'].notna() | 
                              valid_ligands['Affinity_Value'].isna()].drop_duplicates(subset=['Protein_ID', 'Lig_Name'])

# Uniprot mapping csv
for i, u in valid_ligands[['Protein_ID']].itertuples():
    try:
        valid_ligands.at[i, 'UniProt_ID'] = pdb_uni[u]
    except:
        valid_ligands.at[i, 'UniProt_ID'] = np.nan

# Drop if UniProtID not available and reorder columns
confirmed_uniprots = valid_ligands.dropna(subset=['UniProt_ID']).reset_index(drop=True)
order = ['UniProt_ID', 'Protein_ID', 'Lig_Name', 'Chain', 'Resn', 'Smiles_String',
         'Affinity_Measure', 'Relation', 'Affinity_Value', 'Affinity_Unit']
confirmed_uniprots = confirmed_uniprots[order]

# Remove known buffer components and cofactors, then select only ligands with a valid 3 letters HET
confirmed_uniprots = confirmed_uniprots[~confirmed_uniprots['Lig_Name'].isin(cof['code'])]
confirmed_uniprots = confirmed_uniprots[~confirmed_uniprots['Lig_Name'].isin(cbc['code'])]
confirmed_uniprots = confirmed_uniprots[confirmed_uniprots['Lig_Name'].str.len() == 3].reset_index(drop=True)

print(confirmed_uniprots.info(), flush=True)

# 1st retrieve binding affinities from ChEMBL 34 usign our database
uniprots = confirmed_uniprots.UniProt_ID.drop_duplicates().to_list()

# Connect to the database
print('Starting ChEMBL binding affinities search', flush=True)

load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), '../../.env'))

db = pymysql.connect(
    host=os.getenv("DB_HOST"),
    user=os.getenv("DB_USER"),
    password=os.getenv("DB_PASSWORD"),
    database=os.getenv("DB_NAME")



cursor = db.cursor(cursor=pymysql.cursors.DictCursor)

activities_df = pd.DataFrame(columns=['accession', 'molecule_chembl_id', 'standard_type', 'standard_relation',
                                      'standard_value', 'standard_units', 'standard_inchi_key', 'canonical_smiles'])

uni_w_affinities_in_chembl = {}
for i, uni in enumerate(uniprots):
    print(f'{i+1}/{len(uniprots)}', flush=True)
    params = (uni,)
    sql = """
            SELECT
                chembl_34.component_sequences.accession  ,
                chembl_34.target_dictionary.chembl_id target_chembl_id , 
                chembl_34.assays.chembl_id assay_chembl_id,
                chembl_34.assays.variant_id,
                chembl_34.molecule_dictionary.chembl_id molecule_chembl_id,
                chembl_34.molecule_dictionary.chirality,
                chembl_34.activities.activity_id,
                chembl_34.activities.standard_type,
                chembl_34.activities.standard_relation,
                chembl_34.activities.standard_value,
                chembl_34.activities.pchembl_value,
                chembl_34.activities.standard_units,
                chembl_34.compound_structures.standard_inchi_key,
                chembl_34.compound_structures.canonical_smiles

            FROM 
                chembl_34.component_sequences

            JOIN 
                chembl_34.target_components ON
                chembl_34.component_sequences.component_id = chembl_34.target_components.component_id
            JOIN 
                chembl_34.target_dictionary ON
                chembl_34.target_dictionary.tid = chembl_34.target_components.tid
            JOIN 
                chembl_34.assays ON
                chembl_34.target_dictionary.tid = chembl_34.assays.tid
            JOIN 
                chembl_34.activities ON
                chembl_34.activities.assay_id = chembl_34.assays.assay_id
            JOIN 
                chembl_34.molecule_dictionary ON
                chembl_34.molecule_dictionary.molregno = chembl_34.activities.molregno
            JOIN 
                chembl_34.molecule_hierarchy ON
                chembl_34.molecule_dictionary.molregno = chembl_34.molecule_hierarchy.molregno
            JOIN
                chembl_34.compound_structures ON
                chembl_34.molecule_hierarchy.active_molregno = chembl_34.compound_structures.molregno

            WHERE 
                chembl_34.component_sequences.accession = %s AND
                chembl_34.assays.confidence_score >= 8 AND
                ((chembl_34.activities.data_validity_comment is null) OR
                 (chembl_34.activities.data_validity_comment = 'Manually validated')) AND
                chembl_34.molecule_dictionary.chirality != 0 AND
                chembl_34.activities.pchembl_value is not null AND
                ((chembl_34.activities.standard_type = 'IC50') OR 
                 (chembl_34.activities.standard_type = 'Kd') OR 
                 (chembl_34.activities.standard_type = 'Ki'));
          """

    cursor.execute(sql, params)
    rows = cursor.fetchall()

    temp_df = pd.DataFrame(rows)
    #print(f'{uni} has: {len(temp_df)} activities')
    if len(temp_df) != 0:
        uni_w_affinities_in_chembl[uni] = len(temp_df)
    activities_df = activities_df.append(temp_df, ignore_index=True)
print(f'{len(activities_df)} activities have been found', flush=True)
activities_df.to_csv('../data/ChEMBL_activities_for bindingMOAD.csv', index=False)
print('ChEMBL search done', flush=True)

print('Starting PDB HET code ChEMBL ID mapping', flush=True)
# Drop ligands without SMILES and then map the PDB HEt code to the ChEMBL ID (if available)
confirmed_uniprots = confirmed_uniprots.dropna(subset=['Smiles_String']).reset_index(drop=True)
failed = []
for i, lig in confirmed_uniprots[['Lig_Name']].itertuples():
#     if i < 10:
    print(f'{i+1}/{len(confirmed_uniprots)}', flush=True)
    try:
        response = requests.get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{lig}')
        results = response.json()
        for res in results['rcsb_chem_comp_related']:
            if res['resource_name'] == 'ChEMBL':
                confirmed_uniprots.at[i, 'CHEMBL_ID'] = res['resource_accession_code']
    except:
        print(f'ChEMBL ID not found for {lig}', flush=True)
        failed.append(lig)

# Now merge BindingMOAD and ChEMBL using the ChEMBL ID
chembl_no_duplicates = activities_df.drop_duplicates(subset=['accession', 'molecule_chembl_id'], keep='first').reset_index(drop=True)
no_dup_renamed = chembl_no_duplicates.rename(columns={'accession': 'UniProt_ID', 
                                                      'molecule_chembl_id': 'CHEMBL_ID', 
                                                      'standard_type': 'Affinity_Measure', 
                                                      'standard_relation': 'Relation',
                                                      'standard_value': 'Affinity_Value',
                                                      'standard_units': 'Affinity_Unit'})
                                                      
merged = pd.merge(confirmed_uniprots, no_dup_renamed, on=['UniProt_ID', 'CHEMBL_ID'], how='left')
cols_to_fill = ['Affinity_Measure','Relation', 'Affinity_Value', 'Affinity_Unit']

for col in cols_to_fill:
    merged[f'{col}_x'] = merged[f'{col}_x'].fillna(merged[f'{col}_y'])
    merged.drop(columns=[f'{col}_y'], inplace=True)

    
merged.rename(columns={col+'_x': col for col in cols_to_fill}, inplace=True)
merged.drop(columns=['canonical_smiles', 'variant_id', 'chirality', 'activity_id', 'pchembl_value'], inplace=True)

# Now drop all the ligands without any binding affinity, this is the final BindingMOAD dataframe (more filtering will be done later)
selected = merged.dropna(subset=['Affinity_Measure']).reset_index(drop=True)
selected = selected.sort_values(by=['UniProt_ID', 'Protein_ID', 'Lig_Name']).reset_index(drop=True)

# Save BindingMOAD
selected.to_csv('../data/clean_BindingMOAD.csv', index=False)

# Now download PDB files form RCSB web services
print('DOWNLOADING PDB FILES...', flush=True)
for i, uni, pdb in selected[['UniProt_ID', 'Protein_ID']].itertuples():

    os.makedirs('../pdb_files_not_aligned', exist_ok=True)
    ur(f'https://files.rcsb.org/download/{pdb}.pdb', filename=f'../pdb_files_not_aligned/{uni}/{pdb}.pdb')
print('PDB download done successfully', flush=True)







