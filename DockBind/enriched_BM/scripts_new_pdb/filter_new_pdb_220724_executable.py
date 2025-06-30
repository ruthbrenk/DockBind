import os
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from tqdm import tqdm
from urllib.request import urlretrieve as ur
import requests
from pymol import cmd

files = [i for i in os.listdir('pdb_csv_files') if 'csv' in i]
df = pd.DataFrame()
for i, file in enumerate(files):
    temp = pd.read_csv(f'pdb_csv_files/{file}', usecols=['PDB ID', 'Refinement Resolution (Å)', 'Ligand ID', 'Ligand SMILES'])
    df = df.append(temp).reset_index(drop=True)
    
df.info()

cof = pd.read_csv('../known_cofactors.csv')
buf = pd.read_csv('../crystal_buffer_components.csv')
ama = pd.read_csv('../amino_acids.csv')

df['PDB ID'] = df['PDB ID'].fillna(method='ffill')
df['Refinement Resolution (Å)'] = df['Refinement Resolution (Å)'].fillna(method='ffill')

df = df.dropna(subset=['Ligand SMILES', 'Ligand ID'])
df = df[~df['Ligand ID'].isin(buf.code.unique())]
df = df[~df['Ligand ID'].isin(cof.code.unique())]
df = df[~df['Ligand ID'].isin(ama.HET.unique())]
df = df[df['Ligand ID'].str.len() == 3].reset_index(drop=True)

uniprot = pd.read_csv(f'idmapping_2024_07_02.tsv', sep='\t')
uniprot = uniprot.drop_duplicates(subset='From').reset_index(drop=True)
pdb_uni = dict(zip(uniprot.From, uniprot.Entry))
for i, u in df[['PDB ID']].itertuples():
    try:
        df.at[i, 'UniProt_ID'] = pdb_uni[u]
    except Exception as e:
        print(e)
        df.at[i, 'UniProt_ID'] = 'Not found'
        
df = df[df['UniProt_ID'] != 'Not found'].reset_index(drop=True)
df.to_csv('new_pdb_since_2020.csv', index=False)

bdb_full = pd.read_csv('../../../binding_moad_2020_feb2024/data/BindingDB_All_202403.tsv', sep='\t', usecols=['BindingDB Reactant_set_id',
 'Ligand SMILES',
 'BindingDB MonomerID',
 'Ki (nM)',
 'IC50 (nM)',
 'Kd (nM)',                
 'Curation/DataSource',
 'Ligand HET ID in PDB',
 'UniProt (SwissProt) Primary ID of Target Chain'])

bdb_full = bdb_full.dropna(subset=['Ki (nM)', 'Kd (nM)', 'IC50 (nM)'], how='all').reset_index(drop=True)

uniprots = df.UniProt_ID.unique()
bdb_uni = bdb_full['UniProt (SwissProt) Primary ID of Target Chain'].unique()

pdb_in_bdb = df[df['UniProt_ID'].isin(bdb_uni)].reset_index(drop=True)
bdb_in_pdb = bdb_full[bdb_full['UniProt (SwissProt) Primary ID of Target Chain'].isin(uniprots)].reset_index(drop=True)

bad_id = []
for i, ch_id, smi in bdb_in_pdb[['BindingDB MonomerID', 'Ligand SMILES']].itertuples():
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad_id.append(smi)
            
    except:
        continue
        
bdb_in_pdb = bdb_in_pdb[~bdb_in_pdb['Ligand SMILES'].isin(bad_id)]

bad_id = []
for i, ch_id, smi in pdb_in_bdb[['Ligand ID', 'Ligand SMILES']].itertuples():
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad_id.append(smi)
            
    except:
        continue
        
pdb_in_bdb = pdb_in_bdb[~pdb_in_bdb['Ligand SMILES'].isin(bad_id)]

bindingmoad = pd.read_csv('../BindingMOAD_selected_templates.csv')
bm_uni = bindingmoad.UniProt_ID.unique()

new_pdb = pdb_in_bdb[pdb_in_bdb['UniProt_ID'].isin(bm_uni)]
new_pdb_unique = new_pdb['PDB ID'].unique()

new_unis = pdb_in_bdb[~pdb_in_bdb['UniProt_ID'].isin(bm_uni)]
new_uniprots = new_unis.UniProt_ID.unique()

print(f'There are {len(new_uniprots)} new targets')
print(f'There are {len(new_pdb_unique)} new pdb for old targets')

# Now look for similar ligands in BindingDB, to find the binding site
found_activities = []
for job, uni in enumerate(new_uniprots):
    print(f'{job+1}/{len(new_uniprots)}')
    
    temp_pdb = pdb_in_bdb[pdb_in_bdb.UniProt_ID == uni]
    temp_bdb = bdb_in_pdb[bdb_in_pdb['UniProt (SwissProt) Primary ID of Target Chain'] == uni]
    temp_bdb = temp_bdb.drop_duplicates(subset=['BindingDB MonomerID'], keep='first')
    
    print(f'PDB ligands: {len(temp_pdb)}\nBindingDB ligands: {len(temp_bdb)}')
    for i, pdb_id, lig_id, lig_smi in tqdm(temp_pdb[['PDB ID', 'Ligand ID', 'Ligand SMILES']].itertuples(), total=len(temp_pdb)):
        
        ref_mol = Chem.MolFromSmiles(lig_smi)
        
        ref_fps = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, 2048, useFeatures=True)

        for j, ch_id, ch_smi in temp_bdb[['BindingDB MonomerID', 'Ligand SMILES']].itertuples():
            ch_mol = Chem.MolFromSmiles(ch_smi)
            ch_fps = AllChem.GetMorganFingerprintAsBitVect(ch_mol, 2, 2048, useFeatures=True)

            ts = round(DataStructs.TanimotoSimilarity(ref_fps, ch_fps), 2)

            if ts >= 0.7:
                found_activities.append({'UniProt_ID': uni,
                             'PDB_ID': pdb_id,
                             'Ligand_ID': lig_id,
                             'Ligand_SMILES': lig_smi,
                             'BindingDB_ID': ch_id,
                             'BindingDB_Smiles': ch_smi
                            })
                
activities = pd.DataFrame(found_activities)

new_targets_pdb = activities.PDB_ID.unique()
old_targets_pdb = new_pdb['PDB ID'].unique()

pdb_not_found = []

TEMP = f'pdb_files_not_aligned'
print('DOWNLOAD NEW TARGETS PDB FILES')
for i, pdb in enumerate(new_targets_pdb):
    print(f'Staritng {i+1}/{len(new_targets_pdb)}')
    try:
        uni = pdb_uni[pdb]
        os.makedirs(f'{TEMP}/{uni}', exist_ok=True)
        ur(f'https://files.rcsb.org/download/{pdb}.pdb', filename=f'{TEMP}/{uni}/{pdb}.pdb')
    except:
        pdb_not_found.append(pdb)
        
print('DOWNLOAD OLD TARGETS PDB FILES')
for i, pdb in enumerate(old_targets_pdb):
    print(f'Staritng {i+1}/{len(old_targets_pdb)}')
    try:
        uni = pdb_uni[pdb]
        os.makedirs(f'{TEMP}/{uni}', exist_ok=True)
        ur(f'https://files.rcsb.org/download/{pdb}.pdb', filename=f'{TEMP}/{uni}/{pdb}.pdb')
    except:
        pdb_not_found.append(pdb)
        

pdb_not_found = []
for i, uni, pdb in activities[['UniProt_ID', 'PDB_ID']].itertuples():
    if not os.path.exists(f'pdb_files_not_aligned/{uni}/{pdb}.pdb'):
        pdb_not_found.append(pdb)
        

for i, uni, pdb in new_pdb[['UniProt_ID', 'PDB ID']].itertuples():
    if not os.path.exists(f'pdb_files_not_aligned/{uni}/{pdb}.pdb'):
        pdb_not_found.append(pdb)

print(f'PDB not found: {len(pdb_not_found)}')

confirmed_new_targets = activities[~activities['PDB_ID'].isin(pdb_not_found)].reset_index(drop=True)
confirmed_new_pdbs = new_pdb[~new_pdb['PDB ID'].isin(pdb_not_found)].reset_index(drop=True)

confirmed_new_targets.to_csv('tmp/checkpoint_new_targets.csv', index=False)
confirmed_new_pdbs.to_csv('tmp/checkpoint_new_pdbs.csv', index=False)

BASE = '/Home/siv33/ped023/Desktop/new_bm_june24'
OUTPUT = f'{BASE}/data/new_pdb'
INPUT = f'{BASE}/data/new_pdb/pdb_files_not_aligned'
REFERENCE = f'{BASE}/pdb_files_aligned'
REFERENCE_LIGANDS = f'{BASE}/reference_ligands'
URL = 'https://files.rcsb.org/ligands/download'
MODELS = f'{OUTPUT}/models'
TMP = f'{OUTPUT}/tmp'

cmd.reinitialize()

alignment_out = []

new_pdb = pd.read_csv('tmp/checkpoint_new_pdbs.csv')
old_ref = pd.read_csv('../../logs/alignment.csv')

pdb_uni = dict(zip(old_ref.UniProt_ID, old_ref.ref_pdb))

uni_to_new = []

uniprots = new_pdb.UniProt_ID.unique()

for i, uni in enumerate(uniprots):
    print(f'{i+1}/{len(uniprots)}: {uni}')
    pdbs = new_pdb[new_pdb['UniProt_ID'] == uni]['PDB ID'].unique()
    print(f'{len(pdbs)} PDB files to process')
    ref_pdb = pdb_uni[uni]
    files = os.listdir(f'{REFERENCE_LIGANDS}/{uni}')

    found = False
    if files:
        for i in files:
            if ref_pdb in i:
                ref_lig = i
                found = True
                break
    if not found:
        ref_lig = files[0]
        
    for pdb in pdbs:
        print(pdb)

        cmd.load(f'{REFERENCE}/{uni}/{ref_pdb}.pdb', object='reference')
        cmd.load(f'{REFERENCE_LIGANDS}/{uni}/{ref_lig}', 'ref_lig')
        cmd.extract(selection='bc. reference within 2 of ref_lig', name='ref')
        cmd.load(f'{INPUT}/{uni}/{pdb}.pdb', object='new')
            
        rmsd = cmd.align('new', 'ref')
        
        os.makedirs(f'{OUTPUT}/pdb_files_aligned/{uni}', exist_ok=True)
        cmd.save(filename=f'{OUTPUT}/pdb_files_aligned/{uni}/{pdb}.pdb', selection='new')
        cmd.reinitialize()
        
        alignment_out.append({'UniProt_ID': uni,
                              'reference_pdb': ref_pdb,
                              'PDB': pdb,
                              'rmsd': rmsd[0]
                             })     


df_no_dups = confirmed_new_targets.drop_duplicates(subset=['UniProt_ID', 'PDB_ID', 'Ligand_ID']).reset_index(drop=True)

cmd.reinitialize()

MODELS = 'models'
ligands_in_pdb = pd.DataFrame()

for i, uni, pdb, smi, name in df_no_dups[['UniProt_ID', 'PDB_ID', 'Ligand_SMILES', 'Ligand_ID']].itertuples():
    print(f'{i+1}/{len(df_no_dups)}')

    os.makedirs(f'tmp/lig_copies/{uni}', exist_ok=True)
    cmd.load(filename=f'pdb_files_not_aligned/{uni}/{pdb}.pdb', object=pdb)
    cmd.extract(selection=f'r. {name}', name='ligs')
    cmd.save(filename=f'tmp/lig_copies/{uni}/{pdb}_{name}.pdb', selection='ligs')
    cmd.reinitialize()
    ligands_in_pdb.at[i, 'UniProt_ID'] = uni
    ligands_in_pdb.at[i, 'PDB'] = pdb
    ligands_in_pdb.at[i, 'Lig'] = name

    atm = []
    try:
        with open(f'tmp/lig_copies/{uni}/{pdb}_{name}.pdb', 'r') as f:
            for line in f:
                if 'HETATM' in line:
                    if line[16] != ' ':
                        atm.append(line[:16] + ' ' + line[17:])
                    else:
                        atm.append(line)
        keep = ''.join(atm)
        mols = Chem.MolFromPDBBlock(keep, sanitize=False)
        mols = Chem.GetMolFrags(mols, asMols=True, sanitizeFrags=False)

        ligands_in_pdb.at[i, 'Lig_copies'] = len(mols)
    except Exception as e:
#         print(e)
        ligands_in_pdb.at[i, 'Lig_copies'] = -1

ligands_in_pdb.to_csv('tmp/ligands_in_pdb.csv', index=False)

OUT = 'pdb_files_aligned'
cmd.reinitialize()

no_single_ligand = []
for uni in ligands_in_pdb.UniProt_ID.unique():
    
    os.makedirs(f'{OUT}/{uni}', exist_ok=True)
    
    temp_df = ligands_in_pdb[ligands_in_pdb.UniProt_ID == uni].reset_index(drop=True)
    if len(temp_df[temp_df.Lig_copies == 1]) != 0:
        ref = temp_df[temp_df.Lig_copies == 1].iloc[0]
        ref_pdb = ref.PDB
        ref_lig = ref.Lig
        
        for i, pdb, lig in temp_df[['PDB', 'Lig']].itertuples():
            
            cmd.load(f'pdb_files_not_aligned/{uni}/{ref_pdb}.pdb', object='reference')
            if pdb == ref_pdb:
                cmd.save(f'{OUT}/{uni}/{ref_pdb}.pdb', selection='reference')
                alignment_out.append({'UniProt_ID': uni,
                                      'reference_pdb': ref_pdb,
                                      'PDB': pdb,
                                      'rmsd': 0
                                     })
            else:
                cmd.load(f'pdb_files_not_aligned/{uni}/{pdb}.pdb', object=pdb)
                cmd.extract(selection=f'bc. reference within 2 of resn {ref_lig}', name='ref_pdb')
                rmsd = cmd.align(pdb, 'ref_pdb')
                cmd.save(f'{OUT}/{uni}/{pdb}.pdb', selection=pdb)
                cmd.reinitialize()
                
                alignment_out.append({'UniProt_ID': uni,
                                      'reference_pdb': ref_pdb,
                                      'PDB': pdb,
                                      'rmsd': rmsd[0]
                                     })
                
                       
    else:
        no_single_ligand.append(uni)


cmd.reinitialize()

for uni in no_single_ligand:
    
    os.makedirs(f'{OUT}/uni', exist_ok=True)
    
    temp_df = ligands_in_pdb[ligands_in_pdb.UniProt_ID == uni].reset_index(drop=True)
    if len(temp_df[temp_df.Lig_copies == 1]) == 0: 
        print(f'{uni}: {len(temp_df)} PDB to process')
        ref_pdb = temp_df.loc[temp_df.Lig_copies.idxmin(), 'PDB']
        ref_lig = temp_df.loc[temp_df.Lig_copies.idxmin(), 'Lig']
        
        with open(f'tmp/lig_copies/{uni}/{ref_pdb}_{ref_lig}.pdb', 'r') as f:
            chain = f.readlines()[0].strip()[21]
      
        for i, pdb, lig in temp_df[['PDB', 'Lig']].itertuples():
            
            cmd.load(f'pdb_files_not_aligned/{uni}/{ref_pdb}.pdb', object='reference')
            if pdb == ref_pdb:
                cmd.save(f'{OUT}/{uni}/{ref_pdb}.pdb', selection='reference')
                alignment_out.append({'UniProt_ID': uni,
                                      'reference_pdb': ref_pdb,
                                      'PDB': pdb,
                                      'rmsd': 0
                                     })
            else:
                cmd.load(f'pdb_files_not_aligned/{uni}/{pdb}.pdb', object=pdb)
                cmd.extract(selection=f'bc. reference within 2 of c. {chain} and r. {ref_lig}', name='ref_pdb')
                rmsd = cmd.align(pdb, 'ref_pdb')
                
                cmd.save(f'{OUT}/{uni}/{pdb}.pdb', selection=pdb)
                cmd.reinitialize()
                alignment_out.append({'UniProt_ID': uni,
                                      'reference_pdb': ref_pdb,
                                      'PDB': pdb,
                                      'rmsd': rmsd[0]
                                     })



alignment_df = pd.DataFrame(alignment_out)
alignment_df = alignment_df[alignment_df.rmsd <= 2.5]

new_uni_aligned = activities[activities['PDB_ID'].isin(alignment_df.PDB.unique())]
new_uni_aligned = new_uni_aligned.drop_duplicates(subset=['PDB_ID', 'Ligand_ID']).reset_index(drop=True)

new_pdb_aligned = new_pdb[new_pdb['PDB ID'].isin(alignment_df.PDB.unique())]
print(len(new_uni_aligned))
print(len(new_pdb_aligned))

pdb_lig_new = dict(zip(new_uni_aligned.PDB_ID, new_uni_aligned.Ligand_ID))

lig_smiles = dict(zip(new_pdb_aligned['Ligand ID'], new_pdb_aligned['Ligand SMILES']))


new_df = []

for i, uni in enumerate(new_pdb_aligned.UniProt_ID.unique()):
    
    print(f'{i+1}/{len(new_pdb_aligned.UniProt_ID.unique())}: {uni}')
    pdbs = new_pdb_aligned[new_pdb_aligned['UniProt_ID'] == uni]['PDB ID'].unique()
    print(f'{len(pdbs)} PDB files to process')
    ref_pdb = pdb_uni[uni]
    files = os.listdir(f'{REFERENCE_LIGANDS}/{uni}')
    found = False
    if files:
        for i in files:
            if ref_pdb in i:
                ref_lig = i
                found = True
                break
    if not found:
        ref_lig = files[0]    
    for pdb in tqdm(pdbs):
        cmd.load(filename=f'{REFERENCE_LIGANDS}/{uni}/{ref_lig}', object='ref_lig')
        cmd.load(filename=f'pdb_files_aligned/{uni}/{pdb}.pdb', object='new')
        
        cmd.extract(selection='br. new within 5 of ref_lig and not ino.', name='lig')
        cmd.remove('resn HOH')
        cmd.remove('resn NA+K+CA+MG+CL+F+BR+I+SO4+PO4+FE+ZN+CU+MN+CO+NI+CR+MO+V+CD+HG+PB+PT+AU')
        cmd.remove('resn AG+LA+CE+PR+ND+SM+EU+GD+TB+DY+HO+ER+TM+YB+LU+TH+U+NO3+CLO4+ACT+CO3')

        os.makedirs(f'{OUTPUT}/reference_ligands/{uni}', exist_ok=True)
        os.makedirs(f'{TMP}/{uni}', exist_ok=True)
        
        cmd.save(f'{TMP}/{uni}/temp_{pdb}.pdb', selection='lig')
        cmd.reinitialize()
        
        atm = []
        with open(f'{TMP}/{uni}/temp_{pdb}.pdb', 'r') as f:
            for line in f:
                if 'HET' in line:
                    atm.append(line)
        if len(atm) == 0:
            continue

        mol_name = atm[0][17:20]

        try:
            keep = ''.join(atm)
            mols = Chem.MolFromPDBBlock(keep, sanitize=False)
            mols = Chem.GetMolFrags(mols, asMols=True, sanitizeFrags=False)[0]
            mols.SetProp('_Name', mol_name)

            writer = Chem.SDWriter(f'{OUTPUT}/reference_ligands/{uni}/{pdb}_{mol_name}.sdf')
            writer.write(mols)
            writer.close()

            new_df.append({'UniProt_ID': uni,
                           'PDB_ID': pdb,
                           'Lig_Name': mol_name,
                           'Lig_Smiles': lig_smiles[mol_name]}
                         )
        except:
            print(f'{pdb} failed')



uniprots = new_uni_aligned.UniProt_ID.unique()

cmd.reinitialize()

for i, uni in enumerate(uniprots):
    print(f'Starting {i+1}/{len(uniprots)}')
    temp_df = new_uni_aligned[new_uni_aligned.UniProt_ID == uni]
    
    ref_pdb = alignment_df[alignment_df.UniProt_ID == uni].iloc[0].reference_pdb
    ref_lig = pdb_lig_new[ref_pdb]
    
    if ligands_in_pdb[ligands_in_pdb.PDB == ref_pdb].Lig_copies.iloc[0] == 1:
            
        for i, pdb, name, smi in temp_df[['PDB_ID', 'Ligand_ID', 'Ligand_SMILES']].itertuples():
            print(pdb)
            if pdb == ref_pdb:
                cmd.load(filename=f'pdb_files_aligned/{uni}/{ref_pdb}.pdb', object='reference')
                cmd.extract(selection=f'resn {ref_lig}.pdb', name='ref_lig')

                os.makedirs(f'{OUTPUT}/reference_ligands/{uni}', exist_ok=True)
                os.makedirs(f'{TMP}/{uni}', exist_ok=True)
        
                cmd.save(f'{TMP}/{uni}/temp_{ref_pdb}.pdb', selection='ref_lig')
                cmd.reinitialize()
                
            else:
                cmd.load(filename=f'pdb_files_aligned/{uni}/{ref_pdb}.pdb', object='reference')
                cmd.load(filename=f'pdb_files_aligned/{uni}/{pdb}.pdb', object='new')

                cmd.extract(selection=f'resn {ref_lig}', name='ref_lig')
                cmd.extract(selection=f'(br. new within 4 of ref_lig) and r. {name}', name='new_lig')
                
                os.makedirs(f'{OUTPUT}/reference_ligands/{uni}', exist_ok=True)
                os.makedirs(f'{TMP}/{uni}', exist_ok=True)
        
                cmd.save(f'{TMP}/{uni}/temp_{pdb}.pdb', selection='new_lig')
                cmd.reinitialize()                
            
    else:    
        with open(f'tmp/lig_copies/{uni}/{ref_pdb}_{ref_lig}.pdb', 'r') as f:
            chain = f.readlines()[0].strip()[21]
            
        for i, pdb, name, smi in temp_df[['PDB_ID', 'Ligand_ID', 'Ligand_SMILES']].itertuples():
            print(pdb)
            if pdb == ref_pdb:
                cmd.load(filename=f'pdb_files_aligned/{uni}/{ref_pdb}.pdb', object='reference')
                cmd.extract(selection=f'resn {ref_lig} and c. {chain}', name='ref_lig')

                os.makedirs(f'{OUTPUT}/reference_ligands/{uni}', exist_ok=True)
                os.makedirs(f'{TMP}/{uni}', exist_ok=True)
        
                cmd.save(f'{TMP}/{uni}/temp_{ref_pdb}.pdb', selection='ref_lig')
                cmd.reinitialize()
                
            else:
                cmd.load(filename=f'pdb_files_aligned/{uni}/{ref_pdb}.pdb', object='reference')
                cmd.load(filename=f'pdb_files_aligned/{uni}/{pdb}.pdb', object='new')

                cmd.extract(selection=f'resn {ref_lig}', name='ref_lig')
                cmd.extract(selection=f'(br. new within 4 of ref_lig) and r. {name}', name='new_lig')
                
                os.makedirs(f'{OUTPUT}/reference_ligands/{uni}', exist_ok=True)
                os.makedirs(f'{TMP}/{uni}', exist_ok=True)
        
                cmd.save(f'{TMP}/{uni}/temp_{pdb}.pdb', selection='new_lig')
                cmd.reinitialize()
          
    atm = []
    with open(f'{TMP}/{uni}/temp_{pdb}.pdb', 'r') as f:
        for line in f:
            if name in line:
                atm.append(line)
    if len(atm) == 0:
        continue

    try:
        keep = ''.join(atm)
        mols = Chem.MolFromPDBBlock(keep, sanitize=False)
        mols = Chem.GetMolFrags(mols, asMols=True, sanitizeFrags=False)[0]

        writer = Chem.SDWriter(f'{OUTPUT}/reference_ligands/{uni}/{pdb}_{name}.sdf')
        writer.write(mols)
        writer.close()

        new_df.append({'UniProt_ID': uni,
                       'PDB_ID': pdb,
                       'Lig_Name': name,
                       'Lig_Smiles': lig_smiles[name]}
                     )
    except:
        print(f'{pdb} failed')


df_3 = pd.DataFrame(new_df)
df_3.to_csv('new_pdb_confirmed.csv', index=False)
