import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os
from pymol import cmd
import requests
import math
import statistics
from rdkit import Chem
from rdkit.Chem import AllChem

BASE = '..'
DATA = f'{BASE}/data'
PDB = f'{BASE}/pdb_files_not_aligned'
ALIGNED = f'{BASE}/pdb_files_aligned'
MODELS = f'{BASE}/model_ligands'
LIGANDS_OUT = f'{BASE}/reference_ligands'
URL = 'https://files.rcsb.org/ligands/download'
LOGS = f'{BASE}/logs'
TEMP = f'{BASE}/tmp/PDB_alignment_and_ligand_extraction'

log = open(f'{LOGS}/pdb_alignment_and_lig_extraction.log', 'w')

def convert_value(value, unit):
    if unit == 'M':
        return value*(10**9)
    elif unit in ('dM', 'M^-1'):
        return value*(10**8)
    elif unit == 'cM':
        return value*(10**7)
    elif unit == 'mM':
        return value*(10**6)
    elif unit == 'uM':
        return value*(10**3)
    elif unit == 'nM':
        return value
    elif unit == 'pM':
        return value*(10**-3)
    elif unit == 'fM':
        return value*(10**-6)

# Load clean BindingMOAD
df = pd.read_csv('../data/clean_BindingMOAD.csv')

# Drop missing affinity units, select only = relation and convert binding affinities to nM
df = df[df['Relation'] == '=']
df = df[df['Affinity_Measure'].isin(['IC50', 'ic50', 'Kd', 'Ki'])].reset_index(drop=True)
df = df.drop_duplicates(subset=['UniProt_ID', 'Protein_ID'], keep=False).reset_index(drop=True)

for i, measure, value, unit in df[['Affinity_Measure', 'Affinity_Value', 'Affinity_Unit']].itertuples():
    value = convert_value(value, unit)
    df.loc[i, measure] = value

df['IC50'] = df['IC50'].fillna(df['ic50'])
df = df.drop(columns=['ic50', 'Affinity_Measure', 'Affinity_Value', 'Affinity_Unit', 'Relation'])

# Now download ligand model from RCSB resources, this will be used for stereochemistry correction
model_not_available = []
log.write('Starting models download')

for i, uni, pdb, lig in df[['UniProt_ID','Protein_ID', 'Lig_Name']].itertuples():
    print(f'{i+1}/{len(df)}')
    os.makedirs(f'{MODELS}/{uni}', exist_ok=True)
    try:
        response = requests.get(f'{URL}/{lig}_ideal.sdf')
        open(f'{MODELS}/{uni}/{lig}.sdf', 'wb').write(response.content)

        # Use the downloaded SDF to get a clean SMILES with explicit stereochemistry
        model_mol = Chem.SDMolSupplier(f'{MODELS}/{uni}/{lig}.sdf')[0]
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(model_mol)
    except Exception as e:
        log.write(f'{lig}\n')
        print(f'{lig} {e}')
        model_not_available.append(f'{lig}')
        
# Remove ligands without valid model
df = df[~df['Lig_Name'].isin(model_not_available)].reset_index(drop=True)

log.write('Starting ligands check')
# Now check how many ligands are present in each PDB file
ligands_in_pdb =pd.DataFrame()
for i, uni, pdb, lig, chain, resn in df[['UniProt_ID', 'Protein_ID', 'Lig_Name', 'Chain', 'Resn']].itertuples():
    print(f'{i+1}/{len(df)}')
    os.makedirs(f'{TEMP}/lig_copies/{uni}', exist_ok=True)
    cmd.load(filename=f'{PDB}/{uni}/{pdb}.pdb', object=pdb)
    cmd.extract(selection=f'r. {lig}', name='ligs')
    cmd.save(filename=f'{TEMP}/lig_copies/{uni}/{pdb}_{lig}.pdb', selection='ligs')
    cmd.reinitialize()
    ligands_in_pdb.at[i, 'UniProt_ID'] = uni
    ligands_in_pdb.at[i, 'PDB'] = pdb
    ligands_in_pdb.at[i, 'Lig'] = lig
    
    atm = []
    try:
        with open(f'{TEMP}/lig_copies/{uni}/{pdb}_{lig}.pdb', 'r') as f:
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

ligands_in_pdb.to_csv(f'{DATA}/number_of_ligands_in_pdb_file.csv', index=None)

# Remove problematic ligands (Lig_copies = -1)
ligands_to_remove = ligands_in_pdb[ligands_in_pdb['Lig_copies'] == -1]['PDB'].to_list()
df = df[~df['Protein_ID'].isin(ligands_to_remove)].reset_index(drop=True)

ligands_in_pdb_clean = ligands_in_pdb[ligands_in_pdb['Lig_copies'] != -1].reset_index(drop=True)

# Align PDB files
# The idea here is to use PDB files with only one ligand copy as reference
# If it is not possible use another PDB file and takes one chain as reference

uniprots = df.UniProt_ID.drop_duplicates().to_list()
rmsds = open(f'{LOGS}/alignment.csv', 'w')
rmsds.write('UniProt_ID,ref_pdb,aligned_pdb,rmsd\n')
for num, uni in enumerate(uniprots):
    os.makedirs(f'{ALIGNED}/{uni}', exist_ok=True)
    temp_df = df[df['UniProt_ID'] == uni].reset_index(drop=True)
    print(f'Starting {num+1}/{len(uniprots)}: {len(temp_df)} PDB')
    # Get reference PDB file, first PDB with only one ligand
    if len(ligands_in_pdb_clean[(ligands_in_pdb_clean['UniProt_ID'] == uni) &
                                (ligands_in_pdb_clean['Lig_copies'] == 1)] != 0):

        ref_df = ligands_in_pdb_clean[(ligands_in_pdb_clean['UniProt_ID'] == uni) & 
                                      (ligands_in_pdb_clean['Lig_copies'] == 1)].iloc[0]
        ref_pdb = ref_df.PDB
        ref_chain = temp_df[temp_df['Protein_ID'] == ref_pdb].iloc[0].Chain

        for i, pdb in temp_df[['Protein_ID']].itertuples():
            cmd.load(filename=f'{PDB}/{uni}/{ref_pdb}.pdb', object=ref_pdb)
            if pdb == ref_pdb:
                cmd.save(filename=f'{ALIGNED}/{uni}/{pdb}.pdb', selection=ref_pdb)
                cmd.reinitialize()
            else:
                cmd.extract(selection=f'chain {ref_chain} and {ref_pdb}', name='selected')
                cmd.load(filename=f'{PDB}/{uni}/{pdb}.pdb', object=pdb)
                rmsd = cmd.align(pdb, 'selected')
                cmd.save(filename=f'{ALIGNED}/{uni}/{pdb}.pdb', selection=pdb)
                rmsds.write(f'{uni},{ref_pdb},{pdb},{rmsd[0]}\n')
                cmd.reinitialize()
#                 cmd.load(filename=f'{PDB}/{uni}/{pdb}.pdb', object=pdb)
#                 rmsd = cmd.align(pdb, ref_pdb)
#                 cmd.save(filename=f'{ALIGNED}/{uni}/{pdb}.pdb', selection=pdb)
#                 rmsds.write(f'{uni},{ref_pdb},{pdb},{rmsd[0]}\n')
#                 cmd.reinitialize()
    else:
        ref_df = ligands_in_pdb_clean[(ligands_in_pdb_clean['UniProt_ID'] == uni)]
        ref_pdb = ref_df[ref_df['Lig_copies'] == ref_df['Lig_copies'].min()].iloc[0].PDB
        ref_chain = temp_df[temp_df['Protein_ID'] == ref_pdb].iloc[0].Chain
        
        for i, pdb in temp_df[['Protein_ID']].itertuples():
            print(f'Starting PDB {i+1}/{len(temp_df)}')
            cmd.load(filename=f'{PDB}/{uni}/{ref_pdb}.pdb', object=ref_pdb)
            if pdb == ref_pdb:
                cmd.save(filename=f'{ALIGNED}/{uni}/{pdb}.pdb', selection=ref_pdb)
                cmd.reinitialize()
            else:
                cmd.extract(selection=f'chain {ref_chain} and {ref_pdb}', name='selected')
                cmd.load(filename=f'{PDB}/{uni}/{pdb}.pdb', object=pdb)
                rmsd = cmd.align(pdb, 'selected')
                cmd.save(filename=f'{ALIGNED}/{uni}/{pdb}.pdb', selection=pdb)
                rmsds.write(f'{uni},{ref_pdb},{pdb},{rmsd[0]}\n')
                cmd.reinitialize()
print('Alignment done, checking RMSD')            
rmsds.close()

alignment_results = pd.read_csv(f'{LOGS}/alignment.csv')


# Selecet ongly PDB with RMSD within 2.5 Å, to make sure the ligands are actually in the same place
confirmed_aligned = alignment_results[alignment_results['rmsd'] < 2.5].reset_index(drop=True)
confirmed_pdb = list(set(confirmed_aligned.ref_pdb.to_list() + confirmed_aligned.aligned_pdb.to_list()))
final_bm = df[df['Protein_ID'].isin(confirmed_pdb)].reset_index(drop=True)
final_bm = final_bm.rename(columns={'Protein_ID': 'PDB', 
                                    'UniProt_ID': 'UniProt_ID',
                                    'Lig_Name': 'Lig'})


# If there are PDB entries with multiple ligands with binding affinity, remove them.
# For example an enzyme that has a substrate and NADH, there might have ligands that compete for the substrate
# while other ligands compete for NADH. This is a problem for cross-dockign
x = len(final_bm)
final_bm = final_bm.drop_duplicates(subset=['UniProt_ID', 'PDB'], keep=False)
y = len(final_bm)
print(f'Done. There were {x-y} duplicated PDBs')

pdb_lig_dict = dict(zip(final_bm.PDB, final_bm.Lig))
reference_pdb = confirmed_aligned.drop_duplicates(subset=['UniProt_ID', 'ref_pdb']).reset_index(drop=True)
reference_pdb = dict(zip(reference_pdb.UniProt_ID, reference_pdb.ref_pdb))

# Finally extract ligands from the PDB files and save them as reference ligands
extract_log = open(f'{LOGS}/extract_ligand_log.tsv', 'w')
extract_log.write('UniProt_ID\tTotal_ligands\tLigands_not_found\tLigands\tFailed_ligands_to_sdf\n')
missing_ligands = open(f'{LOGS}/missing_ligands.csv', 'w')
uniprots = final_bm.UniProt_ID.drop_duplicates().to_list()
for uni in uniprots:
    os.makedirs(f'{LIGANDS_OUT}/{uni}', exist_ok=True)
    os.makedirs(f'{TEMP}/ligand_extract/{uni}', exist_ok=True)
    os.makedirs(f'{MODELS}/{uni}', exist_ok=True)

    print(f'Starting {uni}')
    ref_pdb = reference_pdb[uni]
    ref_lig = pdb_lig_dict[ref_pdb]

    temp_df = final_bm[final_bm['UniProt_ID'] == uni].reset_index(drop=True)
    lig_list = temp_df.Lig.to_list()
    failed = []
    not_found = []

    for i, pdb, lig in temp_df[['PDB', 'Lig']].itertuples():
        if lig in model_not_available:
            continue
        else:
            model_mol = Chem.SDMolSupplier(f'{MODELS}/{uni}/{lig}.sdf')[0]
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(model_mol)

        cmd.load(filename=f'{ALIGNED}/{uni}/{ref_pdb}.pdb', object=ref_pdb)
        ref_chain = temp_df[temp_df['PDB'] == ref_pdb].Chain.iloc[0]
        cmd.extract(selection=f'r. {ref_lig} and c. {ref_chain}', name='ref_lig')
        
        cmd.load(filename=f'{ALIGNED}/{uni}/{pdb}.pdb', object=pdb)
        cmd.select(selection=f'r. {lig} within 5 of ref_lig', name=f'{pdb}_{lig}')
        cmd.extract(selection=f'br. {pdb}_{lig}', name=f'{pdb}_{lig}')
        if cmd.count_atoms(selection=f'{pdb}_{lig}') != 0:
            cmd.save(filename=f'{TEMP}/ligand_extract/{uni}/{pdb}_{lig}.pdb', selection=f'{pdb}_{lig}')
            cmd.reinitialize()
            try:
                mol = Chem.MolFromPDBFile(f'{TEMP}/ligand_extract/{uni}/{pdb}_{lig}.pdb')
                mol = AllChem.AssignBondOrdersFromTemplate(model_mol, mol)
                mol.SetProp('_Name', f'{pdb}_{lig}')
                writer = Chem.SDWriter(f'{LIGANDS_OUT}/{uni}/{pdb}_{lig}.sdf')
                writer.write(mol)
            except Exception as e:
    #                 print(e)
                failed.append(f'{pdb}_{lig}')
        else:
            not_found.append(f'{pdb}_{lig}')
            missing_ligands.write(f'{uni}, {pdb}, {lig}\n')
            cmd.reinitialize()

    log.write(f'{uni}\t{len(temp_df)}\t{len(not_found)}\t{str(not_found)}\t{str(failed)}\n')

log.close()
missing_ligands.close()

# Remove problematic ligands
log_extraction = pd.read_csv(f'{LOGS}/extract_ligand_log.tsv', sep='\t')
log_extraction['Ligands'] = log_extraction['Ligands'].apply(lambda x: eval(x))
log_extraction['Failed_ligands_to_sdf'] = log_extraction['Failed_ligands_to_sdf'].apply(lambda x: eval(x))
failed_ligands_extraction = log_extraction.apply(lambda row: row['Ligands'] + row['Failed_ligands_to_sdf'], axis=1).tolist()

combined_list = []
for sublist in failed_ligands_extraction:
    combined_list.extend(sublist)

failed_ligands = [i.split(sep='_')[1] for i in combined_list]
final_bm = final_bm[~final_bm['Lig'].isin(failed_ligands)].reset_index(drop=True)
final_bm

final_lig_number = pd.DataFrame(columns=['UniProt_ID', 'Ligands'])
for i, uni in enumerate(log_extraction.UniProt_ID.to_list()):
    final_lig_number.loc[i, 'UniProt_ID'] = uni
    final_lig_number.loc[i, 'Ligands'] = len(final_bm[final_bm['UniProt_ID'] == uni])
    
final_lig_number = final_lig_number.sort_values(by='Ligands')

inal_lig_number_clean = final_lig_number[final_lig_number.Ligands != 0]
final_lig_number_clean['Ligands'] = final_lig_number_clean['Ligands'].astype(int)

# Save final BindingMOAD dataframe
final_bm.to_csv(f'{DATA}/BindingMOAD_selected_templates.csv', index=None)







log.close()


