import os
import pandas as pd
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem
import requests

BASE = '..'
OUTPUT = f'{BASE}/data/new_pdb'
INPUT = f'{BASE}/data/new_pdb/pdb_files_not_aligned'
REFERENCE = f'{BASE}/pdb_files_aligned'
REFERENCE_LIGANDS = f'{BASE}/reference_ligands'
URL = 'https://files.rcsb.org/ligands/download'
MODELS = f'{OUTPUT}/models'
TMP = f'{OUTPUT}/tmp'

cmd.reinitialize()

new_pdb = pd.read_csv('new_pdb_to_align.csv')
old_ref = pd.read_csv('../data/alignment.csv')

pdb_uni = dict(zip(old_ref.UniProt_ID, old_ref.ref_pdb))

new_df = []

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

        if not os.path.exists(f'{INPUT}/{uni}/{pdb}.pdb'):
            continue
       
        cmd.load(f'{REFERENCE}/{uni}/{ref_pdb}.pdb', object='reference')
        cmd.load(f'{REFERENCE_LIGANDS}/{uni}/{ref_lig}', 'ref_lig')
        cmd.extract(selection='bc. reference within 2 of ref_lig', name='ref')
        cmd.load(f'{INPUT}/{uni}/{pdb}.pdb', object='new')
            
        rmsd = cmd.align('new', 'ref')
        if rmsd[0] > 2.5: continue

        os.makedirs(f'{OUTPUT}/pdb_files_aligned/{uni}', exist_ok=True)
        cmd.save(filename=f'{OUTPUT}/pdb_files_aligned/{uni}/{pdb}.pdb', selection='new')

        
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

        model_not_available = []
        os.makedirs(f'{MODELS}/{uni}', exist_ok=True)
        try:
            response = requests.get(f'{URL}/{mol_name}_ideal.sdf')
            open(f'{MODELS}/{uni}/{mol_name}.sdf', 'wb').write(response.content)

            # Use the downloaded SDF to get a clean SMILES with explicit stereochemistry
            model_mol = Chem.SDMolSupplier(f'{MODELS}/{uni}/{mol_name}.sdf')[0]
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(model_mol)
            smiles = Chem.MolToSmiles(model_mol)
        except Exception as e:
            print(f'{mol_name} {e}')
            model_not_available.append(f'{mol_name}')
        if mol_name in model_not_available:
            continue
        else:
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
                               'Lig_Smiles': smiles}
                             )
            except:
                print(f'{pdb} failed')             
df = pd.DataFrame(new_df)
df.to_csv('new_pdb_confirmed.csv', index=False)
