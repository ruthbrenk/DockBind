import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from tqdm import tqdm
from rdkit.Chem import rdFMCS
from spyrmsd.molecule import Molecule
from spyrmsd.rmsd import rmsdwrapper
import os
from tqdm import tqdm
from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize

if not os.path.exists('free_docking/docked'):
    print('To run this you need docked and rescored poses')


# Load Astex

astex = pd.read_csv('../data/astex_selected_for_docking.csv')
astex['Template'] = astex.PDB + '_' + astex.Lig

astex = astex[astex.UniProt_ID != 'P04818']

# Calculate simiarity matrix
temp_smi_dict = dict(zip(astex.Template, astex.Smiles_String))
mol_dict = {name: Chem.MolFromSmiles(smiles) for name, smiles in temp_smi_dict.items()}
fp_dict = {name: AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol) for name, mol in mol_dict.items() if mol is not None}

ligand_names = list(fp_dict.keys())
astex = astex[astex['Template'].isin(ligand_names)]


uniprots = astex.UniProt_ID.unique()
fd_flexx_rmsd = []
fd_hyde_rmsd = []

print('Starting free docking...')
for i, uni in enumerate(uniprots):
    template_df = astex[astex.UniProt_ID == uni]
    templates = template_df.Template.unique()
    print(f'Starting {i+1}/{len(uniprots)}')
    print(f'Templates: {len(templates)}')
    
    for template in templates:
        path_to_docked = f'free_docking/docked/{uni}/docked_in_{template}.sdf'
        path_to_scored = f'free_docking/scored/{uni}/scored_docked_in_{template}.sdf'

        if os.path.exists(path_to_docked) and os.path.getsize(path_to_docked) > 0:

            docked_mols = Chem.SDMolSupplier(path_to_docked)

            for mol in tqdm(docked_mols, total=len(docked_mols)):

                docked = '_'.join(mol.GetProp('_Name').split('_')[:2])
                pose = str(mol.GetProp('_Name').split('_')[2])

                ref_ligand = ref_dict[docked]
                
                rmsd = None
                try:
                    d = Molecule.from_rdkit(ref_ligand)
                    t = Molecule.from_rdkit(mol)

                    rmsd = round(float(rmsdwrapper(d, t)[0]), 3)

                except Exception as e:
                    print(e)

                fd_flexx_rmsd.append({'UniProt_ID': uni,
                                 'Template': template,
                                 'Docked': docked,
                                 'Pose': pose,
                                 'flexx_rmsd': rmsd,
                                 'flexx_score': mol.GetPropsAsDict()['docking-score']}
                                )

        if os.path.exists(path_to_scored) and os.path.getsize(path_to_docked) > 0:

            scored_mols = Chem.SDMolSupplier(path_to_scored)
            for mol in tqdm(scored_mols, total=len(scored_mols)):

                docked = '_'.join(mol.GetProp('_Name').split('_')[:2])
                pose = str(mol.GetProp('_Name').split('_')[2])

                ref_ligand = ref_dict[docked]

                rmsd = None
                try:
                    d = Molecule.from_rdkit(ref_ligand)
                    t = Molecule.from_rdkit(mol)

                    rmsd = round(float(rmsdwrapper(d, t)[0]), 3)

                except Exception as e:
                    print(e)

                fd_hyde_rmsd.append({'UniProt_ID': uni,
                                 'Template': template,
                                 'Docked': docked,
                                 'Pose': pose,    
                                 'hyde_rmsd': rmsd,
                                 'hyde_score': mol.GetPropsAsDict()['BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]']
                                 })

fd_flexx_df = pd.DataFrame(fd_flexx_rmsd)
fd_hyde_df = pd.DataFrame(fd_hyde_rmsd)

fd_merged_df = pd.merge(fd_flexx_df, fd_hyde_df, on=['UniProt_ID', 'Template', 'Docked', 'Pose'])
fd_merged_df = fd_merged_df.drop_duplicates(subset=['Template', 'Docked', 'Pose'], keep='first')
    
fd_uni_temp = fd_merged_df.drop_duplicates(subset=['UniProt_ID', 'Template'])[['UniProt_ID', 'Template']].reset_index(drop=True)
fd_mcs_rmsd_dfs = []
for i, uni, temp in fd_uni_temp.itertuples():
    
    path_to_df = f'free_docking/filtered_1.0Å/{uni}/csv/{temp}_poses_mcs_rmsd.csv'
    
    if os.path.exists(path_to_df):
        temp_df = pd.read_csv(path_to_df, usecols=['UniProt_ID', 'Template', 'Docked', 'Pose_ID', 'mcs_rmsd'])
        fd_mcs_rmsd_dfs.append(temp_df)
        
fd_mcs_df = pd.concat(fd_mcs_rmsd_dfs, ignore_index=True)
fd_mcs_df = fd_mcs_df.rename(columns={'Pose_ID': 'Pose'})

if fd_mcs_df['Pose'].max() <= 10:
    fd_mcs_df['Pose'] = fd_mcs_df['Pose'].apply(lambda x: f"{int(x):02}")
else:
    fd_mcs_df['Pose'] = fd_mcs_df['Pose'].astype(str)
    
fd_flexx_hyde_mcs = pd.merge(fd_merged_df,
                             fd_mcs_df,
                             on=['UniProt_ID', 'Template', 'Docked', 'Pose'],
                             how='left'
                             )


fd_flexx_hyde_mcs.to_csv('../data/free_docking_results.csv', index=False)          

self_fd_flexx_hyde_mcs = fd_flexx_hyde_mcs[fd_flexx_hyde_mcs.Template == fd_flexx_hyde_mcs.Docked]
cros_fd_flexx_hyde_mcs = fd_flexx_hyde_mcs[fd_flexx_hyde_mcs.Template != fd_flexx_hyde_mcs.Docked]

self_fd_flexx_hyde_mcs.to_csv('../data/self_free_docking_results.csv', index=False)
cros_fd_flexx_hyde_mcs.to_csv('../data/cross_free_docking_results.csv', index=False)



print('Starting template-based docking...')
td_flexx_rmsd = []
td_hyde_rmsd = []

for i, uni in enumerate(uniprots):
    template_df = astex[astex.UniProt_ID == uni]
    templates = template_df.Template.unique()
    
    print(f'Starting {i+1}/{len(uniprots)}')
    print(f'Templates: {len(templates)}')
    
    for template in templates:
        path_to_docked = f'template_based_docking/docked/{uni}/docked_in_{template}.sdf'
        path_to_scored = f'template_based_docking/scored/{uni}/scored_docked_in_{template}.sdf'

        if os.path.exists(path_to_docked) and os.path.getsize(path_to_docked) > 0:

            docked_mols = Chem.SDMolSupplier(path_to_docked)

            for mol in tqdm(docked_mols, total=len(docked_mols)):

                docked = '_'.join(mol.GetProp('_Name').split('_')[:2])
                pose = str(mol.GetProp('_Name').split('_')[2])

                ref_ligand = ref_dict[docked]

                rmsd = None
                try:
                    d = Molecule.from_rdkit(ref_ligand)
                    t = Molecule.from_rdkit(mol)

                    rmsd = round(float(rmsdwrapper(d, t)[0]), 3)

                except Exception as e:
                    print(e)

                td_flexx_rmsd.append({'UniProt_ID': uni,
                                 'Template': template,
                                 'Docked': docked,
                                 'Pose': pose,
                                 'flexx_rmsd': rmsd,
                                 'flexx_score': mol.GetPropsAsDict()['docking-score']}
                                )

        if os.path.exists(path_to_scored) and os.path.getsize(path_to_scored) > 0:

            scored_mols = Chem.SDMolSupplier(path_to_scored)
            for mol in tqdm(scored_mols, total=len(scored_mols)):

                docked = '_'.join(mol.GetProp('_Name').split('_')[:2])
                pose = str(mol.GetProp('_Name').split('_')[2])

                ref_ligand = ref_dict[docked]

                rmsd = None
                try:
                    d = Molecule.from_rdkit(ref_ligand)
                    t = Molecule.from_rdkit(mol)

                    rmsd = round(float(rmsdwrapper(d, t)[0]), 3)

                except Exception as e:
                    print(e)

                td_hyde_rmsd.append({'UniProt_ID': uni,
                                 'Template': template,
                                 'Docked': docked,
                                 'Pose': pose,    
                                 'hyde_rmsd': rmsd,
                                 'hyde_score': mol.GetPropsAsDict()['BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]']
                                 })


td_flexx_df = pd.DataFrame(td_flexx_rmsd)
td_hyde_df = pd.DataFrame(td_hyde_rmsd)
                
td_merged_df = pd.merge(td_flexx_df, td_hyde_df, on=['UniProt_ID', 'Template', 'Docked', 'Pose'])

td_uni_temp = td_merged_df.drop_duplicates(subset=['UniProt_ID', 'Template'])[['UniProt_ID', 'Template']].reset_index(drop=True)
td_mcs_rmsd_dfs = []
for i, uni, temp in td_uni_temp.itertuples():
    
    path_to_df = f'template_based_docking/filtered_1.0Å_old/{uni}/csv/{temp}_poses_mcs_rmsd.csv'
    
    if os.path.exists(path_to_df):
        temp_df = pd.read_csv(path_to_df, usecols=['UniProt_ID', 'Template', 'Docked', 'Pose_ID', 'mcs_rmsd'])
        td_mcs_rmsd_dfs.append(temp_df)
        
td_mcs_df = pd.concat(td_mcs_rmsd_dfs, ignore_index=True)
td_mcs_df = td_mcs_df.rename(columns={'Pose_ID': 'Pose'})

if td_mcs_df['Pose'].max() <= 10:
    td_mcs_df['Pose'] = td_mcs_df['Pose'].apply(lambda x: f"{int(x):02}")
else:
    td_mcs_df['Pose'] = td_mcs_df['Pose'].astype(str)
    
td_flexx_hyde_mcs = pd.merge(td_merged_df,
                             td_mcs_df,
                             on=['UniProt_ID', 'Template', 'Docked', 'Pose'],
                             how='left'
                             )
                
td_flexx_hyde_mcs.to_csv('final_data/template_docking_results.csv')


self_td_flexx_hyde_mcs = td_flexx_hyde_mcs[td_flexx_hyde_mcs.Template == td_flexx_hyde_mcs.Docked]
cros_td_flexx_hyde_mcs = td_flexx_hyde_mcs[td_flexx_hyde_mcs.Template != td_flexx_hyde_mcs.Docked]

self_td_flexx_hyde_mcs.to_csv('final_data/self_template_docking_results.csv', index=False)
cros_td_flexx_hyde_mcs.to_csv('final_data/cross_template_docking_results.csv', index=False)


































                
                