import argparse
import pandas as pd
import rdkit
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdFMCS, DataStructs, Draw
import os
import numpy as np
from spyrmsd.molecule import Molecule
from spyrmsd.rmsd import rmsdwrapper
import multiprocessing
from tqdm import tqdm


parser = argparse.ArgumentParser(description='Docking results filtering based on the MCS RMSD of the template-ligand pair')
parser.add_argument('--input_df', '-i', required=True)
parser.add_argument('--docking_poses', '-d', required=True)
parser.add_argument('--reference_ligands', '-r', required=True)
parser.add_argument('--output', '-o', help='Provide output directory location', required=True)
parser.add_argument('--cutoff', '-c', help='Choose MCS RMSD cutoff', required=True)
args=parser.parse_args()


DF = args.input_df
REFERENCE_LIGANDS = args.reference_ligands
OUTPUT = args.output
POSES = args.docking_poses
CUTOFF = float(args.cutoff)
SELECTED_POSES = f'{OUTPUT}/filtered_{CUTOFF}Å'


def main():
    df = pd.read_csv(DF)
    df['Template'] = df.PDB + '_' + df.Lig
    final_rmsd = []

    smi_dict = dict(zip(df.Template, df.Smiles_String))

    for i, uni, temp, smi in df[['UniProt_ID','Template','Smiles_String']].itertuples():
        print(f'Starting {i+1}/{len(df)}')
        path_to_poses = f'{POSES}/{uni}/scored_docked_in_{temp}.sdf'
            
        temp_mol = Chem.SDMolSupplier(f'{REFERENCE_LIGANDS}/{uni}/{temp}.sdf', removeHs=False)[0]
        temp_2d = Chem.MolFromSmiles(smi)
        try:
            temp_mol = AllChem.AssignBondOrdersFromTemplate(mol=temp_mol, refmol=temp_2d)
        except: pass

        if os.path.exists(path_to_poses):
            poses = Chem.SDMolSupplier(path_to_poses, removeHs=False)
            
            mcs_rmsd_list = []
            print('Starting MCS RMSD')
            for pose in tqdm(poses, total=len(poses)):
                dock_smi = smi_dict[pose.GetProp('_Name')[:8]]
                dock_2d = Chem.MolFromSmiles(dock_smi)

                soft_mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([temp_mol,dock_2d],
                                          atomCompare=rdFMCS.AtomCompare.CompareElements,
                                          bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                          ringMatchesRingOnly=False,
                                          timeout=1).smartsString)
                try:
                    pose = AllChem.AssignBondOrdersFromTemplate(mol=pose, refmol=dock_2d)
                except: pass

                rmsd = None
                try:
                    rmsd = get_sub_rmsd(temp_mol, pose, soft_mcs)
                except Exception as e:
                    if 'Graphs are not isomorphic' in str(e):
                        try:
                            updated_mcs = fix_connectivity(temp_mol, dock_2d, soft_mcs)
                            if len(Chem.GetMolFrags(updated_mcs, asMols=True, sanitizeFrags=False)) >=2:
                                updated_mcs = check_fragments(updated_mcs)

                            rmsd = get_sub_rmsd(temp_mol, pose, updated_mcs)
                        except Exception as e2:
                            print(f'Problematic: {uni} {temp} {lig} {e2}')

                mcs_rmsd_list.append({
                        'UniProt_ID': uni, 
                        'Template': temp,
                        'Docked': pose.GetProp('_Name')[:8],
                        'Pose_ID': pose.GetProp('_Name')[9:],
                        'mcs_rmsd': rmsd, 
                        'Hyde_score': pose.GetProp('BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]')
                                    })
                    
            df_out = pd.DataFrame(mcs_rmsd_list)
            if not df_out.empty:
                ligands = df_out.Docked.unique()
                for lig in ligands:

                    os.makedirs(f'{SELECTED_POSES}/{uni}/csv', exist_ok=True)
                    os.makedirs(f'{SELECTED_POSES}/{uni}/poses', exist_ok=True)

                    good_docked = df_out[(df_out.Docked == lig) & (df_out.mcs_rmsd < CUTOFF)]
                    good_docked = good_docked.copy()

                    good_docked.loc[:, 'Hyde_score'] = good_docked['Hyde_score'].astype(float)

                    if len(good_docked) > 0:
                        pose_id = good_docked.loc[good_docked.Hyde_score.idxmin()].Pose_ID
                        os.system(f"sed -n '/{lig}_{pose_id}/,/$$$$/p' {path_to_poses} >> {SELECTED_POSES}/{uni}/poses/filtered_scored_docked_in_{temp}.sdf")
                df_out.to_csv(f'{SELECTED_POSES}/{uni}/csv/{temp}_poses_mcs_rmsd.csv', index=False)

                if os.path.exists(f'{SELECTED_POSES}/{uni}/poses/filtered_scored_docked_in_{temp}.sdf'):
                    selected_poses = Chem.SDMolSupplier(f'{SELECTED_POSES}/{uni}/poses/filtered_scored_docked_in_{temp}.sdf')
                    print('Starting reference RMSD...')
                    for pose in tqdm(selected_poses, total=len(selected_poses)):
                        dock = pose.GetProp('_Name')[:8]
                        dock_mol = Chem.MolFromSmiles(smi_dict[dock])

                        ref_mol = Chem.SDMolSupplier(f'{REFERENCE_LIGANDS}/{uni}/{dock}.sdf')[0]
                        
                        try:
                            pose = AllChem.AssignBondOrdersFromTemplate(mol=pose, refmol=ref_mol)
                        except: pass
                        rmsd = None
                        try:
                            rmsd = get_rmsd(ref_mol, pose)
                        except Exception as e1: print(e1) 
                        final_rmsd.append({'UniProt_ID': uni, 'Template': temp, 'Docked': dock, 
                                           'Pose_ID': pose.GetProp('_Name')[9:],
                                           'spyrmsd': rmsd
                                           })
                        

#    final_df = pd.DataFrame(final_rmsd)
#    final_df.to_csv(f'{OUTPUT}/mcs_rmsd_dataset_{CUTOFF}.csv', index=False)


def get_mcs_conf(mol, matches, mcs):
    # Create a new molecule with selected atoms and their bonds
    new_mol = Chem.RWMol()

#     matches = mol.GetSubstructMatch(mcs)
    # Add atoms and bonds to the editable molecule
    atom_mapping = {}  # To keep track of the mapping between old and new atom indices
    for atom_idx in matches:
        old_atom = mol.GetAtomWithIdx(atom_idx)
        new_atom_idx = new_mol.AddAtom(old_atom)
        atom_mapping[atom_idx] = new_atom_idx


    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in matches and bond.GetEndAtomIdx() in matches:
            # Add the bond if both atoms are in the selected list
            new_mol.AddBond(atom_mapping[bond.GetBeginAtomIdx()],
                            atom_mapping[bond.GetEndAtomIdx()],
                            bond.GetBondType())


    # Convert the editable molecule to a new molecule object
    new_mol = new_mol.GetMol()
#     Chem.SanitizeMol(new_mol)

    conf = mol.GetConformer()

    coords = []
    for s in matches:
        coords.append(tuple(conf.GetAtomPosition(s)))

    new_conf = Chem.Conformer(new_mol.GetNumAtoms())
    for i,(x,y,z) in enumerate(coords):
        new_conf.SetAtomPosition(i,Geometry.Point3D(x,y,z))

    new_conf.Set3D(False)
    new_mol.AddConformer(new_conf)

    return new_mol


def get_rmsd(ref, mol):
    d = Molecule.from_rdkit(ref)
    t = Molecule.from_rdkit(mol)

    rmsd = round(float(rmsdwrapper(d, t)[0]), 3)

    return rmsd


def get_sub_rmsd(dock, temp, mcs):

    dock_matches = dock.GetSubstructMatches(mcs)
    temp_matches = temp.GetSubstructMatches(mcs)

    dock_confs = []
    for dock_match in dock_matches:
        dock_confs.append(get_mcs_conf(dock, dock_match, mcs))

    temp_confs = []
    for temp_match in temp_matches:
        temp_confs.append(get_mcs_conf(temp, temp_match, mcs))

    dock_spy = [Molecule.from_rdkit(dock_conf) for dock_conf in dock_confs]
    temp_spy = [Molecule.from_rdkit(temp_conf) for temp_conf in temp_confs]

    rmsds = []
    for d in dock_spy:
        for t in temp_spy:
            rmsds.append(round(float(rmsdwrapper(d, t)[0]), 3))

    return min(rmsds)


def fix_connectivity(mol1, mol2, mcs):
    
    mol1_matches = mol1.GetSubstructMatches(mcs)
    mol2_matches = mol2.GetSubstructMatches(mcs)

    mol1_confs = []
    for mol1_match in mol1_matches:
        mol1_confs.append(get_mcs_mol(mol1, mol1_match, mcs))

    mol2_confs = []
    for mol2_match in mol2_matches:
        mol2_confs.append(get_mcs_mol(mol2, mol2_match, mcs))
    
    mol1_conf = mol1_confs[0]
    mol2_conf = mol2_confs[0]
    
    to_remove = []
    for i in range(mcs.GetNumAtoms()):

        mol1_atom = mol1_conf.GetAtomWithIdx(i)
        mol2_atom = mol2_conf.GetAtomWithIdx(i)

        if mol1_atom.GetSymbol() == mol2_atom.GetSymbol():
            if not len(mol1_atom.GetNeighbors()) == len(mol2_atom.GetNeighbors()):
                to_remove.append(i)
    if to_remove:
        new_mcs = Chem.EditableMol(mcs)
        for i in reversed(to_remove):
            new_mcs.RemoveAtom(i)
            
        return new_mcs.GetMol()

    
def check_fragments(mol):
    
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    size = []
    for frag in frags:
        size.append(frag.GetNumAtoms())
    updated_mcs = frags[size.index(max(size))]
    
    return updated_mcs


def get_mcs_mol(mol, matches, mcs):
    # Create a new molecule with selected atoms and their bonds
    new_mol = Chem.RWMol()

#     matches = mol.GetSubstructMatch(mcs)
    # Add atoms and bonds to the editable molecule
    atom_mapping = {}  # To keep track of the mapping between old and new atom indices
    for atom_idx in matches:
        old_atom = mol.GetAtomWithIdx(atom_idx)
        new_atom_idx = new_mol.AddAtom(old_atom)
        atom_mapping[atom_idx] = new_atom_idx


    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in matches and bond.GetEndAtomIdx() in matches:
            # Add the bond if both atoms are in the selected list
            new_mol.AddBond(atom_mapping[bond.GetBeginAtomIdx()],
                            atom_mapping[bond.GetEndAtomIdx()],
                            bond.GetBondType())

    new_mol = new_mol.GetMol()

    return new_mol


if __name__ == "__main__":
    main()
