import pandas as pd
import os

# Extract poses through simple text operations =)

df = pd.read_csv('../data/DockBind_w_affinities.csv')

print(f'{len(df.BM_Template.unique())} tmeplates to process')

PATH_TO_SCORED = '../scored'
PATH_TO_PDB = '../docking/pdb_files_aligned_sym'
PATH_TO_REF_LIG = '../docking/reference_ligands_2024/'

path_to_dockbind = '../complexes'

os.makedirs(path_to_dockbind, exist_ok=True)

# Poses with < 10 are problematic: if < 10 add "0" in front (accounts for docking processes that generated less than 10 poses)
for i, uni, temp, dock, pose in df[['UniProt_ID', 'BM_Template', 'Docked', 'Pose_ID']].itertuples():

    if pose < 10:
        pose = f'0{str(pose)}'
    
    path_to_dockbind = f'../complexes/{temp}'
    path_to_poses = f'{PATH_TO_SCORED}/{uni}/scored_docked_in_{temp}.sdf'
    os.makedirs(path_to_dockbind, exist_ok=True)
    
    os.system(f"sed -n '/{dock}_{pose}/,/$$$$/p' {path_to_poses} > \
                {path_to_dockbind}/{dock}.sdf")
    os.system(f'cp {PATH_TO_PDB}/{uni}/{temp.split("_")[0]}.pdb ../complexes/{temp}/')
    
    os.system(f'cp {PATH_TO_REF_LIG}/{uni}/{temp}.sdf ../complexes/{temp}/reference_ligand.sdf')
    
    
# PRevious loop will generate an empty file then pose id doesn't have 0 in front, check for empty files
subdirs = [d for d in os.listdir(path_to_dockbind) if os.path.isdir(os.path.join(path_to_dockbind, d))]

empty_files = []

for sub in tqdm(subdirs):
    sub_path = os.path.join(path_to_dockbind, sub)
    for file_name in os.listdir(sub_path):
        file_path = os.path.join(sub_path, file_name)
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            empty_files.append((sub, file_name))

# Report
print(f"Scanned {len(subdirs)} directories.")
if empty_files:
    print(f"Found {len(empty_files)} empty files:\n")
    
    
# Use empty files list to fix this! Rerun without the "0"
for temp, dock in tqdm(empty_files):
    dock = dock.split('.')[0]

    uni, pose = df[(df['BM_Template'] == temp) & (df['Docked'] == dock)][['UniProt_ID', 'Pose_ID']].iloc[0]


    path_to_poses = f'docking/dataset/tmp/free_docking/scored/{uni}/scored_docked_in_{temp}.sdf'
    # dock
    os.system(f"sed -n '/{dock}_{pose}/,/$$$$/p' {path_to_poses} > \
                 {path_to_dockbind}/{temp}/{dock}.sdf")
                 
                 

# Recheck empty files
empty_files_2 = []

for sub in tqdm(subdirs):
    sub_path = os.path.join(path_to_dockbind, sub)
    for file_name in os.listdir(sub_path):
        file_path = os.path.join(sub_path, file_name)
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            empty_files_2.append((sub, file_name))

# Report
print(f"Scanned {len(subdirs)} directories.")

print(f"Found {len(empty_files_2)} empty files:\n")

