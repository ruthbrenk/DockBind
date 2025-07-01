import pandas as pd

dfs = []
for i in range(1, 4):
    path = f"../data/BindingDB_ligands_selected_for_docking{i}.csv"
    df = pd.read_csv(path)
    dfs.append(df)

# Concatenate them
combined_df = pd.concat(dfs, ignore_index=True)

# Save to a new CSV
combined_df.to_csv("../data_DockBind/BindingDB_ligands_selected_all.csv", index=False)

print("Concatenation complete and file saved.")
