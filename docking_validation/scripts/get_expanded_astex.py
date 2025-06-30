import pandas as pd

bm = pd.read_csv('../data/clean_BindingMOAD.csv')
print(len(bm))
astex = pd.read_csv('../data/astex_diverse_set.csv')

uniprot = pd.read_csv(f'../data/pdb_chain_uniprot.csv',
                      skiprows=[0], usecols=['PDB', 'SP_PRIMARY']).drop_duplicates(subset=['PDB'], keep='first')
uniprot['PDB'] = uniprot['PDB'].str.upper()
pdb_uni = dict(zip(uniprot['PDB'], uniprot['SP_PRIMARY']))

for i, pdb in astex[['PDB']].itertuples():
    try:
        astex.at[i, 'UniProt_ID'] = pdb_uni[pdb]
    except:
        print(f'{pdb} not found')
        
confirmed_uni = astex.UniProt_ID.to_list()
print(len(confirmed_uni))
expanded_astex = bm[bm['UniProt_ID'].isin(confirmed_uni)].reset_index(drop=True)

expanded_astex.to_csv('../data/astex_selected_for_docking.csv', index=False)
