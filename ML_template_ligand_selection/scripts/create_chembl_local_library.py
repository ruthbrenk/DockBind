import pickle
import pymysql
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
rdkit.RDLogger.DisableLog("rdApp.warning")
import time
from dotenv import load_dotenv
from rdkit import rdBase
#rdBase.rdkitVersion
t0=time.asctime(time.localtime(time.time()))

t1=time.time()
print(f'Starting script {t0}')
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), '../../.env'))

db = pymysql.connect(
    host=os.getenv("DB_HOST"),
    user=os.getenv("DB_USER"),
    password=os.getenv("DB_PASSWORD"),
    database=os.getenv("DB_NAME")
)

cursor = db.cursor(cursor=pymysql.cursors.DictCursor)

query = '''SELECT  molecule_dictionary.chembl_id, compound_structures.canonical_smiles
           FROM molecule_dictionary
           JOIN compound_structures ON molecule_dictionary.molregno = compound_structures.molregno;
            
            '''
cursor.execute(query)
rows = cursor.fetchall()

t2 = time.time()
print(f'Retrieving chembl from the databese took {t2-t1:.2f} seconds, starting library generation')
print(f'{len(rows)} molecules have to be processed, starting now...')
data = []
for index, i in enumerate(rows):
    if not ((index+1)%50000):
        print(f"Processed {index+1} molecules in {(time.time()-t2):.1f} seconds")
    try:
        mol_id = i['chembl_id']
        smi = i['canonical_smiles']

        mol = Chem.MolFromSmiles(smi)
        mol.SetProp('_Name', mol_id)

        fp = Chem.PatternFingerprint(mol,fpSize=1024,tautomerFingerprints=True)
        data.append((smi,fp))
    except Exception as e:
        print(f'{mol_id} failed because {e}')

pickle.dump(data,open('data/chembl34_sssdata.pkl','wb+'))
t3=time.time()
print(f"First checkpoint took {t3-t2:.2f} seconds, starting part 2")


mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
fps = rdSubstructLibrary.TautomerPatternHolder(1024)
for smi,fp in data:
    mols.AddSmiles(smi)
    fps.AddFingerprint(fp)
library = rdSubstructLibrary.SubstructLibrary(mols,fps)
t4=time.time()
print(f"Part two took {t4-t3:.2f} seconds. The library has {len(library)} molecules.")
pickle.dump(library,open('data/chembl31_ssslib.pkl','wb+'))
#print(data)
t5 = time.time()
print('Loading library...')

with open('data/chembl34_ssslib.pkl','rb') as f:
    sslib = pickle.load(inf)
t6 = time.time()
print(f'SubstructLibrary loaded in {t6-t5:.2f} seconds with {len(sslib)} molecules')


