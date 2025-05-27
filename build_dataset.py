import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import time

# Initialize clients
activity = new_client.activity
molecule = new_client.molecule

# Targets for hERG and CYP3A4
targets = {
    'herg': 'CHEMBL203',
    'cyp3a4': 'CHEMBL1569'
}

def fetch_activities(target_id):
    activities_list = []
    count = 0
    activities = activity.filter(target__target_components__accession=target_id)
    for act in activities:
        try:
            chembl_id = act['molecule_chembl_id']
            mol_struct = molecule.get(chembl_id).get('molecule_structures', {}).get('canonical_smiles', None)
            if mol_struct is None:
                continue
            std_type = act.get('standard_type')
            std_value = act.get('standard_value')
            if std_value is None:
                continue
            value = float(std_value)
            if std_type and 'p' in std_type:
                value_muM = 10 ** (-value)
            else:
                value_muM = value
            activities_list.append({
                'chembl_id': chembl_id,
                'canonical_smiles': mol_struct,
                'endpoint': target_id,
                'value': value_muM
            })
            count += 1
            if count % 100 == 0:
                print(f"Fetched {count} records for {target_id}...")
        except:
            continue
    print(f"Finished fetching {count} records for {target_id}.")
    return activities_list

# Fetch data for each target
data_records = []
for endpoint, target_id in targets.items():
    print(f"Starting data fetch for {endpoint}...")
    acts = fetch_activities(target_id)
    data_records.extend(acts)
    time.sleep(1)

# Prepare DataFrame
df = pd.DataFrame(data_records)

# Strip salts
df['canonical_smiles'] = df['canonical_smiles'].apply(lambda s: s.split('.')[0] if s else None)

# Remove invalid SMILES
def is_valid_smile(s):
    try:
        return Chem.MolFromSmiles(s) is not None
    except:
        return False

df = df[df['canonical_smiles'].apply(is_valid_smile)]

# Keep only molecules with both endpoints
grouped = df.groupby('chembl_id')
common_ids = [chembl_id for chembl_id, grp in grouped
              if set(grp['endpoint']) == {'herg', 'cyp3a4'}]

# Limit to 300 compounds
selected_ids = common_ids[:300]
final_df = df[df['chembl_id'].isin(selected_ids)]

# Save
final_df.to_csv('dataset.csv', index=False)
print("Dataset saved to 'dataset.csv'.")
