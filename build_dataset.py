import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import time

# Initialize clients
activity = new_client.activity
molecule = new_client.molecule

# ChEMBL target IDs
targets = {
    'herg': 'CHEMBL203',      # hERG
    'cyp3a4': 'CHEMBL1569'    # CYP3A4
}

def fetch_activities(target_id, max_records=1500):
    activities_list = []
    print(f"Fetching activities for target: {target_id}")
    all_activities = list(activity.filter())  # fetch all activities
    count = 0
    for act in all_activities:
        try:
            # Check the 'target' attribute
            target_info = act.get('target')
            # 'target' may be a dict with 'chembl_id'
            target_chembl_id = None
            if isinstance(target_info, dict):
                target_chembl_id = target_info.get('chembl_id')
            # Filter manually: process only if target matches
            if target_chembl_id != target_id:
                continue

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
                print(f"Collected {count} activities for {target_id}")
            if count >= max_records:
                break
        except Exception:
            continue
    print(f"Finished fetching {count} activities for {target_id}")
    return activities_list

# Fetch data for both targets with filter
data = []

for target_name, target_id in targets.items():
    acts = fetch_activities(target_id)
    data.extend(acts)
    time.sleep(1)

# Convert to DataFrame
df = pd.DataFrame(data)

# Cleanup SMILES
df['canonical_smiles'] = df['canonical_smiles'].apply(lambda s: s.split('.')[0] if s else None)

# Validate SMILES
def is_valid_smile(s):
    try:
        return Chem.MolFromSmiles(s) is not None
    except:
        return False

df = df[df['canonical_smiles'].apply(is_valid_smile)]

# Find molecules with both targets
grouped = df.groupby('chembl_id')
common_ids = [
    chembl_id for chembl_id, grp in grouped
    if set(grp['endpoint']) == {'herg', 'cyp3a4'}
]

# Limit to 300 compounds
selected_ids = common_ids[:300]
final_df = df[df['chembl_id'].isin(selected_ids)]

# Save the dataset
final_df.to_csv('dataset.csv', index=False)
print("dataset.csv has been saved.")
