import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import time

# Initialize Clients
activity = new_client.activity
molecule = new_client.molecule

# Target IDs for hERG and CYP3A4
targets = {
    'herg': 'CHEMBL203',      # hERG
    'cyp3a4': 'CHEMBL1569'    # CYP3A4
}

def fetch_activities_for_target(target_chembl_id, max_records=1500):
    """
    Fetch all activities, then filter by target_chembl_id in Python.
    """
    print(f"Fetching all activities for target CHEMBL ID: {target_chembl_id} ...")
    all_activities = list(activity.filter())  # fetch all activities
    filtered = []
    count = 0
    for act in all_activities:
        try:
            target_info = act.get('target')
            # 'target' can be a dict containing 'chembl_id'
            if isinstance(target_info, dict):
                if target_info.get('chembl_id') != target_chembl_id:
                    continue
            else:
                continue  # skip if no target info

            chembl_id = act['molecule_chembl_id']
            smi = molecule.get(chembl_id).get('molecule_structures', {}).get('canonical_smiles', None)
            if not smi:
                continue
            std_type = act.get('standard_type')
            std_value = act.get('standard_value')
            if std_value is None:
                continue
            val = float(std_value)
            # Convert pX to μM if needed
            if std_type and 'p' in std_type:
                val = 10 ** (-val)
            filtered.append({
                'chembl_id': chembl_id,
                'canonical_smiles': smi,
                'endpoint': target_chembl_id,
                'value': val
            })
            count += 1
            if count % 100 == 0:
                print(f"Collected {count} activities for target {target_chembl_id} ...")
            if count >= max_records:
                break
        except:
            continue
    print(f"Total for target {target_chembl_id}: {count} activities.")
    return filtered

# Fetch data for both targets
all_data = []

for target_name, target_id in targets.items():
    data = fetch_activities_for_target(target_id)
    all_data.extend(data)
    time.sleep(1)  # pause between targets

# Convert to DataFrame
df = pd.DataFrame(all_data)

# Remove salts: keep main structure
df['canonical_smiles'] = df['canonical_smiles'].apply(lambda s: s.split('.')[0] if s else None)

# Validate SMILES
def is_valid_smile(s):
    try:
        return Chem.MolFromSmiles(s) is not None
    except:
        return False

df = df[df['canonical_smiles'].apply(is_valid_smile)]

# Find molecules with both targets
groups = df.groupby('chembl_id')
common_ids = [cid for cid, g in groups if set(g['endpoint']) == {'herg', 'cyp3a4'}]

# Limit to 300 compounds
selected_ids = common_ids[:300]
final_df = df[df['chembl_id'].isin(selected_ids)]

# Save dataset
final_df.to_csv('dataset.csv', index=False)
print("dataset.csv has been generated.")
