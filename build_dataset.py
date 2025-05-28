import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import time

# Initialize clients
activity = new_client.activity
molecule = new_client.molecule

# Targets for hERG and CYP3A4
targets = {
    'herg': 'CHEMBL203',      # hERG
    'cyp3a4': 'CHEMBL1569'    # CYP3A4
}

def fetch_activities(target_chembl_id, max_records=1500):
    """
    Fetch activities from ChEMBL and filter in Python by target.
    """
    activities = list(activity.filter())  # fetch all activities
    filtered = []

    count = 0
    for act in activities:
        try:
            target_info = act.get('target')
            # 'target' can be a dict; check its 'chembl_id'
            if isinstance(target_info, dict):
                if target_info.get('chembl_id') != target_chembl_id:
                    continue
            else:
                continue  # skip if no target info or mismatch

            chembl_id = act['molecule_chembl_id']
            smi = molecule.get(chembl_id).get('molecule_structures', {}).get('canonical_smiles', None)
            if not smi:
                continue
            std_type = act.get('standard_type')
            std_value = act.get('standard_value')
            if std_value is None:
                continue
            val = float(std_value)
            if std_type and 'p' in std_type:
                val = 10 ** (-val)  # convert pIC50 etc. to μM
            filtered.append({
                'chembl_id': chembl_id,
                'canonical_smiles': smi,
                'endpoint': target_chembl_id,
                'value': val
            })
            count += 1
            if count % 100 == 0:
                print(f"Collected {count} records for target {target_chembl_id}...")
            if count >= max_records:
                break
        except:
            continue
    print(f"Total for {target_chembl_id}: {count}")
    return filtered

# Fetch data for both targets
data = []
for name, t_id in targets.items():
    data.extend(fetch_activities(t_id))
    time.sleep(1)

# Create DataFrame
df = pd.DataFrame(data)

# Clean SMILES: keep only the first part (remove salts)
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
    cid for cid, g in grouped
    if set(g['endpoint']) == {'herg', 'cyp3a4'}
]

# Limit to 300 compounds
selected_ids = common_ids[:300]
final_df = df[df['chembl_id'].isin(selected_ids)]

# Save result
final_df.to_csv('dataset.csv', index=False)
print("dataset.csv has been generated.")
