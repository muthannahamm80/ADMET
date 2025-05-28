import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import time

# Initialize clients
activity = new_client.activity
molecule = new_client.molecule

# Targets for hERG and CYP3A4
targets = {
    'herg': 'CHEMBL203',      # hERG target ChEMBL ID
    'cyp3a4': 'CHEMBL1569'    # CYP3A4 target ChEMBL ID
}

def fetch_activities(target_id, max_records=1500):
    """
    Fetch activities from ChEMBL. Since filtering may be unsupported,
    fetch a large number and filter afterwards.
    """
    activities_list = []
    print(f"Fetching activities from ChEMBL (target: {target_id}) ...")
    activities = list(activity.filter())  # fetch all activities
    count = 0
    for act in activities:
        try:
            # Check if the target matches by comparing target ChEMBL id
            if 'target' in act:
                target_info = act['target']
                # 'target' info might be a dict
                target_chembl_id = target_info.get('chembl_id') if isinstance(target_info, dict) else None
                if target_chembl_id != target_id:
                    continue  # skip if not matching target

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
                # Convert pX to μM
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
                print(f"Fetched {count} records for target {target_id}...")

            if count >= max_records:
                break
        except Exception as e:
            continue
    print(f"Total activities fetched for {target_id}: {count}")
    return activities_list

# Fetch activities for both targets
data_records = []

for target_name, target_id in targets.items():
    acts = fetch_activities(target_id, max_records=1500)
    data_records.extend(acts)
    time.sleep(1)  # pause between targets

# Convert to DataFrame
df = pd.DataFrame(data_records)

# Remove salts: keep only the first part of SMILES
df['canonical_smiles'] = df['canonical_smiles'].apply(lambda s: s.split('.')[0] if s else None)

# Validate SMILES
def is_valid_smile(s):
    try:
        return Chem.MolFromSmiles(s) is not None
    except:
        return False

df = df[df['canonical_smiles'].apply(is_valid_smile)]

# Find molecules with both endpoints
grouped = df.groupby('chembl_id')
common_ids = [chembl_id for chembl_id, grp in grouped
              if set(grp['endpoint']) == {'herg', 'cyp3a4'}]

# Limit to 300 compounds
selected_ids = common_ids[:300]
final_df = df[df['chembl_id'].isin(selected_ids)]

# Save dataset
final_df.to_csv('dataset.csv', index=False)
print("dataset.csv has been generated.")
