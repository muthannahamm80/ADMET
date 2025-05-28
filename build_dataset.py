import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import SaltRemover

def strip_salt(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    remover = SaltRemover.SaltRemover()
    stripped = remover.StripMol(mol)
    return Chem.MolToSmiles(stripped, canonical=True)

def to_uM(value, units):
    if units == "nM":
        return value / 1000
    if units == "µM" or units == "uM":
        return value
    if units == "mM":
        return value * 1000
    return None

def get_target_id(target_keyword):
    target = new_client.target
    targets = target.filter(target_components__component_synonyms__icontains=target_keyword)
    return targets[0]['target_chembl_id'] if targets else None

herg_target = get_target_id("hERG")
cyp_target = get_target_id("CYP3A4")

herg_df = new_client.activity.filter(target_chembl_id=herg_target, standard_type="IC50").only(
    ["molecule_chembl_id", "canonical_smiles", "standard_value", "standard_units", "pchembl_value"]).to_dataframe()
herg_df = herg_df.dropna(subset=["canonical_smiles", "standard_value", "standard_units"])
herg_df["value"] = herg_df.apply(lambda x: to_uM(float(x["standard_value"]), x["standard_units"]), axis=1)
herg_df["endpoint"] = "herg"
herg_df = herg_df[["molecule_chembl_id", "canonical_smiles", "endpoint", "value", "pchembl_value"]]
herg_df = herg_df.rename(columns={"molecule_chembl_id": "chembl_id", "pchembl_value": "pchembl"})

cyp_df = new_client.activity.filter(target_chembl_id=cyp_target, standard_type="IC50").only(
    ["molecule_chembl_id", "canonical_smiles", "standard_value", "standard_units", "pchembl_value"]).to_dataframe()
cyp_df = cyp_df.dropna(subset=["canonical_smiles", "standard_value", "standard_units"])
cyp_df["value"] = cyp_df.apply(lambda x: to_uM(float(x["standard_value"]), x["standard_units"]), axis=1)
cyp_df["endpoint"] = "cyp3a4"
cyp_df = cyp_df[["molecule_chembl_id", "canonical_smiles", "endpoint", "value", "pchembl_value"]]
cyp_df = cyp_df.rename(columns={"molecule_chembl_id": "chembl_id", "pchembl_value": "pchembl"})

merged = pd.concat([herg_df, cyp_df])
merged["canonical_smiles"] = merged["canonical_smiles"].apply(strip_salt)
merged = merged.dropna(subset=["canonical_smiles", "value"])
pivot = merged.pivot_table(index="chembl_id", columns="endpoint", values="value", aggfunc="first").dropna()
pivot = pivot.reset_index()
valid_ids = pivot["chembl_id"].tolist()
final_df = merged[merged["chembl_id"].isin(valid_ids)]
limited_ids = final_df["chembl_id"].drop_duplicates().head(300).tolist()
final_df = final_df[final_df["chembl_id"].isin(limited_ids)]

assert (final_df["value"] > 0).all(), "Non-positive values detected"
assert not final_df.duplicated(subset=["chembl_id", "endpoint"]).any(), "Duplicate chembl_id-endpoint entries found"

final_df.to_csv("dataset.csv", index=False)
print("Generated dataset.csv with", len(final_df), "rows")
