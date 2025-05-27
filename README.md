# ChEMBL Compound Dataset for hERG and CYP3A4

This script retrieves activity data from the ChEMBL database for the hERG and CYP3A4 targets, converts all potency measures to micromolar (μM), removes salts and invalid SMILES, and filters to include only compounds with data for both targets. The dataset is limited to 300 compounds and saved as `dataset.csv`.

## Usage

### Environment

- Python 3.8+
- Required Python packages:
  - chembl_webresource_client
  - pandas
  - rdkit-pypi

### Installation

Install the necessary packages with:

```bash
pip install chembl_webresource_client pandas rdkit-pypi
