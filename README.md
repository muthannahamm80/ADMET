# ADMET Dataset Builder

This project creates a dataset of ≤300 compounds from ChEMBL that have both hERG and CYP3A4 data.

## Environment
- Python 3.8+
- Packages:
  - pandas
  - rdkit
  - chembl_webresource_client

## Setup
```bash
pip install pandas rdkit-pypi chembl_webresource_client
```

## Run
```bash
python build_dataset.py
```
This will generate `dataset.csv` in the current directory.
