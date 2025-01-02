'''Counts number of unique molecules'''

import pandas as pd
from rdkit import Chem

# Load your CSV file into a DataFrame
df = pd.read_csv('fragment_list_full.csv', header=None)  # Assuming no header in the CSV

# Remove any rows with missing or invalid SMILES
df = df.dropna(subset=[1])  # Assuming SMILES is in the second column (0-based index)

# Convert SMILES to canonical SMILES to ensure consistency
df[1] = df[1].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) is not None else None)

# Count the number of unique molecules based on canonical SMILES
unique_molecules_count = df.nunique()

print(f"Number of unique molecules: {unique_molecules_count}")
