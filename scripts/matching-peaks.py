'''Matches fragments to mass spectral peaks.'''

import csv
from rdkit import Chem
from rdkit.Chem import Descriptors
import ast
import json

# Define a function to calculate molecular weight
def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None

def read_token_mapping(mapping_csv):
    token_mapping = {}
    with open(mapping_csv, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row
        for row in reader:
            smiles, token = row
            token_mapping[token] = smiles
    return token_mapping

def find_smiles_by_number(token_mapping, target_number):
    return token_mapping.get(target_number)

# Set tolerance level for matching algorithm
tolerance = 2.0

# Read the CSV files
input_csv = 'full_output3-token-frag.csv'  
mapping_csv = 'token-frag-mapping.csv'

# Initialize a dictionary to store fragment data
json_data = []

# Read token mapping into memory
token_mapping = read_token_mapping(mapping_csv)

# Process the input CSV
with open(input_csv, 'r', encoding='utf-8-sig') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)  # Skip the header row
    for i, column in enumerate(csvreader, start=1):
        if 10> i >0 :
            spectrum = column[0]
            inchi = column[1]
            original_molecule = column[2]
            fragment_list = column[2:]
            print(i)

            original_molecule_smiles = find_smiles_by_number(token_mapping, original_molecule)
            original_molecule_mw = calculate_molecular_weight(original_molecule_smiles)
            mw_list = []
            matching_fragments = []

            for fragment in fragment_list:
                fragment_smiles = find_smiles_by_number(token_mapping, fragment)
                mw = calculate_molecular_weight(fragment_smiles)
                if mw is not None:
                    mw_list.append({
                        'Token': fragment,
                        'Frag_smiles': fragment_smiles,
                        'MW': mw
                    })

            for peak in spectrum.split():
                mz, intensity = peak.split(":")  

                closest_match = None
                min_difference = float('inf')

                for entry in mw_list:
                    mw_data = entry['MW']
                    fraggy = entry['Frag_smiles']
                    fraggy_token = entry['Token']
                    difference = abs(float(mz) - float(mw_data))

                    if difference <= tolerance and difference < min_difference:
                        closest_match = {
                            'm/z value': mz,
                            'Rel. Intensity': intensity,
                            'Fragment': fraggy,
                            'Token': fraggy_token,
                            'Calculated MW': round(mw_data, 2)
                        } 
                        min_difference = difference

                if closest_match is not None:
                    matching_fragments.append(closest_match)
                else:
                    matching_fragments.append({
                        'm/z value': mz,
                        'Rel. Intensity': intensity,
                        'Fragment': 'not identified',
                    })

            json_data.append({ 
                'OriginalMolecule': original_molecule,
                'Calculated MW': round(original_molecule_mw, 2),
                'Inchi': inchi,
                'Spectrum': matching_fragments
            })

# Write the JSON object to a file
output_json = 'tokenised-examp-27-2.json'  
with open(output_json, 'w') as jsonfile:
    json.dump(json_data, jsonfile, indent=4)
