from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

import itertools
import re
import ast

import json
import csv

# Define a function to calculate molecular weight
def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None
csv.field_size_limit(100000000)

#set tolerance level for matching algorithm
tolerance = 0.4

# Read the CSV file
input_csv = 'fragment_test_slightlybigger.csv'  # Replace with the actual file name

# Initialize a dictionary to store fragment data
molecule_data = []
json_data = []
with open(input_csv, 'r', encoding='utf-8-sig') as csvfile:
    csvreader = csv.reader(csvfile)
    for i, row in enumerate(csvreader):
        if  700 > i > 0:
            fragment_list = row[0]
            original_molecule = row[1]
            spectrum = row[2]
            inchi = row[3]
            fragments = ast.literal_eval(fragment_list)
            original_molecule_mw = calculate_molecular_weight(original_molecule)
            mw_list = []
            peak_list = []
            matching_fragments = []
            
            for fragment in fragments:
                mw = calculate_molecular_weight(fragment.strip("'"))
                if mw is not None:
                    mw_list.append({
                        'Fragment': fragment,
                        'MW': mw
                    })        
                
            for peak in spectrum.split():
                mz, intensity = peak.split(":")  
                matched = False
                for entry in  mw_list:
                    mw_data = entry['MW']
                    fraggy = entry['Fragment']
                    if abs(float(mz) - float(mw_data)) <= tolerance:
                        matching_fragments.append({
                        'm/z value': mz,
                        'Rel. Intensity': intensity,
                        'Fragment': fraggy,
                        'Calculated MW': round(mw_data,2)
                        }) 
                        matched = True
                        break
                if not matched:
                    matching_fragments.append({
                    'm/z value': mz,
                    'Rel. Intensity': intensity,
                    'Fragment': 'not identified',

                    }) 
                        
                    
                                       
                '''
                peak_list.append({
                    'mz': float(mz),
                    'intensity': float(intensity)
                })
                '''
                             
                #sort the list of fragments based on mw       
            sorted_mw_list = sorted(matching_fragments, key=lambda x: x['m/z value'])     
            
            
            json_data.append({ 
                'OriginalMolecule': original_molecule,
                'Calculated MW': round(original_molecule_mw,2),
                'Inchi': inchi,
                'Spectrum': matching_fragments
            })

    # Write the JSON object to a file
    output_json = 'molecule3.json'  # Replace with desired output file name

    with open(output_json, 'a') as jsonfile:
        json.dump(json_data, jsonfile, indent=4)
        '''             
            if mw_list:  # Check if mw_list is not empty
                molecule_data.append({
                    'OriginalMolecule': original_molecule,
                    'MW': original_molecule_mw,
                    'Inchi': inchi,
                    'Fragments': mw_list
                })
            '''










def generate_fragments(parent_molecule_smiles, max_cleaves=3):
    # Create the parent molecule
    parent_molecule = Chem.MolFromSmiles(parent_molecule_smiles)

    if parent_molecule is None:
        print("Invalid SMILES notation for the parent molecule.")
        return

    # Generate fragments
    all_fragments = []
    for cleaves in range(1, max_cleaves + 1):
        # Generate all possible combinations of cleavage sites
        cleavage_combinations = [list(indices) for indices in itertools.combinations(range(parent_molecule.GetNumBonds()), cleaves)]

        for cleavage_sites in cleavage_combinations:
            # Skip cleavage sites containing non-single bonds or aromatic bonds
            if any(
                parent_molecule.GetBondWithIdx(idx).GetBondTypeAsDouble() != 1.0 or
                parent_molecule.GetBondWithIdx(idx).GetIsAromatic()
                for idx in cleavage_sites
            ):
                continue

            # Generate fragment by cleaving the molecule
            cleaved_mol = AllChem.FragmentOnBonds(parent_molecule, cleavage_sites)

            # Convert fragment to SMILES and add to the list
            fragment_smiles = Chem.MolToSmiles(cleaved_mol)
            all_fragments.append(fragment_smiles)

    return all_fragments

# Example usage with the provided SMILES notation
parent_molecule_smiles = "COC(=N3)n(n1)c(C(F)=C3)nc1S(=O)(=O)Nc(c(F)2)c(F)ccc2"
fragments = generate_fragments(parent_molecule_smiles, max_cleaves=3)

# Print generated fragments
for i, fragment in enumerate(fragments, 1):
    print(f"{i}. {fragment}")

# Extract and combine fragments from each line
combined_fragments = []
for line in fragments:
    line1 = re.sub(r'\(\[\d+\*\]\)|\[\d+\*\]|\*', '', line)
    
    line_fragments = line1.split('.')
    # Remove any leading or trailing whitespaces
    line_fragments2 = [fragment.strip() for fragment in line_fragments]
    # Add fragments to the combined list
    combined_fragments.extend(line_fragments2)

# Remove duplicates from the combined fragments
combined_fragments = list(set(combined_fragments))

# Print combined fragments
print("\nCombined Fragments:")
for i, fragment in enumerate(combined_fragments, 1):
    print(f"{i}. {fragment}")
