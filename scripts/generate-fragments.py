'''Generates all possible fragments from parent molecules.'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import itertools
import re
import csv

# This code generates the fragments of the parent molecule, and outputs it in a csv
# It generates fragments by iteratively breaking 1-3 single, non-aromatic bonds in the molecule at a time, until all fragments are found.

# Read CSV file containing SMILES strings
input_csv = 'inchi_sorted4.csv'  # Replace with the actual file name
output_csv = 'compound_list_fragments3.csv'  # Replace with desired output file name
csv.field_size_limit(100000000)
# Function that generates fragments
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

    combined_fragments = []
    for line in all_fragments:
        line1 = re.sub(r'\(\[\d+\*\]\)|\[\d+\*\]|\*', '', line)
        #remove consecutive parentheses ()
        line2 = re.sub(r'\(\)', '', line1)
        #delimit by .
        line_fragments = line2.split('.')
        # Remove any leading or trailing whitespaces
        line_fragments2 = [fragment.strip() for fragment in line_fragments]
        # Add fragments to the combined list
        combined_fragments.extend(line_fragments2)

    # Remove duplicates from the combined fragments
    combined_fragments = list(set(combined_fragments))

    return combined_fragments
#where to start iterating
start_row = 10829
#end_row = 10000

 
with open(input_csv, 'r', encoding='utf-8-sig') as csvfile:
    csvreader = csv.reader(csvfile)
    for i, row in enumerate(csvreader):
        if start_row <= i :
            inchi = row[0]  #assuming inchi is in first column of the csv
            parent_molecules_smiles = row[1]  
            numbers = row[2]
            new_fragment_list = []
            print(i)
            if parent_molecules_smiles is not None:            
                new_fragment_list = generate_fragments(parent_molecules_smiles) #this is the magic maker, generate fragments from parent molecule
                new_fragment_list.insert(0, inchi)
                new_fragment_list.insert(1,numbers)
                new_fragment_list.insert(2, parent_molecules_smiles)
                #write into a csv
                with open(output_csv, mode='a', newline='') as outfile:
                    csvwriter = csv.writer(outfile)
                    csvwriter.writerow(new_fragment_list)
                
