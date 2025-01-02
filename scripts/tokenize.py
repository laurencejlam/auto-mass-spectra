'''Tokenizes all molecules and fragments.'''

import csv

# Define dictionaries to store mappings
compound_mapping = {}
token_counter = 1

# Open the input CSV file for reading and the output CSV file for writing
with open('full_output3.csv', 'r') as infile, open('full_output3-token-frag.csv', 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)

    # Read the header and write it to the output CSV file
    header = next(reader)
    writer.writerow(header)

    # Iterate over each row in the input CSV file
    for row in reader:
        # Extract the spectrum, inchi, parent smiles, and fragments
        spectrum = row[0]
        inchi = row[1]
        parent_smiles = row[2]
        fragments = [frag for frag in row[3:] if frag]  # Remove empty fragments

        # Tokenize parent smiles if not already tokenized
        if parent_smiles not in compound_mapping:
            compound_mapping[parent_smiles] = token_counter
            token_counter += 1
        
        # Tokenize fragments if not already tokenized
        for fragment in fragments:
            if fragment not in compound_mapping:
                compound_mapping[fragment] = token_counter
                token_counter += 1
        
        # Write the extracted data to the output CSV file with tokenized fragments
        writer.writerow([spectrum, inchi, compound_mapping[parent_smiles]] + [compound_mapping[frag] for frag in fragments])

# Write the compound mapping to a separate CSV file
with open('token-frag-mapping.csv', 'w', newline='') as mapping_file:
    writer = csv.writer(mapping_file)
    writer.writerow(['Compound', 'Token'])
    for compound, token in compound_mapping.items():
        writer.writerow([compound, token])
