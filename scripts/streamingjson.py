'''Streams JSON to csv'''

import json
import csv

# Step 1: Read the JSON file
with open('data.json', 'r', encoding="utf8") as json_file:
    data = json.load(json_file)

#Step 2: Open the csv file, a for append
with open('book1.csv', 'a', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    
    # Add column headers
    csv_writer.writerow(['Inchi', 'SMILES', 'Spectrum'])

    #Step 3: Run the iterations

    #this is the number of iterations
    y = 98152
    
    for x in range(y):
        spectrum_data = data[x]['spectrum']

        # this looks up the inchi, if not available, looks up inchikey
        metadata_data = data[x]["compound"][0]["metaData"]
        for input in metadata_data:
            if input.get("name") =="InChI":
                inchi_data = input.get("value")
                break
        
        #in case no inchikey, this looks for inchi
        if inchi_data == None: 
            for input2 in metadata_data:
                if input2.get("name") =="InChIkey":
                    inchi_data = input.get("value")
                    break 
        
            # this looks up the SMILES value
        for entry in metadata_data:
            if entry.get("name") =="SMILES":
                smiles_data = entry.get("value")
                break
            else: smiles_data = None
        
            #this forgoes the data if it is longer than the character limit for a cell, 
            #then writes the data in a csv file
        
        csv_writer.writerow([inchi_data, smiles_data, spectrum_data])
        print(x)
        x+1
       
            
