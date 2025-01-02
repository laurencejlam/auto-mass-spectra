import pandas as pd

# Load your original CSV file into a pandas DataFrame
original_df = pd.read_csv('book1.csv')

# Create an empty DataFrame to store the InChI and corresponding rows
result_df = pd.DataFrame(columns=['InChI', 'SMILES','Row Numbers'])

# Iterate through unique compounds in the first column of the original DataFrame
for inchi_str in original_df.iloc[:, 0].unique():
    # Find the rows where the current compound appears
    rows = original_df.index[original_df.iloc[:, 0] == inchi_str].tolist()
    # Get the corresponding SMILES from the second column
    smiles = original_df.loc[original_df.iloc[:, 0] == inchi_str, original_df.columns[1]].values[0]
    # Adjust row numbers by adding 1 to each row value
    rows = [row + 1 for row in rows]
    # Append the InChI and corresponding rows to the result DataFrame
    result_df = pd.concat([result_df, pd.DataFrame({'InChI': [inchi_str], 'SMILES':[smiles], 'Row Numbers': [rows]})], ignore_index=True)

# Save the final result DataFrame to a new CSV file
result_df.to_csv('inchi_sorted5.csv', index=False, header = False)
