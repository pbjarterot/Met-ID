import sys
import sqlite3
from rdkit import Chem

def substruct_matches(smiles_strings, smarts_pattern):
    smarts_mol = Chem.MolFromSmarts(smarts_pattern)
    if smarts_mol is None:
        return "Invalid SMARTS pattern"
    
    results = []
    for smiles in smiles_strings:
        input_mol = Chem.MolFromSmiles(smiles)
        if input_mol is None:
            results.append(0)
        else:
            matches = len(input_mol.GetSubstructMatches(smarts_mol))
            results.append(matches)
    return results

def chunk_list(input_list, n):
    """Break a list into chunks of size n."""
    return [input_list[i:i+n] for i in range(0, len(input_list), n)]

def get_smiles_from_db(db_path, table_name):
    conn = sqlite3.connect(db_path)
    print(conn)
    cursor = conn.cursor()

    # Select the 'smiles' column from the table
    cursor.execute(f'SELECT smiles FROM {table_name}')

    # Fetch all rows and print them
    rows = cursor.fetchall()
    smiles_values = [row[0] for row in rows]

    # Close the connection
    cursor.close()
    conn.close()
    return smiles_values


if __name__ == "__main__":
    # Check if the number of arguments is valid
    if len(sys.argv) < 3:
        sys.exit(1)
    print("hi from the exe")
    smarts_pattern = sys.argv[1]
    db_path = sys.argv[2]
    table_name = sys.argv[3]
    chunk_size = sys.argv[4]

    smiles_strings = get_smiles_from_db(db_path, table_name)
    chunked_smiles = chunk_list(smiles_strings, int(chunk_size))

    for smiless in chunked_smiles:
        results = substruct_matches(smiless, smarts_pattern)

        print(results)