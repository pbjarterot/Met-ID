import sys
import json
import struct
import sqlite3
import numpy as np
from argparse import ArgumentParser
from rdkit import Chem
import time

from concurrent.futures import ProcessPoolExecutor
from functools import partial

def substruct_match_worker(smiles_chunk, smarts_mol):
    """ Substructure matching for a chunk of SMILES. """
    results = []
    for smiles in smiles_chunk:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            results.append(len(mol.GetSubstructMatches(smarts_mol)))
        else:
            results.append(0)
    return results

def optimized_transactions(conn, data_to_insert, smarts, table_name, name):
    try:
        # Preprocess SMARTS once
        smarts_mol = Chem.MolFromSmarts(smarts)
        if not smarts_mol:
            raise ValueError("Invalid SMARTS pattern")

        # Split data into chunks
        batch_size = 500  # Larger batch size for efficiency
        batches = [data_to_insert[i:i+batch_size] for i in range(0, len(data_to_insert), batch_size)]

        # Process in parallel
        print("Starting parallel RDKit processing...")
        with ProcessPoolExecutor() as executor:
            for i, batch in enumerate(batches):
                indices = [i[0] for i in batch]
                smiles_list = [i[1] for i in batch]
                print(i/len(batches))
                # Parallel computation of matches
                worker_func = partial(substruct_match_worker, smarts_mol=smarts_mol)
                results = executor.submit(worker_func, smiles_list).result()

                # Update database sequentially to avoid contention
                updates = [(value, index) for value, index in zip(results, indices)]
                conn.executemany(f"UPDATE {table_name} SET '{name}' = ? WHERE rowid = ?", updates)
        
        conn.commit()
        print("Transaction committed successfully.")

    except sqlite3.Error as e:
        print(f"SQLite error: {e}")
        conn.rollback()
    finally:
        conn.close()

def main():
    path = sys.argv[1]
    smarts = sys.argv[2]
    name = sys.argv[3]

    print(path)
    conn = sqlite3.connect(path)
    cursor = conn.cursor()

    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    print(cursor.fetchall())

    cursor.execute("SELECT id, smiles FROM metabolites")
    results = cursor.fetchall()

    optimized_transactions(cursor, results, smarts, "functional_groups", name)


if __name__ == "__main__":
    main()
