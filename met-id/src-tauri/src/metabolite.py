import sys
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

if __name__ == "__main__":
    # Check if the number of arguments is valid
    if len(sys.argv) < 3:
        sys.exit(1)
    #print("hi from the exe")
    smarts_pattern = sys.argv[1]
    smiles_strings = sys.argv[2:]


    #for smiles in smiles_strings:
    results = substruct_matches(smiles_strings, smarts_pattern)

    print(results)