import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re

def is_valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def is_valid_smarts(smarts):
    mol = Chem.MolFromSmarts(smarts)
    return mol is not None

def get_mass(formula):
    parts = re.findall("[A-Z][a-z]?|[0-9]+", formula)
    mass = 0

    for index in range(len(parts)):
        if parts[index].isnumeric():
            continue

        atom = Chem.MolFromSmiles("[" + parts[index] + "]")
        multiplier = int(parts[index + 1]) if len(parts) > index + 1 and parts[index + 1].isnumeric() else 1
        mass += rdMolDescriptors.CalcExactMolWt(atom) * multiplier

    return mass

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py function_name input_data")
        sys.exit(1)

    function_name = sys.argv[1]
    input_data = sys.argv[2]

    try:
        if function_name == "is_valid_smiles":
            result = is_valid_smiles(input_data)
            print(f"Is valid SMILES: {result}")
        elif function_name == "is_valid_smarts":
            result = is_valid_smarts(input_data)
            print(f"Is valid SMARTS: {result}")
        elif function_name == "get_mass":
            result = get_mass(input_data)
            print(f"Molecular weight: {result:.2f}")
        else:
            print("Invalid function name. Available functions: is_valid_smiles, is_valid_smarts, get_mass")
    except Exception as e:
        print(e)