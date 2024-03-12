import sys
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt


def calculate_mz_if_smiles(inp):
    # Convert the SMILES string to a molecule object
    molecule = Chem.MolFromSmiles(inp)

    # Calculate the monoisotopic mass
    monoisotopic_mass = CalcExactMolWt(molecule)

    return monoisotopic_mass

if __name__ == "__main__":
    # Check if the number of arguments is valid
    if len(sys.argv) < 2:
        sys.exit(1)
    #print("hi from the exe")
    input_string = sys.argv[1]

    try:
        
        res = calculate_mz_if_smiles(input_string)
    except:
        try:
            res = float(input_string)
        except: res = "[]"
        
    print(res)