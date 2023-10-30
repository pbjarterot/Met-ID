import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def return_formula(mol):
    return rdMolDescriptors.CalcMolFormula(mol)

def return_mz(mol): 
    return Descriptors.ExactMolWt(mol)


def substruct_matches(smiles, smarts_patterns):
    input_mol = Chem.MolFromSmiles(smiles)
    
    if input_mol is None:
        return "Invalid SMILES pattern"
    
    results = []
    for smarts in smarts_patterns:
        smarts_mol = Chem.MolFromSmarts(smarts)
        if input_mol is None:
            results.append(0)
        elif smarts_mol is None:
            results.append(0)
        else:
            matches = len(input_mol.GetSubstructMatches(smarts_mol))
            results.append(matches)
    return return_formula(input_mol), return_mz(input_mol), results

if __name__ == "__main__":
    # Check if the number of arguments is valid
    if len(sys.argv) < 3:
        sys.exit(1)
    #print("hi from the exe")
    smiles = sys.argv[1]
    smarts = sys.argv[2:]


    #for smiles in smiles_strings:
    formula, mz, results = substruct_matches(smiles, smarts)

    print(formula, mz, results)


#let (mz, formula) = self.get_mz_and_formula_from_smiles().unwrap();
#vec![self.name.clone().replace("\"", "'"), self.smiles.clone(), formula, mz]