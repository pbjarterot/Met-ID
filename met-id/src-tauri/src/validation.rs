use pyo3::{prelude::*, prepare_freethreaded_python};
use pyo3::types::{PyDict, PyString};


#[tauri::command]
pub fn check_smiles(smiles: String) -> bool {
    prepare_freethreaded_python();
    let mut val: bool = true;
    Python::with_gil(|py| {

        let locals: &PyDict = PyDict::new(py);
        locals.set_item("code", PyString::new(py, r#"
from rdkit import Chem

def is_valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

result = is_valid_smiles(smiles)    
"#)).unwrap();
    locals.set_item("smiles", smiles).unwrap();

    //execute the function
    if let Err(e) = py.run("exec(code)", Some(locals), None) {
        e.print(py);
    }

    // Retrieve the result from the locals dictionary
    let result: &PyAny = locals.get_item("result").unwrap();
    val = result.extract().unwrap();
    });


    val
}

#[tauri::command]
pub fn check_smarts(smarts: String) -> bool {
    prepare_freethreaded_python();
    let mut val: bool = true;
    Python::with_gil(|py| {

        let locals: &PyDict = PyDict::new(py);
        locals.set_item("code", PyString::new(py, r#"
from rdkit import Chem

def is_valid_smarts(smarts):
    mol = Chem.MolFromSmarts(smarts)
    return mol is not None

result = is_valid_smarts(smarts)    
"#)).unwrap();
    locals.set_item("smarts", smarts).unwrap();

    //execute the function
    if let Err(e) = py.run("exec(code)", Some(locals), None) {
        e.print(py);
    }

    // Retrieve the result from the locals dictionary
    let result: &PyAny = locals.get_item("result").unwrap();
    val = result.extract().unwrap();
    });

    val
}

/* 
pub fn mass_from_smiles(smiles: &String) -> (String, f64) {
    let mut monoisotopic_mass: f64 = 0.0;
    let mut formula: String = "".to_string();

    Python::with_gil(|py| {
        let code = r#"
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

molecule = Chem.MolFromSmiles(smiles)
if molecule is not None:
    monoisotopic_mass = Descriptors.ExactMolWt(molecule)
    formula = CalcMolFormula(molecule)
else:
    monoisotopic_mass = 0.0
    formula = ""
"#;
        let locals = PyDict::new(py);
        locals.set_item("smiles", smiles).unwrap();

        // Execute the Python code
        if let Err(e) = py.run(code, Some(locals), None) {
            e.print(py);
            formula = "".to_string();
            monoisotopic_mass = 0.0;
        } else {
            // Retrieve the monoisotopic mass from the locals dictionary
            if let Ok(result) = locals.get_item("monoisotopic_mass").unwrap().extract::<f64>() {
                monoisotopic_mass = result;
            }

            // Retrieve the formula from the locals dictionary
            if let Ok(result) = locals.get_item("formula").unwrap().extract::<String>() {
                formula = result;
            }
        }
    });

    (formula, monoisotopic_mass)
}
*/

pub fn mass_from_formula(formula: &String) -> f64 {
    let mut molecular_weight: f64 = 0.0;

    Python::with_gil(|py| {
        let code = r#"
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


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

try:
    molecular_weight = get_mass(formula)
except Exception as e:
    print(e)
    molecular_weight = 0.0
"#;
        let locals = PyDict::new(py);
        locals.set_item("formula", formula).unwrap();

        // Execute the Python code
        if let Err(e) = py.run(code, Some(locals), None) {
            e.print(py);
            molecular_weight = 0.0;
        } else {
            // Retrieve the molecular weight from the locals dictionary
            if let Ok(result) = locals.get_item("molecular_weight").unwrap().extract::<f64>() {
                molecular_weight = result;
            }
        }
    });

    molecular_weight
}