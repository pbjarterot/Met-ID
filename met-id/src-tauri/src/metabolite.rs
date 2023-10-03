use pyo3::{Python, PyErr, prepare_freethreaded_python, types::PyModule, PyResult};
//use pyo3::types::{PyDict, PyString};
//use pyo3::prelude::*;
//use pyo3::types::IntoPyDict;
use std::collections::HashMap;
use rayon::prelude::*;


#[derive(Debug, Clone)]
pub struct Metabolite {
    pub name: String,
    pub smiles: String,
}

impl Metabolite {
    pub fn get(&self) -> Vec<String> {
        let (mz, formula) = self.get_mz_and_formula_from_smiles().unwrap();
        vec![self.name.clone().replace("\"", "'"), self.smiles.clone(), formula, mz]
    }

    fn get_mz_and_formula_from_smiles(&self) -> Result<(String, String), PyErr> {
        let smiles = &self.smiles;

        Python::with_gil(|py| {
            let rdkit = PyModule::import(py, "rdkit.Chem")?;
            let descriptors = PyModule::import(py, "rdkit.Chem.Descriptors")?;
            let formula = PyModule::import(py, "rdkit.Chem.rdMolDescriptors")?;

            let mol = rdkit.call_method1("MolFromSmiles", (smiles,))?;

            if mol.is_none() {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid SMILES string"));
            }

            let mz: f64 = descriptors.call_method1("ExactMolWt", (mol.clone(),))?.extract()?;
            let mol_formula: String = formula.call_method1("CalcMolFormula", (mol.clone(),))?.extract()?;


            Ok((mz.to_string(), mol_formula))
        })
    }

    pub fn functional_group(&self, functional_smarts: &HashMap<String, String>) -> Result<HashMap<String, String>, PyErr> {
        prepare_freethreaded_python();
        let smils = &self.smiles;

        let result: Result<Vec<Option<(usize, String)>>, PyErr> = functional_smarts
            .par_iter()
            .map(|(key, value)| {
                Python::with_gil(|py| {
                    let rdkit = py.import("rdkit.Chem")?;
                    let mol = rdkit.getattr("MolFromSmiles")?.call1((smils,))?;
                    if !mol.is_none() {
                        let func = rdkit.getattr("MolFromSmarts")?.call1((value,))?;
                        let a = mol.call_method1("GetSubstructMatches", (func,))?;
                        let aa: Vec<Vec<usize>> = a.extract()?;
                        let tuple_len = aa.len();
                        Ok(Some((tuple_len, key.clone())))
                    } else {
                        Ok(None)
                    }
                })
            })
            .collect();

        match result {
            Ok(mut results) => {
                let mut fgh: HashMap<String, String> = HashMap::new();
                results.drain(..).filter_map(|x| x).for_each(|(u, s)| {
                    fgh.insert(s, u.to_string());
                });
                Ok(fgh)
            }
            Err(e) => Err(e),
        }
    }

}


pub fn single_functional_group(smiles_list: &Vec<String>, functional_smarts: &str) -> Vec<usize> {
    prepare_freethreaded_python();

    const BATCH_SIZE: usize = 10; 

    Python::with_gil(|py| -> PyResult<Vec<usize>> {
        let rdkit = py.import("rdkit.Chem")?;
        let smarts_to_mol = rdkit.getattr("MolFromSmarts")?;
        let func = smarts_to_mol.call1((functional_smarts,))?;  // Created only once

        let mol_from_smiles = rdkit.getattr("MolFromSmiles")?;  // Created only once

        let mut results = Vec::with_capacity(smiles_list.len()); // Pre-allocate the vector

        for chunk in smiles_list.chunks(BATCH_SIZE) {
            for smiles in chunk {
                match mol_from_smiles.call1((smiles,)) {
                    Ok(ref actual_mol) if !actual_mol.is_none() => {
                        let matches = actual_mol.call_method1("GetSubstructMatches", (func,))?;
                        let matches_count: Vec<Vec<usize>> = matches.extract()?;
                        results.push(matches_count.len());
                    },
                    _ => results.push(0),
                }
            }
        }
        Ok(results)
    })
    .unwrap_or_else(|_| vec![]) // Handle errors by returning an empty vec
}
