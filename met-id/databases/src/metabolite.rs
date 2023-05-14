use pyo3::{Python, prepare_freethreaded_python};
use pyo3::types::PyModule;
use std::collections::HashMap;
use rayon::prelude::*;


#[derive(Debug, Clone)]
pub struct Metabolite {
    pub accession: String,
    pub name: String,
    pub smiles: String,
    pub formula: String,
    pub mz: f64,
    pub endo_exo: Vec<String>
}

impl Metabolite {
    pub fn get(&self) -> Vec<String> {
        vec![self.accession.clone(), self.name.clone(), self.smiles.clone(),
             self.formula.clone(), self.mz.to_string()]
    }

    pub fn functional_group(&self, functional_smarts: &HashMap<String, String>) -> HashMap<String, String> {
        let smils: String = self.smiles.clone();
        let smarts_: HashMap<String, String> = functional_smarts.clone();

        let result: Vec<Option<(usize, String)>> = smarts_
            .par_iter()
            .map(|(key, value)| {
                let value_clone = value.clone();
                let key_clone = key.clone();
                let smiles_clone = smils.clone();
                prepare_freethreaded_python();
                Python::with_gil(|py| {
                    let rdkit: &PyModule = py.import("rdkit.Chem").unwrap();
                    let mol: &pyo3::PyAny = rdkit.getattr("MolFromSmiles").unwrap().call1((&smiles_clone,)).unwrap();
                    if !mol.is_none() {
                        let func: &pyo3::PyAny = rdkit.getattr("MolFromSmarts").unwrap().call1((&value_clone,)).unwrap();
                        let a: &pyo3::PyAny = mol.call_method1("GetSubstructMatches", (func,)).unwrap();
                        let aa: Vec<Vec<usize>> = a.extract().unwrap();
                        let tuple_len = aa.len();

                        Some((tuple_len, key_clone))
                    } else {
                        None
                    }
                })
            })
            .collect();

        let mut fgh: HashMap<String, String> = HashMap::new();
        for index in 0..functional_smarts.values().len() {
            //let (length, group): (usize, String) = &result[index];
            match &result[index] {
                Some((u, s)) => fgh.insert(s.clone(), u.to_string()),
                None => None,
            };

        }

        fgh
    }
    
}