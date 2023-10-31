use pyo3::{Python, PyErr, prepare_freethreaded_python};
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
        vec![self.accession.clone(), self.name.clone().replace("\"", "'"), self.smiles.clone(),
             self.formula.clone(), self.mz.to_string()]
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