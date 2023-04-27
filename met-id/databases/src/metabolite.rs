use pyo3::{Python, PyAny};
use pyo3::types::PyModule;
use std::collections::HashMap;
use std::sync::mpsc;
use pyo3::types::PyTuple;

#[derive(Debug, Clone)]
pub struct Metabolite {
    pub accession: String,
    pub name: String,
    pub smiles: String,
    pub formula: String,
    pub mz: f64
}

impl Metabolite {
    pub fn get(&self) -> Vec<String> {
        vec![self.accession.clone(), self.name.clone(), self.smiles.clone(),
             self.formula.clone(), self.mz.to_string()]
    }

    pub fn functional_group(&self, functional_smarts: &HashMap<String, String>) -> HashMap<String, String> {
        let (tx, rx) = mpsc::channel();
        let smils = self.smiles.clone();
        let smarts_ = functional_smarts.clone();

        for (key, value) in &smarts_ {
            let tx1 = tx.clone();
            let value_clone = value.clone();
            let key_clone = key.clone();
            let smiles = smils.clone();
            
            pyo3::prepare_freethreaded_python();
            Python::with_gil(|py: Python| {
            let rdkit: &PyModule = py.import("rdkit.Chem").unwrap();
            let mol: &PyAny = rdkit.getattr("MolFromSmiles").unwrap().call((smiles, ), None).unwrap();
            
            let tuple_len: usize;
            if !mol.is_none() {
                let func: &PyAny = rdkit.getattr("MolFromSmarts").unwrap().call((value_clone, ), None).unwrap();
                let a: &PyAny = mol.getattr("GetSubstructMatches").unwrap().call((func, ), None).unwrap();
                let aa: &PyTuple = a.extract().unwrap();
                tuple_len = aa.len();
            } else {
                tuple_len = 0 as usize;
            }
            
            tx1.send((tuple_len, key_clone)).unwrap();
            });
        };
        
        let mut fgh: HashMap<String, String> = HashMap::new();
        for _ in functional_smarts.values() {
            let (length, group): (usize, String) = rx.recv().unwrap();
            fgh.insert(group, length.to_string());
        };
    
        fgh
    }
}