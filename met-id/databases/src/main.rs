
mod parse_hmdb;
mod metabolite;
use std::thread;

use std::collections::HashMap;
use std::io::{BufReader, BufRead};
use std::fs::File;
use std::sync::{Arc, Mutex};


use rusqlite::{Connection, Result};
use std::error::Error as StdError;
use std::fmt::{self, Display};

use crossbeam_channel::{unbounded, Sender, Receiver};
use pyo3::{Python, PyAny};
//use pyo3::types::PyModule;
//use pyo3::types::PyTuple;
//use pyo3::PyResult;

#[derive(Debug)]
struct MyError(String);

impl Display for MyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl StdError for MyError {}


fn import_tsv(path: &str) -> HashMap<String, String> {
    let file: File = File::open(path).expect("Failed to open file");
    let reader: BufReader<File> = BufReader::new(file);

    let mut my_map: HashMap<String, String> = HashMap::new();

    for line in reader.lines() {
        let line: String = line.expect("Failed to read line");
        let parts: Vec<_> = line.split("\t").collect();

        if parts.len() == 2 {
            my_map.insert(parts[0].to_string(), parts[1].to_string());
        }
    }
    my_map
}

fn import_tsv2(path: &str) -> HashMap<String, Vec<String>> {
    let file = File::open(path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut my_map: HashMap<String, Vec<String>> = HashMap::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let parts: Vec<_> = line.split("\t").collect();

        if parts.len() == 3 {

            my_map.insert(parts[0].to_string(), vec![parts[1].to_string(), parts[2].to_string()]);
        }

    }
    my_map
}

fn remove_db(db_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    if std::fs::metadata(db_path).is_ok() {
        // If it does, remove it
        std::fs::remove_file(db_path)?;
    }
    else {
        println!("Database not deleted");
    }
    Ok(())
}

fn make_metabolite_db(conn: &Connection) -> Result<(), Box<dyn std::error::Error + Send>> {
    let functional_smarts: HashMap<String, String> = import_tsv("database/smarts.tsv");
    let metabolites: Vec<metabolite::Metabolite> = parse_hmdb::parse_xml("database/hmdb_metabolites_reduced.xml".to_string()).unwrap();
    let arc_metabolites: Arc<Vec<Arc<metabolite::Metabolite>>> = Arc::new(metabolites.into_iter().map(Arc::new).collect::<Vec<_>>());

    conn.execute(
        "CREATE TABLE metabolites (id INTEGER PRIMARY KEY, accession TEXT NOT NULL, name TEXT NOT NULL, mz TEXT NOT NULL, 
                                       chemicalformula TEXT NOT NULL, smiles TEXT NOT NULL, phenols TEXT NOT NULL, catechols TEXT NOT NULL, 
                                       carbonyls TEXT NOT NULL, aldehydes TEXT NOT NULL, carboxylicacids TEXT NOT NULL, primaryamines TEXT NOT NULL, 
                                       triols TEXT NOT NULL, diols TEXT NOT NULL, hydroxyls TEXT NOT NULL)",
         [],
    ).unwrap();
     

    let num_threads = 8;
    let mut mets_vec: Vec<Vec<String>> = Vec::new();
    let mut fgh_vec: Vec<HashMap<String, String>> = Vec::new();

    let chunk_size = (arc_metabolites.len() + num_threads - 1) / num_threads;
    let (tx, rx): (Sender<(Vec<String>, HashMap<String, String>)>, Receiver<(Vec<String>, HashMap<String, String>)>) = unbounded();

    let mut thread_handles = Vec::new();


    for i in 0..num_threads{
        let chunk_start: usize = i * chunk_size;
        let chunk_end: usize = std::cmp::min((i + 1) * chunk_size, arc_metabolites.len());
        let chunk: Arc<Vec<Arc<metabolite::Metabolite>>> = Arc::clone(&arc_metabolites);
        let chunk_mutex: Arc<Mutex<Arc<Vec<Arc<metabolite::Metabolite>>>>> = Arc::new(Mutex::new(chunk));
        let fnsmarts: HashMap<String, String> = functional_smarts.clone();
        let tx_clone = tx.clone();
        
        
        let thread_handle = thread::spawn(move || {
            let chunk: std::sync::MutexGuard<Arc<Vec<Arc<metabolite::Metabolite>>>> = chunk_mutex.lock().unwrap();
            for item in &chunk[chunk_start..chunk_end] {
                let mets: Vec<String> = item.get();
                let fs_: HashMap<String, String> = fnsmarts.clone();
                let fgh: HashMap<String, String> = item.functional_group(&fs_);
                tx_clone.send((mets, fgh)).unwrap();
            }
        });
        thread_handles.push(thread_handle);
    }
    for thread_handle in thread_handles {
        thread_handle.join().unwrap();
    }

    while mets_vec.len() < arc_metabolites.len(){
        let mets: Vec<String>;
        let fgh: HashMap<String, String>;
        let rx_clone = rx.clone();
        
        match rx_clone.recv() {
            Ok(message) => {
                (mets, fgh) = message;
                mets_vec.push(mets);
                fgh_vec.push(fgh);
            }
            Err(_e) => {},
        }
    }

    for i in 0..arc_metabolites.len() {
        let mets = mets_vec[i].clone();
        let fgh = fgh_vec[i].clone();
        let a: [&String; 14] = [&mets[0], &mets[1], &mets[2], &mets[3], &mets[4],
                                    &fgh.get(&"phenols".to_string()).unwrap(), &fgh.get(&"catechols".to_string()).unwrap(), 
                                    &fgh.get(&"carbonyls".to_string()).unwrap(), &fgh.get(&"aldehydes".to_string()).unwrap(), 
                                    &fgh.get(&"carboxylic acids".to_string()).unwrap(), &fgh.get(&"primary amines".to_string()).unwrap(), 
                                    &fgh.get(&"triols".to_string()).unwrap(), &fgh.get(&"diols".to_string()).unwrap(), 
                                    &fgh.get(&"hydroxyls".to_string()).unwrap()];
        
        let mut stmt: rusqlite::Statement = conn.prepare("INSERT INTO metabolites (accession, name, mz, chemicalformula, smiles, 
            phenols, catechols, carbonyls, aldehydes, carboxylicacids, primaryamines, triols, diols, hydroxyls) 
            VALUES (:accession, :name, :mz, :chemicalformula, :smiles,
            :phenols, :catechols, :carbonyls, :aldehydes, :carboxylicacids, :primaryamines, :triols, :diols, :hydroxyls)")
            .map_err(|err| MyError(format!("Error: {}", err))).unwrap();

        stmt.execute(&a).unwrap();
    }

    Ok(())
}

fn make_lipid_db(conn: &Connection) -> Result<()> {
    conn.execute(
        "CREATE TABLE lipids (
             id INTEGER PRIMARY KEY,
             name TEXT NOT NULL,
             mz TEXT NOT NULL,
             formula TEXT NOT NULL,
             smiles TEXT NOT NULL
         )",
         [],
    )?;

    pyo3::prepare_freethreaded_python();
    Python::with_gil(|py: Python| {
        let filename = "database/structures.sdf";
        let rdkit = py.import("rdkit.Chem").unwrap();
        let desc = py.import("rdkit.Chem.Descriptors").unwrap();
        let supplier = rdkit.call_method1("SDMolSupplier", (filename,)).unwrap().cast_as::<PyAny>().unwrap();

        for _ in 0..supplier.len().unwrap() {
            let mol = supplier.call_method1("__next__", ()).unwrap().cast_as::<PyAny>().unwrap();
            let name: String;
            let formula: String;
            let smiles: String;
            let mz: String;

            let result_name = mol.call_method1("GetProp", ("NAME", ));
            if result_name.is_ok() {name = result_name.unwrap().extract::<String>().unwrap();} else {name = "".to_string();}

            let result_formula = mol.call_method1("GetProp", ("FORMULA", ));
            if result_formula.is_ok() {formula = result_formula.unwrap().extract::<String>().unwrap();} else {formula = "".to_string();}

            let result_smiles = mol.call_method1("GetProp", ("SMILES", ));
            if result_smiles.is_ok() {smiles = result_smiles.unwrap().extract::<String>().unwrap();} else {smiles = "".to_string()}

            let mass = desc.call_method1("ExactMolWt", (mol, ));
            if mass.is_ok() {mz = mass.unwrap().extract::<f64>().unwrap().to_string();} else {mz = "".to_string();}


            let a: [&String; 4] = [&name, &mz, &formula, &smiles];
            
            let mut stmt: rusqlite::Statement = conn.prepare("INSERT INTO lipids (name, mz, formula, smiles) 
                                                                            VALUES (:name, :mz, :formula, :smiles)")
                .map_err(|err| MyError(format!("Error: {}", err))).unwrap();

            stmt.execute(&a).unwrap();
        
        }
        
    });
    
    Ok(())
}

fn make_matrix_db(conn: &Connection) -> Result<()> {
    let matrices = import_tsv2("database/matrix.tsv");
    conn.execute(
        "CREATE TABLE matrices (
             id INTEGER PRIMARY KEY,
             name TEXT NOT NULL,
             formula TEXT NOT NULL,
             deltamass TEXT NOT NULL
         )",
         [],
    )?;

    for (key, value) in matrices {
        let a: [&String; 3] = [&key, &value[0], &value[1]];    
        let mut stmt: rusqlite::Statement = conn.prepare("INSERT INTO matrices (name, formula, deltamass) 
                                                               VALUES (:name, :formula, :deltamass)")
                                .map_err(|err| MyError(format!("Error: {}", err))).unwrap();
        stmt.execute(&a).unwrap();
    }

    Ok(())
}

fn main() -> () {
    let db_path = "db.db";
    remove_db(db_path).unwrap();
    let conn = Connection::open(db_path).map_err(|err| MyError(format!("Error: {}", err))).unwrap();
    make_metabolite_db(&conn).unwrap();
    make_lipid_db(&conn).unwrap();
    make_matrix_db(&conn).unwrap();
}
