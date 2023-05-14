
mod parse_hmdb;
mod metabolite;

use std::collections::HashMap;
use std::io::{BufReader, BufRead};
use std::fs::File;
use rusqlite::{Connection, Result};
use pyo3::Python;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;


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
    let file: File = File::open(path).expect("Failed to open file");
    let reader: BufReader<File> = BufReader::new(file);
    let mut my_map: HashMap<String, Vec<String>> = HashMap::new();

    for line in reader.lines() {
        let line: String = line.expect("Failed to read line");
        let parts: Vec<_> = line.split("\t").collect();

        if parts.len() == 3 {

            my_map.insert(parts[0].to_string(), vec![parts[1].to_string(), parts[2].to_string()]);
        }

    }
    my_map
}

fn get_derivs(deriv: &str, fgh: &HashMap<String, String>) -> usize{
    let ret: usize = match fgh.get(&deriv.to_string()) {
        Some(num) => num.parse::<usize>().unwrap(),
        None => 0 as usize
    };
    ret
}

fn get_fgh(input: &str, fgh: &HashMap<String, String>) -> String {
    match fgh.get(&input.to_string()) {
        Some(x) => x.to_string(),
        None => "0".to_string()
    }
}

fn make_metabolite_db(conn: &Connection) -> Result<(), Box<dyn std::error::Error + Send>> {
    let functional_smarts: HashMap<String, String> = import_tsv("database/smarts.tsv");
    let metabolites: Vec<metabolite::Metabolite> = parse_hmdb::parse_xml("database/hmdb_metabolites.xml".to_string()).unwrap();
    let csf_metabolites: Vec<String> = parse_hmdb::parse_only_accession(&"database/csf_metabolites.xml".to_string()).unwrap();
    let urine_metabolites: Vec<String> = parse_hmdb::parse_only_accession(&"database/urine_metabolites.xml".to_string()).unwrap();
    let serum_metabolites: Vec<String> = parse_hmdb::parse_only_accession(&"database/serum_metabolites.xml".to_string()).unwrap();

    let (vecs, maps): (Vec<Vec<String>>, Vec<HashMap<String, String>>) = metabolites.par_iter()
        .map(|metabolite: &metabolite::Metabolite| {
            let mets: Vec<String> = metabolite.get();
            let fgh: HashMap<String, String> = metabolite.functional_group(&functional_smarts);
            (mets, fgh)
        }).unzip();

    conn.execute(
        "CREATE TABLE metabolites (accession TEXT NOT NULL PRIMARY KEY, name TEXT NOT NULL, smiles TEXT NOT NULL, 
                                        chemicalformula TEXT NOT NULL, mz TEXT NOT NULL, csf TEXT NOT NULL, urine TEXT NOT NULL, serum TEXT NOT NULL,
                                        endogenous TEXT NOT NULL, exogenous TEXT NOT NULL, unspecified TEXT NOT NULL,
                                        phenols TEXT NOT NULL, catechols TEXT NOT NULL, 
                                        carbonyls TEXT NOT NULL, aldehydes TEXT NOT NULL, carboxylicacids TEXT NOT NULL, primaryamines TEXT NOT NULL, 
                                        triols TEXT NOT NULL, diols TEXT NOT NULL, hydroxyls TEXT NOT NULL,
                                        fmp integer, dpp integer, tahs integer, ca integer, girardp integer, girardt integer, ampp integer, boronicacid integer )",
            [],
    ).unwrap();

    for i in 0..metabolites.len() {
        let mets: Vec<String> = vecs[i].clone();
        let fgh: HashMap<String, String> = maps[i].clone();

        let in_csf: i32 = if csf_metabolites.contains(&mets[0]) { 1 } else { 0 };
        let in_urine: i32 = if urine_metabolites.contains(&mets[1]) { 1 } else { 0 };
        let in_serum: i32 = if serum_metabolites.contains(&mets[2]) { 1 } else { 0 };
        let endo: i32 = if metabolites[i].endo_exo.contains(&"Endogenous".to_string()) { 1 } else { 0 };
        let exo: i32 = if metabolites[i].endo_exo.contains(&"Exogenous".to_string()) { 1 } else { 0 };
        let unspec: i32 = if metabolites[i].endo_exo.len() == 0 { 1 } else { 0 };

        let fmp10: usize = get_derivs(&"phenols", &fgh) + get_derivs(&"primary amines", &fgh);
        let dpp_tahs_ca: usize = get_derivs(&"primary amines", &fgh);
        let girard: usize = get_derivs(&"carbonyls", &fgh);
        let ampp: usize = get_derivs(&"carboxylic acids", &fgh) + get_derivs(&"aldehydes", &fgh);
        let ba: usize = get_derivs(&"diols", &fgh) + get_derivs(&"catechols", &fgh);
        let a: [&String; 28] = [&mets[0], &mets[1], &mets[2], &mets[3], &mets[4],
                                &in_csf.to_string(), &in_urine.to_string(), &in_serum.to_string(), 
                                &endo.to_string(), &exo.to_string(), &unspec.to_string(),
                                &get_fgh(&"phenols", &fgh), &get_fgh(&"catechols", &fgh), &get_fgh(&"carbonyls", &fgh), 
                                &get_fgh(&"aldehydes", &fgh), &get_fgh(&"carboxylic acids", &fgh), &get_fgh(&"primary amines", &fgh), 
                                &get_fgh(&"triols", &fgh), &get_fgh(&"diols", &fgh), &get_fgh(&"hydroxyls", &fgh),
                                &fmp10.to_string(), &dpp_tahs_ca.to_string(), &dpp_tahs_ca.to_string(), &dpp_tahs_ca.to_string(), 
                                &girard.to_string(), &girard.to_string(), &ampp.to_string(), &ba.to_string(),
                                ];    
        let mut stmt: rusqlite::Statement = conn.prepare("INSERT INTO metabolites (accession, name, smiles, chemicalformula, mz, csf, urine, serum, endogenous, exogenous, unspecified,
            phenols, catechols, carbonyls, aldehydes, carboxylicacids, primaryamines, triols, diols, hydroxyls,
            fmp, dpp, tahs, ca, girardp, girardt, ampp, boronicacid) 
            VALUES (:accession, :name, :smiles, :chemicalformula, :mz, :csf, :urine, :serum, :endogenous, :exogenous, :unspecified,
            :phenols, :catechols, :carbonyls, :aldehydes, :carboxylicacids, :primaryamines, :triols, :diols, :hydroxyls,
            :fmp, :dpp, :tahs, :ca, :girardp, :girardt, :ampp, :boronicacid)").unwrap();

        stmt.execute(&a).unwrap();
    }
    Ok(())
}

fn get_lipid_attr(mol: &pyo3::PyAny, arg: &str) -> String{
    let a: String;
    let res: std::result::Result<&pyo3::PyAny, pyo3::PyErr> = mol.call_method1("GetProp", (arg, ));
    if res.is_ok() {a = res.unwrap().extract::<String>().unwrap();} else {a = "".to_string();}
    a
}

fn make_lipid_db(conn: &Connection) -> Result<()> {
    conn.execute(
        "CREATE TABLE lipids (name TEXT NOT NULL, mz TEXT NOT NULL, formula TEXT NOT NULL, smiles TEXT NOT NULL)",
         [],
    )?;

    pyo3::prepare_freethreaded_python();
    Python::with_gil(|py: Python| {
        let filename: &str = "database/structures.sdf";
        let rdkit: &pyo3::types::PyModule = py.import("rdkit.Chem").unwrap();
        let desc: &pyo3::types::PyModule = py.import("rdkit.Chem.Descriptors").unwrap();
        let supplier: &pyo3::PyAny = rdkit.call_method1("SDMolSupplier", (filename,)).unwrap();

        for _ in 0..supplier.len().unwrap() {
            let mol: &pyo3::PyAny = supplier.call_method1("__next__", ()).unwrap();
            let mz: String;

            let name: String = get_lipid_attr(mol, "NAME");
            let formula: String = get_lipid_attr(mol, "FORMULA");
            let smiles: String = get_lipid_attr(mol, "SMILES");

            let mass: std::result::Result<&pyo3::PyAny, pyo3::PyErr> = desc.call_method1("ExactMolWt", (mol, ));
            if mass.is_ok() {mz = mass.unwrap().extract::<f64>().unwrap().to_string();} else {mz = "".to_string();}

            let a: [&String; 4] = [&name, &mz, &formula, &smiles];
            
            let mut stmt: rusqlite::Statement = conn.prepare("INSERT INTO lipids (name, mz, formula, smiles) 
                                                                            VALUES (:name, :mz, :formula, :smiles)").unwrap();
            stmt.execute(&a).unwrap();
        
        }
        
    });
    
    Ok(())
    
}

fn make_matrix_db(conn: &Connection) -> Result<()> {
    let matrices: HashMap<String, Vec<String>> = import_tsv2("database/matrix.tsv");
    conn.execute(
        "CREATE TABLE matrices (mname TEXT NOT NULL, formula TEXT NOT NULL, deltamass TEXT NOT NULL)",
        [],
    )?;

    for (key, value) in matrices {
        let a: [&String; 3] = [&key, &value[0], &value[1]];    
        let mut stmt: rusqlite::Statement = conn.prepare(
            "INSERT INTO matrices (mname, formula, deltamass) VALUES (:mname, :formula, :deltamass)").unwrap();
        stmt.execute(&a).unwrap();
    }

    Ok(())
}

fn main() -> () {
    let db_path = "db.db";
    //remove the database if it already exists
    std::fs::remove_file(db_path).expect("Database could not be removed");

    let conn: Connection = Connection::open(db_path).expect("Could not open db");
    make_metabolite_db(&conn).unwrap();
    make_lipid_db(&conn).unwrap();
    make_matrix_db(&conn).unwrap();
}
