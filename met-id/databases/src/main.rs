
mod parse_hmdb;
mod metabolite;
//mod draw_molecules;


use rusqlite::{ Batch, Connection };
use std::collections::{ HashMap, HashSet};
use std::io::{BufReader, BufRead};
use std::fs::File;
//use std::io::Write;

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

        if parts.len() == 5 {
            my_map.insert(parts[0].to_string(), vec![parts[1].to_string(), parts[2].to_string(), parts[3].to_string(), parts[4].to_string()]);
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

fn get_fgh(input: &str, fgh: &HashMap<String, String>) -> usize {
    match fgh.get(&input.to_string()) {
        Some(x) => x.parse::<usize>().unwrap(),
        None => 0
    }
}

fn create_tables(conn: &Connection) -> () {
    conn.execute("CREATE TABLE metabolites (id INTEGER PRIMARY KEY AUTOINCREMENT, 
        name TEXT, 
        smiles TEXT, 
        chemicalformula TEXT, 
        mz TEXT)", 
    []
    ).unwrap();

    /* 
    conn.execute("CREATE TABLE functional_groups (id INTEGER PRIMARY KEY AUTOINCREMENT, 
                phenols INTEGER NOT NULL, 
                catechols INTEGER NOT NULL, 
                carbonyls INTEGER NOT NULL, 
                aldehydes INTEGER NOT NULL, 
                carboxylicacids INTEGER NOT NULL, 
                primaryamines INTEGER NOT NULL,
                triols INTEGER NOT NULL, 
                diols INTEGER NOT NULL, 
                hydroxyls INTEGER NOT NULL)",
                []
    ).unwrap();
    */

    conn.execute("CREATE TABLE functional_groups (id INTEGER PRIMARY KEY AUTOINCREMENT, 
        phenols INTEGER NOT NULL,  
        aldehydes INTEGER NOT NULL, 
        carboxylicacids INTEGER NOT NULL, 
        primaryamines INTEGER NOT NULL)",
        []
    ).unwrap();

    conn.execute("CREATE TABLE derivatized_by (id INTEGER PRIMARY KEY AUTOINCREMENT,
            fmp INTEGER NOT NULL, 
            ampp INTEGER NOT NULL)",
    []
    ).unwrap();
    /* 
    conn.execute("CREATE TABLE derivatized_by (id INTEGER PRIMARY KEY AUTOINCREMENT,
        fmp INTEGER NOT NULL, 
        dpp INTEGER NOT NULL, 
        tahs INTEGER NOT NULL, 
        ca INTEGER NOT NULL, 
        girardp INTEGER NOT NULL, 
        girardt INTEGER NOT NULL, 
        ampp INTEGER NOT NULL, 
        boronicacid INTEGER NOT NULL)",
    []
    ).unwrap();
    */

    conn.execute("CREATE TABLE in_tissue (id INTEGER PRIMARY KEY AUTOINCREMENT,
        csf INTEGER NOT NULL, 
        urine INTEGER NOT NULL, 
        serum INTEGER NOT NULL)", 
    []
    ).unwrap();

    conn.execute("CREATE TABLE endogeneity (id INTEGER PRIMARY KEY AUTOINCREMENT,
        endogenous INTEGER NOT NULL, 
        exogenous INTEGER NOT NULL, 
        unspecified INTEGER NOT NULL)", 
    []
    ).unwrap();

    conn.execute("CREATE TABLE db_accessions (id INTEGER PRIMARY KEY AUTOINCREMENT,
                    hmdb TEXT)", 
    []
    ).unwrap();

    ()
}


fn make_metabolite_db(conn: &Connection) -> Result<(), Box<dyn std::error::Error + Send>> {
    create_tables(conn);


    let functional_smarts: HashMap<String, String> = import_tsv("database/smarts2.tsv");
    let metabolites: Vec<metabolite::Metabolite> = parse_hmdb::parse_xml("database/hmdb_metabolites.xml".to_string()).unwrap();
    let csf_metabolites: Vec<String> = parse_hmdb::parse_only_accession(&"database/csf_metabolites.xml".to_string()).unwrap();
    let urine_metabolites: Vec<String> = parse_hmdb::parse_only_accession(&"database/urine_metabolites.xml".to_string()).unwrap();
    let serum_metabolites: Vec<String> = parse_hmdb::parse_only_accession(&"database/serum_metabolites.xml".to_string()).unwrap();

    let (vecs, maps): (Vec<Vec<String>>, Vec<HashMap<String, String>>) = metabolites.par_iter()
        .map(|metabolite: &metabolite::Metabolite| {
            let mets: Vec<String> = metabolite.get();
            let fgh: HashMap<String, String> = metabolite.functional_group(&functional_smarts).unwrap();
            (mets, fgh)
        }).unzip();

    let csf_metabolites_set: HashSet<_> = csf_metabolites.into_iter().collect();
    let urine_metabolites_set: HashSet<_> = urine_metabolites.into_iter().collect();
    let serum_metabolites_set: HashSet<_> = serum_metabolites.into_iter().collect();

    let endogenous_str = "Endogenous".to_string();
    let exogenous_str = "Exogenous".to_string();

    let metabolites_len = metabolites.len();

    // Building the SQL statement dynamically
    let mut sql = String::from("BEGIN;\n");

    let total: u64 = metabolites_len.try_into().unwrap();

    for i in 0..metabolites_len {
        //println!("{:?}, {:?}", i, metabolites_len);
        let mets: &Vec<String> = &vecs[i];
        let fgh: &HashMap<String, String> = &maps[i];

        let in_csf: i32 = csf_metabolites_set.contains(&mets[0]) as i32;
        let in_urine: i32 = urine_metabolites_set.contains(&mets[0]) as i32;
        let in_serum: i32 = serum_metabolites_set.contains(&mets[0]) as i32;
        let endo: i32 = metabolites[i].endo_exo.contains(&endogenous_str) as i32;
        let exo: i32 = metabolites[i].endo_exo.contains(&exogenous_str) as i32;
        let unspec: i32 = metabolites[i].endo_exo.is_empty() as i32;
        if endo > 0 {
            println!("{:?}, {:?}, {:?}", endo, exo, unspec);
        }
        

        let fmp10: i32 = (get_derivs(&"phenols", fgh) + get_derivs(&"primary amines", fgh)).try_into().unwrap();
        //let dpp_tahs_ca:i32 = get_derivs(&"primary amines", fgh).try_into().unwrap();
        //let girard: i32 = get_derivs(&"carbonyls", fgh).try_into().unwrap();
        let ampp: i32 = (get_derivs(&"carboxylic acids", fgh) + get_derivs(&"aldehydes", fgh)).try_into().unwrap();
        //let ba: i32 = (get_derivs(&"diols", fgh) + get_derivs(&"catechols", fgh)).try_into().unwrap();

        // Append SQL statements to the string
        sql.push_str(&format!(
            "INSERT INTO metabolites (name, smiles, chemicalformula, mz) VALUES (\"{}\", \"{}\", \"{}\", \"{}\");\n",
            mets[1], mets[2], mets[3], mets[4]
        ));

        sql.push_str(&format!(
            "INSERT INTO functional_groups (phenols, aldehydes, carboxylicacids, primaryamines) VALUES ({}, {}, {}, {});\n",
            get_fgh("phenols", fgh),
            get_fgh("aldehydes", fgh),
            get_fgh("carboxylic acids", fgh),
            get_fgh("primary amines", fgh),
        ));

        /* 
        sql.push_str(&format!(
            "INSERT INTO functional_groups (phenols, catechols, carbonyls, aldehydes, carboxylicacids, primaryamines, triols, diols, hydroxyls) VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {});\n",
            get_fgh("phenols", fgh),
            //get_fgh("catechols", fgh),
            //get_fgh("carbonyls", fgh),
            get_fgh("aldehydes", fgh),
            get_fgh("carboxylic acids", fgh),
            get_fgh("primary amines", fgh),
            //get_fgh("triols", fgh),
            //get_fgh("diols", fgh),
            //get_fgh("hydroxyls", fgh)
        ));
        */
        sql.push_str(&format!(
            "INSERT INTO derivatized_by (fmp, ampp) VALUES ({}, {});\n",
            fmp10, ampp
        ));

        sql.push_str(&format!(
            "INSERT INTO in_tissue (csf, urine, serum) VALUES ({}, {}, {});\n",
            in_csf, in_urine, in_serum
        ));

        sql.push_str(&format!(
            "INSERT INTO endogeneity (endogenous, exogenous, unspecified) VALUES ({}, {}, {});\n",
            endo, exo, unspec
        ));

        sql.push_str(&format!(
            "INSERT INTO db_accessions (hmdb) VALUES ('{}');\n",
            mets[0]
        ));

    }

    sql.push_str("COMMIT;\n");


    // Prepare and execute the batch
    let mut batch: Batch<'_, '_> = Batch::new(&conn, &sql);

    while let Some(mut stmt) = batch.next().unwrap() {
        //println!("{:?}\n{:?}", &stmt, &sql);
        //println!("new batch");
        stmt.execute([]).unwrap();
    }
    println!("done");
    /* 
    //batch.execute()?;
    conn.execute("CREATE INDEX IF NOT EXISTS index_metabolites_name ON metabolites(name)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_metabolites_smiles ON metabolites(smiles)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_metabolites_chemicalformula ON metabolites(chemicalformula)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_metabolites_mz ON metabolites(mz)", []).unwrap();

    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_phenols ON functional_groups(phenols)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_catechols ON functional_groups(catechols)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_carbonyls ON functional_groups(carbonyls)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_aldehydes ON functional_groups(aldehydes)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_carboxylicacids ON functional_groups(carboxylicacids)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_primaryamines ON functional_groups(primaryamines)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_triols ON functional_groups(triols)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_diols ON functional_groups(diols)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_functional_groups_hydroxyls ON functional_groups(hydroxyls)", []).unwrap();

    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_fmp ON derivatized_by(fmp)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_dpp ON derivatized_by(dpp)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_tahs ON derivatized_by(tahs)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_ca ON derivatized_by(ca)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_girardp ON derivatized_by(girardp)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_girardt ON derivatized_by(girardt)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_ampp ON derivatized_by(ampp)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_derivatized_by_boronicacid ON derivatized_by(boronicacid)", []).unwrap();

    conn.execute("CREATE INDEX IF NOT EXISTS index_in_tissue_csf ON in_tissue(csf)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_in_tissue_urine ON in_tissue(urine)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_in_tissue_serum ON in_tissue(serum)", []).unwrap();

    conn.execute("CREATE INDEX IF NOT EXISTS index_endogeneity_endogenous ON endogeneity(endogenous)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_endogeneity_exogenous ON endogeneity(exogenous)", []).unwrap();
    conn.execute("CREATE INDEX IF NOT EXISTS index_endogeneity_unspecified ON endogeneity(unspecified)", []).unwrap();
    */

    Ok(())
}

fn get_lipid_attr(mol: &pyo3::PyAny, arg: &str) -> String{
    let a: String;
    let res: std::result::Result<&pyo3::PyAny, pyo3::PyErr> = mol.call_method1("GetProp", (arg, ));
    if res.is_ok() {a = res.unwrap().extract::<String>().unwrap();} else {a = "".to_string();}
    a
}

fn make_lipid_db(conn: &Connection) -> Result<(), rusqlite::Error> {
    println!("doing lipids now");
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

fn make_adduct_db(conn: &Connection) -> Result<(), rusqlite::Error> {
    let matrices: HashMap<String, Vec<String>> = import_tsv2("database/adducts.tsv");
    println!("{:?}", matrices);
    conn.execute(
        "CREATE TABLE adducts (adduct TEXT NOT NULL, mname TEXT NOT NULL, numfunctionalgroups TEXT NOT NULL, formula TEXT NOT NULL, deltamass TEXT NOT NULL)",
        [],
    )?;

    for (key, value) in matrices {
        println!("{}, {:?}", key, value);
        let a: [&String; 5] = [&key, &value[0], &value[1], &value[2], &value[3]];    
        let mut stmt: rusqlite::Statement = conn.prepare(
            "INSERT INTO adducts (adduct, mname, numfunctionalgroups, formula, deltamass) VALUES (:adduct, :mname, :numfunctionalgroups, :formula, :deltamass)").unwrap();
        stmt.execute(&a).unwrap();
    }

    Ok(())
}

fn make_matrices_db(conn: &Connection) -> Result<(), rusqlite::Error> {
    let file: File = File::open("database/matrices.tsv").expect("Failed to open file");
    let reader: BufReader<File> = BufReader::new(file);

    let mut my_map: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line: String = line.expect("Failed to read line");
        let parts: Vec<_> = line.split("\t").collect();
        my_map.push(parts[0].to_string());
    }
    conn.execute(
        "CREATE TABLE matrices (matrix TEXT NOT NULL, phenols TEXT NOT NULL, 'primary amines' TEXT NOT NULL, aldehydes TEXT NOT NULL, 'carboxylic acids' TEXT NOT NULL)",
        [],
    )?;

    for matrix in my_map {
        if matrix == "FMP10".to_string() {
            let a: [&String; 5] = [&matrix, &"1".to_string(), &"1".to_string(), &"0".to_string(), &"0".to_string()];    
            let mut stmt: rusqlite::Statement = conn.prepare(
                "INSERT INTO matrices (matrix, phenols, 'primary amines', aldehydes, 'carboxylic acids') VALUES (:matrix, :phenols, :primary_amines, :aldehydes, :carboxylic_acids)").unwrap();
            stmt.execute(&a).unwrap();

        } else if matrix=="AMPP"{
            let a: [&String; 5] = [&matrix, &"0".to_string(), &"0".to_string(), &"1".to_string(), &"1".to_string()];    
            let mut stmt: rusqlite::Statement = conn.prepare(
                "INSERT INTO matrices (matrix, phenols, 'primary amines', aldehydes, 'carboxylic acids') VALUES (:matrix, :phenols, :primary_amines, :aldehydes, :carboxylic_acids)").unwrap();
            stmt.execute(&a).unwrap();
        } else {
            let a: [&String; 5] = [&matrix, &"0".to_string(), &"0".to_string(), &"0".to_string(), &"0".to_string()];    
            let mut stmt: rusqlite::Statement = conn.prepare(
                "INSERT INTO matrices (matrix, phenols, 'primary amines', aldehydes, 'carboxylic acids') VALUES (:matrix, :phenols, :primary_amines, :aldehydes, :carboxylic_acids)").unwrap();
            stmt.execute(&a).unwrap();
        }
            
    }

    Ok(())
}

fn make_functional_group_smarts_table(conn: &Connection) -> Result<(), rusqlite::Error> {

    let functional_smarts: HashMap<String, String> = import_tsv("database/smarts2.tsv");
    // Create the table
    conn.execute(
        "CREATE TABLE IF NOT EXISTS functional_group_smarts (
                  name   TEXT NOT NULL,
                  smarts TEXT NOT NULL
                  )",
        [],
    )?;

    
    // Insert data
    for (key, value) in functional_smarts.iter() {
        let a: [&String; 2] = [&key, &value];
        let mut stmt = conn.prepare("INSERT OR REPLACE INTO functional_group_smarts (name, smarts) VALUES (?1, ?2)").unwrap();
        stmt.execute(&a)?;
    }
    
    Ok(())
}


fn drop_columns_except(conn: &Connection, table_name: &str, columns_to_keep: Vec<&str>) -> rusqlite::Result<()> {
    // Generate comma-separated list of columns to keep
    let columns_str = columns_to_keep.join(", ");

    // Create a new temporary table
    let create_tmp_table_sql = format!("CREATE TABLE tmp AS SELECT {} FROM {}", columns_str, table_name);
    conn.execute(&create_tmp_table_sql, [])?;

    // Drop the original table
    let drop_original_table_sql = format!("DROP TABLE {}", table_name);
    conn.execute(&drop_original_table_sql, [])?;

    // Rename the temporary table to the original table name
    let rename_tmp_table_sql = format!("ALTER TABLE tmp RENAME TO {}", table_name);
    conn.execute(&rename_tmp_table_sql, [])?;

    Ok(())
}


fn main() -> () {
    
    let db_path = "db.db";
    //let db_path = ":memory:";
    
    //remove the database if it already exists
    std::fs::remove_file(db_path).expect("Database could not be removed");

    let conn: Connection = Connection::open(db_path).expect("Could not open db");
    make_metabolite_db(&conn).unwrap();
    make_functional_group_smarts_table(&conn).unwrap();
    make_lipid_db(&conn).unwrap();
    //conn.execute("DROP TABLE adducts", []).unwrap();
    make_adduct_db(&conn).unwrap();
    //conn.execute("DROP TABLE matrices", []).unwrap();
    make_matrices_db(&conn).unwrap();
    //let columns_to_keep = vec!["id", "phenols", "aldehydes", "carboxylicacids", "primaryamines"];
    //drop_columns_except(&conn, "functional_groups", columns_to_keep).unwrap();

    //conn.execute("DROP TABLE user_metabolites", []).unwrap();
    //conn.execute("DROP TABLE user_derivatized_by", []).unwrap();
    //conn.execute("DROP TABLE user_endogeneity", []).unwrap();
    //conn.execute("DROP TABLE user_in_tissue", []).unwrap();
    //conn.execute("DROP TABLE user_db_accessions", []).unwrap();
    //conn.execute("DROP TABLE user_functional_groups", []).unwrap();
    //conn.execute("DROP TABLE user_matrices", []).unwrap();
    //conn.execute("DROP TABLE user_adducts", []).unwrap();
    //conn.execute("DROP TABLE user_functional_group_smarts", []).unwrap();

}
