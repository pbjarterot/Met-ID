use std::collections::HashMap;
use serde::{Serialize, Serializer, Deserialize};
use crate::mass_match::mass_matcher;
use maplit::hashmap;
use thiserror::Error;
use r2d2_sqlite::SqliteConnectionManager;
use std::mem;
use rusqlite::{params, Result};//, TransactionBehavior};//, Connection};
use byteorder::{ByteOrder, LittleEndian};
use std::cmp::Ordering;
//use crate::validation::mass_from_smiles;
use crate::sql_build_query::{build_query, build_condition_query};
use crate::add_to_db::add_to_db_functions::*;
use crate::{get_connection, get_msms_connection};

use std::sync::mpsc;



#[derive(Debug)]
struct MS1DbRow {
    mz: f64,
    name: String,
    mname: String,
    accession: String,
    smiles: String,
    formula: String
}

#[derive(Debug)]
struct MS2DbRow {
    identifier: String,
    adduct: String,
    cid: String,
    window: String,
    tof: String,
    mz: String,
    data: (Vec<f64>, Vec<i64>),
}

#[derive(Deserialize, Debug)]
pub struct Args {
    pub metabolome: String,
    pub matrix: String,
    pub met_type: Vec<String>,
    pub adducts: Vec<String>
}

#[derive(Debug, Error)]
#[error("Serde error: {0}")]
pub struct SerdeError(#[from] serde_json::Error);

impl Serialize for SerdeError {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

#[derive(Serialize)]
pub enum MyError {
    SerdeJsonError(SerdeError),
}

impl From<serde_json::Error> for MyError {
    fn from(err: serde_json::Error) -> MyError {
        MyError::SerdeJsonError(SerdeError::from(err))
    }
}


fn sql_query(query: &String) -> (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>) {
    //connect to db
    let conn = get_connection().unwrap();
    //query
    let mut stmt: rusqlite::Statement = conn.prepare(query).expect("Query cannot be run");
    //collect results from the database query
    
    let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
        Ok(MS1DbRow {
            mz: row.get(0).unwrap_or(0.0 as f64),
            name: row.get(1).unwrap_or("".to_string()),
            mname: row.get(2).unwrap_or("".to_string()),
            accession: row.get(3).unwrap_or("".to_string()),
            smiles: row.get(4).unwrap_or("".to_string()),
            formula: row.get(5).unwrap_or("".to_string())
        })
    }).unwrap();

    let mut mzs: Vec<f64> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut mnames: Vec<String> = Vec::new();
    let mut accessions: Vec<String> = Vec::new();
    let mut smiless: Vec<String> = Vec::new();
    let mut formulas: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {

        let row: MS1DbRow = item.unwrap();
        mzs.insert(index, row.mz);
        names.insert(index, row.name);
        mnames.insert(index, row.mname);
        accessions.insert(index, row.accession);
        smiless.insert(index, row.smiles);
        formulas.insert(index, row.formula);

    }
    (mzs, names, mnames, accessions, smiless, formulas)
}

#[tauri::command]
pub fn sql_handler(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>, _mass_error: String, masses: Vec<String>, mzwindow: String) -> Vec<Vec<HashMap<String, String>>> {
    let count: bool = false;

    let args: Args = Args {
        metabolome: met,
        matrix: mat,
        met_type: typ,
        adducts
    };
    // Use my_class instance as needed
    let min_peak: f64 = 0.0;
    let max_peak: f64 = 1200.0;

    let query_str: String = build_query(args, min_peak, max_peak, count);
    println!("{:?}", query_str);

    let db_input: (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>) = sql_query(&query_str);

    let mut input_masses: Vec<f64> = Vec::<f64>::new();
    let mut latent_variables: Vec<String> = Vec::<String>::new();
    for i in masses {
        match i.parse::<f64>() {
            Ok(num) => input_masses.push(num),
            Err(_e ) => {
                println!("Treating {:?} as a latent variable", i);
                latent_variables.push(i);
            },
        }
    }
    let db_input_array: Box<[f64]> = db_input.0.into_boxed_slice();
    let db_input_array_ref: &[f64] = &*db_input_array;
                                                            //                                     names       mnames      accessions  smiless      formulas
    let output: Vec<Vec<HashMap<String, String>>> = mass_matcher(input_masses, db_input_array_ref, db_input.1, db_input.2, db_input.3, db_input.4, db_input.5, _mass_error, mzwindow);    
    output
}



fn build_count_query(met: String, matrix: String, typ: Vec<String>, adducts: Vec<String>) -> String {
    let mut query: String = "SELECT COUNT(*) FROM ".to_string();

    let functional_groups_map: HashMap<&str, &str> = hashmap!{
        "Phenols" => "functional_groups.phenols",
        "Aldehydes" => "functional_groups.aldehydes",
        "Carboxylic" => "functional_groups.carboxylicacids",
        "Primary" => "functional_groups.primaryamines",
    };

    let met_type_map: HashMap<&str, &str> = hashmap!{
        "Endogenous" => "endogeneity.endogenous",
        "Exogenous" => "endogeneity.exogenous",
        "Unspecified" => "endogeneity.unspecified"
    };

    if met.starts_with("HMDB") {
        query += "metabolites INNER JOIN endogeneity ON metabolites.id = endogeneity.id INNER JOIN functional_groups ON metabolites.id = functional_groups.id";
    } else {
        query += "lipids";
        return query;
    };
    if matrix == "Positive mode".to_string() || matrix == "Negative mode".to_string() {
        query += &format!(" WHERE {met_type} ", 
            met_type= build_condition_query(&typ, &met_type_map, false).unwrap())[..];
        return query;
    } else {
        query += &format!(" WHERE {met_type} {functional_group}", 
            met_type = build_condition_query(&typ, &met_type_map, true).unwrap(),
            functional_group = build_condition_query(&adducts, &functional_groups_map, false).unwrap())[..];
        return query;
    }
}

#[tauri::command]
pub fn sql_counter(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>) -> i64 {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    if typ.len() == 0 || adducts.len() == 0 {
        return 0;
    }
    let query_str: String = build_count_query(met, mat, typ, adducts);
    println!("{:?}", query_str);

    let count: i64 = conn.query_row(
        &query_str[..], 
        rusqlite::params![],
        |row: &rusqlite::Row<'_>| row.get(0),
    ).unwrap();

    count
}

fn collect_msms_spectra(binary_data: BinaryData) -> (Vec<f64>, Vec<i64>) {
    // Ensure that the binary data is properly structured
    assert!(binary_data.binary.len() % (mem::size_of::<f64>() + mem::size_of::<i64>()) == 0);

    let mut f64_values = Vec::new();
    let mut i64_values = Vec::new();

    // Iterate through the binary data, taking slices of the appropriate sizes
    for i in (0..binary_data.binary.len()).step_by(mem::size_of::<f64>() + mem::size_of::<i64>()) {
        let f64_bytes = &binary_data.binary[i..i + mem::size_of::<f64>()];
        let i64_bytes = &binary_data.binary[i + mem::size_of::<f64>()..i + mem::size_of::<f64>() + mem::size_of::<i64>()];

        // Use the byteorder crate to convert the slices to numbers
        let f64_value = f64::from_le_bytes(f64_bytes.try_into().unwrap());
        let i64_value = i64::from_le_bytes(i64_bytes.try_into().unwrap());

        f64_values.push(f64_value);
        i64_values.push(i64_value);
    }

    let result = (f64_values, i64_values);
    result
}

#[derive(Debug)]
struct IDentifierMap {
    name: String,
    hmdb: String,
}

fn get_names_from_identifiers(input_data: &Vec<String>) -> rusqlite::Result<HashMap<String, String>> {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    let mut name_id_conversion: HashMap<String, String> = HashMap::new();

    let params: String = input_data.iter().map(|s| format!("'{}'", s)).collect::<Vec<_>>().join(", ");

    let sql: String = format!("SELECT metabolites.name, db_accessions.hmdb FROM metabolites JOIN db_accessions ON metabolites.id = db_accessions.id WHERE db_accessions.hmdb IN ({})", params);
    // Execute the query
    let mut stmt = conn.prepare(&sql)?;
    let rows = stmt.query_map([], |row| {
        Ok(IDentifierMap {
            name: row.get(0).unwrap_or("".to_string()), 
            hmdb: row.get(1).unwrap_or("".to_string()), 
        })
    }).unwrap();

    // Iterate through the rows and do something with the results
    for item in rows {
        let row: IDentifierMap = item.unwrap();
        name_id_conversion.insert(row.hmdb, row.name);
    }
    name_id_conversion.insert("HMDBM000001".to_string(), "FMP10".to_string());

    Ok(name_id_conversion)
}

#[tauri::command]
pub fn get_msms() -> Vec<Vec<String>>{
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_msms").unwrap();
    let mut stmt: rusqlite::Statement = conn.prepare("SELECT identifier, adduct, cid, window, tof, mz, spectra FROM MSMS").expect("Query cannot be run");
    let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
        let adduct_value = row.get(1).unwrap_or("".to_string());
        let adduct_final_value = if adduct_value.is_empty() {
            row.get(5).unwrap_or("".to_string())
        } else {
            adduct_value
        };
        Ok(MS2DbRow {
            identifier: row.get(0).unwrap_or("".to_string()),
            adduct: adduct_final_value,
            cid: row.get(2).unwrap_or("".to_string()),
            window: row.get(3).unwrap_or("".to_string()),
            tof: row.get(4).unwrap_or("".to_string()),
            mz: row.get(5).unwrap_or("".to_string()),
            data: (Vec::new(), Vec::new()),
            //data: row.get(6).unwrap_or(Vec::new()),
        })
    }).unwrap();

    let mut identifiers: Vec<String> = Vec::new();
    let mut adducts: Vec<String> = Vec::new();
    let mut cids: Vec<String> = Vec::new();
    let mut windows: Vec<String> = Vec::new();
    let mut tofs: Vec<String> = Vec::new();
    let mut mzs: Vec<String> = Vec::new();
    let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        let row: MS2DbRow = item.unwrap();
        identifiers.insert(index, row.identifier);
        adducts.insert(index, row.adduct);
        cids.insert(index, row.cid);
        windows.insert(index, row.window);
        tofs.insert(index, row.tof);
        mzs.insert(index, row.mz);
        datas.insert(index, row.data);
    }
    let mut identifiers2: Vec<String> = identifiers.clone();
    identifiers2.sort();
    identifiers2.dedup();
    let conversions: HashMap<String, String> = get_names_from_identifiers(&identifiers2).unwrap();

    let mut names: Vec<String> = Vec::new();
    for (index, id) in identifiers.iter().enumerate() {
        let name = conversions.get(id).unwrap_or(&"".to_string()).to_owned();
        names.insert(index, name);
    }

    let return_data: Vec<Vec<String>> = vec![names, identifiers, adducts, cids, windows, tofs, mzs];
    return_data
}

#[derive(Serialize, Deserialize)]
struct SpectrumPoint {
    x: f64,
    y: i64,
}

struct BinaryData {
    binary: Vec<u8>
}

struct MSMSBinary {
    binary: Vec<u8>,
    cid: String
}

fn get_single_msms_spectra(identifier: String, adduct: String) -> (Vec<(Vec<f64>, Vec<i64>)>, Vec<String>) {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_single_msms_spectra").unwrap();

    let sql: String = format!("SELECT spectra, cid FROM MSMS WHERE identifier = '{}' AND adduct = '{}'", identifier, adduct);

    let mut stmt = conn.prepare(&sql).unwrap();
    let rows = stmt.query_map([], |row| {
        Ok(MSMSBinary {
            binary: row.get(0).unwrap_or(Vec::new()), 
            cid: row.get(1).unwrap_or(String::new()), 
        })
    }).unwrap();
    

    let mut spectras: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();
    let mut cids: Vec<String> = Vec::new();

    for (index, item) in rows.enumerate() {
        let row: MSMSBinary = item.unwrap();
        let bin: BinaryData = BinaryData{ binary: row.binary};
        let data: (Vec<f64>, Vec<i64>) = collect_msms_spectra(bin);
        spectras.insert(index, data);
        cids.insert(index, row.cid);

    }


    (spectras, cids)
}

#[tauri::command]
pub fn get_msms_spectra(identifier: String, adduct: String) -> String {
    let spectra: (Vec<(Vec<f64>, Vec<i64>)>, Vec<String>) = get_single_msms_spectra(identifier, adduct);

    let mut data: HashMap<String, Vec<SpectrumPoint>> = HashMap::new();
    for i in 0..spectra.0.len() {
        let spectrum: &(Vec<f64>, Vec<i64>) = &spectra.0[i];
        let cid: &String = &spectra.1[i];
        let mut dat: Vec<SpectrumPoint> = Vec::new();
        for idx in 0..spectrum.0.len() {
            dat.push(SpectrumPoint {x: spectrum.0[idx], y: spectrum.1[idx]});
        }
        data.insert(cid.to_string(), dat);
    }
    serde_json::to_string(&data).expect("Failed to serialize data")
}

#[tauri::command]
pub fn get_name_from_identifier_msms(identifier: String) -> String {

    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    let sql: String = format!("SELECT metabolites.name, db_accessions.hmdb FROM metabolites JOIN db_accessions ON metabolites.id = db_accessions.id WHERE db_accessions.hmdb IN (?)");
    let name: Result<String, rusqlite::Error> = conn.query_row(&sql, &[&identifier], |row| {
        row.get(0)
    });

    match name {
        Ok(name) => return name,
        Err(e) => {
            // Handle error or return None if no row found
            if e == rusqlite::Error::QueryReturnedNoRows {
                return "".to_string();
            } else {
                panic!("Database error: {:?}", e);
            }
        }
    }

}


fn match_fragment(blob: &[u8], min_mz: f64, max_mz: f64) -> Option<(f64, i64)> {
    let float_bytes = &blob[0..8]; // First 8 bytes for f64
    let int_bytes = &blob[8..16]; // Next 8 bytes for i64

    let float_val = LittleEndian::read_f64(float_bytes);
    let int_val = LittleEndian::read_i64(int_bytes);

    if min_mz < float_val && float_val < max_mz {
        return Some((float_val, int_val));
    }
    
    None
}


#[tauri::command]
pub fn find_msms_fragments(mass: f64, window: f64) -> () {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("find_msms_fragments").unwrap();

    let mut stmt: rusqlite::Statement<'_> = conn.prepare("SELECT name, cid, adduct, spectra FROM MSMS").unwrap();
    let mut rows: rusqlite::Rows<'_> = stmt.query(params![]).unwrap();


    let min_mz: f64 = mass - (window/2000.0);
    let max_mz: f64 = mass + (window/2000.0);



    while let Some(row) = rows.next().unwrap() {
        let name: String = row.get(0).unwrap();
        let cid: i32 = row.get(1).unwrap();
        let adduct: String = row.get(2).unwrap();
        let spectra: Vec<u8> = row.get(3).unwrap();

        if let Some((float_val, int_val)) = match_fragment(&spectra, min_mz, max_mz) {
            println!("Column name: {}, cid: {}, adduct: {}, Float Value: {}, Int Value: {}", name, cid, adduct, float_val, int_val);
        }
    }
}



#[allow(unused_assignments)]
fn contains_between_a_and_b(vec: &[f64], min_mz: f64, max_mz: f64) -> bool {
    let mut lower_bound: usize = 0;
    let mut upper_bound: usize = 0;

    // Find the lower bound
    match vec.binary_search_by(|&x| x.partial_cmp(&min_mz).unwrap_or(Ordering::Greater)) {
        Ok(index) => lower_bound = index + 1,  // Skip exact matches
        Err(index) => lower_bound = index,      // Err gives us the point of insertion.
    }

    // Find the upper bound
    match vec.binary_search_by(|&x| x.partial_cmp(&max_mz).unwrap_or(Ordering::Less)) {
        Ok(index) => upper_bound = index,       // Skip exact matches
        Err(index) => upper_bound = index,      // Err gives us the point of insertion.
    }

    lower_bound < upper_bound
}


#[tauri::command]
pub fn ms2_search_spectra(name: String, fragment: String, ms1mass: String, fragmentslider: String, _ms1massslider: String) -> Vec<Vec<String>> {
    //let conn: r2d2::PooledConnection<SqliteConnectionManager> = MSMS_POOL.get().unwrap();
    let table: &str = "MSMS";
    let name_condition: &str = &name.as_str();
    let sort_field: &str = "name";
    
    let mut query = format!("SELECT identifier, adduct, cid, window, tof, mz, spectra FROM {}", table);
    
    if name_condition != "" {
        query.push_str(" WHERE name LIKE '%");
        query.push_str(name_condition);
        query.push_str("%' ")
    }
    

    if let Ok(ms1) = ms1mass.parse::<f64>() {
        if query.contains(" WHERE ") {
            query.push_str(" AND CAST(mz AS REAL) < ");
        } else {
            query.push_str(" WHERE CAST(mz AS REAL) < ");
        }
        query.push_str(&(ms1 + 0.5).to_string());
        query.push_str(" AND CAST(mz AS REAL) > ");
        query.push_str(&(ms1 - 0.5).to_string());
    }

    if !sort_field.is_empty() {
        query.push_str(" ORDER BY ");
        query.push_str(sort_field);
    }

    println!("{}", &query);

    if let Ok(mut stmt) = get_msms_connection("ms2_search_spectra").unwrap().prepare(&query) {
        let mut identifiers: Vec<String> = Vec::new();
        let mut adducts: Vec<String> = Vec::new();
        let mut cids: Vec<String> = Vec::new();
        let mut windows: Vec<String> = Vec::new();
        let mut tofs: Vec<String> = Vec::new();
        let mut mzs: Vec<String> = Vec::new();
        let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();



        let db_iter: Box<dyn Iterator<Item = rusqlite::Result<Option<MS2DbRow>>>> = if let Ok(fragment_mz) = fragment.parse::<f64>() {
            let min_mz = fragment_mz - (fragmentslider.parse::<f64>().unwrap_or(0.0) /2000.0);
            let max_mz = fragment_mz + (fragmentslider.parse::<f64>().unwrap_or(0.0) /2000.0);
            println!("{:?}, {:?}, {:?}, {:?}", fragment_mz, min_mz, max_mz, fragmentslider);
            Box::new(stmt.query_map([], move |row: &rusqlite::Row<'_>| {
                let spectrum = collect_msms_spectra(BinaryData {binary: row.get(6).unwrap_or(Vec::new())});
                if contains_between_a_and_b(&spectrum.0, min_mz, max_mz) {
                    let adduct_value = row.get(1).unwrap_or("".to_string());
                    let adduct_final_value = if adduct_value.is_empty() {
                        row.get(5).unwrap_or("".to_string())
                    } else {
                        adduct_value
                    };
                    Ok(Some(MS2DbRow {
                        identifier: row.get(0).unwrap_or("".to_string()),
                        adduct: adduct_final_value,
                        cid: row.get(2).unwrap_or("".to_string()),
                        window: row.get(3).unwrap_or("".to_string()),
                        tof: row.get(4).unwrap_or("".to_string()),
                        mz: row.get(5).unwrap_or("".to_string()),
                        data: spectrum,

                    }))
                } else {
                    Ok(None)
                }
            }).unwrap())
        } else {
            Box::new(stmt.query_map([], |row: &rusqlite::Row<'_>| {
                let adduct_value = row.get(1).unwrap_or("".to_string());
                let adduct_final_value = if adduct_value.is_empty() {
                    row.get(5).unwrap_or("".to_string())
                } else {
                    adduct_value
                };
                Ok(Some(MS2DbRow {
                    identifier: row.get(0).unwrap_or("".to_string()),
                    adduct: adduct_final_value,
                    cid: row.get(2).unwrap_or("".to_string()),
                    window: row.get(3).unwrap_or("".to_string()),
                    tof: row.get(4).unwrap_or("".to_string()),
                    mz: row.get(5).unwrap_or("".to_string()),
                    data: (Vec::new(), Vec::new()),
                }))
            }).unwrap())
        };

        /* 
        // Execute query and process results here
        let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
            let spectrum = collect_msms_spectra(BinaryData {binary: row.get(6).unwrap_or(Vec::new())});
            if contains_between_a_and_b(&spectrum.0, 0.0, 1.0) {
                Ok(Some(MS2DbRow {
                    identifier: row.get(0).unwrap_or("".to_string()),
                    adduct: row.get(1).unwrap_or("".to_string()),
                    cid: row.get(2).unwrap_or("".to_string()),
                    window: row.get(3).unwrap_or("".to_string()),
                    tof: row.get(4).unwrap_or("".to_string()),
                    mz: row.get(5).unwrap_or("".to_string()),
                    data: spectrum,
                }))
            } else {
                Ok(None)
            }

        }).unwrap();

        */
        let filtered_iter = db_iter.filter_map(|result| {
            match result {
                Ok(Some(db_row)) => Some(Ok(db_row)),
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            }
        });
        //let filtered_iter = db_iter.filter_map(Result::ok);
        //parse results for passing back to the parent function
        for (index, item) in filtered_iter.enumerate() {

            let row: MS2DbRow = item.unwrap();
            

            identifiers.insert(index, row.identifier);
            adducts.insert(index, row.adduct);
            cids.insert(index, row.cid);
            windows.insert(index, row.window);
            tofs.insert(index, row.tof);
            mzs.insert(index, row.mz);
            datas.insert(index, row.data);//collect_msms_spectra(BinaryData { binary: row.data }));
            
        }


        let mut identifiers2: Vec<String> = identifiers.clone();
        identifiers2.sort();
        identifiers2.dedup();
        let conversions: HashMap<String, String> = get_names_from_identifiers(&identifiers2).unwrap();

        let mut names: Vec<String> = Vec::new();
        for (index, id) in identifiers.iter().enumerate() {
            let name = conversions.get(id).unwrap_or(&"".to_string()).to_owned();
            names.insert(index, name);
        }
        println!("{:?}", names);
        let return_data: Vec<Vec<String>> = vec![names, identifiers, adducts, cids, windows, tofs, mzs];
        return return_data
    } else {
        return Vec::new()
    }
}







pub fn add_to_db_rust(name: String, smiles_smarts_mz: String, met_type: String, 
                    endo_exo_or_other: HashMap<String, bool>, in_tissue: HashMap<String, bool>, adducts: HashMap<String, String>, progress_sender: &mpsc::Sender<f32>) -> bool {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    match &met_type[..] {
        "metabolite" => add_metabolite_to_db(conn, name, smiles_smarts_mz, in_tissue, endo_exo_or_other),
        "matrix" => add_matrix_to_db(conn, name, smiles_smarts_mz, endo_exo_or_other, adducts),
        "fg" => add_fg_to_db(conn, name, smiles_smarts_mz, endo_exo_or_other, progress_sender),
        _ => ()
    }
    true
}






fn get_fg() -> Result<Vec<String>> {
    let conn = get_connection().unwrap();
    let mut stmt = conn.prepare("PRAGMA table_info(matrices);")?;
    let column_info = stmt.query_map(params![], |row| {
        let name: String = row.get(1)?; // 'name' is the second column, thus index 1
        Ok(name)
    })?;

    let mut column_names = Vec::new();
    for name in column_info {
        if let Ok(name) = name {
            column_names.push(name);
        }
    }
    let value_to_remove = "matrix".to_string();
    
    if let Some(pos) = column_names.iter().position(|x| *x == value_to_remove) {
        column_names.remove(pos);
    }
    Ok(column_names)
}

#[tauri::command]
pub fn get_functional_groups() -> Vec<String> {
    match get_fg() {  // added parentheses here
        Ok(result) => return result,
        Err(_) => return Vec::new(),
    }
}


#[tauri::command]
pub fn get_tissues() -> Vec<String> {
    let conn = get_connection().unwrap();
    let mut stmt = conn.prepare("PRAGMA table_info(in_tissue);").unwrap();
    let column_info = stmt.query_map(params![], |row| {
        let name: String = row.get(1)?; // 'name' is the second column, thus index 1
        Ok(name)
    }).unwrap();

    let mut column_names = Vec::new();
    for name in column_info {
        if let Ok(name) = name {
            column_names.push(name);
        }
    }
    let value_to_remove = "id".to_string();
    
    if let Some(pos) = column_names.iter().position(|x| *x == value_to_remove) {
        column_names.remove(pos);
    }

    column_names

}






fn get_mtx() -> Result<Vec<String>> {
    let conn = get_connection().unwrap();
    let mut stmt = conn.prepare("SELECT matrix FROM matrices WHERE matrix NOT IN ('POSITIVE_MODE', 'NEGATIVE_MODE')").unwrap();
    let name_iter = stmt.query_map(params![], |row| {
        Ok(row.get::<_, String>(0).unwrap())
    })?;
    
    // Collect names into a vector, handling errors as they occur
    let names = name_iter.collect();
    names
}

#[tauri::command]
pub fn get_matrices() -> Vec<String> {
    match get_mtx() {  // added parentheses here
        Ok(result) => return result,
        Err(_) => return Vec::new(),
    }
}


pub fn check_if_table_exists(source_table: &str, new_table: &str) -> Result<()>{
    extern crate rusqlite;
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    let mut stmt = conn.prepare(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?1",
    )?;
    let mut rows = stmt.query(params![new_table])?;
    if rows.next()?.is_some() {
        return Ok(());
    }

    let mut stmt = conn.prepare(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name=?1",
    )?;
    let mut rows = stmt.query(params![source_table])?;
    if let Some(row) = rows.next()? {
        let schema: String = row.get(0)?;

        let new_schema = schema.replace(source_table, new_table);

        conn.execute(&new_schema, [])?;

        println!("Table {} has been created.", new_table);
    } else {
        println!("Source table {} does not exist.", source_table);
    }

    Ok(())
}

