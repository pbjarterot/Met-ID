use std::collections::HashMap;
use rusqlite::Connection;
use serde::{Serialize, Serializer, Deserialize};
use crate::mass_match::mass_matcher;
use maplit::hashmap;
use thiserror::Error;

#[derive(Debug)]
struct MS1DbRow {
    mz: f64,
    name: String,
    mname: String,
    accession: String,
    smiles: String
}

#[derive(Deserialize, Debug)]
pub struct Args {
    metabolome: String,
    matrix: String,
    met_type: Vec<String>,
    adducts: Vec<String>
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

fn build_q_string(types: &Vec<String>, hash_map: &HashMap<&str, &str>) -> Option<String> {
    let mut string: String = Default::default();
    match types.len() {
        0 => return None,
        x if x > 0 => {
            for (index, type_) in types.iter().enumerate() {
                match hash_map.get(&type_[..]) {
                    Some(a) => {
                        string += "CAST(a.";
                        string += a;
                        string += " AS REAL)=1 ";
                    },
                    None => ()
                }
                //string += hash_map.get(&type_[..]).unwrap();
                if index < types.len() - 1 {
                     string += " OR ";
                }
            }
        },
        _ => (),
    }
    Some(string)
}

fn build_q_string2(types: &Vec<String>, hash_map: &HashMap<&str, &str>) -> Option<String> {
    let mut string: String = Default::default();
    match types.len() {
        0 => return None,
        x if x > 0 => {
            for (index, type_) in types.iter().enumerate() {
                match hash_map.get(&type_[..]) {
                    Some(a) => {
                        string += &format!("'{b}'", b=a)[..];
                        //string += &"\"".to_string()[..];
                        //string += a;
                        //string += &"\"".to_string()[..];
                    },
                    None => ()
                }
                //string += hash_map.get(&type_[..]).unwrap();
                if index < types.len() - 1 {
                     string += ", ";
                }
            }
        },
        _ => (),
    }
    Some(string)
}

fn build_query(args: Args, min_mz: f64, max_mz: f64) -> String {
    //println!("adducts: {:?}", &args.adducts);
    let met_type_map: HashMap<&str, &str> = hashmap!{
        "Endogenous" => "endogenous",
        "Exogenous" => "exogenous",
        "Unspecified" => "unspecified"
    };

    let _matrix_map: HashMap<&str, &str> = hashmap!{
        "FMP-10" => "fmp",
        "DPP" => "dpp",
        "TAHS" => "tahs",
        "CA" => "ca",
        "Girard P" => "girardp",
        "Girard T" => "girardt",
        "AMPP" => "ampp",
        "Boronic Acid" => "boronicacid"
    };

    let functional_groups_map: HashMap<&str, &str> = hashmap!{
        "Phenols" => "phenols",
        "Catechols" => "catechols",
        "Carbonyls" => "carbonyls",
        "Aldehydes" => "aldehydes",
        "Carboxylic" => "carboxylicacids",
        "Primary" => "primaryamines",
        "Triols" => "triols",
        "Diols" => "diols",
        "Hydroxyls" => "hydroxyls"
    };

    let adduct_map: HashMap<&str, &str> = hashmap!{
        "[M+K]+" => "M+K",
        "[M+H]+" => "M+H",
        "[M+Na]+" => "M+Na",
        "[M+K-H2]-" => "'M+K-H2'",
        "[M-H]-" => "'M-H'",
        "[M+Na-H2]-" => "'M+Na-H2'",
        "[M+Cl]-" => "'M+Cl'",
    };

    let tissue_map: HashMap<&str, &str> = hashmap!{
        "HMDB (CSF)" => "csf",
        "HMDB (Urine)" => "urine",
        "HMDB (Serum)" =>"serum",
    };

    let mut query: String = String::new();
    query += "SELECT CAST(a.mz AS REAL) + CAST(b.deltamass AS REAL), name, b.mname, accession, smiles FROM";

    let mut db_string: String;

    match args.metabolome.starts_with("HMDB") {
        true => {
            query += " metabolites ";
            match tissue_map.get(&args.metabolome[..]){
                Some(x) => {
                    db_string = "AND (".to_string();
                    db_string += x;
                    db_string += &String::from(" >0 )");
                },
                None => db_string = "".to_string(),
            }

        },
        false => {
            query += " lipids ";
            db_string = String::from("");
        },
    }
    query += "a INNER JOIN matrices b ";

    //println!("{:?}", args.matrix);
    //println!("{:?}", args.adducts);
    let if_fmp10: bool = args.matrix == "FMP-10".to_string();
    let pos_neg: bool = args.matrix == "Positive mode".to_string() || args.matrix == "Negative mode".to_string();
    let met_type_string: String = build_q_string(&args.met_type, &met_type_map).unwrap();
    let chosen_deriv_string: &str = &args.matrix;

    let a: String = if if_fmp10 {
        format!(r#"ON CASE WHEN b.mname='M+Anth' AND fmp = 1 AND a.mz + b.deltamass < {max} AND a.mz + b.deltamass > {min} {db} AND ({met_type}) AND ({func_group}) THEN 1 WHEN b.mname IN ('M+Anth', 'M+2Anth', 'M+2Anth_2') AND fmp = 2 AND a.mz + b.deltamass < {max} AND a.mz + b.deltamass > {min} {db} AND ({met_type}) AND ({func_group}) THEN 1 WHEN b.mname IN ('M+Anth', 'M+2Anth', 'M+2Anth_2', 'M+3Anth1', 'M+3Anth2', 'M+3Anth3') AND fmp>2 AND a.mz + b.deltamass < {max} AND a.mz + b.deltamass > {min} {db} AND ({met_type}) AND ({func_group}) THEN 1 ELSE 0 END=1"#, 
        max=&max_mz.to_string(), min=&min_mz.to_string(), db=db_string, met_type = met_type_string, func_group= build_q_string(&args.adducts, &functional_groups_map).unwrap())
    } else if pos_neg {
        format!(r#"ON (CAST(a.mz AS REAL) + CAST(b.deltamass AS REAL) < {max}) AND (CAST(a.mz AS REAL) + CAST(b.deltamass AS REAL) > {min}) {db} AND ({met_type}) AND (b.mname IN ({adduct}))"#, 
        max=&max_mz.to_string(), min=&min_mz.to_string(), db=db_string, met_type=met_type_string, adduct=build_q_string2(&args.adducts, &adduct_map).unwrap().replace("\\", ""))

    } else {
        format!(r#"ON b.mname ='M+{deriv}' AND a.mz + b.deltamass < {max} AND a.mz + b.deltamass > {min} {db} AND ({met_type}) AND ({func_group})"#, 
        max=&max_mz.to_string(), 
        min=&min_mz.to_string(),
        db=db_string, 
        met_type = met_type_string, 
        func_group=build_q_string(&args.adducts, &functional_groups_map).unwrap(),
        deriv = chosen_deriv_string
        )
    };
    query += &a;
    query += " ORDER BY CAST(a.mz AS REAL) + CAST(b.deltamass AS REAL)";
    //println!("{:?}", query);
    query
}



fn sql_query(query: &String) -> (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>) {
    //connect to db
    let conn: Connection = Connection::open("../databases/db.db").expect("Database file not found");

    //query
    let mut stmt: rusqlite::Statement = conn.prepare(query).expect("Query cannot be run");
    //collect results from the database query
    let db_iter = stmt.query_map([], |row| {
        Ok(MS1DbRow {
            mz: row.get(0)?,
            name: row.get(1)?,
            mname: row.get(2)?,
            accession: row.get(3)?,
            smiles: row.get(4)?
        })
    }).unwrap();

    let mut mzs: Vec<f64> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut mnames: Vec<String> = Vec::new();
    let mut accessions: Vec<String> = Vec::new();
    let mut smiless: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        //println!("{:?}", item);
        let row: MS1DbRow = item.unwrap();
        mzs.insert(index, row.mz);
        names.insert(index, row.name);
        mnames.insert(index, row.mname);
        accessions.insert(index, row.accession);
        smiless.insert(index, row.smiles)

    }
    (mzs, names, mnames, accessions, smiless)
}




#[tauri::command]
pub fn sql_handler(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>, _mass_error: String, masses: Vec<String>) -> Vec<Vec<HashMap<String, String>>> {
    println!("masses: {:?}", masses);
    let args: Args = Args {
        metabolome: met,
        matrix: mat,
        met_type: typ,
        adducts
    };
    // Use my_class instance as needed
    let min_peak: f64 = 110.0;
    let max_peak: f64 = 420.0;
    //println!("{:?}", args);
    let query_str: String = build_query(args, min_peak, max_peak);
    let db_input: (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>) = sql_query(&query_str);
    //println!("{:?}", db_input);


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

    let output: Vec<Vec<HashMap<String, String>>> = mass_matcher(input_masses, db_input_array_ref, db_input.1, db_input.2, db_input.3, db_input.4);
    println!("{:?}", output);


    
    output
}
