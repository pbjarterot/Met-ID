use crate::add_to_db::{
    db_accessions::fill_user_db_accessions, derivatized_by::fill_user_derivatized_by,
};
use crate::sidecar::sidecar_function;
use crate::sql_mod::table::check_if_table_exists;
use r2d2_sqlite::SqliteConnectionManager;
use regex::Regex;
use rusqlite::{Result, ToSql};
use std::collections::HashMap;
use std::iter::repeat;

use super::{
    endogeneity::fill_user_endogeneity,
    functional_group::{fill_user_functional_groups, get_matrices_fgs},
    get_hashmap_from_table,
    in_tissue::fill_user_in_tissue,
};

fn metabolite_for_db_sidecar(
    smiles: &String,
    smarts_map: &mut HashMap<String, String>,
) -> (String, String, HashMap<String, String>) {
    let mut smarts_vec: Vec<_> = smarts_map.values().cloned().collect();
    smarts_vec.insert(0, smiles.to_owned());

    let sidecar_output = sidecar_function("metabolite_for_db".to_string(), smarts_vec.to_owned());

    parse_input(&sidecar_output.unwrap()[..], smarts_map).unwrap()
}

fn parse_input(
    input: &str,
    smarts_map: &mut HashMap<String, String>,
) -> Result<(String, String, HashMap<String, String>), &'static str> {
    //let input = "C3H8 44.062600255999996 [0, 0, 0, 0]\r";

    // Extract molecule formula
    let molecule_re: Regex = Regex::new(r"([A-Z][a-z]*\d*)+").unwrap();
    let molecule: &str = molecule_re.find(input).unwrap().as_str();

    // Extract float number
    let float_re: Regex = Regex::new(r"\d+\.\d+").unwrap();
    let float_number: f64 = float_re.find(input).unwrap().as_str().parse().unwrap();

    // Extract numbers from the array
    let array_start: usize = input.find('[').unwrap();
    let array_end: usize = input.find(']').unwrap();
    let array_str: &str = &input[array_start + 1..array_end];
    let array: Vec<i32> = array_str
        .split(", ")
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    // Assuming keys and values have the same length
    let keys: Vec<String> = smarts_map.keys().cloned().collect();
    for (key, value) in keys.iter().zip(array.iter()) {
        smarts_map.insert(key.clone(), value.clone().to_string());
    }

    Ok((
        molecule.to_string(),
        float_number.to_string(),
        smarts_map.to_owned(),
    ))
}

pub fn fill_user_metabolites(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
    name: &String,
    smiles: &String,
    formula: &String,
    mz: &String,
) {
    let data = vec![
        ("name".to_string(), name),
        ("smiles".to_string(), smiles),
        ("chemicalformula".to_string(), formula),
        ("mz".to_string(), mz),
    ];

    let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
    let placeholders: Vec<&str> = repeat("?").take(data.len()).collect();

    let sql = format!(
        "INSERT INTO user_metabolites ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();
    conn.execute(&sql, values.as_slice()).unwrap();
}

pub fn add_metabolite_to_db(
    conn: r2d2::PooledConnection<SqliteConnectionManager>,
    name: String,
    smiles: String,
    in_tissue: HashMap<String, bool>,
    endo_exo: HashMap<String, bool>,
) -> () {
    check_if_table_exists("metabolites", "user_metabolites").unwrap();
    check_if_table_exists("derivatized_by", "user_derivatized_by").unwrap();
    check_if_table_exists("endogeneity", "user_endogeneity").unwrap();
    check_if_table_exists("in_tissue", "user_in_tissue").unwrap();
    check_if_table_exists("db_accessions", "user_db_accessions").unwrap();
    check_if_table_exists("functional_groups", "user_functional_groups").unwrap();

    let mut functional_smarts: HashMap<String, String> = get_hashmap_from_table(&conn);

    let (formula, mz, fgh) = metabolite_for_db_sidecar(&smiles, &mut functional_smarts);
    //let metabolite: Metabolite = Metabolite{ name: name.clone(), smiles: smiles.clone() };
    //let fgh: HashMap<String, String> = metabolite.functional_group(&functional_smarts).unwrap();
    let matrices_fgs: HashMap<String, Vec<String>> = get_matrices_fgs(&conn).unwrap();

    fill_user_derivatized_by(&conn, &fgh, &matrices_fgs);
    fill_user_metabolites(&conn, &name, &smiles, &formula, &mz);
    fill_user_endogeneity(&conn, &endo_exo);
    fill_user_in_tissue(&conn, &in_tissue);
    fill_user_functional_groups(&conn, &fgh);
    fill_user_db_accessions(&conn);
}
