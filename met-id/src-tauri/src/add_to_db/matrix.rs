use regex::Regex;
use crate::sidecar::sidecar_function;
use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, ToSql};
use std::iter::repeat;
use std::collections::HashMap;
use crate::sql_mod::table::check_if_table_exists;



fn matrix_for_db_sidecar(smiles: String) -> String {
    let mut smiles_vec: Vec<_> = Vec::new();
    smiles_vec.insert(0, smiles.clone());
    let sidecar_output: std::prelude::v1::Result<String, crate::sidecar::CommandError> = sidecar_function("matrix_for_db".to_string(), smiles_vec);
    
    let re: Regex = Regex::new(r"[\r\n]").unwrap();
    let sc_o: String = re.replace_all(&sidecar_output.unwrap_or("".to_string()), "").to_string();

    sc_o
}
fn split_key(key: &str) -> Option<(String, i32)> {
    // Reverse the string and find the first non-numeric character
    let reversed: String = key.chars().rev().collect();
    let idx = reversed.chars().position(|ch| !ch.is_numeric())?;

    // Get the prefix and number using the identified index
    let number: i32 = reversed[0..idx].chars().rev().collect::<String>().parse().ok()?;
    let prefix: String = reversed[idx..].chars().rev().collect();

    Some((prefix, number))
}


pub fn add_matrix_to_db(conn: r2d2::PooledConnection<SqliteConnectionManager>, name: String, _smiles_mz: String, checkboxes: HashMap<String, bool>, adducts: HashMap<String, String>) -> () {
    check_if_table_exists("matrices", "user_matrices").unwrap();
    check_if_table_exists("adducts", "user_adducts").unwrap();    

    let array_length: usize = adducts.len()/3;
    
    //filling user_adducts
    let mut adduct_arr: [&str; 20] = [""; 20];
    let mut mx_formula_arr: [&str; 20] = [""; 20];
    let mut fg_arr: [&str; 20] = [""; 20];
    let mut deltamass_vec: Vec<String> = Vec::new();

    for (key, value) in adducts.iter() {
        let (prefix, number) = if let Some(result) = split_key(key) {
            result
        } else {
            continue; // Skip if key doesn't match expected format
        };

        match prefix.as_str() {
            "add-adduct-name-" => adduct_arr[number as usize] = value,
            "add-matrix-formula-" => {
                deltamass_vec.insert(number as usize, matrix_for_db_sidecar(value.clone()));
                mx_formula_arr[number as usize] = &value[..];
            },
            "add-fg-" => {
                fg_arr[number as usize] = if value.is_empty() { "1" } else { value };
            },
            _ => {}, // Other cases
        }
    }
    // Iterate over the vectors and insert the data into the table
    for i in 0..array_length {
        conn.execute(
            "INSERT INTO user_adducts (adduct, mname, numfunctionalgroups, formula, deltamass, maxcoverage) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![adduct_arr[i], name, fg_arr[i], mx_formula_arr[i], deltamass_vec[i], array_length],
        ).unwrap();
    }

    // filling user_matrices
    let mut column_names: Vec<String> = checkboxes.iter()
        .map(|(col, _)| format!("'{}'", col))
        .collect();
    let mut placeholders: Vec<&str> = repeat("?").take(checkboxes.len()).collect();

    column_names.push("matrix".to_string());
    placeholders.push("?");


    let sql = format!(
        "INSERT INTO user_matrices ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    let mut values: Vec<&dyn ToSql> = checkboxes.iter().map(|(_, value)| value as &dyn ToSql).collect();
    values.push(&name as &dyn ToSql);

    conn.execute(&sql, values.as_slice()).unwrap();

}