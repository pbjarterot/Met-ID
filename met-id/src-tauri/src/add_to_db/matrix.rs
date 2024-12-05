use crate::database_setup::get_connection;
use crate::sql_mod::table::check_if_table_exists;
use rusqlite::{params, ToSql};
use std::collections::HashMap;
use std::iter::repeat;

#[derive(Debug)]
struct MatrixAdduct {
    name: String,
    mz: String,
    derivs: String,
}

#[tauri::command]
pub fn add_matrix_to_db_rust(
    name: String,
    checkboxes: HashMap<String, bool>,
    adducts: Vec<String>,
) -> () {
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();
    check_if_table_exists("matrices", "user_matrices").unwrap();
    check_if_table_exists("adducts", "user_adducts").unwrap();

    let array_length: usize = adducts.len() / 3;

    let adducts: Vec<MatrixAdduct> = adducts
        .chunks(3)
        .map(|chunk| {
            let cnk = chunk.to_vec();
            MatrixAdduct {
                name: cnk[0].clone(),
                mz: cnk[1].clone(),
                derivs: cnk[2].clone(),
            }
        })
        .collect();

    // Iterate over the vectors and insert the data into the table
    for i in 0..array_length {
        conn.execute(
            "INSERT INTO user_adducts (adduct, mname, numfunctionalgroups, formula, deltamass, maxcoverage) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![adducts[i].name, name, adducts[i].derivs, "".to_string(), adducts[i].mz, array_length],
        ).unwrap();
    }

    // filling user_matrices
    let mut column_names: Vec<String> = checkboxes
        .iter()
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

    let mut values: Vec<&dyn ToSql> = checkboxes
        .iter()
        .map(|(_, value)| value as &dyn ToSql)
        .collect();
    values.push(&name as &dyn ToSql);

    conn.execute(&sql, values.as_slice()).unwrap();
}
