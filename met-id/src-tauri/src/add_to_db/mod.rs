pub mod db_accessions;
pub mod derivatized_by;
pub mod endogeneity;
pub mod functional_group;
pub mod functional_group_macos;
mod functional_group_server;
pub mod in_tissue;
pub mod lipids_functional_groups;
pub mod matrix;
pub mod metabolite;
pub mod msms;

use crate::add_to_db::functional_group::add_fg_to_db;
use crate::add_to_db::metabolite::add_metabolite_to_db;

use crate::database_setup::get_connection;
use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::Result;
use std::collections::HashMap;

fn get_hashmap_from_table(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
) -> HashMap<String, String> {
    let mut hashmap: HashMap<String, String> = HashMap::new();

    for table in &["functional_group_smarts", "user_functional_groups_smarts"] {
        match conn.prepare(&format!("SELECT name, smarts FROM {}", table)) {
            Ok(mut stmt) => {
                let mut rows = stmt.query([]).unwrap();

                while let Some(row) = rows.next().unwrap() {
                    let key: String = row.get(0).unwrap();
                    let value: String = row.get(1).unwrap();
                    hashmap.insert(key, value);
                }
            }
            Err(_) => (),
        }
    }
    hashmap
}

pub fn add_to_db_rust(
    name: String,
    smiles_smarts_mz: String,
    met_type: String,
    endo_exo_or_other: HashMap<String, bool>,
    in_tissue: HashMap<String, bool>,
    _adducts: Vec<String>,
    progress_sender: std::sync::mpsc::Sender<f32>,
) -> bool {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    match &met_type[..] {
        "metabolite" => {
            add_metabolite_to_db(conn, name, smiles_smarts_mz, in_tissue, endo_exo_or_other)
        }
        "fg" => add_fg_to_db(
            conn,
            name,
            smiles_smarts_mz,
            endo_exo_or_other,
            progress_sender,
        ),
        _ => (),
    }

    true
}

fn get_table_column_names(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
    table_name: &str,
) -> Result<Vec<String>> {
    // Prepare a query, doesn't matter if we're not going to actually execute it.
    let stmt = conn.prepare(&format!("SELECT * FROM {table_name} LIMIT 0"))?;

    // Get the column count.
    let col_count = stmt.column_count();

    // Get the column names.
    let mut col_names = Vec::new();
    for i in 0..col_count {
        let col_name = stmt.column_name(i)?;
        col_names.push(col_name.to_string());
    }

    Ok(col_names)
}
