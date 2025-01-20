use crate::add_to_db::functional_group_server::functional_group_server;
use crate::database_setup::get_connection;
use crate::sidecar::sidecar_function3;
use crate::sql_mod::table::check_if_table_exists;
use crate::{get_app_handle, MSDBPATH};

use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, Result, ToSql};
use std::collections::HashMap;
use std::iter::repeat;
use std::sync::mpsc;


use super::{functional_group_server, get_table_column_names};
/* 
fn single_functional_group(smiles: &mut Vec<String>, smarts: &String) -> Vec<usize> {
    smiles.insert(0, smarts.to_owned());
    let sidecar_output = sidecar_function3(get_app_handle().unwrap(), "metabolite".to_string(), smiles.to_owned());
    let result: Vec<usize> = match sidecar_output {
        Ok(s) => s.iter().map(|n| *n as usize).collect(),
        Err(_) => Vec::new(),
    };

    result
}
*/
/* 
pub fn update_functional_groups(
    table_name: &str,
    table2_name: &str,
    name: &String,
    smarts: &String,
    progress_sender: &mpsc::Sender<f32>,
) {
    check_if_table_exists("metabolites", "user_metabolites").unwrap();
    const BATCH_SIZE: usize = 100;

    let mut conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    let conn2: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    let mode: String = conn
        .query_row("PRAGMA journal_mode=WAL;", [], |row| row.get(0))
        .unwrap();

    assert_eq!(mode, "wal");

    let mut select_stmt = conn2
        .prepare(&format!("SELECT id, smiles FROM {}", table2_name))
        .unwrap();
    let total_rows = select_stmt.query_map([], |_row| Ok(())).unwrap().count();
    let mut processed_rows = 0;

    let mut batch_smiles = Vec::with_capacity(BATCH_SIZE);
    let mut batch_rowids = Vec::with_capacity(BATCH_SIZE);

    for row in select_stmt
        .query_map([], |row| {
            let rowid: usize = row.get(0)?;
            let smiles: String = row.get(1)?;
            Ok((rowid, smiles))
        })
        .unwrap()
    {
        let (rowid, smiles) = row.unwrap();
        batch_smiles.push(smiles);
        batch_rowids.push(rowid);

        if batch_smiles.len() == BATCH_SIZE {
            let match_counts = single_functional_group(&mut batch_smiles, &smarts);

            let transaction = conn.transaction().unwrap();

            for (i, &count) in match_counts.iter().enumerate() {
                transaction
                    .execute(
                        &format!("UPDATE {} SET '{}' = ?1 WHERE rowid = ?2", table_name, name),
                        params![count, batch_rowids[i]],
                    )
                    .unwrap();
            }

            transaction.commit().unwrap();

            processed_rows += BATCH_SIZE;
            let progress = (processed_rows as f32 / total_rows as f32) * 100.0;
            progress_sender.send(progress).unwrap();

            batch_smiles.clear();
            batch_rowids.clear();
        }
    }

    // Process any remaining items in the batch
    if !batch_smiles.is_empty() {
        let match_counts = single_functional_group(&mut batch_smiles, &smarts);

        let transaction = conn.transaction().unwrap();

        for (i, &count) in match_counts.iter().enumerate() {
            transaction
                .execute(
                    &format!("UPDATE {} SET '{}' = ?1 WHERE rowid = ?2", table_name, name),
                    params![count, batch_rowids[i]],
                )
                .unwrap();
        }

        transaction.commit().unwrap();

        processed_rows += batch_smiles.len();
        let progress = (processed_rows as f32 / total_rows as f32) * 100.0;
        progress_sender.send(progress).unwrap();
    }
}
*/

pub fn add_fg_to_db(
    conn: r2d2::PooledConnection<SqliteConnectionManager>,
    name: String,
    smarts: String,
    matrices: HashMap<String, bool>,
    progress_sender: std::sync::mpsc::Sender<f32>,
) -> () {
    check_if_table_exists("functional_group_smarts", "user_functional_group_smarts").unwrap();
    check_if_table_exists("functional_groups", "user_functional_groups").unwrap();
    check_if_table_exists("matrices", "user_matrices").unwrap();
    //let conn2: r2d2::PooledConnection<SqliteConnectionManager> = POOL.get().unwrap();

    conn.execute(
        "INSERT INTO user_functional_group_smarts (name, smarts) VALUES (?1, ?2)",
        params![name, smarts],
    )
    .unwrap();

    conn.execute(
        &format!("ALTER TABLE functional_groups ADD COLUMN '{}' INTEGER", name),
        [],
    )
    .unwrap();

    conn.execute(
        &format!(
            "ALTER TABLE user_functional_groups ADD COLUMN '{}' INTEGER",
            name
        ),
        [],
    )
    .unwrap();

    functional_group_server(progress_sender, smarts.clone(), "functional_groups".to_string(), name.clone()).unwrap();

    //update functional_groups & user_functional_groups
    // Add the new column (as before)
    
    conn.execute(&format!("INSERT INTO functional_group_smarts (name, smarts) VALUES ('{}', '{}')", name, smarts), []).unwrap();
    //if any matrix is pressed, update derivatized_by and user_matrices
    update_matrix_table_with_functional_group("matrices", &name, &matrices);
    update_matrix_table_with_functional_group("user_matrices", &name, &matrices);
}
/* 
fn update_functional_groups2(
    table_name: &str,
    table2_name: &str,
    name: &String,
    smarts: &String,
    progress_sender: &mpsc::Sender<f32>
) -> () {
    sidecar_function3(get_app_handle().unwrap(), progress_sender);
}
*/


/*
fn get_hashmap_from_table(conn: &r2d2::PooledConnection<SqliteConnectionManager>) -> HashMap<String, String> {
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
            Err(_) => {
                ()
            }
        }
    }
    hashmap
}
*/

fn update_matrix_table_with_functional_group(
    table_name: &str,
    name: &str,
    matrices: &HashMap<String, bool>,
) {
    let mut conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    conn.execute(
        &format!("ALTER TABLE {} ADD COLUMN '{}' TEXT", table_name, name),
        [],
    )
    .unwrap();

    let tx = conn.transaction().unwrap();
    for (key, value) in matrices {
        tx.execute(
            &format!("UPDATE {table_name} SET '{name}' = ?1 WHERE matrix = ?2"),
            params![value, key],
        )
        .unwrap();
    }
    tx.commit().unwrap();
}

pub fn get_matrices_fgs(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
) -> Result<HashMap<String, Vec<String>>, rusqlite::Error> {
    let mut stmt = conn.prepare("SELECT * FROM matrices").unwrap();
    let col_count = stmt.column_count();

    // Fetch column names beforehand
    let mut col_names = Vec::new();
    for i in 0..col_count {
        let name = stmt.column_name(i)?.to_string();
        col_names.push(name);
    }

    let mut matrices_fgs: HashMap<String, Vec<String>> = HashMap::new();

    // Execute the query and map rows
    let mut rows = stmt.query([])?;

    while let Some(row) = rows.next()? {
        let mut fgs: Vec<String> = Vec::new();
        for i in 0..col_count {
            let value: String = row.get(i)?;
            if value == "1" {
                fgs.push(col_names[i].clone())
            }
        }

        let name_: String = row.get(0)?;
        matrices_fgs.insert(name_, fgs);
    }
    matrices_fgs.remove("Negative Mode");
    matrices_fgs.remove("Positive Mode");

    Ok(matrices_fgs)
}

pub fn fill_user_functional_groups(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
    fgh: &HashMap<String, String>,
) {
    let mut functional_groups: Vec<String> =
        get_table_column_names(conn, "functional_groups").unwrap();
    functional_groups.remove(0);

    let data: Vec<_> = fgh.iter().map(|(k, v)| (format!("'{}'", k), v)).collect();

    let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
    let placeholders: Vec<&str> = repeat("?").take(data.len()).collect();

    let sql = format!(
        "INSERT INTO user_functional_groups ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();

    conn.execute(&sql, values.as_slice()).unwrap();
}
