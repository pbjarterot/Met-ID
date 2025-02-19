use crate::add_to_db::functional_group_server::functional_group_server;
use crate::database_setup::get_connection;
use crate::sql_mod::table::check_if_table_exists;

use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, Result, ToSql};
use std::collections::HashMap;
use std::iter::repeat;

use super::get_table_column_names;


//functional_groups adder for MacOS
#[cfg(target_os = "macos")]
fn functional_group_target(progress_sender: std::sync::mpsc::Sender<f32>, smarts: String, table: String, name: String) -> Result<()> {
    use crate::add_to_db::functional_group_macos::functional_group_macos;
    functional_group_macos(progress_sender, smarts, table, name.clone())
}

//functional_groups adder for windows, linux
#[cfg(any(target_os = "windows", target_os = "linux"))]
fn functional_group_target(progress_sender: std::sync::mpsc::Sender<f32>, smarts: String, table: String, name: String) -> Result<()>{
    functional_group_server(progress_sender, smarts, table, name.clone()).unwrap();
    Ok(())
}



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

    functional_group_target(progress_sender, smarts.clone(), "functional_groups".to_string(), name.clone()).unwrap();

    //update functional_groups & user_functional_groups
    // Add the new column (as before)
    
    conn.execute(&format!("INSERT INTO functional_group_smarts (name, smarts) VALUES ('{}', '{}')", name, smarts), []).unwrap();
    //if any matrix is pressed, update derivatized_by and user_matrices
    update_matrix_table_with_functional_group("matrices", &name, &matrices);
    update_matrix_table_with_functional_group("user_matrices", &name, &matrices);
}

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

pub fn get_smiles_from_db() -> Vec<Vec<String>> {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    const BATCH_SIZE: usize = 1_000;

    let mode: String = conn
        .query_row("PRAGMA journal_mode=WAL;", [], |row| row.get(0))
        .unwrap();

    assert_eq!(mode, "wal");

    let mut select_stmt = conn
        .prepare(&format!("SELECT smiles FROM metabolites", ))
        .unwrap();

    //let total_rows = select_stmt.query_map([], |_row| Ok(())).unwrap().count();
    //let mut processed_rows = 0;

    let mut smiles_vec = Vec::new();

    for row in select_stmt
        .query_map([], |row| {
            let smiles: String = row.get(0)?;
            Ok(smiles)
        })
        .unwrap()
        {
            let smiles = row.unwrap();
            smiles_vec.push(smiles);
        };

    smiles_vec.chunks(BATCH_SIZE).map(|chunk| chunk.to_vec()).collect()

}