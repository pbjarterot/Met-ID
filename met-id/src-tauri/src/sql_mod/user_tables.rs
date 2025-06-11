use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, Connection, Result};

use crate::database_setup::get_connection;

use super::table::check_if_table_exists;

struct MS1UserMetabolitesRow {
    id: usize,
    name: String,
    smiles: String,
    formula: String,
    mz: String,
}

pub fn update_user_metabolites() -> Vec<Vec<String>> {
    check_if_table_exists("metabolites", "user_metabolites").unwrap();
    check_if_table_exists("derivatized_by", "user_derivatized_by").unwrap();
    check_if_table_exists("endogeneity", "user_endogeneity").unwrap();
    check_if_table_exists("in_tissue", "user_in_tissue").unwrap();
    check_if_table_exists("db_accessions", "user_db_accessions").unwrap();
    check_if_table_exists("functional_groups", "user_functional_groups").unwrap();
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    let mut stmt: rusqlite::Statement = conn
        .prepare("SELECT id, name, smiles, chemicalformula, mz FROM user_metabolites")
        .expect("Query cannot be run");
    let db_iter = stmt
        .query_map([], |row: &rusqlite::Row<'_>| {
            Ok(MS1UserMetabolitesRow {
                id: row.get(0).unwrap(),
                name: row.get(1).unwrap_or("".to_string()),
                smiles: row.get(2).unwrap_or("".to_string()),
                formula: row.get(3).unwrap_or("".to_string()),
                mz: row.get(4).unwrap_or("".to_string()),
            })
        })
        .unwrap();

    let mut ids: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut smiles: Vec<String> = Vec::new();
    let mut formulas: Vec<String> = Vec::new();
    let mut mzs: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        let row: MS1UserMetabolitesRow = item.unwrap();
        ids.insert(index, row.id.to_string());
        names.insert(index, row.name);
        smiles.insert(index, row.smiles);
        formulas.insert(index, row.formula);
        mzs.insert(index, row.mz)
    }

    let return_data: Vec<Vec<String>> = vec![ids, names, smiles, formulas, mzs];
    return_data
}

pub fn remove_row_from_user_metabolites(rowid: usize) -> usize {
    remove_row_from_table(rowid, "user_metabolites");
    remove_row_from_table(rowid, "user_derivatized_by");
    remove_row_from_table(rowid, "user_endogeneity");
    remove_row_from_table(rowid, "user_in_tissue");
    remove_row_from_table(rowid, "user_functional_groups");
    remove_row_from_table(rowid, "user_db_accessions");

    return 1;
}

struct MS1UserMatricesRow {
    id: usize,
    name: String,
}

pub fn update_user_matrices() -> Vec<Vec<String>> {
    check_if_table_exists("matrices", "user_matrices").unwrap();
    check_if_table_exists("adducts", "user_adducts").unwrap();

    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    ensure_id_column(&conn, "user_matrices").unwrap();

    let mut stmt: rusqlite::Statement = conn
        .prepare("SELECT id, matrix FROM user_matrices")
        .expect("Query cannot be run");
    let db_iter = stmt
        .query_map([], |row: &rusqlite::Row<'_>| {
            let id_value = row.get(0).unwrap();
            Ok(MS1UserMatricesRow {
                id: id_value,
                name: row.get(1).unwrap_or("".to_string()),
            })
        })
        .unwrap();

    let mut ids: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut derivatizes: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        let row: MS1UserMatricesRow = item.unwrap();
        ids.insert(index, row.id.to_string());
        names.insert(index, row.name);
        derivatizes.insert(index, "".to_string());
    }

    let return_data: Vec<Vec<String>> = vec![ids, names, derivatizes];
    return_data
}

pub fn remove_row_from_user_matrices(rowid: usize) -> usize {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    let mut matrix_name = String::from("");
    let mut stmt = conn
        .prepare(&format!(
            "SELECT matrix FROM user_matrices WHERE id = {}",
            rowid
        ))
        .expect("Query cannot be run");
    let db_iter = stmt
        .query_map([], |row| Ok(row.get(0).unwrap_or("".to_string())))
        .unwrap();
    for (_index, item) in db_iter.enumerate() {
        matrix_name = item.unwrap()
    }
    let sql: String = format!("DELETE FROM user_adducts WHERE mname = '{}'", matrix_name);
    match conn.execute(&sql[..], params![]) {
        Ok(_) => println!("Columns from {:?} have been deleted", matrix_name),
        Err(e) => println!("Column could not be deleted: {}", e),
    };

    remove_row_from_table(rowid, "user_matrices");

    return 1;
}

struct MS1UserFgsRow {
    id: usize,
    name: String,
    smarts: String,
}

pub fn update_user_fgs() -> Vec<Vec<String>> {
    check_if_table_exists("functional_group_smarts", "user_functional_group_smarts").unwrap();
    check_if_table_exists("functional_groups", "user_functional_groups").unwrap();
    check_if_table_exists("matrices", "user_matrices").unwrap();

    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    ensure_id_column(&conn, "user_functional_group_smarts").unwrap();

    let mut stmt: rusqlite::Statement = conn
        .prepare("SELECT id, name, smarts FROM user_functional_group_smarts")
        .expect("Query cannot be run");

    let db_iter = stmt
        .query_map([], |row: &rusqlite::Row<'_>| {
            Ok(MS1UserFgsRow {
                id: row.get(0).unwrap_or(1),
                name: row.get(1).unwrap_or("".to_string()),
                smarts: row.get(2).unwrap_or("".to_string()),
            })
        })
        .unwrap();

    let mut ids: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut smarts: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        let row: MS1UserFgsRow = item.unwrap();
        ids.insert(index, row.id.to_string());
        names.insert(index, row.name);
        smarts.insert(index, row.smarts);
    }
    println!("{:?}", smarts);

    let return_data: Vec<Vec<String>> = vec![ids, names, smarts];
    return_data
}

pub fn remove_user_fgs(rowid: usize, column_to_remove: &str) -> usize {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    let sql: &str = "DELETE FROM user_functional_group_smarts WHERE id = ?1";
    conn.execute(sql, params![rowid]).unwrap();

    println!("column_to_remove: {}", column_to_remove);
    match remove_column_from_table("matrices", column_to_remove) {
        Ok(_) => (),
        Err(e) => eprintln!("{:?}", e),
    };

    match remove_column_from_table("user_matrices", column_to_remove) {
        Ok(_) => (),
        Err(e) => eprintln!("{:?}", e),
    };

    match remove_column_from_table("functional_groups", column_to_remove) {
        Ok(_) => (),
        Err(e) => eprintln!("{:?}", e),
    };
    match remove_column_from_table("user_functional_groups", column_to_remove) {
        Ok(_) => (),
        Err(e) => eprintln!("{:?}", e),
    };

    match conn.execute(
        &format!(
            "DELETE FROM functional_group_smarts WHERE name = '{}'",
            column_to_remove
        ),
        [],
    ) {
        Ok(_) => (),
        Err(e) => eprintln!("{:?}", e),
    };

    return 1;
}

fn remove_row_from_table(rowid: usize, table_name: &str) -> usize {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    let sql: String = format!("DELETE FROM {} WHERE id = {}", table_name, rowid);
    conn.execute(&sql[..], params![]).unwrap();
    return 1;
}

fn remove_column_from_table(table_name: &str, column_to_drop: &str) -> Result<()> {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    // Step 1: Fetch the column names
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({});", table_name))?;
    let mut columns: Vec<String> = stmt
        .query_map([], |row| row.get::<_, String>(1))? // Get column names
        .filter_map(|res| res.ok()) // Ignore errors
        .filter(|col| col != column_to_drop) // Exclude the column to drop
        .collect();

    columns = columns.iter().map(|item| format!("\"{}\"", item)).collect();

    // Step 2: Build a new table without the dropped column
    let new_table_name = format!("{}_new", table_name);
    let columns_list = columns.join(", "); // List of remaining columns
    let create_new_table_sql = format!(
        "CREATE TABLE {new_table} AS SELECT {columns} FROM {old_table};",
        new_table = new_table_name,
        columns = columns_list,
        old_table = table_name
    );

    conn.execute(&create_new_table_sql, [])?;

    // Step 3: Drop the old table
    conn.execute(&format!("DROP TABLE {};", table_name), [])?;

    // Step 4: Rename the new table to the original table name
    conn.execute(
        &format!("ALTER TABLE {} RENAME TO {};", new_table_name, table_name),
        [],
    )?;

    println!(
        "Successfully dropped column '{}' from table '{}'.",
        column_to_drop, table_name
    );
    Ok(())
}

fn ensure_id_column(conn: &Connection, table_name: &str) -> Result<()> {
    println!("ensuring id column for {:?}", table_name);

    let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table_name))?;
    let columns = stmt.query_map([], |row| Ok(row.get::<_, String>(1)?))?;

    let mut column_names = Vec::new();
    let mut has_id_column = false;
    for column in columns {
        let name = column?;
        if name == "id" {
            has_id_column = true;
        } else {
            let mut name_to_push: String = "'".to_string();
            name_to_push += &name[..];
            name_to_push += "'";

            column_names.push(name_to_push);
        }
    }

    if !has_id_column {
        let new_table_name = format!("{}_new", table_name);
        let new_table_columns = format!(
            "id INTEGER PRIMARY KEY AUTOINCREMENT, {}",
            column_names.join(", ")
        );
        conn.execute(
            &format!("CREATE TABLE {} ({})", new_table_name, new_table_columns),
            [],
        )?;

        let old_table_columns = column_names.join(", ");
        conn.execute(
            &format!(
                "INSERT INTO {} ({}) SELECT {} FROM {}",
                new_table_name, old_table_columns, old_table_columns, table_name
            ),
            [],
        )?;

        conn.execute(&format!("DROP TABLE {}", table_name), [])?;
        conn.execute(
            &format!("ALTER TABLE {} RENAME TO {}", new_table_name, table_name),
            [],
        )?;
    }

    Ok(())
}

pub fn check_fg_duplicate(name: String) -> bool {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    // Function to check if a column exists in a table
    fn column_exists(conn: &Connection, table_name: &str, column_name: String) -> Result<bool> {
        let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table_name))?;
        let column_names: Vec<String> = stmt
            .query_map([], |row| row.get::<_, String>(1))? // Column 1 contains the column name
            .filter_map(Result::ok)
            .collect();
        Ok(column_names.contains(&column_name.to_string()))
    }

    // Example usage
    if column_exists(&conn, "functional_groups", name).unwrap() {
        return true;
    }

    return false;
}
