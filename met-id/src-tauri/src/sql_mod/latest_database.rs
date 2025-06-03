use crate::add_to_db::lipids_functional_groups::complete_lipids_functional_groups;
use crate::database_setup::get_connection;
use rusqlite::Result;
use std::collections::HashMap;

use super::table::check_if_table_exists;

pub fn check_latest_database() -> Result<()> {
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();
    match change_lipidmaps_column_type(&conn) {
        Ok(_) => println!("lipid table updated"),
        Err(e) => println!("lipid table up to date: {:?}", e),
    };
    match change_metabolite_column_type(&conn) {
        Ok(_) => println!("metabolite table updated"),
        Err(e) => println!("metabolite table up to date: {:?}", e),
    };

    match lipid_functional_groups_table(&conn) {
        Ok(_) => println!("lipid functional groups updated"),
        Err(e) => println!("lipid functional groups up to date: {:?}", e),
    }

    Ok(())
}

fn change_lipidmaps_column_type(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
) -> Result<()> {
    // Step 1: Create a new table with the updated column type

    //conn.execute("CREATE INDEX id")
    conn.execute(
        "CREATE TABLE new_lipids (
            id INTEGER PRIMARY KEY,
            name TEXT,
            mz REAL,
            formula TEXT,
            smiles TEXT
        )",
        [],
    )?;

    // Step 2: Copy data from the old table to the new table
    conn.execute(
        "INSERT INTO new_lipids (name, mz, formula, smiles)
        SELECT name, mz, formula, smiles FROM lipids",
        [],
    )?;

    // Step 3: Drop the old table
    conn.execute("DROP TABLE lipids", [])?;

    // Step 4: Rename the new table to the old name
    conn.execute("ALTER TABLE new_lipids RENAME TO lipids", [])?;

    Ok(())
}

fn change_metabolite_column_type(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
) -> Result<()> {
    // Step 1: Create a new table with the updated column type

    //conn.execute("CREATE INDEX id")
    conn.execute(
        "CREATE TABLE new_metabolites (
            id INTEGER PRIMARY KEY,
            name TEXT,
            mz REAL,
            formula TEXT,
            smiles TEXT
        )",
        [],
    )?;

    // Step 2: Copy data from the old table to the new table
    conn.execute(
        "INSERT INTO new_metabolites (name, mz, formula, smiles)
        SELECT name, mz, formula, smiles FROM metabolites",
        [],
    )?;

    // Step 3: Drop the old table
    conn.execute("DROP TABLE metabolites", [])?;

    // Step 4: Rename the new table to the old name
    conn.execute("ALTER TABLE new_metabolites RENAME TO metabolites", [])?;

    Ok(())
}

fn table_exists(conn: &rusqlite::Connection, table_name: &str) -> rusqlite::Result<bool> {
    let mut stmt =
        conn.prepare("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?1")?;
    let count: i64 = stmt.query_row([table_name], |row| row.get(0))?;
    Ok(count > 0)
}

fn column_count(conn: &rusqlite::Connection, table_name: &str) -> rusqlite::Result<usize> {
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table_name))?;
    let count = stmt.query_map([], |_| Ok(()))?.count();
    Ok(count)
}

fn has_same_column_count(
    conn: &rusqlite::Connection,
    table1: &str,
    table2: &str,
) -> rusqlite::Result<bool> {
    let count1 = column_count(conn, table1)?;
    let count2 = column_count(conn, table2)?;

    Ok(count1 == count2)
}

fn columns_with_empty_values(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
    table_name: &str,
) -> Result<Vec<String>> {
    // Step 1: Get column names from PRAGMA
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table_name))?;
    let column_names: Vec<String> = stmt
        .query_map([], |row| row.get(1))? // Column 1 is the name
        .collect::<Result<Vec<_>, _>>()?;

    let mut columns_with_empty = Vec::new();

    // Step 2: For each column, check if it has NULL or empty string
    for col in &column_names {
        let query = format!(
            "SELECT EXISTS(SELECT 1 FROM {} WHERE \"{}\" IS NULL OR \"{}\" = '')",
            table_name, col, col
        );

        let has_empty: bool = conn.query_row(&query, [], |row| row.get(0))?;
        if has_empty {
            columns_with_empty.push(col.clone());
        }
    }

    Ok(columns_with_empty)
}

fn add_missing_columns(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
    from_table: &str,
    to_table: &str,
) -> Result<Vec<String>> {
    // Get column info from both tables
    let get_columns = |table: &str| -> Result<Vec<(String, String)>> {
        let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table))?;
        let cols = stmt
            .query_map([], |row| {
                let name: String = row.get(1)?; // column name
                let col_type: String = row.get(2)?; // type
                Ok((name, col_type))
            })?
            .collect::<Result<Vec<_>, _>>()?;
        Ok(cols)
    };

    let from_columns = get_columns(from_table)?;
    let to_columns: Vec<String> = get_columns(to_table)?
        .into_iter()
        .map(|(name, _)| name)
        .collect();

    // Add any missing columns
    let mut altered_columns = Vec::new();
    for (col_name, col_type) in from_columns {
        if !to_columns.contains(&col_name) {
            let alter_sql = format!(
                "ALTER TABLE {} ADD COLUMN '{}' {}",
                to_table, col_name, col_type
            );
            altered_columns.push(col_name);
            conn.execute(&alter_sql, [])?;
        }
    }

    Ok(altered_columns)
}

fn get_columns(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
    table_name: &str,
) -> Result<Vec<String>> {
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table_name))?;
    let mut columns = stmt
        .query_map([], |row| row.get(1))? // column index 1 = name
        .collect::<Result<Vec<_>, _>>()?;
    columns.retain(|col| col != "id");
    Ok(columns)
}

fn get_smarts_for_names(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
    names: &Vec<String>,
) -> Result<HashMap<String, String>> {
    if names.is_empty() {
        return Ok(HashMap::new());
    }
    check_if_table_exists("functional_group_smarts", "user_functional_group_smarts").unwrap();

    let csv_quoted = names
        .iter()
        .map(|s| format!("'{}'", s)) // add single quotes
        .collect::<Vec<_>>()
        .join(", ");

    // Combined SQL query for both tables
    let sql = format!(
        "SELECT name, smarts FROM functional_group_smarts WHERE name IN ({})
        UNION
        SELECT name, smarts FROM user_functional_group_smarts WHERE name IN ({})",
        csv_quoted, csv_quoted
    );

    let mut stmt = conn.prepare(&sql)?;

    // We use params_from_iter twice because placeholders appear twice
    let mut params = names.clone();
    params.extend(names.clone());

    let results = stmt
        .query_map([], |row| {
            let name: String = row.get(0)?;
            let smarts: String = row.get(1)?;
            Ok((name, smarts))
        })?
        .collect::<Result<HashMap<_, _>, _>>()?;
    Ok(results)
}

fn lipid_functional_groups_table(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
) -> Result<()> {
    if table_exists(conn, "lipids_functional_groups").unwrap() {
        if has_same_column_count(conn, "lipids_functional_groups", "functional_groups").unwrap() {
            let missing_columns: Vec<String> =
                columns_with_empty_values(conn, "lipids_functional_groups").unwrap();
            let missing_smarts: HashMap<String, String> =
                get_smarts_for_names(conn, &missing_columns).unwrap();
            if !missing_columns.is_empty() {
                complete_lipids_functional_groups(missing_smarts).unwrap();
            } else {
                return Ok(());
            }
        } else {
            let missing_columns: Vec<String> =
                add_missing_columns(conn, "functional_groups", "lipids_functional_groups").unwrap();
            let missing_smarts: HashMap<String, String> =
                get_smarts_for_names(conn, &missing_columns).unwrap();
            complete_lipids_functional_groups(missing_smarts).unwrap();
        }
    } else {
        conn.execute(
            "CREATE TABLE lipids_functional_groups AS SELECT * FROM functional_groups WHERE 0",
            [],
        )
        .unwrap();
        conn.execute("INSERT INTO lipids_functional_groups (id) SELECT row_number FROM (SELECT ROW_NUMBER() OVER () AS row_number FROM lipids)", []).unwrap();
        let missing_columns: Vec<String> = get_columns(conn, "lipids_functional_groups").unwrap();
        let missing_smarts: HashMap<String, String> =
            get_smarts_for_names(conn, &missing_columns).unwrap();
        complete_lipids_functional_groups(missing_smarts).unwrap();
    }
    Ok(())
}
