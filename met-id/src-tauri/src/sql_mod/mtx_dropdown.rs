use crate::database_setup::get_connection;
use rusqlite::Result;
use std::collections::HashMap;

fn create_hashmap(
    user_mtx: &str,
) -> Result<HashMap<String, Vec<String>>, Box<dyn std::error::Error>> {
    // Connect to the SQLite database
    let conn = get_connection()?; // Assuming `get_connection` returns a Result with a connection

    // Get column names, excluding 'id' and 'matrix'
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", user_mtx)[..])?;
    let columns_iter = stmt.query_map([], |row| Ok(row.get::<_, String>(1)?))?;

    let mut columns: Vec<String> = Vec::new();
    for column in columns_iter {
        let col_name: String = column?;
        if col_name != "id" && col_name != "matrix" {
            columns.push(col_name);
        }
    }

    // Construct and execute the query to fetch all relevant columns
    let query = format!(
        "SELECT matrix, {} FROM {}",
        columns
            .iter()
            .map(|col| format!("\"{}\"", col))
            .collect::<Vec<String>>()
            .join(", "),
        user_mtx
    );
    let mut stmt = conn.prepare(&query)?;
    let rows = stmt.query_map([], |row| {
        let matrix: String = row.get(0)?;
        let mut col_values: Vec<String> = Vec::new();
        for (i, col) in columns.iter().enumerate() {
            // Attempt to retrieve the value as an Option<i32> first
            let int_value_result = row.get::<_, Option<i32>>(i + 1);

            // Check if we got an i32 or need to try parsing a String
            let should_add_col = if let Ok(Some(int_value)) = int_value_result {
                // We successfully got an i32, now check if it's > 0
                int_value > 0
            } else {
                // If we didn't get an i32, try getting a String and parsing it
                if let Ok(Some(string_value)) = row.get::<_, Option<String>>(i + 1) {
                    // Attempt to parse the string as i32 and check if it's > 0
                    string_value
                        .parse::<i32>()
                        .map_or(false, |parsed| parsed > 0)
                } else {
                    // If neither an i32 nor a String, or any error occurs, do not add the column
                    false
                }
            };

            // If the value is an i32 and > 0, or a String that can be parsed to an i32 > 0, add the column name
            if should_add_col {
                col_values.push(col.clone());
            }
        }
        Ok((matrix, col_values))
    })?;

    // Build the hashmap
    let mut hashmap: HashMap<String, Vec<String>> = HashMap::new();
    for row in rows {
        let (matrix, col_values) = row?;
        hashmap.insert(matrix, col_values);
    }

    Ok(hashmap)
}

fn create_hashmap2() -> Result<HashMap<String, Vec<String>>> {
    // Connect to the SQLite database
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> =
        get_connection().unwrap();

    // Update the SQL query to include a WHERE clause that filters the rows
    // where "foo" is either "foo1" or "foo2"
    let mut stmt = conn.prepare(
        "SELECT mname, adduct FROM adducts WHERE mname IN ('Positive Mode', 'Negative Mode')",
    )?;
    let rows = stmt.query_map([], |row| {
        Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
    })?;

    let mut data_map: HashMap<String, Vec<String>> = HashMap::new();

    for row in rows {
        let (foo, bar) = row?;
        data_map.entry(foo).or_insert_with(Vec::new).push(bar);
    }
    Ok(data_map)
}

pub fn matrix_dropdown() -> HashMap<String, Vec<String>> {
    let a: HashMap<String, Vec<String>> = match create_hashmap("matrices") {
        Ok(hashmap) => hashmap,
        Err(e) => {
            eprintln!("Error creating hashmap: {}", e);
            HashMap::new()
        }
    };

    let b: HashMap<String, Vec<String>> = match create_hashmap("user_matrices") {
        Ok(hashmap) => hashmap,
        Err(e) => {
            eprintln!("Error creating hashmap: {}", e);
            HashMap::new()
        }
    };

    let c: HashMap<String, Vec<String>> = match create_hashmap2() {
        Ok(hashmap) => hashmap,
        Err(e) => {
            eprintln!("Error creating hashmap: {}", e);
            HashMap::new()
        }
    };

    let d: HashMap<String, Vec<String>> = a.into_iter().chain(b).chain(c).collect();
    d
}
