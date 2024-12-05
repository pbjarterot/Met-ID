use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::ToSql;

use crate::{
    database_setup::get_msms_connection, files::read_mzml_for_msms_to_add_to_db,
    sql_mod::table::check_if_table_exists_msms,
};

pub fn fill_user_msms(bin_data: Vec<u8>) {
    // Convert &[u8] to Vec<u8>
    let bin_data_vec = bin_data.to_vec();

    // Change data to Vec of tuples with String and Box<dyn ToSql>
    let data: Vec<(String, Box<dyn ToSql>)> = vec![
        ("name".to_string(), Box::new("User Input".to_string())),
        ("identifier".to_string(), Box::new("User Input".to_string())),
        ("adduct".to_string(), Box::new("1A".to_string())),
        ("cid".to_string(), Box::new("0eV".to_string())),
        ("window".to_string(), Box::new("1Da".to_string())),
        ("tof".to_string(), Box::new("0.0".to_string())),
        ("mz".to_string(), Box::new("0.0".to_string())),
        ("spectra".to_string(), Box::new(bin_data_vec)),
    ];

    let conn: r2d2::PooledConnection<SqliteConnectionManager> =
        get_msms_connection("connecting to MSMS db").unwrap();

    // Start a transaction
    conn.execute("BEGIN", []).unwrap();

    // Prepare an update statement
    let update_sql = format!(
        "UPDATE MSMS SET {} WHERE identifier = ?",
        data.iter()
            .map(|(col, _)| format!("{} = ?", col))
            .collect::<Vec<_>>()
            .join(", ")
    );
    let binding = "User Input".to_string();
    let update_values: Vec<&dyn ToSql> = data
        .iter()
        .map(|(_, value)| &**value)
        .chain(std::iter::once(&binding as &dyn ToSql))
        .collect();
    let updated_rows = conn.execute(&update_sql, &*update_values).unwrap();

    if updated_rows == 0 {
        // If no rows were updated, do an insert
        let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
        let placeholders: Vec<&str> = std::iter::repeat("?").take(data.len()).collect();

        let insert_sql = format!(
            "INSERT INTO MSMS ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        // Map the data to a vector of references to the trait object
        let insert_values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| &**value).collect();
        conn.execute(&insert_sql, &*insert_values).unwrap();
    }

    // Commit the transaction
    conn.execute("COMMIT", []).unwrap();
}

pub fn add_to_usermsms(
    name: String,
    adduct: String,
    mz: String,
    cid: String,
    tof: String,
    mzwindow: String,
    identifier: String,
    path: String,
    matrix: String,
) {
    check_if_table_exists_msms("MSMS", "user_MSMS").unwrap();
    // Convert &[u8] to Vec<u8>
    //let bin_data_vec = bin_data.to_vec();
    let bin_data_vec = read_mzml_for_msms_to_add_to_db(path);

    // Change data to Vec of tuples with String and Box<dyn ToSql>
    let data: Vec<(String, Box<dyn ToSql>)> = vec![
        ("name".to_string(), Box::new(name)),
        ("identifier".to_string(), Box::new(identifier)),
        ("adduct".to_string(), Box::new(adduct)),
        ("cid".to_string(), Box::new(cid)),
        ("window".to_string(), Box::new(mzwindow)),
        ("tof".to_string(), Box::new(tof)),
        ("mz".to_string(), Box::new(mz)),
        ("spectra".to_string(), Box::new(bin_data_vec)),
        ("matrix".to_string(), Box::new(matrix)),
    ];

    let conn: r2d2::PooledConnection<SqliteConnectionManager> =
        get_msms_connection("connecting to MSMS db").unwrap();

    // If no rows were updated, do an insert
    let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
    let placeholders: Vec<&str> = std::iter::repeat("?").take(data.len()).collect();

    let sql = format!(
        "INSERT INTO user_MSMS ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    // Map the data to a vector of references to the trait object
    let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();
    conn.execute(&sql, values.as_slice()).unwrap();
}
