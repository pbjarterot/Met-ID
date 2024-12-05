use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::ToSql;
use std::collections::HashMap;
use std::iter::repeat;

pub fn fill_user_in_tissue(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
    in_tissue: &HashMap<String, bool>,
) {
    let column_names: Vec<String> = in_tissue.iter().map(|(col, _)| col.to_string()).collect();
    let placeholders: Vec<&str> = repeat("?").take(in_tissue.len()).collect();

    let sql = format!(
        "INSERT INTO user_in_tissue ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    let values: Vec<&dyn ToSql> = in_tissue
        .iter()
        .map(|(_, value)| value as &dyn ToSql)
        .collect();

    conn.execute(&sql, values.as_slice()).unwrap();
}
