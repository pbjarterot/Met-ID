use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::ToSql;
use std::collections::HashMap;
use std::iter::repeat;

pub fn fill_user_endogeneity(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
    endo_exo: &HashMap<String, bool>,
) {
    let column_names: Vec<String> = endo_exo.iter().map(|(col, _)| col.to_string()).collect();
    let placeholders: Vec<&str> = repeat("?").take(endo_exo.len()).collect();

    let sql = format!(
        "INSERT INTO user_endogeneity ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    let values: Vec<&dyn ToSql> = endo_exo
        .iter()
        .map(|(_, value)| value as &dyn ToSql)
        .collect();

    conn.execute(&sql, values.as_slice()).unwrap();
}
