use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::ToSql;
use std::collections::HashMap;
use std::iter::repeat;

use super::get_table_column_names;

pub fn fill_user_derivatized_by(
    conn: &r2d2::PooledConnection<SqliteConnectionManager>,
    fgh: &HashMap<String, String>,
    matrices_fgs: &HashMap<String, Vec<String>>,
) {
    let mut matrices = get_table_column_names(conn, "derivatized_by").unwrap();
    matrices.remove(0);

    let data: Vec<(String, String)> = calculate_category_sums(fgh, matrices_fgs);

    // Generate the column names and placeholder strings dynamically
    let column_names: Vec<String> = data
        .iter()
        .map(|(col, _)| format!("'{}'", col.to_string()))
        .collect();
    let placeholders: Vec<&str> = repeat("?").take(data.len()).collect();

    let sql: String = format!(
        "INSERT INTO user_derivatized_by ({}) VALUES ({})",
        column_names.join(", "),
        placeholders.join(", ")
    );

    let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();

    conn.execute(&sql, values.as_slice()).unwrap();
}

fn calculate_category_sums(
    foo: &HashMap<String, String>,
    bar: &HashMap<String, Vec<String>>,
) -> Vec<(String, String)> {
    let mut result: Vec<(String, String)> = Vec::new();

    for (category, names) in bar {
        let mut sum: i32 = 0;
        for name in names {
            if let Some(value_str) = foo.get(name) {
                if let Ok(value) = value_str.parse::<i32>() {
                    sum += value;
                }
            }
        }
        result.push((category.clone(), sum.to_string()));
    }

    result
}
