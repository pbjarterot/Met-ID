use crate::database_setup::get_connection;
use rusqlite::params;

pub fn get_tissues() -> Vec<String> {
    let conn = get_connection().unwrap();
    let mut stmt = conn.prepare("PRAGMA table_info(in_tissue);").unwrap();
    let column_info = stmt
        .query_map(params![], |row| {
            let name: String = row.get(1)?; // 'name' is the second column, thus index 1
            Ok(name)
        })
        .unwrap();

    let mut column_names = Vec::new();
    for name in column_info {
        if let Ok(name) = name {
            column_names.push(name);
        }
    }
    let value_to_remove = "id".to_string();

    if let Some(pos) = column_names.iter().position(|x| *x == value_to_remove) {
        column_names.remove(pos);
    }

    column_names
}
