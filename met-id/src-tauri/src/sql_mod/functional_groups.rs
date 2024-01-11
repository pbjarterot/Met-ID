use rusqlite::{params, Result};

use crate::database_setup::get_connection;

pub fn get_functional_groups() -> Vec<String> {
  match get_fg() {  // added parentheses here
      Ok(result) => return result,
      Err(_) => return Vec::new(),
  }
}

fn get_fg() -> Result<Vec<String>> {
  let conn = get_connection().unwrap();
  let mut stmt = conn.prepare("PRAGMA table_info(matrices);")?;
  let column_info = stmt.query_map(params![], |row| {
      let name: String = row.get(1)?; // 'name' is the second column, thus index 1
      Ok(name)
  })?;

  let mut column_names = Vec::new();
  for name in column_info {
      if let Ok(name) = name {
          column_names.push(name);
      }
  }
  let value_to_remove = "matrix".to_string();
  
  if let Some(pos) = column_names.iter().position(|x| *x == value_to_remove) {
      column_names.remove(pos);
  }
  Ok(column_names)
}