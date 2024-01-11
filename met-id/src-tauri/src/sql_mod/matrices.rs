use crate::database_setup::get_connection;
use rusqlite::{params, Result};


fn get_mtx() -> Result<Vec<String>> {
  let conn = get_connection().unwrap();
  let mut stmt = conn.prepare("SELECT matrix FROM matrices WHERE matrix NOT IN ('POSITIVE_MODE', 'NEGATIVE_MODE')").unwrap();
  let name_iter = stmt.query_map(params![], |row| {
      Ok(row.get::<_, String>(0).unwrap())
  })?;
  
  // Collect names into a vector, handling errors as they occur
  let names = name_iter.collect();
  names
}


pub fn get_matrices() -> Vec<String> {
  match get_mtx() {  // added parentheses here
      Ok(result) => return result,
      Err(_) => return Vec::new(),
  }
}