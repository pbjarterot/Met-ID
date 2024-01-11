use crate::database_setup::get_connection;
use super::MS1DbRow;



pub fn sql_query(query: &String) -> (Vec<f64>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>) {
  //connect to db
  let conn = get_connection().unwrap();
  //query
  let mut stmt: rusqlite::Statement = conn.prepare(query).expect("Query cannot be run");
  
  //collect results from the database query
  
  let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
      Ok(MS1DbRow {
          mz: row.get(0).unwrap_or(0.0 as f64),
          name: row.get(1).unwrap_or("".to_string()),
          mname: row.get(2).unwrap_or("".to_string()),
          accession: row.get(3).unwrap_or("".to_string()),
          smiles: row.get(4).unwrap_or("".to_string()),
          formula: row.get(5).unwrap_or("".to_string())
      })
  }).unwrap();

  let mut mzs: Vec<f64> = Vec::new();
  let mut names: Vec<String> = Vec::new();
  let mut mnames: Vec<String> = Vec::new();
  let mut accessions: Vec<String> = Vec::new();
  let mut smiless: Vec<String> = Vec::new();
  let mut formulas: Vec<String> = Vec::new();

  //parse results for passing back to the parent function
  for (index, item) in db_iter.enumerate() {

      let row: MS1DbRow = item.unwrap();
      mzs.insert(index, row.mz);
      names.insert(index, row.name);
      mnames.insert(index, row.mname);
      accessions.insert(index, row.accession);
      smiless.insert(index, row.smiles);
      formulas.insert(index, row.formula);

  }
  (mzs, names, mnames, accessions, smiless, formulas)
}