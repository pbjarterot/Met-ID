use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, Connection, Result};

use crate::database_setup::get_connection;

use super::table::check_if_table_exists;

struct MS1UserMetabolitesRow {
  id: usize,
  name: String, 
  smiles: String,
  formula: String,
  mz: String
}


pub fn update_user_metabolites() -> Vec<Vec<String>> {
  check_if_table_exists("metabolites",       "user_metabolites").unwrap();
  check_if_table_exists("derivatized_by",    "user_derivatized_by").unwrap();
  check_if_table_exists("endogeneity",       "user_endogeneity").unwrap();
  check_if_table_exists("in_tissue",         "user_in_tissue").unwrap();
  check_if_table_exists("db_accessions",     "user_db_accessions").unwrap();
  check_if_table_exists("functional_groups", "user_functional_groups").unwrap();
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
  let mut stmt: rusqlite::Statement = conn.prepare("SELECT id, name, smiles, chemicalformula, mz FROM user_metabolites").expect("Query cannot be run");
  let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
      Ok(MS1UserMetabolitesRow {
          id: row.get(0).unwrap(),
          name: row.get(1).unwrap_or("".to_string()),
          smiles: row.get(2).unwrap_or("".to_string()),
          formula: row.get(3).unwrap_or("".to_string()),
          mz: row.get(4).unwrap_or("".to_string()),
      })
  }).unwrap();

  let mut ids: Vec<String> = Vec::new();
  let mut names: Vec<String> = Vec::new();
  let mut smiles: Vec<String> = Vec::new();
  let mut formulas: Vec<String> = Vec::new();
  let mut mzs: Vec<String> = Vec::new();

  //parse results for passing back to the parent function
  for (index, item) in db_iter.enumerate() {
      let row: MS1UserMetabolitesRow = item.unwrap();
      ids.insert(index, row.id.to_string());
      names.insert(index, row.name);
      smiles.insert(index, row.smiles);
      formulas.insert(index, row.formula);
      mzs.insert(index, row.mz)
  }
  
  let return_data: Vec<Vec<String>> = vec![ids, names, smiles, formulas, mzs];
  return_data
  
}

pub fn remove_row_from_user_metabolites(rowid: usize) -> usize {
  remove_row_from_table(rowid, "user_metabolites");
  remove_row_from_table(rowid, "user_derivatized_by");
  remove_row_from_table(rowid, "user_endogeneity");
  remove_row_from_table(rowid, "user_in_tissue");
  remove_row_from_table(rowid, "user_functional_groups");
  remove_row_from_table(rowid, "user_db_accessions");

  return 1
}


struct MS1UserMatricesRow {
  id: usize,
  name: String, 
}


pub fn update_user_matrices() -> Vec<Vec<String>> {
  check_if_table_exists("matrices", "user_matrices").unwrap();
  check_if_table_exists("adducts", "user_adducts").unwrap();  
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
  ensure_id_column(&conn, "user_matrices").unwrap();
  let mut stmt: rusqlite::Statement = conn.prepare("SELECT id, matrix FROM user_matrices").expect("Query cannot be run");
  let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {

      let id_value = row.get(0).unwrap();
      Ok(MS1UserMatricesRow {
          id: id_value,
          name: row.get(1).unwrap_or("".to_string()),
      })
  }).unwrap();

  let mut ids: Vec<String> = Vec::new();
  let mut names: Vec<String> = Vec::new();
  let mut derivatizes: Vec<String> = Vec::new();

  //parse results for passing back to the parent function
  for (index, item) in db_iter.enumerate() {
      let row: MS1UserMatricesRow = item.unwrap();
      ids.insert(index, row.id.to_string());
      names.insert(index, row.name);
      derivatizes.insert(index, "".to_string());
  }
  
  let return_data: Vec<Vec<String>> = vec![ids, names, derivatizes];
  return_data
  
}



pub fn remove_row_from_user_matrices(rowid: usize) -> usize {
	let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
	let mut matrix_name = String::from("");
	let mut stmt = conn.prepare(&format!("SELECT matrix FROM user_matrices WHERE id = {}", rowid)).expect("Query cannot be run");
	let db_iter = stmt.query_map([], |row| {
		Ok(row.get(0).unwrap_or("".to_string()))
	}).unwrap();
	for (_index, item) in db_iter.enumerate() {
		matrix_name = item.unwrap()
	};
	let sql: String  = format!("DELETE FROM user_adducts WHERE mname = '{}'", matrix_name);
	match conn.execute(&sql[..], params![]) {
		Ok(_) => println!("Columns from {:?} have been deleted", matrix_name),
		Err(e) => println!("Column could not be deleted: {}", e)
	};

	remove_row_from_table(rowid, "user_matrices");

	return 1
}

struct MS1UserFgsRow {
  id: usize,
  name: String, 
  smarts: String,
}

pub fn update_user_fgs() -> Vec<Vec<String>> {
  check_if_table_exists("functional_group_smarts", "user_functional_group_smarts").unwrap();
  check_if_table_exists("functional_groups",       "user_functional_groups").unwrap();
  check_if_table_exists("matrices",                "user_matrices").unwrap();
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
  ensure_id_column(&conn, "user_functional_group_smarts").unwrap();
  let mut stmt: rusqlite::Statement = conn.prepare("SELECT id, name, smarts FROM user_functional_group_smarts").expect("Query cannot be run");
  let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
      Ok(MS1UserFgsRow {
          id: row.get(0).unwrap_or(1),
          name: row.get(1).unwrap_or("".to_string()),
          smarts: row.get(2).unwrap_or("".to_string()),
      })
  }).unwrap();

  let mut ids: Vec<String> = Vec::new();
  let mut names: Vec<String> = Vec::new();
  let mut smarts: Vec<String> = Vec::new();
  let mut derivatizes2: Vec<String> = Vec::new();

  //parse results for passing back to the parent function
  for (index, item) in db_iter.enumerate() {
      let row: MS1UserFgsRow = item.unwrap();
      ids.insert(index, row.id.to_string());
      names.insert(index, row.name);
      smarts.insert(index, row.smarts);
      derivatizes2.insert(index, "".to_string());
  }
  println!("{:?}", smarts);
  
  let return_data: Vec<Vec<String>> = vec![ids, names, smarts, derivatizes2];
  return_data
  
}




pub fn remove_user_fgs(rowid: usize, column_to_remove: &str) -> usize {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
  let sql: &str = "DELETE FROM user_functional_group_smarts WHERE id = ?1";
  conn.execute(sql, params![rowid]).unwrap();

  remove_column_from_table("matrices", column_to_remove).unwrap();
  remove_column_from_table("matrices", column_to_remove).unwrap();

  remove_column_from_table("functional_groups", column_to_remove).unwrap();
  remove_column_from_table("user_functional_groups", column_to_remove).unwrap();

  return 1
}










fn remove_row_from_table(rowid: usize, table_name: &str) -> usize {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
  let sql: String  = format!("DELETE FROM {} WHERE id = {}", table_name, rowid);
  conn.execute(&sql[..], params![]).unwrap();
  return 1
}


fn remove_column_from_table(table_name: &str, column_to_remove: &str) -> Result<()> {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    // Get the column names
    let mut stmt = conn.prepare("PRAGMA table_info(?1);")?;
    let columns = stmt.query_map(params![table_name], |row| {
        Ok(row.get::<_, String>(1)?)
    })?;

    // Filter out the column to remove and join the remaining columns
    let filtered_columns: Vec<String> = columns
        .filter_map(Result::ok)
        .filter(|col| col != column_to_remove)
        .collect();
    let column_list = filtered_columns.join(", ");

    // Start a transaction
    conn.execute("BEGIN TRANSACTION;", params![])?;

    // Create a temporary table
    conn.execute(
        &format!("CREATE TEMPORARY TABLE ?1_backup AS SELECT * FROM ?2 WHERE 1=0;"),
        params![table_name, table_name],
    )?;

    // Copy data to temporary table
    conn.execute(
        &format!("INSERT INTO bar_backup SELECT {} FROM {};", column_list, table_name),
        params![],
    )?;

    // Drop the original table
    conn.execute("DROP TABLE ?1;", params![table_name])?;

    // Rename the temporary table
    conn.execute("ALTER TABLE ?1_backup RENAME TO ?2;", params![table_name, table_name])?;

    // Commit the transaction
    conn.execute("COMMIT;", params![])?;

    Ok(())
}
























fn ensure_id_column(conn: &Connection, table_name: &str) -> Result<()> {
  let mut stmt = conn.prepare(&format!("PRAGMA table_info({})", table_name))?;
  let columns = stmt.query_map([], |row| Ok(row.get::<_, String>(1)?))?;

  let mut column_names = Vec::new();
  let mut has_id_column = false;
  for column in columns {
      let name = column?;
      if name == "id" {
          has_id_column = true;
      } else {
          let mut name_to_push: String = "'".to_string();
          name_to_push += &name[..];
          name_to_push += "'";

          column_names.push(name_to_push);
      }
  }

  if !has_id_column {
      let new_table_name = format!("{}_new", table_name);
      let new_table_columns = format!("id INTEGER PRIMARY KEY AUTOINCREMENT, {}", column_names.join(", "));
      conn.execute(&format!("CREATE TABLE {} ({})", new_table_name, new_table_columns), [])?;

      let old_table_columns = column_names.join(", ");
      conn.execute(&format!("INSERT INTO {} ({}) SELECT {} FROM {}", new_table_name, old_table_columns, old_table_columns, table_name), [])?;

      conn.execute(&format!("DROP TABLE {}", table_name), [])?;
      conn.execute(&format!("ALTER TABLE {} RENAME TO {}", new_table_name, table_name), [])?;
  }

  Ok(())
}