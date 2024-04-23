use std::collections::HashMap;
use crate::database_setup::get_connection;


struct DBData {
  adjusted_mz: f64,
  name: String,
  adduct: String,
  db_accession: String,
  smiles: String,
  chemicalformula: String,
  mname: String
}
struct DBNamesIDs {
  name: String,
  id: usize
}

pub fn get_db_data(index:usize) -> (String, String, String, String, HashMap<String, HashMap<String, f64>>) {
  let conn = get_connection().unwrap();
  let mut query: String = "SELECT CAST(metabolites.mz AS REAL) + CAST(adducts.deltamass AS REAL) AS adjusted_mz, ".to_string();
  query += "metabolites.name, adducts.adduct, db_accessions.hmdb, metabolites.smiles, metabolites.chemicalformula, adducts.mname FROM metabolites ";
  query += "INNER JOIN db_accessions ON metabolites.id = db_accessions.id INNER JOIN endogeneity ON metabolites.id = endogeneity.id ";
  query += "INNER JOIN functional_groups ON metabolites.id = functional_groups.id INNER JOIN derivatized_by ON metabolites.id = derivatized_by.id ";
  query += "INNER JOIN in_tissue ON metabolites.id = in_tissue.id ";
  query += "CROSS JOIN (SELECT adduct FROM adducts) as m LEFT JOIN adducts on m.adduct = adducts.adduct ";
  query += &format!("WHERE adducts.numfunctionalgroups <= 100 AND metabolites.id = {} ORDER BY adjusted_mz", index);
  //println!("{:?}", query);


  let mut stmt: rusqlite::Statement = conn.prepare(&query).expect("Query cannot be run");
  let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
      Ok(DBData {
          adjusted_mz: row.get(0).unwrap_or(0.0),
          name: row.get(1).unwrap_or("".to_string()),
          adduct: row.get(2).unwrap_or("".to_string()),
          db_accession: row.get(3).unwrap_or("".to_string()),
          smiles: row.get(4).unwrap_or("".to_string()),
          chemicalformula: row.get(5).unwrap_or("".to_string()),
          mname: row.get(6).unwrap_or("".to_string()),
      })
  }).unwrap();

  let mut name_str: String = String::new();
  let mut db_accession_str: String = String::new();
  let mut smiles_str: String = String::new();
  let mut formula_str: String = String::new();
  let mut mname_vec: Vec<String> = Vec::new();
  let mut adduct_vec: Vec<String> = Vec::new();
  let mut adj_mz_vec: Vec<f64> = Vec::new();


  for (index, item) in db_iter.enumerate() {
    let row: DBData = item.unwrap();
    if index == 0 {
      name_str = row.name;
      db_accession_str = row.db_accession;
      smiles_str = row.smiles;
      formula_str = row.chemicalformula;
    }
    mname_vec.insert(index, row.mname);
    adduct_vec.insert(index, row.adduct.clone());
    adj_mz_vec.insert(index, row.adjusted_mz);
    //println!("{:?}, {:?}, {:?}, {:?}", row.adduct, row.adjusted_mz, row.adduct, row.mname);
  }

  // Create an empty HashMap to store your category and sub-HashMaps.
  let mut main_hashmap: HashMap<String, HashMap<String, f64>> = HashMap::new();

   // Iterate over the indices of the categories vector.
   for i in 0..mname_vec.len() {
    // Extract the current category, key, and value.
    let category = &mname_vec[i];
    let key = &adduct_vec[i];
    let value = adj_mz_vec[i];

    // Use the entry API to handle the insertion more cleanly.
    main_hashmap.entry(category.to_string())
        .or_insert_with(HashMap::new) // This will create a new sub-HashMap if the category does not exist.
        .insert(key.to_string(), value);         // Inserts the key-value pair into the sub-HashMap.
}

  // Print the resulting HashMap of HashMaps to verify the structure.
  for (category, sub_map) in &main_hashmap {
      println!("Category: {}", category);
      for (key, value) in sub_map {
          println!("  {} -> {}", key, value);
      }
  }


  (name_str, db_accession_str, smiles_str, formula_str, main_hashmap)
}


pub fn db_ids_and_names(inputvalue: String) -> HashMap<String, usize> {
  let conn = get_connection().unwrap();
  let mut stmt: rusqlite::Statement = conn.prepare("SELECT name, id FROM metabolites WHERE name LIKE ?1 ORDER BY id LIMIT 1000").expect("Query cannot be run");
  let db_iter = stmt.query_map([format!("%{}%", inputvalue)], |row: &rusqlite::Row<'_>| {
      Ok(DBNamesIDs {
          name: row.get(0).unwrap_or("".to_string()),
          id: row.get(1).unwrap_or(0),
      })
  }).unwrap();
  let mut map = HashMap::new();

  for (index, item) in db_iter.enumerate() {
    let row: DBNamesIDs = item.unwrap();
    map.insert(row.name, row.id);
  }
  map
}