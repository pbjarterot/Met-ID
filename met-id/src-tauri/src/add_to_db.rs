
pub mod add_to_db_functions {
    use r2d2_sqlite::SqliteConnectionManager;
    use rusqlite::{params, Result, ToSql};
    //use crate::validation::mass_from_formula;
    use crate::sql_mod::table::{check_if_table_exists, check_if_table_exists_msms};
    use crate::sidecar::sidecar_function;
    //use crate::metabolite::{ Metabolite, single_functional_group };
    use crate::database_setup::{get_connection, get_msms_connection};
    use std::collections::HashMap;
    use std::iter::repeat;
    use std::sync::mpsc;
    use crate::files::read_mzml_for_msms_to_add_to_db;

    fn single_functional_group(smiles: &mut Vec<String>, smarts: &String) -> Vec<usize> {
        smiles.insert(0, smarts.to_owned());
        let sidecar_output: std::prelude::v1::Result<String, crate::sidecar::CommandError> = sidecar_function("metabolite".to_string(), smiles.to_owned());
        let result: Vec<usize> = match sidecar_output {
            Ok(s) => {
                s.trim_end_matches('\r')
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .split(", ")
                    .filter_map(|n| n.parse::<usize>().ok())
                    .collect()
            },
            Err(_) => Vec::new(),
        };
        
        result
    }

    fn metabolite_for_db_sidecar(smiles: &String, smarts_map: &mut HashMap<String, String>) -> (String, String, HashMap<String, String>) {
        let mut smarts_vec: Vec<_> = smarts_map.values().cloned().collect();
        smarts_vec.insert(0, smiles.to_owned());

        let sidecar_output = sidecar_function("metabolite_for_db".to_string(), smarts_vec.to_owned());

        parse_input(&sidecar_output.unwrap()[..], smarts_map).unwrap()
    }

    fn parse_input(input: &str, smarts_map: &mut HashMap<String, String>) -> Result<(String, String, HashMap<String, String>), &'static str> {
        use regex::Regex;
        //let input = "C3H8 44.062600255999996 [0, 0, 0, 0]\r";

        // Extract molecule formula
        let molecule_re = Regex::new(r"([A-Z][a-z]*\d*)+").unwrap();
        let molecule = molecule_re.find(input).unwrap().as_str();

        // Extract float number
        let float_re = Regex::new(r"\d+\.\d+").unwrap();
        let float_number: f64 = float_re.find(input).unwrap().as_str().parse().unwrap();

        // Extract numbers from the array
        let array_start = input.find('[').unwrap();
        let array_end = input.find(']').unwrap();
        let array_str = &input[array_start + 1..array_end];
        let array: Vec<i32> = array_str
            .split(", ")
            .filter_map(|s| s.trim().parse().ok())
            .collect();

        // Assuming keys and values have the same length
        let keys: Vec<String> = smarts_map.keys().cloned().collect();
        for (key, value) in keys.iter().zip(array.iter()) {
            smarts_map.insert(key.clone(), value.clone().to_string());
        }
    
        Ok((molecule.to_string(), float_number.to_string(), smarts_map.to_owned()))
    }

    /* 
    fn update_fg(smarts: &String) {
        const BATCH_SIZE: usize = 50;
        let match_counts: Vec<usize> = single_functional_group(smarts);
    }
    */
     
    pub fn update_functional_groups(table_name: &str, table2_name: &str, name: &String, smarts: &String, progress_sender: &mpsc::Sender<f32>) {
        const BATCH_SIZE: usize = 50;
        
    
        let mut conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
        let conn2: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
        
        conn.execute(
            &format!("ALTER TABLE {} ADD COLUMN {} NUMBER", table_name, name),
            [],
        ).unwrap();
        let mode: String = conn.query_row("PRAGMA journal_mode=WAL;", [], |row| {
            row.get(0)
        }).unwrap();
        
        assert_eq!(mode, "wal");
    
        let mut select_stmt = conn2.prepare(&format!("SELECT id, smiles FROM {}", table2_name)).unwrap();
        let total_rows = select_stmt.query_map([], |_row| Ok(())).unwrap().count();
        let mut processed_rows = 0;
    
        let mut batch_smiles = Vec::with_capacity(BATCH_SIZE);
        let mut batch_rowids = Vec::with_capacity(BATCH_SIZE);
    
        for row in select_stmt.query_map([], |row| {
            let rowid: usize = row.get(0)?;
            let smiles: String = row.get(1)?;
            Ok((rowid, smiles))
        }).unwrap() {
            let (rowid, smiles) = row.unwrap();
            batch_smiles.push(smiles);
            batch_rowids.push(rowid);
    
            if batch_smiles.len() == BATCH_SIZE {
                let match_counts = single_functional_group(&mut batch_smiles, &smarts);
                
                let transaction = conn.transaction().unwrap();
    
                for (i, &count) in match_counts.iter().enumerate() {
                    transaction.execute(&format!("UPDATE {} SET '{}' = ?1 WHERE rowid = ?2", table_name, name), params![count, batch_rowids[i]]).unwrap();
                }
    
                transaction.commit().unwrap();
    
                processed_rows += BATCH_SIZE;
                let progress = (processed_rows as f32 / total_rows as f32) * 100.0;
                progress_sender.send(progress).unwrap();
    
                batch_smiles.clear();
                batch_rowids.clear();
            }
        }
    
        // Process any remaining items in the batch
        if !batch_smiles.is_empty() {
            let match_counts = single_functional_group(&mut batch_smiles, &smarts);
    
            let transaction = conn.transaction().unwrap();
            
            for (i, &count) in match_counts.iter().enumerate() {
                transaction.execute(&format!("UPDATE {} SET '{}' = ?1 WHERE rowid = ?2", table_name, name), params![count, batch_rowids[i]]).unwrap();
            }
    
            transaction.commit().unwrap();
    
            processed_rows += batch_smiles.len();
            let progress = (processed_rows as f32 / total_rows as f32) * 100.0;
            progress_sender.send(progress).unwrap();
        }
    }
    
    fn update_matrix_table_with_functional_group(table_name: &str, name: &str, matrices: &HashMap<String, bool>) {
        let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
        conn.execute(
            &format!("ALTER TABLE {} ADD COLUMN {} TEXT", table_name, name),
            [],
        ).unwrap();


        for (key, value) in matrices {
            conn.execute(
                &format!("INSERT OR REPLACE INTO {} (matrix, {}) VALUES (?1, ?2)", table_name, name),
                params![key, *value as i32],
            ).unwrap();
        }
    }

    pub fn add_fg_to_db(conn: r2d2::PooledConnection<SqliteConnectionManager>, name: String, smarts: String, matrices: HashMap<String, bool>, progress_sender: &mpsc::Sender<f32>) -> () {
        check_if_table_exists("functional_group_smarts", "user_functional_group_smarts").unwrap();
        check_if_table_exists("functional_groups",       "user_functional_groups").unwrap();
        check_if_table_exists("matrices",                "user_matrices").unwrap();
        //let conn2: r2d2::PooledConnection<SqliteConnectionManager> = POOL.get().unwrap();

        conn.execute(
            "INSERT INTO user_functional_group_smarts (name, smarts) VALUES (?1, ?2)",
            params![name, smarts],
        ).unwrap();

        //update functional_groups & user_functional_groups
        // Add the new column (as before)
        
        update_functional_groups("functional_groups", "metabolites", &name, &smarts, progress_sender);
        update_functional_groups("user_functional_groups", "user_metabolites", &name, &smarts, progress_sender);

        //update_fg(&smarts);
        
        //if any matrix is pressed, update derivatized_by and user_matrices
        update_matrix_table_with_functional_group("matrices", &name, &matrices);
        update_matrix_table_with_functional_group("user_matrices", &name, &matrices);


    }

    fn get_hashmap_from_table(conn: &r2d2::PooledConnection<SqliteConnectionManager>) -> HashMap<String, String> {
        let mut hashmap: HashMap<String, String> = HashMap::new();

        for table in &["functional_group_smarts", "user_functional_groups_smarts"] {
            match conn.prepare(&format!("SELECT name, smarts FROM {}", table)) {
                Ok(mut stmt) => {
                    let mut rows = stmt.query([]).unwrap();

                    while let Some(row) = rows.next().unwrap() {
                        let key: String = row.get(0).unwrap();
                        let value: String = row.get(1).unwrap();
                        hashmap.insert(key, value);
                    }
                }
                Err(_) => {
                    ()
                }
            }
        }
        hashmap
    }

    fn get_table_column_names(conn: &r2d2::PooledConnection<SqliteConnectionManager>, table_name: &str) -> Result<Vec<String>> {
        // Prepare a query, doesn't matter if we're not going to actually execute it.
        let stmt = conn.prepare(&format!("SELECT * FROM {table_name} LIMIT 0"))?;

        // Get the column count.
        let col_count = stmt.column_count();

        // Get the column names.
        let mut col_names = Vec::new();
        for i in 0..col_count {
            let col_name = stmt.column_name(i)?;
            col_names.push(col_name.to_string());
        }

        Ok(col_names)
    }

    fn get_matrices_fgs(conn: &r2d2::PooledConnection<SqliteConnectionManager>) -> Result<HashMap<String, Vec<String>>, rusqlite::Error> {
        let mut stmt = conn.prepare("SELECT * FROM matrices").unwrap();
        let col_count = stmt.column_count();
    
        // Fetch column names beforehand
        let mut col_names = Vec::new();
        for i in 0..col_count {
            let name = stmt.column_name(i)?.to_string();
            col_names.push(name);
        }
    
        let mut matrices_fgs: HashMap<String, Vec<String>> = HashMap::new();
    
        // Execute the query and map rows
        let mut rows = stmt.query([])?;
        
        while let Some(row) = rows.next()? {

            let mut fgs: Vec<String> = Vec::new();
            for i in 0..col_count {
                let value: String = row.get(i)?;
                if value == "1" {
                    fgs.push(col_names[i].clone())
                }
            }
            
            let name_: String = row.get(0)?;
            matrices_fgs.insert(name_, fgs);
        }
        matrices_fgs.insert("fmp".to_string(), vec!["phenols".to_string(), "primary amines".to_string()]);
        matrices_fgs.insert("ampp".to_string(), vec!["aldehydes".to_string(), "carboxylic acids".to_string()]);
    
        Ok(matrices_fgs)
    }
    
    fn fill_user_derivatized_by(conn: &r2d2::PooledConnection<SqliteConnectionManager>, fgh: &HashMap<String, String>, matrices_fgs: &HashMap<String, Vec<String>>) {

        let mut matrices = get_table_column_names(conn, "derivatized_by").unwrap();
        matrices.remove(0);

        let mut data: Vec::<(String, String)> = Vec::new();

        for matrix in matrices {
            let functional_groups = matrices_fgs.get(&matrix);
            let result = match functional_groups {
                Some(func_groups) => func_groups.iter().map(|fg| fgh.get(fg).unwrap().parse::<usize>().unwrap()).sum(),
                None => 0,
            };

            data.push((matrix, result.to_string()));
        }

        // Generate the column names and placeholder strings dynamically
        let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
        let placeholders: Vec<&str> = repeat("?").take(data.len()).collect();

        let sql = format!(
            "INSERT INTO user_derivatized_by ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();

        conn.execute(&sql, values.as_slice()).unwrap();

    }
    
    fn fill_user_metabolites(conn: &r2d2::PooledConnection<SqliteConnectionManager>, name: &String, smiles: &String, formula: &String, mz: &String) {
        let data = vec![
            ("name".to_string(), name), 
            ("smiles".to_string(), smiles),
            ("chemicalformula".to_string(), formula),
            ("mz".to_string(), mz)
        ];

        let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
        let placeholders: Vec<&str> = repeat("?").take(data.len()).collect();

        let sql = format!(
            "INSERT INTO user_metabolites ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();
        conn.execute(&sql, values.as_slice()).unwrap();
    }

    pub fn fill_user_msms(bin_data: Vec<u8>) {
        // Convert &[u8] to Vec<u8>
        let bin_data_vec = bin_data.to_vec();
    
        // Change data to Vec of tuples with String and Box<dyn ToSql>
        let data: Vec<(String, Box<dyn ToSql>)> = vec![
            ("name".to_string(), Box::new("User Input".to_string())),
            ("identifier".to_string(), Box::new("USER2".to_string())),
            ("adduct".to_string(), Box::new("1".to_string())),
            ("cid".to_string(), Box::new("0eV".to_string())),
            ("window".to_string(), Box::new("1Da".to_string())),
            ("tof".to_string(), Box::new("0.0".to_string())),
            ("mz".to_string(), Box::new("0.0".to_string())),
            ("spectra".to_string(), Box::new(bin_data_vec)),
        ];
    
        let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("connecting to MSMS db").unwrap();
    
        // Start a transaction
        conn.execute("BEGIN", []).unwrap();
    
        // Prepare an update statement
        let update_sql = format!(
            "UPDATE MSMS SET {} WHERE identifier = ?",
            data.iter().map(|(col, _)| format!("{} = ?", col)).collect::<Vec<_>>().join(", ")
        );
        let binding = "USER2".to_string();
        let update_values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| &**value).chain(std::iter::once(&binding as &dyn ToSql)).collect();
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
    
    pub fn add_to_usermsms(name: String, adduct: String, mz: String, cid: String, tof: String, mzwindow: String, identifier: String, path: String, matrix: String) {
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
            ("matrix".to_string(), Box::new(matrix))

        ];
    
        let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("connecting to MSMS db").unwrap();
    
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

    fn fill_user_endogeneity(conn: &r2d2::PooledConnection<SqliteConnectionManager>, endo_exo: &HashMap<String, bool>) {
        let column_names: Vec<String> = endo_exo.iter().map(|(col, _)| col.to_string()).collect();
        let placeholders: Vec<&str> = repeat("?").take(endo_exo.len()).collect();

        let sql = format!(
            "INSERT INTO user_endogeneity ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        let values: Vec<&dyn ToSql> = endo_exo.iter().map(|(_, value)| value as &dyn ToSql).collect();

        conn.execute(&sql, values.as_slice()).unwrap();
    }

    fn fill_user_in_tissue(conn: &r2d2::PooledConnection<SqliteConnectionManager>, in_tissue: &HashMap<String, bool>) {
        let column_names: Vec<String> = in_tissue.iter().map(|(col, _)| col.to_string()).collect();
        let placeholders: Vec<&str> = repeat("?").take(in_tissue.len()).collect();

        let sql = format!(
            "INSERT INTO user_in_tissue ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        let values: Vec<&dyn ToSql> = in_tissue.iter().map(|(_, value)| value as &dyn ToSql).collect();

        conn.execute(&sql, values.as_slice()).unwrap();
    }

    fn fill_user_functional_groups(conn: &r2d2::PooledConnection<SqliteConnectionManager>, fgh: &HashMap<String, String>) {
        let mut functional_groups: Vec<String> = get_table_column_names(conn, "functional_groups").unwrap();
        functional_groups.remove(0);

        let new_fgh: HashMap<String, String> = fgh
            .iter()
            .map(|(k, v)| (k.replace(" ", ""), v.clone()))
            .collect();


        let mut data: HashMap<String, &String> = HashMap::new();
        for fg in functional_groups {
            if let Some(value) = new_fgh.get(&fg) {
                data.insert(fg.clone(), value);
            }
        }



        let column_names: Vec<String> = data.iter().map(|(col, _)| col.to_string()).collect();
        let placeholders: Vec<&str> = repeat("?").take(data.len()).collect();

        let sql = format!(
            "INSERT INTO user_functional_groups ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        let values: Vec<&dyn ToSql> = data.iter().map(|(_, value)| value as &dyn ToSql).collect();

        conn.execute(&sql, values.as_slice()).unwrap();
    }

    fn fill_user_db_accessions(conn: &r2d2::PooledConnection<SqliteConnectionManager>) {
        conn.execute(
            "INSERT INTO user_db_accessions (hmdb) VALUES (?1)",
            params!["placeholder"],  // using 0 as a placeholder
        ).unwrap();
    
        // Get the last inserted row's ID
        let last_id = conn.last_insert_rowid();

        // Create the special_id based on the last inserted row's ID
        let special_id = format!("METID{:05}", last_id);

        // Update the newly inserted row to set the special_id
        conn.execute(
            "UPDATE user_db_accessions SET hmdb = ?1 WHERE id = ?2",
            params![special_id, last_id],
        ).unwrap();
    }

    pub fn add_metabolite_to_db(conn: r2d2::PooledConnection<SqliteConnectionManager>, name: String, smiles: String, in_tissue: HashMap<String, bool>, endo_exo: HashMap<String, bool>) -> () {        
        check_if_table_exists("metabolites",       "user_metabolites").unwrap();
        check_if_table_exists("derivatized_by",    "user_derivatized_by").unwrap();
        check_if_table_exists("endogeneity",       "user_endogeneity").unwrap();
        check_if_table_exists("in_tissue",         "user_in_tissue").unwrap();
        check_if_table_exists("db_accessions",     "user_db_accessions").unwrap();
        check_if_table_exists("functional_groups", "user_functional_groups").unwrap();

        let mut functional_smarts: HashMap<String, String> = get_hashmap_from_table(&conn);

        let (formula, mz, fgh) = metabolite_for_db_sidecar(&smiles, &mut functional_smarts);
        //let metabolite: Metabolite = Metabolite{ name: name.clone(), smiles: smiles.clone() };
        //let fgh: HashMap<String, String> = metabolite.functional_group(&functional_smarts).unwrap();
        let matrices_fgs: HashMap<String, Vec<String>> = get_matrices_fgs(&conn).unwrap();

        fill_user_derivatized_by(&conn, &fgh, &matrices_fgs);
        fill_user_metabolites(&conn, &name, &smiles, &formula, &mz);
        fill_user_endogeneity(&conn, &endo_exo);
        fill_user_in_tissue(&conn, &in_tissue);
        fill_user_functional_groups(&conn, &fgh);
        fill_user_db_accessions(&conn);
    }

    fn split_key(key: &str) -> Option<(String, i32)> {
        // Reverse the string and find the first non-numeric character
        let reversed: String = key.chars().rev().collect();
        let idx = reversed.chars().position(|ch| !ch.is_numeric())?;
    
        // Get the prefix and number using the identified index
        let number: i32 = reversed[0..idx].chars().rev().collect::<String>().parse().ok()?;
        let prefix: String = reversed[idx..].chars().rev().collect();
    
        Some((prefix, number))
    }

    pub fn add_matrix_to_db(conn: r2d2::PooledConnection<SqliteConnectionManager>, name: String, _smiles_mz: String, checkboxes: HashMap<String, bool>, adducts: HashMap<String, String>) -> () {
        check_if_table_exists("matrices", "user_matrices").unwrap();
        check_if_table_exists("adducts", "user_adducts").unwrap();    
        
        //filling user_adducts
        let mut adduct_vec: Vec<String> = Vec::new();
        let mut mx_formula_vec: Vec<String> = Vec::new();
        let mut fg_vec: Vec<String> = Vec::new();
        let mut deltamass_vec: Vec<String> = Vec::new();

        for (key, value) in adducts.iter() {
            if let Some((prefix, number)) = split_key(key) {
                let _column = match prefix.as_str() {
                    "add-adduct-name-" => adduct_vec.insert(number as usize, value.to_string()),
                    "add-matrix-formula-" => {
                        mx_formula_vec.insert(number as usize, value.to_string()); 
                        deltamass_vec.insert(number as usize, "0".to_string());//mass_from_formula(value).to_string());
                    },
                    "add-fg-" => {
                        if value.to_string() == "".to_string() {
                            fg_vec.insert(number as usize, "1".to_string());
                        } else {
                            fg_vec.insert(number as usize, value.to_string());
                        }
                    },
                    _ => continue, // skip if not one of the expected prefixes
                };
            }
        }

        // Iterate over the vectors and insert the data into the table
        for i in 0..adduct_vec.len() {
            conn.execute(
                "INSERT INTO user_adducts (adduct, mname, numfunctionalgroups, formula, deltamass) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![adduct_vec[i], name, fg_vec[i], mx_formula_vec[i], deltamass_vec[i]],
            ).unwrap();
        }

        // filling user_matrices
        let mut column_names: Vec<String> = checkboxes.iter()
            .map(|(col, _)| format!("'{}'", col))
            .collect();
        let mut placeholders: Vec<&str> = repeat("?").take(checkboxes.len()).collect();

        column_names.push("matrix".to_string());
        placeholders.push("?");


        let sql = format!(
            "INSERT INTO user_matrices ({}) VALUES ({})",
            column_names.join(", "),
            placeholders.join(", ")
        );

        let mut values: Vec<&dyn ToSql> = checkboxes.iter().map(|(_, value)| value as &dyn ToSql).collect();
        values.push(&name as &dyn ToSql);

        conn.execute(&sql, values.as_slice()).unwrap();

    }

    pub fn add_to_db_rust(name: String, smiles_smarts_mz: String, met_type: String, endo_exo_or_other: HashMap<String, bool>, in_tissue: HashMap<String, bool>, adducts: HashMap<String, String>, progress_sender: &mpsc::Sender<f32>) -> bool {
        let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

        match &met_type[..] {
        "metabolite" => add_metabolite_to_db(conn, name, smiles_smarts_mz, in_tissue, endo_exo_or_other),
        "matrix" => add_matrix_to_db(conn, name, smiles_smarts_mz, endo_exo_or_other, adducts),
        "fg" => add_fg_to_db(conn, name, smiles_smarts_mz, endo_exo_or_other, progress_sender),
        _ => ()
        }

        true
    }
}
