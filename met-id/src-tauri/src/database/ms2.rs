

/* 
use std::collections::HashMap;
use r2d2_sqlite::SqliteConnectionManager;



////////////////////////////////////////////TAURI FUNCTIONS ///////////////////////////////////////////////////////////////////
/* 
#[tauri::command]
pub fn get_msms_tauri() -> Vec<Vec<String>> {
    //MS2DbRow::get_msms()
}
*/




/////////////////////////////////////// CODE ///////////////////////////////////////////////////////////////////////////////////

#[derive(Debug)]
struct MS2DbRow {
    identifier: String,
    adduct: String,
    cid: String,
    window: String,
    tof: String,
    mz: String,
    data: (Vec<f64>, Vec<i64>),
}
/* 
impl MS2DbRow {

    fn get_names_from_identifiers(input_data: &Vec<String>) -> rusqlite::Result<HashMap<String, String>> {
        let conn: r2d2::PooledConnection<SqliteConnectionManager> = super::POOL.get().unwrap();
    
        let mut name_id_conversion: HashMap<String, String> = HashMap::new();
    
        let params: String = input_data.iter().map(|s| format!("'{}'", s)).collect::<Vec<_>>().join(", ");
    
        let sql: String = format!("SELECT metabolites.name, db_accessions.hmdb FROM metabolites JOIN db_accessions ON metabolites.id = db_accessions.id WHERE db_accessions.hmdb IN ({})", params);
        // Execute the query
        let mut stmt = conn.prepare(&sql)?;
        let rows = stmt.query_map([], |row| {
            Ok(super::IDentifierMap {
                name: row.get(0).unwrap_or("".to_string()), 
                hmdb: row.get(1).unwrap_or("".to_string()), 
            })
        }).unwrap();
    
        // Iterate through the rows and do something with the results
        for item in rows {
            let row: super::IDentifierMap = item.unwrap();
            name_id_conversion.insert(row.hmdb, row.name);
        }
        name_id_conversion.insert("".to_string(), "FMP10".to_string());
    
        Ok(name_id_conversion)
    }







    pub fn get_msms() -> Vec<Vec<String>> {
        let conn = super::MSMS_POOL.get().expect("Couldn't get the connection from the pool");
        let mut stmt = conn.prepare("SELECT identifier, adduct, cid, window, tof, mz, spectra FROM MSMS")
            .expect("Query cannot be run");
    
        let rows: Vec<MS2DbRow> = stmt.query_map([], |row| {
            Ok(MS2DbRow {
                identifier: row.get(0).unwrap_or_default(),
                adduct: row.get(1).unwrap_or_else(|_| row.get(5).unwrap_or_default()),
                cid: row.get(2).unwrap_or_default(),
                window: row.get(3).unwrap_or_default(),
                tof: row.get(4).unwrap_or_default(),
                mz: row.get(5).unwrap_or_default(),
                data: (vec![], vec![]),
            })
        }).unwrap().filter_map(Result::ok).collect();
    
        let identifiers: Vec<String> = rows.iter().map(|row| row.identifier.clone()).collect();
        let mut unique_identifiers = identifiers.clone();
        unique_identifiers.sort();
        unique_identifiers.dedup();
    
        let conversions = Self::get_names_from_identifiers(&unique_identifiers).unwrap();
    
        let names = identifiers.iter()
            .map(|id| conversions.get(id).unwrap_or(&String::new()).clone())
            .collect();
    
        vec![names, rows.iter().map(|row| row.identifier.clone()).collect(), rows.iter().map(|row| row.adduct.clone()).collect(), rows.iter().map(|row| row.cid.clone()).collect(), 
        rows.iter().map(|row| row.window.clone()).collect(), rows.iter().map(|row| row.tof.clone()).collect(), rows.iter().map(|row| row.mz.clone()).collect(),
        ]
    }




    #[tauri::command]
    pub fn ms2_search_spectra(name: String, fragment: String, ms1mass: String, fragmentslider: String, _ms1massslider: String) -> Vec<Vec<String>> {
        //let conn: r2d2::PooledConnection<SqliteConnectionManager> = MSMS_POOL.get().unwrap();
        let table: &str = "MSMS";
        let name_condition: &str = &name.as_str();
        let sort_field: &str = "name";
        
        let mut query = format!("SELECT identifier, adduct, cid, window, tof, mz, spectra FROM {}", table);
        
        if name_condition != "" {
            query.push_str(" WHERE name LIKE '%");
            query.push_str(name_condition);
            query.push_str("%' ")
        }
        

        if let Ok(ms1) = ms1mass.parse::<f64>() {
            if query.contains(" WHERE ") {
                query.push_str(" AND CAST(mz AS REAL) < ");
            } else {
                query.push_str(" WHERE CAST(mz AS REAL) < ");
            }
            query.push_str(&(ms1 + 0.5).to_string());
            query.push_str(" AND CAST(mz AS REAL) > ");
            query.push_str(&(ms1 - 0.5).to_string());
        }

        if !sort_field.is_empty() {
            query.push_str(" ORDER BY ");
            query.push_str(sort_field);
        }

        println!("{}", &query);

        


        if let Ok(mut stmt) = MSMS_POOL.get().unwrap().prepare(&query) {
            let mut identifiers: Vec<String> = Vec::new();
            let mut adducts: Vec<String> = Vec::new();
            let mut cids: Vec<String> = Vec::new();
            let mut windows: Vec<String> = Vec::new();
            let mut tofs: Vec<String> = Vec::new();
            let mut mzs: Vec<String> = Vec::new();
            let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();



            let db_iter: Box<dyn Iterator<Item = rusqlite::Result<Option<MS2DbRow>>>> = if let Ok(fragment_mz) = fragment.parse::<f64>() {
                let min_mz = fragment_mz - (fragmentslider.parse::<f64>().unwrap_or(0.0) /2000.0);
                let max_mz = fragment_mz + (fragmentslider.parse::<f64>().unwrap_or(0.0) /2000.0);
                println!("{:?}, {:?}, {:?}, {:?}", fragment_mz, min_mz, max_mz, fragmentslider);
                Box::new(stmt.query_map([], move |row: &rusqlite::Row<'_>| {
                    let spectrum = collect_msms_spectra(BinaryData {binary: row.get(6).unwrap_or(Vec::new())});
                    if contains_between_a_and_b(&spectrum.0, min_mz, max_mz) {
                        let adduct_value = row.get(1).unwrap_or("".to_string());
                        let adduct_final_value = if adduct_value.is_empty() {
                            row.get(5).unwrap_or("".to_string())
                        } else {
                            adduct_value
                        };
                        Ok(Some(MS2DbRow {
                            identifier: row.get(0).unwrap_or("".to_string()),
                            adduct: adduct_final_value,
                            cid: row.get(2).unwrap_or("".to_string()),
                            window: row.get(3).unwrap_or("".to_string()),
                            tof: row.get(4).unwrap_or("".to_string()),
                            mz: row.get(5).unwrap_or("".to_string()),
                            data: spectrum,

                        }))
                    } else {
                        Ok(None)
                    }
                }).unwrap())
            } else {
                Box::new(stmt.query_map([], |row: &rusqlite::Row<'_>| {
                    let adduct_value = row.get(1).unwrap_or("".to_string());
                    let adduct_final_value = if adduct_value.is_empty() {
                        row.get(5).unwrap_or("".to_string())
                    } else {
                        adduct_value
                    };
                    Ok(Some(MS2DbRow {
                        identifier: row.get(0).unwrap_or("".to_string()),
                        adduct: adduct_final_value,
                        cid: row.get(2).unwrap_or("".to_string()),
                        window: row.get(3).unwrap_or("".to_string()),
                        tof: row.get(4).unwrap_or("".to_string()),
                        mz: row.get(5).unwrap_or("".to_string()),
                        data: (Vec::new(), Vec::new()),
                    }))
                }).unwrap())
            };

            let filtered_iter = db_iter.filter_map(|result| {
                match result {
                    Ok(Some(db_row)) => Some(Ok(db_row)),
                    Ok(None) => None,
                    Err(e) => Some(Err(e)),
                }
            });
            //let filtered_iter = db_iter.filter_map(Result::ok);
            //parse results for passing back to the parent function
            for (index, item) in filtered_iter.enumerate() {

                let row: MS2DbRow = item.unwrap();
                

                identifiers.insert(index, row.identifier);
                adducts.insert(index, row.adduct);
                cids.insert(index, row.cid);
                windows.insert(index, row.window);
                tofs.insert(index, row.tof);
                mzs.insert(index, row.mz);
                datas.insert(index, row.data);//collect_msms_spectra(BinaryData { binary: row.data }));
                
            }


            let mut identifiers2: Vec<String> = identifiers.clone();
            identifiers2.sort();
            identifiers2.dedup();
            let conversions: HashMap<String, String> = get_names_from_identifiers(&identifiers2).unwrap();

            let mut names: Vec<String> = Vec::new();
            for (index, id) in identifiers.iter().enumerate() {
                let name = conversions.get(id).unwrap_or(&"".to_string()).to_owned();
                names.insert(index, name);
            }
            println!("{:?}", names);
            let return_data: Vec<Vec<String>> = vec![names, identifiers, adducts, cids, windows, tofs, mzs];
            return return_data
        } else {
            return Vec::new()
        }
}
}

*/

*/