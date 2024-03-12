use std::cmp::Ordering;
use std::collections::HashMap;
use std::time::Instant;

use r2d2_sqlite::SqliteConnectionManager;
use crate::add_to_db::add_to_db_functions::add_to_usermsms;
use crate::database_setup::{get_msms_connection, get_connection};
use super::table::check_if_table_exists_msms;
use crate::ms2_match::ms2_matcher;
use super::{MS2DbRow, MS2DbRow2};
use super::IDentifierMap;
use super::SpectrumPoint;
use std::mem;
use rusqlite::params;


struct MSMSBinary {
  binary: Vec<u8>,
  cid: String
}

struct BinaryData {
  binary: Vec<u8>
}


pub fn get_msms() -> Vec<Vec<String>> {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_msms").unwrap();
    check_if_table_exists_msms("MSMS", "user_MSMS").unwrap();
    let mut stmt: rusqlite::Statement = conn.prepare("SELECT identifier, adduct, cid, window, tof, mz, spectra, matrix FROM MSMS UNION ALL SELECT identifier, adduct, cid, window, tof, mz, spectra, matrix FROM user_MSMS").expect("Query cannot be run");
    let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
        let adduct_value = row.get(1).unwrap_or("".to_string());
        let adduct_final_value = if adduct_value.is_empty() {
            row.get(5).unwrap_or("".to_string())
        } else {
            adduct_value
        };
        Ok(MS2DbRow {
            name: "".to_string(),
            identifier: row.get(0).unwrap_or("".to_string()),
            adduct: adduct_final_value,
            cid: row.get(2).unwrap_or("".to_string()),
            window: row.get(3).unwrap_or("".to_string()),
            tof: row.get(4).unwrap_or("".to_string()),
            mz: row.get(5).unwrap_or("".to_string()),
            data: (Vec::new(), Vec::new()),
            matrix: row.get(7).unwrap_or("".to_string())
            //data: row.get(6).unwrap_or(Vec::new()),
        })
    }).unwrap();

    let mut identifiers: Vec<String> = Vec::new();
    let mut adducts: Vec<String> = Vec::new();
    let mut cids: Vec<String> = Vec::new();
    let mut windows: Vec<String> = Vec::new();
    let mut tofs: Vec<String> = Vec::new();
    let mut mzs: Vec<String> = Vec::new();
    let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();
    let mut matrices: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        let row: MS2DbRow = item.unwrap();
        identifiers.insert(index, row.identifier);
        adducts.insert(index, row.adduct);
        cids.insert(index, row.cid);
        windows.insert(index, row.window);
        tofs.insert(index, row.tof);
        mzs.insert(index, row.mz);
        datas.insert(index, row.data);
        matrices.insert(index, row.matrix);
    }
    let mut identifiers2: Vec<String> = identifiers.clone();
    identifiers2.sort();
    identifiers2.dedup();
    let conversions: HashMap<String, String> = get_names_from_identifiers(&identifiers2).unwrap();

    let mut names: Vec<String> = Vec::new();
    for (index, id) in identifiers.iter().enumerate() {
        let name = conversions.get(id).unwrap_or(&"U".to_string()).to_owned();
        names.insert(index, name);
    }
    let return_data: Vec<Vec<String>> = vec![names, identifiers, adducts, cids, mzs, windows, tofs, matrices];
    return_data
}

fn get_names_from_identifiers(input_data: &Vec<String>) -> rusqlite::Result<HashMap<String, String>> {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

  let mut name_id_conversion: HashMap<String, String> = HashMap::new();

  let params: String = input_data.iter().map(|s| format!("'{}'", s)).collect::<Vec<_>>().join(", ");

  let sql: String = format!("SELECT metabolites.name, db_accessions.hmdb FROM metabolites JOIN db_accessions ON metabolites.id = db_accessions.id WHERE db_accessions.hmdb IN ({})", params);
  // Execute the query
  let mut stmt = conn.prepare(&sql)?;
  let rows = stmt.query_map([], |row| {
      Ok(IDentifierMap {
          name: row.get(0).unwrap_or("".to_string()), 
          hmdb: row.get(1).unwrap_or("".to_string()), 
      })
  }).unwrap();

  // Iterate through the rows and do something with the results
  for item in rows {
      let row: IDentifierMap = item.unwrap();
      name_id_conversion.insert(row.hmdb, row.name);
  }
  name_id_conversion.insert("HMDBM000001".to_string(), "FMP10".to_string());
  name_id_conversion.insert("USER2".to_string(), "User Input".to_string());

  Ok(name_id_conversion)
}

pub fn get_msms_spectra(identifier: String, adduct: String, cid: String) -> String {
  let spectra: (Vec<(Vec<f64>, Vec<i64>)>, Vec<String>) = get_single_msms_spectra(identifier, adduct, cid);

    let mut data: HashMap<String, Vec<SpectrumPoint>> = HashMap::new();
    for i in 0..spectra.0.len() {
        let spectrum: &(Vec<f64>, Vec<i64>) = &spectra.0[i];
        let cid: &String = &spectra.1[i];
        let mut dat: Vec<SpectrumPoint> = Vec::new();
        for idx in 0..spectrum.0.len() {
            dat.push(SpectrumPoint {x: spectrum.0[idx], y: spectrum.1[idx]});
        }
        data.insert(cid.to_string(), dat);
    }
    serde_json::to_string(&data).expect("Failed to serialize data")
}

fn get_single_msms_spectra(identifier: String, adduct: String, cid:String) -> (Vec<(Vec<f64>, Vec<i64>)>, Vec<String>) {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_single_msms_spectra").unwrap();

  let sql: String = format!("SELECT spectra, cid FROM MSMS WHERE identifier = '{}' AND adduct = '{}' AND cid = '{}'  UNION SELECT spectra, cid FROM user_MSMS WHERE identifier = '{}' AND adduct = '{}' AND cid = '{}'", identifier, adduct, cid, identifier, adduct, cid);

  let mut stmt = conn.prepare(&sql).unwrap();
  let rows = stmt.query_map([], |row| {
      Ok(MSMSBinary {
          binary: row.get(0).unwrap_or(Vec::new()), 
          cid: row.get(1).unwrap_or(String::new()), 
      })
  }).unwrap();
  

  let mut spectras: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();
  let mut cids: Vec<String> = Vec::new();

  for (index, item) in rows.enumerate() {
      let row: MSMSBinary = item.unwrap();
      let bin: BinaryData = BinaryData{ binary: row.binary};
      let data: (Vec<f64>, Vec<i64>) = collect_msms_spectra(bin);
      spectras.insert(index, data);
      cids.insert(index, row.cid);

  }


  (spectras, cids)
}

fn collect_msms_spectra(binary_data: BinaryData) -> (Vec<f64>, Vec<i64>) {
  // Ensure that the binary data is properly structured
  assert!(binary_data.binary.len() % (mem::size_of::<f64>() + mem::size_of::<i64>()) == 0);

  let mut f64_values = Vec::new();
  let mut i64_values = Vec::new();

  // Iterate through the binary data, taking slices of the appropriate sizes
  for i in (0..binary_data.binary.len()).step_by(mem::size_of::<f64>() + mem::size_of::<i64>()) {
      let f64_bytes = &binary_data.binary[i..i + mem::size_of::<f64>()];
      let i64_bytes = &binary_data.binary[i + mem::size_of::<f64>()..i + mem::size_of::<f64>() + mem::size_of::<i64>()];

      // Use the byteorder crate to convert the slices to numbers
      let f64_value = f64::from_le_bytes(f64_bytes.try_into().unwrap());
      let i64_value = i64::from_le_bytes(i64_bytes.try_into().unwrap());

      f64_values.push(f64_value);
      i64_values.push(i64_value);
  }

  let result = (f64_values, i64_values);
  result
}

pub fn get_name_from_identifier_msms(identifier: String) -> String {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

  let sql: String = format!("SELECT metabolites.name, db_accessions.hmdb FROM metabolites JOIN db_accessions ON metabolites.id = db_accessions.id WHERE db_accessions.hmdb IN (?)");
  let name: Result<String, rusqlite::Error> = conn.query_row(&sql, &[&identifier], |row| {
      row.get(0)
  });

  match name {
      Ok(name) => return name,
      Err(e) => {
          // Handle error or return None if no row found
          if e == rusqlite::Error::QueryReturnedNoRows {
              return "".to_string();
          } else {
              panic!("Database error: {:?}", e);
          }
      }
  }
}


pub fn ms2_search_spectra(name: String, fragment: String, ms1mass: String, fragmentslider: String, _ms1massslider: String) -> Vec<Vec<String>> {

    let name_condition: &str = &name.as_str();
    let sort_field: &str = "name";
    
    let mut query: String = "SELECT identifier, adduct, cid, window, tof, mz, spectra, name, matrix FROM MSMS ".to_string();
    
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

    query.push_str(" UNION ALL SELECT identifier, adduct, cid, window, tof, mz, spectra, name, matrix FROM user_MSMS");

    if let Ok(ms1) = ms1mass.parse::<f64>() {

        query.push_str(" WHERE CAST(mz AS REAL) < ");
        
        query.push_str(&(ms1 + 0.5).to_string());
        query.push_str(" AND CAST(mz AS REAL) > ");
        query.push_str(&(ms1 - 0.5).to_string());
    } 

    if !sort_field.is_empty() {
        query.push_str(" ORDER BY ");
        query.push_str(sort_field);
    }

    if let Ok(mut stmt) = get_msms_connection("ms2_search_spectra").unwrap().prepare(&query) {
        let mut identifiers: Vec<String> = Vec::new();
        let mut adducts: Vec<String> = Vec::new();
        let mut cids: Vec<String> = Vec::new();
        let mut windows: Vec<String> = Vec::new();
        let mut tofs: Vec<String> = Vec::new();
        let mut mzs: Vec<String> = Vec::new();
        let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();
        let mut names: Vec<String> = Vec::new();
        let mut matrices: Vec<String> = Vec::new();

        

        let db_iter: Box<dyn Iterator<Item = rusqlite::Result<Option<MS2DbRow>>>> = if let Ok(fragment_mz) = fragment.parse::<f64>() {
            let min_mz = fragment_mz - (fragmentslider.parse::<f64>().unwrap_or(0.0) /2000.0);
            let max_mz = fragment_mz + (fragmentslider.parse::<f64>().unwrap_or(0.0) /2000.0);
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
                        name: row.get(7).unwrap_or("".to_string()),
                        identifier: row.get(0).unwrap_or("".to_string()),
                        adduct:  adduct_final_value,
                        cid: row.get(2).unwrap_or("".to_string()),
                        window: row.get(3).unwrap_or("".to_string()),
                        tof: row.get(4).unwrap_or("".to_string()),
                        mz: row.get(5).unwrap_or("".to_string()),
                        data: spectrum,
                        matrix: row.get(8).unwrap_or("".to_string())

                    }))
                } else {
                    Ok(None)
                }
            }).unwrap())
        } else if let Ok(fragment_mz) = ms1mass.parse::<f64>() {
            let min_mz = fragment_mz - (_ms1massslider.parse::<f64>().unwrap_or(0.0) /2000.0);
            let max_mz = fragment_mz + (_ms1massslider.parse::<f64>().unwrap_or(0.0) /2000.0);
            Box::new(stmt.query_map([], move |row: &rusqlite::Row<'_>| {
                let spectrum = collect_msms_spectra(BinaryData {binary: row.get(6).unwrap_or(Vec::new())});
                let mz: f64 = row.get(5).unwrap_or("".to_string()).parse::<f64>().unwrap_or(0.0);
                
                if mz < max_mz && mz > min_mz {
                    Ok(Some(MS2DbRow {
                        name: row.get(7).unwrap_or("".to_string()),
                        identifier: row.get(0).unwrap_or("".to_string()),
                        adduct: row.get(1).unwrap_or("".to_string()),
                        cid: row.get(2).unwrap_or("".to_string()),
                        window: row.get(3).unwrap_or("".to_string()),
                        tof: row.get(4).unwrap_or("".to_string()),
                        mz: row.get(5).unwrap_or("".to_string()),
                        data: spectrum,
                        matrix: row.get(8).unwrap_or("".to_string())

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
                    name: row.get(7).unwrap_or("".to_string()),
                    identifier: row.get(0).unwrap_or("".to_string()),
                    adduct: adduct_final_value,
                    cid: row.get(2).unwrap_or("".to_string()),
                    window: row.get(3).unwrap_or("".to_string()),
                    tof: row.get(4).unwrap_or("".to_string()),
                    mz: row.get(5).unwrap_or("".to_string()),
                    data: (Vec::new(), Vec::new()),
                    matrix: row.get(8).unwrap_or("".to_string())
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
            names.insert(index, row.name);
            identifiers.insert(index, row.identifier);
            adducts.insert(index, row.adduct);
            cids.insert(index, row.cid);
            windows.insert(index, row.window);
            tofs.insert(index, row.tof);
            mzs.insert(index, row.mz);
            datas.insert(index, row.data);
            matrices.insert(index, row.matrix);
            
        }

        let vector_len = names.len();
        move_to_front(&mut names, &mut identifiers, &mut adducts, &mut cids, &mut windows, &mut tofs, &mut mzs, &mut matrices, "User Input".to_string(), vector_len);

        let return_data: Vec<Vec<String>> = vec![names, identifiers, adducts, cids, windows, tofs, mzs, Vec::new(), matrices];
        return return_data
    } else {
        return Vec::new()
    }
}

#[allow(unused_assignments)]
fn contains_between_a_and_b(vec: &[f64], min_mz: f64, max_mz: f64) -> bool {
    let mut lower_bound: usize = 0;
    let mut upper_bound: usize = 0;

    // Find the lower bound
    match vec.binary_search_by(|&x| x.partial_cmp(&min_mz).unwrap_or(Ordering::Greater)) {
        Ok(index) => lower_bound = index + 1,  // Skip exact matches
        Err(index) => lower_bound = index,      // Err gives us the point of insertion.
    }

    // Find the upper bound
    match vec.binary_search_by(|&x| x.partial_cmp(&max_mz).unwrap_or(Ordering::Less)) {
        Ok(index) => upper_bound = index,       // Skip exact matches
        Err(index) => upper_bound = index,      // Err gives us the point of insertion.
    }

    lower_bound < upper_bound
}

fn move_to_front<T: PartialEq>(names: &mut Vec<T>, identifiers: &mut Vec<T>, adducts: &mut Vec<T>, cids: &mut Vec<T>, windows: &mut Vec<T>, tofs: &mut Vec<T>,  mzs: &mut Vec<T>, matrices: &mut Vec<T>, item: T, vector_len: usize) {
  if let Some(index) = names.iter().position(|x| *x == item) {
    names.rotate_right(vector_len - index);
    identifiers.rotate_right(vector_len - index);
    adducts.rotate_right(vector_len - index);
    cids.rotate_right(vector_len - index);
    windows.rotate_right(vector_len - index);
    tofs.rotate_right(vector_len- index);
    mzs.rotate_right(vector_len - index);
    matrices.rotate_right(vector_len - index);

  }
}


pub fn match_msms_to_ui(binsize: f64, threshold: f64) -> (Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<f64>, Vec<String>){
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_msms").unwrap();
  let mut stmt: rusqlite::Statement = conn.prepare("SELECT name, identifier, adduct, cid, window, tof, mz, spectra, matrix FROM MSMS UNION SELECT name, identifier, adduct, cid, window, tof, mz, spectra, matrix FROM user_MSMS").expect("Query cannot be run");
  let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
      let adduct_value = row.get(2).unwrap_or("".to_string());
      let spectrum = collect_msms_spectra(BinaryData {binary: row.get(7).unwrap_or(Vec::new())});
      let adduct_final_value = if adduct_value.is_empty() {
          row.get(6).unwrap_or("".to_string())
      } else {
          adduct_value
      };
      Ok(MS2DbRow {
          name: row.get(0).unwrap_or("".to_string()),
          identifier: row.get(1).unwrap_or("".to_string()),
          adduct: adduct_final_value,
          cid: row.get(3).unwrap_or("".to_string()),
          window: row.get(4).unwrap_or("".to_string()),
          tof: row.get(5).unwrap_or("".to_string()),
          mz: row.get(6).unwrap_or("".to_string()),
          //data: (Vec::new(), Vec::new()),
          data: spectrum,
          //data: row.get(6).unwrap_or((Vec::new(), Vec::new())),
          matrix: row.get(8).unwrap_or("".to_string())
      })
  }).unwrap();
  let mut names: Vec<String> = Vec::new();
  let mut identifiers: Vec<String> = Vec::new();
  let mut adducts: Vec<String> = Vec::new();
  let mut cids: Vec<String> = Vec::new();
  let mut windows: Vec<String> = Vec::new();
  let mut tofs: Vec<String> = Vec::new();
  let mut mzs: Vec<String> = Vec::new();
  let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();
  let mut matrices: Vec<String> = Vec::new();

  //parse results for passing back to the parent function
  for (index, item) in db_iter.enumerate() {
      
      let row: MS2DbRow = item.unwrap();
      names.insert(index, row.name);
      identifiers.insert(index, row.identifier);
      adducts.insert(index, row.adduct);
      cids.insert(index, row.cid);
      windows.insert(index, row.window);
      tofs.insert(index, row.tof);
      mzs.insert(index, row.mz);
      datas.insert(index, row.data);
      matrices.insert(index, row.matrix);
  }


  let mut index_to_remove = None;

  // Find the index of the item "foo"
  for (index, item) in identifiers.iter().enumerate() {
      if *item == "USER2" {
          index_to_remove = Some(index);
          break;
      }
  }
  let mut cossim: Vec<f64>;
  match index_to_remove {
      Some(index) => {
          let user_spectrum: (Vec<f64>, Vec<i64>) = datas[index].clone();
          let start = Instant::now();
          cossim = ms2_matcher(&user_spectrum, &datas.clone(), binsize, threshold);
          let duration = start.elapsed();
          println!("Time elapsed in ms2 matcher is: {:?}", duration);
      },
      None => {
          cossim = Vec::new();
      },
      
  } 
  sort_by_corresponding_f64(&mut cossim, &mut identifiers, &mut adducts, &mut cids, &mut names, &mut mzs, &mut matrices);


  (names, identifiers, adducts, cids, mzs, cossim, matrices)

}

fn sort_by_corresponding_f64(
  vec_f64: &mut Vec<f64>,
  vec_str1: &mut Vec<String>,
  vec_str2: &mut Vec<String>,
  vec_str3: &mut Vec<String>,
  vec_str4: &mut Vec<String>,
  vec_str5: &mut Vec<String>,
  vec_str6: &mut Vec<String>
  ) {
  let mut combined: Vec<(_, _, _, _, _, _, _)> = vec_f64
      .iter()
      .cloned()
      .zip(vec_str1.iter().cloned())
      .zip(vec_str2.iter().cloned())
      .zip(vec_str3.iter().cloned())
      .zip(vec_str4.iter().cloned())
      .zip(vec_str5.iter().cloned())
      .zip(vec_str6.iter().cloned())
      .map(|((((((a, b), c), d), e), f), g)| (a, b, c, d, e, f, g))
      .collect();

  combined.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

  for (i, (num, text1, text2, text3, text4, text5, text6)) in combined.into_iter().enumerate() {
      vec_f64[i] = num;
      vec_str1[i] = text1;
      vec_str2[i] = text2;
      vec_str3[i] = text3;
      vec_str4[i] = text4;
      vec_str5[i] = text5;
      vec_str6[i] = text6;
  }
}


pub fn add_msms_to_db(name: String, adduct: String, mz: String, cid: String, tof: String, mzwindow: String, identifier: String, path: String, matrix: String) -> String {
  add_to_usermsms(name, adduct, mz, cid, tof, mzwindow, identifier, path, matrix);
  
  String::from("Done")
}


pub fn show_user_msms_db() -> Vec<Vec<String>> {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_msms").unwrap();
    let mut stmt: rusqlite::Statement = conn.prepare("SELECT id, name, identifier, adduct, cid, window, tof, mz, spectra, matrix FROM user_MSMS").expect("Query cannot be run");
    let db_iter = stmt.query_map([], |row: &rusqlite::Row<'_>| {
        let adduct_value = row.get(3).unwrap_or("".to_string());
        let adduct_final_value = if adduct_value.is_empty() {
            row.get(5).unwrap_or("".to_string())
        } else {
            adduct_value
        };
        let id_value = row.get(0).unwrap();
        Ok(MS2DbRow2 {
            id: id_value,
            name: row.get(1).unwrap_or("".to_string()),
            identifier: row.get(2).unwrap_or("".to_string()),
            adduct: adduct_final_value,
            cid: row.get(4).unwrap_or("".to_string()),
            window: row.get(5).unwrap_or("".to_string()),
            tof: row.get(6).unwrap_or("".to_string()),
            mz: row.get(7).unwrap_or("".to_string()),
            data: (Vec::new(), Vec::new()),
            matrix: row.get(9).unwrap_or("".to_string()),
            //data: row.get(6).unwrap_or(Vec::new()),
        })
    }).unwrap();

    let mut ids: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut identifiers: Vec<String> = Vec::new();
    let mut adducts: Vec<String> = Vec::new();
    let mut cids: Vec<String> = Vec::new();
    let mut windows: Vec<String> = Vec::new();
    let mut tofs: Vec<String> = Vec::new();
    let mut mzs: Vec<String> = Vec::new();
    let mut datas: Vec<(Vec<f64>, Vec<i64>)> = Vec::new();
    let mut matrices: Vec<String> = Vec::new();

    //parse results for passing back to the parent function
    for (index, item) in db_iter.enumerate() {
        let row: MS2DbRow2 = item.unwrap();
        ids.insert(index, row.id.to_string());
        names.insert(index, row.name);
        identifiers.insert(index, row.identifier);
        adducts.insert(index, row.adduct);
        cids.insert(index, row.cid);
        windows.insert(index, row.window);
        tofs.insert(index, row.tof);
        mzs.insert(index, row.mz);
        datas.insert(index, row.data);
        matrices.insert(index, row.matrix)
    }
    let mut identifiers2: Vec<String> = identifiers.clone();
    identifiers2.sort();
    identifiers2.dedup();

    
    let return_data: Vec<Vec<String>> = vec![ids, names, identifiers, adducts, cids, windows, tofs, mzs, matrices];
    return_data
    
}

pub fn remove_row_from_msms_user_db(rowid: usize) -> usize {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("get_msms").unwrap();
  let sql = "DELETE FROM user_MSMS WHERE id = ?1";
  conn.execute(sql, params![rowid]).unwrap();
  return 1
}