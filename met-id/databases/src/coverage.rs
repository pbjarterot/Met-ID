use rusqlite::Connection;

pub fn coverage_string(conn: &Connection, fg_string: &str, matrix: &str) -> String {
  let sql: &str = &format!("SELECT numfunctionalgroups, maxcoverage FROM adducts WHERE mname = '{}'", matrix);
  let mut stmt = conn.prepare(sql).unwrap();
  let adducts_iter = stmt.query_map([], |row| {
    Ok((
      row.get::<_,i32>(0).unwrap(),
      row.get::<_,i32>(1).unwrap(),
    ))
  }).unwrap();

  let mut fg_vec: Vec<i32> = Vec::new();
  let mut cov_vec: Vec<i32> = Vec::new();

  for adduct in adducts_iter {
    match adduct {
      Ok((fg, cov)) => {
        if !fg_vec.contains(&fg) {
          fg_vec.push(fg);
          cov_vec.push(cov);
        }
      },
      Err(e) => println!("Error: {}", e)
    }
  }
  println!("fg_vec:  {:?}\ncov_vec: {:?}", fg_vec, cov_vec);

  let max_index = fg_vec.iter()
                         .enumerate()
                         .max_by_key(|&(_, value)| value)
                         .map(|(index, _)| index)
                         .unwrap();


  let mut sql = "CASE ".to_string();

  sql += &format!("WHEN {} > {} THEN {} ", fg_string, fg_vec[max_index], cov_vec[max_index]);

  for index in 0..fg_vec.len() {
    sql += &format!("WHEN {} = {} THEN {} ", fg_string, fg_vec[index], cov_vec[index])
  }
  sql += "END as coverage_value";
  println!("{:?}", sql);
  sql
}