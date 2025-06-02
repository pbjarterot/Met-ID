use super::MS1DbRow;
use crate::database_setup::get_connection;

pub fn sql_query(
    query: &String,
) -> (
    Vec<f64>,
    Vec<String>,
    Vec<String>,
    Vec<String>,
    Vec<String>,
    Vec<String>,
    Vec<i32>,
) {
    // Connect to db
    let conn = get_connection().unwrap();
    // Query
    println!("QUERY:::::::{:?}", query);
    let mut stmt = match conn.prepare(query) {
        Ok(stmt) => stmt,
        Err(e) => {
            println!("Error2: {}", e);
            return (
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            );
        }
    };

    let db_iter = match stmt.query_map([], |row| {
        Ok(MS1DbRow {
            mz: row.get(0).unwrap_or(0.0),
            name: row.get(1).unwrap_or_else(|_| "".to_string()),
            mname: row.get(2).unwrap_or_else(|_| "".to_string()),
            accession: row.get(3).unwrap_or_else(|_| "".to_string()),
            smiles: row.get(4).unwrap_or_else(|_| "".to_string()),
            formula: row.get(5).unwrap_or_else(|_| "".to_string()),
            possible_derivs: row.get(6).unwrap_or(0),
        })
    }) {
        Ok(iter) => iter,
        Err(e) => {
            println!("Error3: {}", e);
            return (
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            );
        }
    };

    let mut mzs = Vec::new();
    let mut names = Vec::new();
    let mut mnames = Vec::new();
    let mut accessions = Vec::new();
    let mut smiless = Vec::new();
    let mut formulas = Vec::new();
    let mut possible_derivs = Vec::new();

    // Parse results for passing back to the parent function
    for item in db_iter {
        match item {
            Ok(row) => {
                mzs.push(row.mz);
                names.push(row.name);
                mnames.push(row.mname);
                accessions.push(row.accession);
                smiless.push(row.smiles);
                formulas.push(row.formula);
                possible_derivs.push(row.possible_derivs);
            }
            Err(e) => println!("Error processing row: {}", e),
        }
    }

    (
        mzs,
        names,
        mnames,
        accessions,
        smiless,
        formulas,
        possible_derivs,
    )
}
