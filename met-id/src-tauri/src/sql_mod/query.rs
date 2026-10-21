use super::MS1DbRow;
use crate::database_setup::get_connection;

#[derive(Debug)]
pub struct QueryResults {
    pub db_masses: Vec<f64>,
    pub db_names: Vec<String>,
    pub db_mnames: Vec<String>,
    pub db_accessions: Vec<String>,
    pub db_smiles: Vec<String>,
    pub db_formulas: Vec<String>,
    pub db_pos_derivs: Vec<i32>, 
    pub db_ccs: Vec<f64>,
}

impl QueryResults {
    fn new() -> Self {
        QueryResults {
            db_masses: Vec::new(),
            db_names: Vec::new(),
            db_mnames: Vec::new(),
            db_accessions: Vec::new(),
            db_smiles: Vec::new(),
            db_formulas: Vec::new(),
            db_pos_derivs: Vec::new(),
            db_ccs: Vec::new(),
        }
    }
}

pub fn sql_query(
    query: &str
) -> QueryResults {
    // Connect to db
    let mut res: QueryResults = QueryResults::new();
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> = get_connection().unwrap();
    
    // Query
    let mut stmt = match conn.prepare(query) {
        Ok(stmt) => stmt,
        Err(e) => {
            eprintln!("Error: {}", e);
            return res;
        }
    };

    let db_iter = match stmt.query_map([], |row| {
        Ok(MS1DbRow {
            mz: row.get(0)?,
            name: row.get(1)?,
            mname: row.get(2)?,
            accession: row.get(3)?,
            smiles: row.get(4)?,
            formula: row.get(5)?,
            possible_derivs: row.get(6)?,
            ccs: row.get::<_, Option<f64>>(7)?.unwrap_or(0.0),
        })
    }) {
        Ok(iter) => iter,
        Err(e) => {
            eprintln!("Error3: {}", e);
            return res;
        }
    };

    // Parse results for passing back to the parent function
    for item in db_iter {
        match item {
            Ok(row) => {
                res.db_masses.push(row.mz);
                res.db_names.push(row.name);
                res.db_mnames.push(row.mname);
                res.db_accessions.push(row.accession);
                res.db_smiles.push(row.smiles);
                res.db_formulas.push(row.formula);
                res.db_pos_derivs.push(row.possible_derivs);
                res.db_ccs.push(row.ccs);
            }
            Err(e) => eprintln!("Error processing row: {}", e),
        }
    }
    res
}
