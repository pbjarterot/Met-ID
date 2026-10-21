use crate::{api::errors::{JsonParseError, APIResult}, database_setup::get_connection};
use axum::Json;
use base64::{Engine as _, engine::general_purpose};
use serde::{Serialize, Deserialize};
use rusqlite::{Result, types::Value};

#[derive(Debug, Deserialize)]
pub struct CustomSqlData (String);

#[derive(Debug, Deserialize)]
pub struct SearchByNameInputData {
    name: String,
    ccs: bool,
}

#[derive(Debug, Serialize)]
pub struct CustomSqlOutputData (Vec<String>);

#[derive(Debug, Serialize)]
pub struct SQLTable {
    col_names: Vec<String>,
    rows: Vec<Vec<JsonCell>>
}

#[derive(Debug, Serialize)]
//#[serde(tag = "kind", content = "data")]
pub enum SQLGet {
    Tables(Vec<String>),
    Table(SQLTable),
}


struct TableNames {
    pub table: String
}

#[derive(Debug, Serialize, Clone)]
#[serde(untagged)]
pub enum JsonCell {
    Null,
    Integer(i64),
    Real(f64),
    Text(String),
    Blob(String), // maybe base64
}

impl From<Value> for JsonCell {
    fn from(v: Value) -> Self {
        match v {
            Value::Null => JsonCell::Null,
            Value::Integer(i) => JsonCell::Integer(i),
            Value::Real(f) => JsonCell::Real(f),
            Value::Text(s) => JsonCell::Text(s),
            Value::Blob(b) => JsonCell::Blob(general_purpose::STANDARD.encode(b)),
        }
    }
}

fn get_table_columns(conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>, table_name: &str) -> Vec<String> {
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({});", table_name)).unwrap();
    let columns: Vec<String> = stmt
        .query_map([], |row| row.get::<_, String>(1)).unwrap() // Get column names
        .filter_map(|res| res.ok()) // Ignore errors
        .filter(|col| col != "id")
        .collect();
    columns.iter().map(|item| format!("\"{}\"", item)).collect()
}

pub fn fetch_table(
    conn: &r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager>,
    sql: &str
) -> Result<SQLTable> {
    let mut stmt = conn.prepare(sql)?;

    let col_count = stmt.column_count();
    let col_names: Vec<String> = stmt.column_names().iter().map(|s| s.to_string()).collect();

    let rows_iter = stmt.query_map([], move |row| {
        let mut out = Vec::with_capacity(col_count);
        for i in 0..col_count {
            let v: Value = row.get(i)?;
            out.push(JsonCell::from(v));
        }
        Ok(out)
    })?;

    let mut rows = Vec::new();
    for r in rows_iter {
        rows.push(r?);
    }

    Ok(SQLTable { col_names, rows })
}


/// TODO Add search by name for metabolites to get info on functional groups, masses etc like the database tab
/// ALSO search by chemical formula

pub async fn search_by_name(Json(sql): Json<SearchByNameInputData>) -> APIResult<Json<SQLGet>> {
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> = get_connection().unwrap();

    let adducts_table: SQLTable = fetch_table(&conn, "SELECT adduct, mname FROM adducts").unwrap();
    let adducts: Vec<JsonCell> = adducts_table.rows.iter().map(|row| row[0].clone()).collect::<_>();

    let adducts_string = if sql.ccs {
        adducts
        .iter()
        .filter_map(|cell| {
            if let JsonCell::Text(s) = cell {
                Some(format!(
                    "MAX(CASE WHEN a.adduct = '{}' THEN m.mz + a.deltamass END) AS '{} m/z',
                     MAX(CASE WHEN metabolites_ccs.adduct = '{}' THEN metabolites_ccs.CCS END) AS '{} CCS'",
                    s.replace('\'', "''"), s,
                    s.replace('\'', "''"), s,
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
    } else{
        adducts
        .iter()
        .filter_map(|cell| {
            if let JsonCell::Text(s) = cell {
                Some(format!(
                    "MAX(CASE WHEN a.adduct = '{}' THEN m.mz + a.deltamass END) AS '{} m/z'",
                    s.replace('\'', "''"), s,
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
    };

    let functional_groups: Vec<String> = get_table_columns(&conn, "functional_groups");
    let functional_groups_string = functional_groups
        .iter()
        .map(|fg| format!("functional_groups.{fg}"))
        .collect::<Vec<_>>()
        .join(", ");
    
    let sql_string: &String = if sql.ccs {
        &format!(
        "SELECT m.name, m.mz, m.chemicalformula, \
        {}, {}\
        FROM metabolites AS m \
        INNER JOIN functional_groups ON m.id = functional_groups.id \
        INNER JOIN metabolites_ccs ON m.id = metabolites_ccs.metabolite_id \
        CROSS JOIN adducts a \
        WHERE name LIKE '%{}%' AND (a.adduct = metabolites_ccs.adduct) \
        GROUP BY m.id, m.name, m.mz, m.chemicalformula \
        ", functional_groups_string, adducts_string, sql.name.as_str())
    } else {
        &format!("SELECT m.name, m.mz, m.chemicalformula, \
        {}, {} \
        FROM metabolites AS m \
        INNER JOIN functional_groups ON m.id = functional_groups.id \
        CROSS JOIN adducts a \
        WHERE name LIKE '%{}%' \
        GROUP BY m.id, m.name, m.mz, m.chemicalformula \
        ", functional_groups_string, adducts_string, sql.name.as_str())
    };

    let res: SQLTable = fetch_table(&conn, sql_string).unwrap();
    return Ok(Json(SQLGet::Table(res)))
}

pub async fn get_database(Json(sql): Json<CustomSqlData>) -> APIResult<Json<SQLGet>> {
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> = get_connection().unwrap();
    let mut res: Vec<String> = Vec::new();

    let sql_string: &str = match sql.0.as_str() {
        "tables" => "SELECT name FROM sqlite_master WHERE type='table';",
        _ => &format!("SELECT * FROM {}", sql.0)
    };

    if sql.0.as_str() == "tables" {
        let mut stmt = match conn.prepare(sql_string) {
            Ok(stmt) => stmt,
            Err(e) => {
                eprintln!("Error: {}", e);
                return Ok(Json(SQLGet::Tables(res)));
            }
        };

        let db_iter = match stmt.query_map([], |row| {
            Ok(TableNames {
                table: row.get(0)?,
            })
        }) {
            Ok(iter) => iter, 
            Err(e) => {
                println!("Error: {}", e);
                return Ok(Json(SQLGet::Tables(res)));
            }
        };

        for item in db_iter {
            match item {
                Ok(row) => {
                    res.push(row.table)
                }
                Err(e) => println!("Error processing row: {}", e),
            }
        }
    }
    
    if sql.0.as_str() != "tables" {
        let res: SQLTable  = fetch_table(&conn, sql_string).unwrap();
        return Ok(Json(SQLGet::Table(res)))
    }

    Ok(Json(SQLGet::Tables(res)))
}

pub async fn post_database() {
    let _conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> = get_connection().unwrap();
    
    //functions to add matrices, functional groups, metabolites


}

pub async fn custom_sql(Json(sql): Json<CustomSqlData>) 
    -> Result<Json<CustomSqlOutputData>, JsonParseError>
{
    let conn: r2d2::PooledConnection<r2d2_sqlite::SqliteConnectionManager> = get_connection().unwrap();

    let _res = conn.execute(&sql.0, []).unwrap();
    let res: Vec<String> = Vec::new();
    Ok(Json(CustomSqlOutputData(res)))
}