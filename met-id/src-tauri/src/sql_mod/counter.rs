use r2d2_sqlite::SqliteConnectionManager;
use crate::database_setup::get_connection;

use super::build_query::build_count_query;


pub fn sql_counter(met: String, mat: String, typ: Vec<String>, adducts: Vec<String>) -> i64 {
  let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    if typ.len() == 0 || adducts.len() == 0 {
        return 0;
    }
    let query_str: String = build_count_query(met, mat, typ, adducts);

    let count: i64 = conn.query_row(
        &query_str[..], 
        rusqlite::params![],
        |row: &rusqlite::Row<'_>| row.get(0),
    ).unwrap();

    count
}