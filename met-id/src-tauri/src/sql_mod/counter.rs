use r2d2_sqlite::SqliteConnectionManager;
use crate::database_setup::get_connection;
use super::build_query::build_count_query;


pub fn sql_counter(met: String, matrix: String, typ: Vec<String>, adducts: Vec<String>) -> i64 {
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
    if typ.len() == 0 || adducts.len() == 0 {
        return 0;
    }
    let query_str: String = build_count_query(met, matrix, typ, adducts);
    println!("{:?}", query_str);
    let count: i64 = match conn.query_row(
        &query_str[..], 
        rusqlite::params![],
        |row: &rusqlite::Row<'_>| row.get(0),
    ) {
        Ok(num) => {
            num
        },
        Err(e) => {
            println!("Error: {:?}", e);
            0 as i64
        }
    };

    count
}