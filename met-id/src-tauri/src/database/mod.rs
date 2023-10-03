/* 
pub mod ms1;
pub mod ms2;


use lazy_static::lazy_static;
use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;


lazy_static! {
    static ref POOL: Pool<SqliteConnectionManager> = {
        let manager = SqliteConnectionManager::file("../databases/db.db");
        Pool::new(manager).expect("Failed to create pool.")
    };

    static ref MSMS_POOL: Pool<SqliteConnectionManager> = {
        let manager = SqliteConnectionManager::file("../MetID_parse/msms_db.db");
        Pool::new(manager).expect("Failed to create pool.")
    };
}

#[derive(Debug)]
struct IDentifierMap {
    name: String,
    hmdb: String,
}

*/