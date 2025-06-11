use log::info;
use once_cell::sync::OnceCell;
use r2d2::{Pool, PooledConnection};
use r2d2_sqlite::SqliteConnectionManager;

pub static POOL: OnceCell<Pool<SqliteConnectionManager>> = OnceCell::new();
pub static MSMS_POOL: OnceCell<Pool<SqliteConnectionManager>> = OnceCell::new();

pub fn get_connection() -> Result<PooledConnection<SqliteConnectionManager>, r2d2::Error> {
    // First, get the Pool instance from the OnceCell.
    info!("making an MS1 connection");
    let pool: &Pool<SqliteConnectionManager> = POOL.get().expect("Pool has not been initialized");
    // Then, get a connection from the pool.
    pool.get()
}

pub fn get_msms_connection(
    info_str: &str,
) -> Result<PooledConnection<SqliteConnectionManager>, r2d2::Error> {
    // First, get the Pool instance from the OnceCell.
    info!("making an MS2 connection - {:?}", info_str);
    let msms_pool: &Pool<SqliteConnectionManager> = MSMS_POOL.get().unwrap(); //.expect("Pool has not been initialized");

    // Then, get a connection from the pool.
    msms_pool.get()
}
