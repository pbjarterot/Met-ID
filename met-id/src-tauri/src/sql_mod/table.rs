use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::{params, Result};

use crate::database_setup::{get_connection, get_msms_connection};

pub fn check_if_table_exists(source_table: &str, new_table: &str) -> Result<()> {
    extern crate rusqlite;
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();

    let mut stmt = conn.prepare("SELECT name FROM sqlite_master WHERE type='table' AND name=?1")?;
    let mut rows = stmt.query(params![new_table])?;
    if rows.next()?.is_some() {
        return Ok(());
    }

    let mut stmt = conn.prepare("SELECT sql FROM sqlite_master WHERE type='table' AND name=?1")?;
    let mut rows = stmt.query(params![source_table])?;
    if let Some(row) = rows.next()? {
        let schema: String = row.get(0)?;

        let new_schema = schema.replace(source_table, new_table);

        conn.execute(&new_schema, [])?;

        println!("Table {} has been created.", new_table);
    } else {
        println!("Source table {} does not exist.", source_table);
    }

    Ok(())
}

pub fn check_if_table_exists_msms(source_table: &str, new_table: &str) -> Result<()> {
    extern crate rusqlite;
    let conn: r2d2::PooledConnection<SqliteConnectionManager> = get_msms_connection("").unwrap();

    let mut stmt = conn.prepare("SELECT name FROM sqlite_master WHERE type='table' AND name=?1")?;
    let mut rows = stmt.query(params![new_table])?;
    if rows.next()?.is_some() {
        return Ok(());
    }

    let mut stmt = conn.prepare("SELECT sql FROM sqlite_master WHERE type='table' AND name=?1")?;
    let mut rows = stmt.query(params![source_table])?;
    if let Some(row) = rows.next()? {
        let schema: String = row.get(0)?;

        let new_schema = schema.replace(source_table, new_table);

        conn.execute(&new_schema, [])?;

        println!("Table {} has been created.", new_table);
    } else {
        println!("Source table {} does not exist.", source_table);
    }

    Ok(())
}
