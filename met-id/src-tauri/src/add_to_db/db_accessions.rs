use r2d2_sqlite::SqliteConnectionManager;
use rusqlite::params;

pub fn fill_user_db_accessions(conn: &r2d2::PooledConnection<SqliteConnectionManager>) {
    conn.execute(
        "INSERT INTO user_db_accessions (hmdb) VALUES (?1)",
        params!["placeholder"], // using 0 as a placeholder
    )
    .unwrap();

    // Get the last inserted row's ID
    let last_id = conn.last_insert_rowid();

    // Create the special_id based on the last inserted row's ID
    let special_id = format!("METID{:05}", last_id);

    // Update the newly inserted row to set the special_id
    conn.execute(
        "UPDATE user_db_accessions SET hmdb = ?1 WHERE id = ?2",
        params![special_id, last_id],
    )
    .unwrap();
}
