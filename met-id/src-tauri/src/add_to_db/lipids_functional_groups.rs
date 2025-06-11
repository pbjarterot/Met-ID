use rusqlite::Result;
use std::collections::HashMap;
use std::sync::mpsc;

use crate::add_to_db::functional_group::functional_group_target;

pub fn complete_lipids_functional_groups(map: HashMap<String, String>) -> Result<()> {
    println!("---------------------- Updating ---------------------------");
    let (tx, _rx) = mpsc::channel();
    for (name, smart) in map.iter() {
        functional_group_target(
            tx.clone(),
            smart.clone(),
            "lipids_functional_groups".to_string(),
            name.clone(),
            "lipids".to_string(),
        )
        .unwrap();
    }
    println!("---------------------- Updating ---------------------------");
    Ok(())
}

/*
    let tx = conn.transaction()?;

    // Get columns (excluding "id")
    let mut stmt = tx.prepare(&format!("PRAGMA table_info({})", table))?;
    let columns: Vec<String> = stmt
        .query_map([], |row| row.get::<_, String>(1))?
        .filter_map(|res| res.ok())
        .filter(|name| name != "id")  // assuming id is the primary key
        .collect();

    for col in &columns {
        let sql = format!(
            "UPDATE {table}
                SET {col} = (
                SELECT b.{col} FROM {ref_table} b WHERE b.id = {table}.id
            )
            WHERE ({col} IS NULL OR {col} = '')
            AND EXISTS (
                SELECT 1 FROM {ref_table} b WHERE b.id = {table}.id AND b.{col} IS NOT NULL AND b.{col} != ''
            )",
            table = table,
            ref_table = ref_table,
            col = col
        );

        tx.execute(&sql, [])?;
    }

    tx.commit()
}
*/
