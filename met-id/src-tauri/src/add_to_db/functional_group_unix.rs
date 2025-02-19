#[cfg(target_family = "unix")]
pub fn functional_group_unix(
    progress_sender: std::sync::mpsc::Sender<f32>, 
    smarts: String, 
    table_name: String, 
    column_name: String) 
    -> Result<()> 
    {
        use r2d2_sqlite::SqliteConnectionManager;
        use rusqlite::{params, Result, ToSql};
        use crate::database_setup::get_connection;
        use super::functional_group::get_smiles_from_db;
        use rdkit::{ROMol, RWMol, substruct_match, SubstructMatchParameters};

        let mut conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
        let transaction = conn.transaction().unwrap();

        //retrieves smiles in batches of 1000
        let smiles: Vec<Vec<String>> = get_smiles_from_db(/*Maybe add argument to set batch size */);
        let smarts_mol: ROMol = RWMol::from_smarts(&smarts).unwrap.to_ro_mol();
        let params: SubstructMatchParameters = SubstructMatchParameters::new();

        for (batch_idx, batch) in smiles.iter().enumerate() {
            for (smiles_idx, smile) in batch.iter().enumerate() {
                let mol: ROMol = ROMol::from_smiles(&smile).unwrap();
                let res: usize = substruct_match(&mol, &smarts_mol, &params).len();
                let new_idx: usize = (batch_idx * 1_000) + smiles_idx;
                transaction.execute(&format!("UPDATE '{}' SET '{}' = ?1 WHERE rowid = ?2", table_name, column_name), params![res, new_idx]).unwrap();
            }
            progress_sender.send(batch_idx as f32 / smiles.len() as f32).unwrap();
        }
        transaction.commit().unwrap();
        Ok(())
    }