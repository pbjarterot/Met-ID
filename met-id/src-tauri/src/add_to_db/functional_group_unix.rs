use rusqlite::{params, Result, ToSql};
use rdkit::ROMol;
use rdkit::RWMol;
use rdkit::substruct_match;
use rdkit::SubstructMatchParameters;

pub fn functional_group_unix(
    progress_sender: std::sync::mpsc::Sender<f32>, 
    smarts: String, 
    table_name: String, 
    column_name: String) 
    -> Result<()> 
    {
        let mut conn: r2d2::PooledConnection<SqliteConnectionManager> = get_connection().unwrap();
        let progress_sender: std::sync::mpsc::Sender<f32> = self.progress_sender.clone();
        let transaction = conn.transaction().unwrap();

        //retrieves smiles in batches of 1000
        let smiles: Vec<Vec<String>> = get_smiles_from_db(/*Maybe add argument to set batch size */);
        let smarts_mol: ROMol = RWMol::from_smarts(&smarts).unwrap.to_ro_mol();
        let params: SubstructMatchParameters = SubstructMatchParameters::new();

        for (batch_idx, batch) in smiles.iter().enumerate() {
            for (smiles_idx, smile) in batch.iter().enumerate() {
                let mol: ROMol = ROMol::from_smiles(&smile).unwrap();
                let res: usize = substruct_march(&mol, &query, &params).len();
                let new_idx: usize = (batch_idx * 1_000) + idx;
                transaction.execute(&format!("UPDATE '{}' SET '{}' = ?1 WHERE rowid = ?2", table_name, column_name), params![res, new_idx]).unwrap();
            }
            progress_sender.send(batch_idx / smiles.len()).unwrap();
        }
        transaction.commit().unwrap();
        Ok(())
    }