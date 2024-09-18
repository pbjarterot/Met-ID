// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod files;
mod mass_match;
mod ms2_match;
mod regression;
pub mod sql_mod;
mod add_to_db;
mod multiprocessing;
mod logging;
mod install_helper_functions;
mod sidecar;
mod database_setup;
mod binary_setup;
mod splashscreen;

#[cfg(test)]
mod testing;

use std::path::PathBuf;
use database_setup::get_connection;
use log::{error, info};
use logging::LOGGER;
use once_cell::sync::OnceCell;
use std::panic;
use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;
use std::env;



pub static MSMSPATH: OnceCell<String> = OnceCell::new(); 


#[tauri::command]
fn is_backend_ready() -> bool {
    // Return true if setup is complete, otherwise return false.
    info!("checking backend ready");
    // Here's a very basic example. You may need more logic based on your actual setup.
    database_setup::MSMS_POOL.get().is_some() && database_setup::POOL.get().is_some()
}


fn main() {
    // Initialize logger
    log::set_logger(&*LOGGER).unwrap();
    log::set_max_level(log::LevelFilter::Info);

    // Set custom panic hook
    let default_hook = panic::take_hook();
    panic::set_hook(Box::new(move |info| {
        error!("Uncaught panic: {}", info);
        default_hook(info);
    }));

    tauri::Builder::default()
        .setup(|app| {
            //let splashscreen_window: Window = app.get_window("splashscreen").unwrap();
            //let main_window: Window = app.get_window("main").unwrap();

            let db_path: PathBuf = install_helper_functions::ensure_database_in_appdata(&app, "db.db");
            let msms_db_path: PathBuf = install_helper_functions::ensure_database_in_appdata(&app, "msms_db.db");

            
            MSMSPATH.set(msms_db_path.to_string_lossy().to_string()).expect("Failed to set MSMS DB PATH");
            //let _python_path: PathBuf = install_helper_functions::ensure_python_in_appdata(&app, "");
            /* 
            let binary_path = match std::env::consts::OS {
                "windows" => "metabolite-x86_64-pc-windows-msvc.exe",
                "macos" => {
                    match std::env::consts::ARCH {
                        "x86_64" => "metabolite_for_db-x86_64-apple-darwin",
                        "aarch64" => "metabolite_for_db-aarch64-apple-darwin",
                        _ => "",
                    }
                },
                "linux" => "metabolite_for_db-x86_64-unknown-linux-gnu",
                _ => "",
            };
            */
            //let metabolite_bin_path: PathBuf = install_helper_functions::ensure_bin_in_appdata(&app, binary_path);

            //binary_setup::METABOLITE_BIN_PATH.set(metabolite_bin_path);

            let pool: Pool<SqliteConnectionManager> = install_helper_functions::create_pool_from_app_path(&app, db_path);
            let msms_pool: Pool<SqliteConnectionManager> = install_helper_functions::create_pool_from_app_path(&app, msms_db_path);

            database_setup::POOL.set(pool).expect("Failed to initialize POOL");
            database_setup::MSMS_POOL.set(msms_pool).expect("Failed to initialize MSMS_POOL");
            info!("setting the database pools \n{:?}\n{:?}", database_setup::POOL, database_setup::MSMS_POOL);


            let mut conn = get_connection().unwrap();

            let tx = conn.transaction().unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_adducts AS
                SELECT * FROM adducts
                UNION ALL
                SELECT * FROM user_adducts", 
                [],
            ).unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_db_accessions AS
                SELECT * FROM db_accessions
                UNION ALL
                SELECT * FROM user_db_accessions", 
                [],
            ).unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_metabolites AS
                SELECT * FROM metabolites
                UNION ALL
                SELECT * FROM user_metabolites", 
                [],
            ).unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_derivatized_by AS
                SELECT * FROM derivatized_by
                UNION ALL
                SELECT * FROM user_derivatized_by", 
                [],
            ).unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_endogeneity AS
                SELECT * FROM endogeneity
                UNION ALL
                SELECT * FROM user_endogeneity", 
                [],
            ).unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_functional_groups AS
                SELECT * FROM functional_groups
                UNION ALL
                SELECT * FROM user_functional_groups", 
                [],
            ).unwrap();
            
            tx.execute(
                "CREATE TEMP TABLE temp_concat_in_tissue AS
                SELECT * FROM in_tissue
                UNION ALL
                SELECT * FROM user_in_tissue", 
                [],
            ).unwrap();
            
            tx.execute("CREATE INDEX idx_temp_concat_db_accessions_id ON temp_concat_db_accessions(id)", []).unwrap();
            
            tx.execute("CREATE INDEX idx_temp_concat_metabolites_id ON temp_concat_metabolites(id)", []).unwrap();
            
            tx.execute("CREATE INDEX idx_temp_concat_derivatized_by_id ON temp_concat_derivatized_by(id)", []).unwrap();
            
            tx.execute("CREATE INDEX idx_temp_concat_endogeneity_id ON temp_concat_endogeneity(id)", []).unwrap();
            tx.execute("CREATE INDEX idx_temp_concat_endogeneity_criteria ON temp_concat_endogeneity(endogenous, exogenous, unspecified)", []).unwrap();
            
            tx.execute("CREATE INDEX idx_temp_concat_functional_groups_id ON temp_concat_functional_groups(id)", []).unwrap();
            tx.execute("CREATE INDEX idx_temp_concat_functional_groups_criteria ON temp_concat_functional_groups('Phenolic Hydroxyls', 'Primary Amines')", []).unwrap();
            
            tx.execute("CREATE INDEX idx_temp_concat_in_tissue_id ON temp_concat_in_tissue(id)", []).unwrap();

            tx.commit().unwrap();
        
            // we perform the initialization code on a new task so the app doesn't freeze
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![splashscreen::close_splashscreen, 
                                                 sidecar::sidecar_function,
                                                 multiprocessing::my_command,
                                                 //files::load_msms, 
                                                 files::parse_ms1_csv, 
                                                 files::read_ms1_csv,
                                                 files::read_ms1_ctrlv,
                                                 files::read_mass_error_csv,
                                                 files::save_csv,
                                                 files::read_mzml_for_msms,
                                                 mass_match::calculate_adjusted_mass,
                                                 sql_mod::matrix_dropdown_tauri,
                                                 sql_mod::sql_handler_tauri,
                                                 sql_mod::sql_counter_tauri,
                                                 sql_mod::get_msms_spectra_tauri,
                                                 sql_mod::get_msms_tauri,
                                                 sql_mod::get_name_from_identifier_msms_tauri,
                                                 sql_mod::ms2_search_spectra_tauri,
                                                 sql_mod::match_msms_to_ui_tauri,
                                                 sql_mod::get_functional_groups_tauri,
                                                 sql_mod::get_matrices_tauri,
                                                 sql_mod::get_tissues_tauri,
                                                 sql_mod::add_msms_to_db_tauri,
                                                 sql_mod::show_user_msms_db_tauri,
                                                 sql_mod::remove_row_from_user_metabolites_tauri,
                                                 sql_mod::remove_row_from_msms_user_db_tauri,
                                                 sql_mod::remove_row_from_user_matrices_tauri,
                                                 sql_mod::remove_from_user_fgs_tauri,
                                                 sql_mod::update_user_metabolites_tauri,
                                                 sql_mod::update_user_matrices_tauri,
                                                 sql_mod::update_user_fgs_tauri,
                                                 sql_mod::db_data_tauri,
                                                 sql_mod::db_ids_and_names_tauri,
                                                 regression::mass_error_regression,
                                                 is_backend_ready
                                                 ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
