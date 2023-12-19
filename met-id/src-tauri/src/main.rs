// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod files;
mod sql;
mod mass_match;
mod ms2_match;
mod regression;
pub mod database;
mod validation;
mod sql_build_query;
mod add_to_db;
mod multiprocessing;
mod logging;
mod install_helper_functions;
mod sidecar;
mod database_setup;
mod splashscreen;

#[cfg(test)]
mod testing;

use std::path::PathBuf;
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

            let pool: Pool<SqliteConnectionManager> = install_helper_functions::create_pool_from_app_path(&app, db_path);
            let msms_pool: Pool<SqliteConnectionManager> = install_helper_functions::create_pool_from_app_path(&app, msms_db_path);

            database_setup::POOL.set(pool).expect("Failed to initialize POOL");
            database_setup::MSMS_POOL.set(msms_pool).expect("Failed to initialize MSMS_POOL");
            info!("setting the database pools \n{:?}\n{:?}", database_setup::POOL, database_setup::MSMS_POOL);
        
            // we perform the initialization code on a new task so the app doesn't freeze
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![splashscreen::close_splashscreen, 
                                                 sidecar::sidecar_function,
                                                 multiprocessing::my_command,
                                                 files::load_msms, 
                                                 files::parse_ms1_csv, 
                                                 files::read_ms1_csv,
                                                 files::read_ms1_ctrlv,
                                                 files::read_mass_error_csv,
                                                 files::save_csv,
                                                 files::read_mzml_for_msms,
                                                 sql::sql_handler,
                                                 sql::sql_counter,
                                                 sql::get_msms,
                                                 sql::get_msms_spectra,
                                                 sql::get_name_from_identifier_msms,
                                                 //sql::find_msms_fragments,
                                                 sql::ms2_search_spectra,
                                                 sql::match_msms_to_ui,
                                                 sql::get_functional_groups,
                                                 sql::get_matrices,
                                                 sql::get_tissues,
                                                 sql::add_msms_to_db,
                                                 sql::show_user_msms_db,
                                                 sql::remove_row_from_msms_user_db,
                                                 regression::mass_error_regression,
                                                 //validation::check_smiles,
                                                 //validation::check_smarts,
                                                 //mass_match::calculate_adjusted_mass,
                                                 is_backend_ready
                                                 ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
