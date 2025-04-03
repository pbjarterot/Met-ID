// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod add_to_db;
mod binary_setup;
mod database_setup;
mod files;
mod install_helper_functions;
mod logging;
mod mass_match;
mod ms2_match;
mod multiprocessing;
mod regression;
mod sidecar;
mod splashscreen;
pub mod sql_mod;
pub mod updater;
#[cfg(test)]
mod testing;

use log::{error, info};
use logging::LOGGER;
use once_cell::sync::OnceCell;
use r2d2::Pool;
use r2d2_sqlite::SqliteConnectionManager;
use sql_mod::build_query::check_temp_tables;
use std::collections::HashMap;
use std::env;
use std::panic;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::Mutex;
use tauri::AppHandle;

pub static MSMSPATH: OnceCell<String> = OnceCell::new();
pub static MSDBPATH: OnceCell<String> = OnceCell::new();

pub static APP_HANDLE: OnceCell<AppHandle> = OnceCell::new();
fn set_app_handle(handle: AppHandle) {
    if APP_HANDLE.set(handle).is_err() {
        eprintln!("AppHandle is already set!");
    }
}
fn get_app_handle() -> Option<&'static AppHandle> {
    APP_HANDLE.get()
}

#[tauri::command]
fn is_backend_ready() -> bool {
    // Return true if setup is complete, otherwise return false.
    info!("checking backend ready");
    // Here's a very basic example. You may need more logic based on your actual setup.
    database_setup::MSMS_POOL.get().is_some() && database_setup::POOL.get().is_some()
}




#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    // Initialize logger
    log::set_logger(&*LOGGER).unwrap();
    log::set_max_level(log::LevelFilter::Info);

    // Set custom panic hook
    let default_hook = panic::take_hook();
    panic::set_hook(Box::new(move |info| {
        error!("Uncaught panic: {}", info);
        default_hook(info);
    }));
    
    let callback_map: updater::CallbackMap = Arc::new(Mutex::new(HashMap::new()));


    tauri::Builder::default()
        .plugin(tauri_plugin_updater::Builder::new().build())
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_fs::init())
        .manage(callback_map.clone())
        .setup(|app| {

            let handle = app.handle().clone();
            tauri::async_runtime::spawn(async move {
                updater::update(handle, callback_map).await.unwrap();               
            });
            
            set_app_handle(app.handle().clone());

            
            //let splashscreen_window: Window = app.get_window("splashscreen").unwrap();
            //let main_window: Window = app.get_window("main").unwrap();

            let db_path: PathBuf =
                install_helper_functions::ensure_database_in_appdata(&app, "db.db");
            let msms_db_path: PathBuf =
                install_helper_functions::ensure_database_in_appdata(&app, "msms_db.db");

            MSMSPATH
                .set(msms_db_path.to_string_lossy().to_string())
                .expect("Failed to set MSMS DB PATH");

            let pool: Pool<SqliteConnectionManager> =
                install_helper_functions::create_pool_from_app_path(&app, db_path);
            let msms_pool: Pool<SqliteConnectionManager> =
                install_helper_functions::create_pool_from_app_path(&app, msms_db_path);

            database_setup::POOL
                .set(pool)
                .expect("Failed to initialize POOL");
            database_setup::MSMS_POOL
                .set(msms_pool)
                .expect("Failed to initialize MSMS_POOL");
            info!(
                "setting the database pools \n{:?}\n{:?}",
                database_setup::POOL,
                database_setup::MSMS_POOL
            );

            check_temp_tables();

            // we perform the initialization code on a new task so the app doesn't freeze
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            splashscreen::close_splashscreen,
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
            add_to_db::matrix::add_matrix_to_db_rust,
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
            sql_mod::check_fg_duplicate_tauri,
            regression::mass_error_regression,
            is_backend_ready,
            updater::frontend_bool_response
        ])
        .plugin(tauri_plugin_dialog::init())
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
