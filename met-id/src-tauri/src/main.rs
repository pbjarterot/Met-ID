// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use logging::LOGGER;
use tauri::{App, Manager, Window};

mod files;
mod sql;
mod mass_match;
mod regression;
pub mod database;
mod validation;
mod sql_build_query;
mod add_to_db;
mod metabolite;
mod multiprocessing;
mod logging;
mod install_helper_functions;

use std::path::PathBuf;
use log::{error, info};
use std::panic;
use r2d2::{Pool, PooledConnection};
use r2d2_sqlite::SqliteConnectionManager;
use once_cell::sync::OnceCell;
use std::env;

pub static POOL: OnceCell<Pool<SqliteConnectionManager>> = OnceCell::new();
pub static MSMS_POOL: OnceCell<Pool<SqliteConnectionManager>> = OnceCell::new();

pub fn get_connection() -> Result<PooledConnection<SqliteConnectionManager>, r2d2::Error> {
    // First, get the Pool instance from the OnceCell.
    info!("making an MS1 connection");
    let pool: &Pool<SqliteConnectionManager> = POOL.get().expect("Pool has not been initialized");
    // Then, get a connection from the pool.
    pool.get()
}

pub fn get_msms_connection(info_str: &str) -> Result<PooledConnection<SqliteConnectionManager>, r2d2::Error> {
    // First, get the Pool instance from the OnceCell.
    info!("making an MS2 connection - {:?}", info_str);
    let msms_pool: &Pool<SqliteConnectionManager> = MSMS_POOL.get().unwrap();//.expect("Pool has not been initialized");
    
    // Then, get a connection from the pool.
    msms_pool.get()
}


#[tauri::command]
async fn close_splashscreen(window: Window) {
    // Close splashscreen
    window.get_window("splashscreen").expect("no window labeled 'splashscreen' found").close().unwrap();
    // Show main window
    window.get_window("main").expect("no window labeled 'main' found").show().unwrap();
}

fn create_pool_from_app_path(_app: &App, resource_path: PathBuf) -> Pool<SqliteConnectionManager> {
    info!("looking in: {:?}", resource_path);
    let manager = SqliteConnectionManager::file(resource_path);
    Pool::new(manager).expect("Failed to create pool.")
}

fn get_executable_path() -> std::path::PathBuf {
    env::current_exe().expect("Failed to get current executable path")
}




#[derive(serde::Deserialize)]
struct IsBackendReady;

#[tauri::command]
fn is_backend_ready(_window: tauri::Window, _: IsBackendReady) -> bool {
    // Return true if setup is complete, otherwise return false.
    info!("checking backend ready");
    // Here's a very basic example. You may need more logic based on your actual setup.
    MSMS_POOL.get().is_some() && POOL.get().is_some()
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
            /* 
            let python_resource_path = app.path_resolver()
                .resolve_resource("./bin/python3.dll")
                .expect("failed to resolve resource");

            info!("looking in: {:?}", python_resource_path);
            let output = Command::new(python_resource_path)
                .output()
                .expect("Failed to execute Python binary.");
            */

            //let splashscreen_window: Window = app.get_window("splashscreen").unwrap();
            //let main_window: Window = app.get_window("main").unwrap();
            

            info!("Initializing");

            let db_path: PathBuf = install_helper_functions::ensure_database_in_appdata(&app, "db.db");
            let msms_db_path: PathBuf = install_helper_functions::ensure_database_in_appdata(&app, "msms_db.db");
            let _python_path: PathBuf = install_helper_functions::ensure_python_in_appdata(&app, "");

            let pool: Pool<SqliteConnectionManager> = create_pool_from_app_path(&app, db_path);
            let msms_pool: Pool<SqliteConnectionManager> = create_pool_from_app_path(&app, msms_db_path);

            POOL.set(pool).expect("Failed to initialize POOL");
            MSMS_POOL.set(msms_pool).expect("Failed to initialize MSMS_POOL");
            info!("setting the database pools \n{:?}\n{:?}", POOL, MSMS_POOL);
        
            // we perform the initialization code on a new task so the app doesn't freeze
            /* 
            tauri::async_runtime::spawn(async move {
                // initialize your app here instead of sleeping :)
                println!("Initializing...");
                std::thread::sleep(std::time::Duration::from_secs(4));
                println!("Done initializing.");

        
                // After it's done, close the splashscreen and display the main window
                splashscreen_window.close().unwrap();
                main_window.show().unwrap();
            });
            */
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![close_splashscreen, 
                                                 multiprocessing::my_command,
                                                 files::load_msms, 
                                                 files::parse_ms1_csv, 
                                                 files::read_ms1_csv,
                                                 files::read_ms1_ctrlv,
                                                 sql::sql_handler,
                                                 sql::sql_counter,
                                                 sql::get_msms,
                                                 sql::get_msms_spectra,
                                                 sql::get_name_from_identifier_msms,
                                                 sql::find_msms_fragments,
                                                 sql::ms2_search_spectra,
                                                 sql::get_functional_groups,
                                                 sql::get_matrices,
                                                 sql::get_tissues,
                                                 regression::mass_error_regression,
                                                 validation::check_smiles,
                                                 validation::check_smarts,
                                                 is_backend_ready
                                                 ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
